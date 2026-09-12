#!/usr/bin/env python3
"""Run AutoDock Vina in parallel over a prepared ligand library."""

from __future__ import annotations

import csv
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from .config import ResolvedConfig
from .utils import (
    DockingError,
    ToolNotFoundError,
    find_tool,
    parse_vina_affinities,
    run_command,
    safe_name,
    tail,
    write_json,
)

RESULT_FIELDS = [
    "id",
    "replicate",
    "seed",
    "smiles",
    "affinity",
    "mode",
    "rmsd_lb",
    "rmsd_ub",
    "pose_file",
    "status",
    "error",
    "wall_seconds",
]
_CSV_LOCK = threading.Lock()


def run_docking(cfg: ResolvedConfig, log):
    receptor = cfg.receptor_output()
    if not receptor.exists():
        raise DockingError(f"receptor PDBQT not found: {receptor}")
    manifest_path = cfg.manifest_path()
    if not manifest_path.exists():
        raise DockingError(
            f"ligand manifest not found: {manifest_path}; run prepare-ligands first"
        )

    docked_dir = cfg.docked_dir()
    docked_dir.mkdir(parents=True, exist_ok=True)
    results_path = cfg.results_path()

    vina = find_tool(cfg.get("docking", "executable", "vina"))
    if not vina:
        raise ToolNotFoundError(
            "AutoDock Vina not found; set docking.executable in the config "
            "or install autodock-vina"
        )

    rows = list(_read_csv(manifest_path))
    prepared = [
        row for row in rows
        if row.get("pdbqt") and Path(row["pdbqt"]).exists()
    ]
    seeds = cfg.get("docking", "seeds") or []
    if not isinstance(seeds, (list, tuple)):
        seeds = [seeds]
    seeds = [int(seed) for seed in seeds if str(seed).strip()]
    if not seeds:
        seeds = [int(cfg.get("docking", "seed", 42))]
    replicates = max(1, int(cfg.get("docking", "replicates", 1)))
    if replicates > len(seeds):
        seeds = [seeds[index % len(seeds)] for index in range(replicates)]
    tasks = []
    for row in prepared:
        for replicate in range(1, replicates + 1):
            task = dict(row)
            task["replicate"] = replicate
            task["seed"] = int(seeds[(replicate - 1) % len(seeds)])
            tasks.append(task)
    if not tasks:
        raise DockingError("no prepared PDBQT ligands found in the manifest")
    if not cfg.get("docking", "resume", True) and results_path.exists():
        results_path.unlink()

    done: set[tuple[str, str]] = set()
    if cfg.get("docking", "resume", True) and results_path.exists():
        done = {
            (
                str(row.get("id", "")),
                str(row.get("replicate", "1") or "1"),
            )
            for row in _read_csv(results_path)
            if row.get("status") == "ok"
        }
        if done:
            log.info("resuming: %s ligands already docked", len(done))

    max_workers = int(cfg.get("docking", "max_workers", 4))
    pending = [
        row for row in tasks
        if (
            str(row.get("id", "")),
            str(row.get("replicate", "1")),
        ) not in done
    ]
    log.info(
        "docking %s ligands with %s workers",
        len(pending),
        max_workers,
    )

    if pending:
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(_dock_one, cfg, vina, row, docked_dir, log): row
                for row in pending
            }
            for future in as_completed(futures):
                row = futures[future]
                try:
                    record = future.result()
                except Exception as exc:
                    record = _record(row, status="error", error=str(exc))
                _append_row(results_path, record)

    positive_control = _run_positive_control(
        cfg,
        vina,
        docked_dir,
        log,
    )

    all_results = list(_read_csv(results_path)) if results_path.exists() else []
    ok_results = [r for r in all_results if r.get("status") == "ok"]
    affinities = [float(r["affinity"]) for r in ok_results if r.get("affinity")]
    summary = {
        "total": len(all_results),
        "ok": len(ok_results),
        "failed": sum(1 for r in all_results if r.get("status") != "ok"),
        "replicates_per_ligand": replicates,
        "seeds": seeds,
        "positive_control": positive_control,
        "best_affinity": min(affinities) if affinities else None,
        "results_csv": str(results_path),
        "docked_dir": str(docked_dir),
    }
    write_json(docked_dir / "summary.json", summary)
    log.info(
        "docking complete: %s ok, best affinity %s",
        len(ok_results),
        summary["best_affinity"],
    )
    return summary


def _run_positive_control(
    cfg: ResolvedConfig,
    vina: str,
    docked_dir: Path,
    log,
) -> dict:
    """Dock an optional known ligand and compare it with a preset threshold."""
    value = cfg.get("docking", "positive_control_pdbqt")
    if not value:
        return {"status": "not_configured"}
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = cfg.workdir / path
    if not path.exists():
        return {"status": "failed", "reason": f"control ligand not found: {path}"}
    row = {
        "id": cfg.get("docking", "positive_control_id", "positive_control"),
        "smiles": "",
        "pdbqt": str(path),
        "replicate": 1,
        "seed": cfg.get("docking", "seed", 42),
    }
    try:
        record = _dock_one(cfg, vina, row, docked_dir, log)
    except Exception as exc:
        record = _record(row, status="error", error=str(exc))
    threshold = cfg.get("docking", "positive_control_max_affinity")
    affinity = record.get("affinity", "")
    passed = None
    if threshold not in (None, "") and affinity not in ("", None):
        try:
            passed = float(affinity) <= float(threshold)
        except (TypeError, ValueError):
            passed = None
    result = {
        **record,
        "threshold": threshold,
        "passed": passed,
    }
    pd.DataFrame([result]).to_csv(
        docked_dir / "positive_control_results.csv",
        index=False,
    )
    if passed is False:
        log.warning(
            "positive control docked at %s, worse than threshold %s",
            affinity,
            threshold,
        )
    else:
        log.info("positive control result: %s", result)
    return result


def _dock_one(cfg: ResolvedConfig, vina: str, row: dict, docked_dir: Path, log):
    lig_id = safe_name(row.get("id"), "ligand")
    replicate = int(row.get("replicate", 1) or 1)
    out_path = docked_dir / f"{lig_id}_rep{replicate}.pdbqt"
    cmd = build_vina_command(
        cfg,
        vina,
        str(cfg.receptor_output()),
        str(Path(row["pdbqt"])),
        str(out_path),
        seed=row.get("seed"),
    )
    started = time.time()
    try:
        result = run_command(
            cmd,
            timeout=int(cfg.get("docking", "timeout_seconds", 600)),
        )
    except DockingError as exc:
        return _record(row, status="error", error=str(exc))
    elapsed = round(time.time() - started, 2)
    if result.returncode != 0:
        return _record(
            row,
            status="error",
            error=tail(result.stderr or result.stdout),
            wall=elapsed,
        )
    modes = parse_vina_affinities(result.stdout)
    if not modes:
        return _record(
            row,
            status="no_pose",
            error=tail(result.stdout or result.stderr),
            wall=elapsed,
        )
    best = min(modes, key=lambda item: item["affinity"])
    return {
        "id": lig_id,
        "replicate": replicate,
        "seed": row.get("seed", ""),
        "smiles": row.get("smiles", ""),
        "affinity": best["affinity"],
        "mode": best.get("mode", 1),
        "rmsd_lb": best.get("rmsd_lb", ""),
        "rmsd_ub": best.get("rmsd_ub", ""),
        "pose_file": str(out_path),
        "status": "ok",
        "error": "",
        "wall_seconds": elapsed,
    }


def _record(row: dict, status: str, error: str, wall: float = 0.0) -> dict:
    return {
        "id": safe_name(row.get("id"), "ligand"),
        "replicate": row.get("replicate", 1),
        "seed": row.get("seed", ""),
        "smiles": row.get("smiles", ""),
        "affinity": "",
        "mode": "",
        "rmsd_lb": "",
        "rmsd_ub": "",
        "pose_file": "",
        "status": status,
        "error": error,
        "wall_seconds": wall,
    }


def build_vina_command(
    cfg: ResolvedConfig,
    vina: str,
    receptor: str,
    ligand: str,
    out: str,
    seed: int | str | None = None,
) -> list[str]:
    if vina.lower().endswith(".py"):
        cmd = [sys.executable, vina]
    else:
        cmd = [vina]
    center = cfg.receptor_center()
    size = cfg.receptor_size()
    cmd += [
        "--receptor", receptor,
        "--ligand", ligand,
        "--out", out,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(int(cfg.get("docking", "exhaustiveness", 8))),
        "--num_modes", str(int(cfg.get("docking", "num_modes", 9))),
        "--energy_range", str(float(cfg.get("docking", "energy_range", 3.0))),
        "--cpu", "1",
        "--seed", str(
            int(seed)
            if seed not in (None, "")
            else int(cfg.get("docking", "seed", 42))
        ),
    ]
    scoring = cfg.get("docking", "scoring", "")
    if scoring:
        cmd += ["--scoring", str(scoring)]
    for flex in cfg.receptor_flexible():
        cmd += ["--flex", str(flex)]
    return cmd


def _read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _append_row(path: Path, row: dict) -> None:
    with _CSV_LOCK:
        is_new = not path.exists()
        with path.open("a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=RESULT_FIELDS)
            if is_new:
                writer.writeheader()
            writer.writerow({key: row.get(key, "") for key in RESULT_FIELDS})
