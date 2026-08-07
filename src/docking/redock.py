#!/usr/bin/env python3
"""Re-dock top-ranked hits with higher exhaustiveness."""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .config import ResolvedConfig
from .docking import _append_row, _dock_one, _read_csv
from .utils import DockingError, ToolNotFoundError, find_tool, write_json


def run_redock(cfg: ResolvedConfig, log) -> dict:
    ranked = cfg.reports_dir() / "ranked_results.csv"
    if not ranked.exists():
        raise DockingError(f"ranked results not found: {ranked}; run analyze first")
    manifest = cfg.manifest_path()
    if not manifest.exists():
        raise DockingError(f"ligand manifest not found: {manifest}")

    import pandas as pd

    top_n = int(cfg.get("redock", "top_n", 20))
    frame = pd.read_csv(ranked, dtype={"id": str}).head(top_n)
    top_ids = set(frame["id"].astype(str))
    manifest_rows = {
        row["id"]: row for row in _read_csv(manifest)
        if row.get("pdbqt") and Path(row["pdbqt"]).exists()
    }
    tasks = [manifest_rows[i] for i in top_ids if i in manifest_rows]
    if not tasks:
        raise DockingError("no prepared PDBQT files found for top hits")

    vina = find_tool(cfg.get("docking", "executable", "vina"))
    if not vina:
        raise ToolNotFoundError("AutoDock Vina not found; run check-env")

    redock_dir = cfg.output_dir / "redock"
    redock_dir.mkdir(parents=True, exist_ok=True)
    results_path = redock_dir / "results.csv"

    done: set[str] = set()
    if cfg.get("redock", "resume", True) and results_path.exists():
        done = {
            row["id"] for row in _read_csv(results_path)
            if row.get("status") == "ok"
        }

    max_workers = int(cfg.get("redock", "max_workers", 4))
    pending = [row for row in tasks if row["id"] not in done]
    log.info(
        "redocking %s top hits with exhaustiveness %s",
        len(pending),
        cfg.get("redock", "exhaustiveness", 32),
    )

    original = {
        key: cfg.data["docking"].get(key)
        for key in ["exhaustiveness", "num_modes", "energy_range"]
    }
    cfg.data["docking"]["exhaustiveness"] = cfg.get("redock", "exhaustiveness", 32)
    cfg.data["docking"]["num_modes"] = cfg.get("redock", "num_modes", 9)
    cfg.data["docking"]["energy_range"] = cfg.get("redock", "energy_range", 3.0)
    try:
        if pending:
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = {
                    pool.submit(_dock_one, cfg, vina, row, redock_dir, log): row
                    for row in pending
                }
                for future in as_completed(futures):
                    row = futures[future]
                    try:
                        record = future.result()
                    except Exception as exc:
                        from .docking import _record

                        record = _record(row, status="error", error=str(exc))
                    _append_row(results_path, record)
    finally:
        for key, value in original.items():
            if value is None:
                cfg.data["docking"].pop(key, None)
            else:
                cfg.data["docking"][key] = value

    all_results = _read_csv(results_path)
    ok_results = [r for r in all_results if r.get("status") == "ok"]
    affinities = [float(r["affinity"]) for r in ok_results if r.get("affinity")]
    summary = {
        "top_n": len(tasks),
        "ok": len(ok_results),
        "best_affinity": min(affinities) if affinities else None,
        "results_csv": str(results_path),
        "redock_dir": str(redock_dir),
    }
    write_json(redock_dir / "summary.json", summary)
    log.info("redock complete: %s ok, best affinity %s", len(ok_results), summary["best_affinity"])
    return summary
