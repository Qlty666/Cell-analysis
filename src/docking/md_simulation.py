#!/usr/bin/env python3
"""GROMACS molecular dynamics simulation for top docking poses.

The module can run in two modes:

- ``prepare``: turn docking poses into a GROMACS-ready complex without
  requiring ligand force-field tools. This is useful for review and for
  running the generated inputs on a separate workstation.
- ``auto``: run minimization, NVT/NPT equilibration and a short production
  simulation with GROMACS. Ligand parameters are generated with ACPYPE
  (GAFF2/AM1-BCC) when available, or taken from ``md_simulation.topology_dir``.
"""

from __future__ import annotations

import html
import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .config import ResolvedConfig
from .utils import (
    DockingError,
    ToolNotFoundError,
    find_tool,
    run_command,
    tail,
)

MD_DIR_NAME = "06_md"


def run_md_simulation(
    cfg: ResolvedConfig,
    log,
    mode: str | None = None,
    force: bool = False,
) -> dict:
    """Prepare or run GROMACS MD for the top ranked docking poses."""

    mode = (mode or cfg.get("md_simulation", "mode", "prepare") or "prepare").lower()
    if mode not in ("prepare", "auto"):
        raise DockingError("md_simulation.mode must be 'prepare' or 'auto'")

    out_dir = cfg.md_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = _pick_hits(cfg)
    top_n = int(
        cfg.get("md_simulation", "top_n")
        or cfg.get("md", "top_n", 1)
        or 1
    )
    selected = hits.head(top_n)
    log.info(
        "MD %s for %s top poses -> %s",
        mode,
        len(selected),
        out_dir,
    )

    rows: list[dict] = []
    ok_count = 0
    for _, row in selected.iterrows():
        lig_id = str(row.get("id", "ligand"))
        run_dir = out_dir / lig_id
        entry = {
            "id": lig_id,
            "affinity": row.get("affinity", ""),
            "smiles": row.get("smiles", ""),
            "mode": mode,
            "status": "prepared",
            "error": "",
            "time_ns": "",
            "rmsd_protein_mean_nm": "",
            "rmsd_ligand_mean_nm": "",
            "rmsf_ligand_mean_nm": "",
        }
        try:
            entry.update(_prepare_hit_dir(cfg, row, run_dir, log))
            if mode == "auto":
                metrics = _run_gromacs_hit(cfg, row, run_dir, log, force=force)
                entry.update(metrics)
                entry["status"] = "completed"
                ok_count += 1
            else:
                log.info("MD inputs prepared for %s -> %s", lig_id, run_dir)
        except Exception as exc:
            entry["status"] = "failed"
            entry["error"] = str(exc)
            log.error("MD simulation failed for %s: %s", lig_id, exc)
        rows.append(entry)

    results_path = out_dir / "md_simulation_results.csv"
    pd.DataFrame(rows).to_csv(results_path, index=False)
    completed = sum(1 for row in rows if row["status"] == "completed")
    prepared = sum(1 for row in rows if row["status"] == "prepared")
    failed = sum(1 for row in rows if row["status"] == "failed")
    summary = {
        "mode": mode,
        "requested": int(len(selected)),
        "completed": completed,
        "prepared": prepared,
        "failed": failed,
        "output_dir": str(out_dir),
        "results_csv": str(results_path),
    }
    write_json(out_dir / "md_simulation_summary.json", summary)
    _write_report(out_dir, rows, summary)

    if failed == len(rows) and len(rows) > 0:
        errors = [row["error"] for row in rows if row.get("error")]
        raise DockingError(
            ("no MD simulation completed: " if mode == "auto" else "no MD inputs prepared: ")
            + (errors[-1] if errors else "unknown error")
        )
    log.info(
        "MD simulation %s complete: %s completed, %s prepared, %s failed",
        mode,
        completed,
        prepared,
        failed,
    )
    return summary


def _pick_hits(cfg: ResolvedConfig) -> pd.DataFrame:
    reports = cfg.reports_dir()
    candidates = [
        reports / "01_analysis" / "data" / "fig_48_diverse_hits.csv",
        reports / "01_analysis" / "data" / "fig_47_top_hits.csv",
        reports / "01_analysis" / "data" / "fig_46_47_ranked_results.csv",
        cfg.results_path(),
    ]
    for path in candidates:
        if path.exists():
            frame = pd.read_csv(path, dtype={"id": str})
            if not frame.empty:
                return frame
    raise DockingError(
        "no docking analysis found; run dock/analyze before MD simulation"
    )


def _prepare_hit_dir(
    cfg: ResolvedConfig,
    row,
    run_dir: Path,
    log,
) -> dict:
    run_dir.mkdir(parents=True, exist_ok=True)
    receptor_pdb = _receptor_pdb(cfg)
    if not receptor_pdb.exists():
        raise DockingError(f"receptor PDB not found for MD: {receptor_pdb}")

    pose_file = _pose_path(cfg, row)
    if not pose_file.exists():
        raise DockingError(f"docking pose not found: {pose_file}")

    protein_pdb = _write_protein_pdb(receptor_pdb, run_dir)
    ligand_pdb = _write_ligand_pdb(pose_file, run_dir)
    _write_complex_pdb(protein_pdb, ligand_pdb, run_dir / "complex.pdb")
    _write_mdp_files(run_dir, cfg)
    _write_prepare_readme(run_dir, cfg)

    return {
        "protein_pdb": str(protein_pdb),
        "ligand_pdb": str(ligand_pdb),
        "complex_pdb": str(run_dir / "complex.pdb"),
    }


def _receptor_pdb(cfg: ResolvedConfig) -> Path:
    configured = cfg.get("md_simulation", "receptor_pdb")
    if configured:
        return cfg._resolve(str(configured), cfg.workdir)
    candidate = cfg.receptor_input()
    if candidate.suffix.lower() not in (".pdb", ".ent"):
        default = cfg.workdir / "data" / "receptors" / "receptor.pdb"
        if default.exists():
            return default
    return candidate


def _pose_path(cfg: ResolvedConfig, row) -> Path:
    pose = row.get("pose_file") if hasattr(row, "get") else None
    if pose and Path(str(pose)).exists():
        return Path(str(pose))
    return cfg.docked_dir() / f"{str(row.get('id', 'ligand'))}.pdbqt"


def _write_protein_pdb(source: Path, run_dir: Path) -> Path:
    """Keep only standard ATOM records so pdb2gmx can build the topology."""

    out = run_dir / "protein.pdb"
    kept: list[str] = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM  "):
            resname = line[17:20].strip().upper()
            if resname in {"HOH", "WAT", "NA", "CL", "K", "CA", "MG", "ZN"}:
                continue
            kept.append(line[:6] + line[6:])
    if not kept:
        raise DockingError(f"no protein ATOM records found in {source}")
    text = "\n".join(kept) + "\nTER\nEND\n"
    out.write_text(text, encoding="utf-8")
    return out


def _write_ligand_pdb(pose_file: Path, run_dir: Path) -> Path:
    out = run_dir / "ligand.pdb"
    if pose_file.suffix.lower() == ".pdb":
        shutil.copyfile(pose_file, out)
        return out
    clean_pdbqt = run_dir / "pose_first_model.pdbqt"
    _write_first_pdbqt_model(pose_file, clean_pdbqt)
    obabel = find_tool("obabel")
    if obabel:
        result = run_command(
            [obabel, str(clean_pdbqt), "-O", str(out), "-p", "7.4"],
            timeout=120,
            cwd=run_dir,
        )
        if result.returncode == 0 and out.exists() and out.stat().st_size > 0:
            return out
    _pdbqt_to_pdb_fallback(clean_pdbqt, out)
    return out


def _write_first_pdbqt_model(pose_file: Path, out: Path) -> None:
    lines = pose_file.read_text(encoding="utf-8", errors="replace").splitlines()
    if not any(line.strip().startswith("MODEL") for line in lines):
        out.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return
    inside = False
    kept: list[str] = []
    for line in lines:
        if line.strip().startswith("MODEL"):
            inside = True
            continue
        if line.strip().startswith("ENDMDL") and inside:
            break
        if inside and line.strip() and not line.strip().startswith("END"):
            kept.append(line)
    if not kept:
        raise DockingError(f"cannot find first MODEL in {pose_file}")
    out.write_text("\n".join(kept) + "\n", encoding="utf-8")


def _pdbqt_to_pdb_fallback(pose_file: Path, out: Path) -> None:
    kept: list[str] = []
    for line in pose_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if line.startswith("MODEL") or line.startswith("ENDMDL"):
            continue
        record = "HETATM" + line[6:]
        record = record[:17] + "LIG" + record[20:]
        kept.append(record)
    if not kept:
        raise DockingError(f"cannot extract ligand coordinates from {pose_file}")
    out.write_text("\n".join(kept) + "\nEND\n", encoding="utf-8")


def _write_complex_pdb(protein: Path, ligand: Path, out: Path) -> None:
    protein_lines = [
        line
        for line in protein.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM", "TER"))
    ]
    ligand_lines = [
        line
        for line in ligand.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM  ", "HETATM", "TER"))
    ]
    text = "\n".join(protein_lines + ligand_lines).replace("\nTER\nTER", "\nTER")
    text = text.rstrip() + "\nEND\n"
    out.write_text(text, encoding="utf-8")


def _run_gromacs_hit(
    cfg: ResolvedConfig,
    row,
    run_dir: Path,
    log,
    force: bool = False,
) -> dict:
    if not force and (run_dir / "md.tpr").exists():
        log.info("reusing existing GROMACS simulation for %s", row.get("id"))
        return _analyze_gromacs_output(cfg, run_dir)

    gmx = _find_gmx(cfg)
    if not gmx:
        raise ToolNotFoundError(
            "GROMACS gmx not found; add it to PATH or set md_simulation.executable"
        )
    if not (run_dir / "complex.pdb").exists():
        _prepare_hit_dir(cfg, row, run_dir, log)

    log.info("GROMACS: building protein topology for %s", row.get("id"))
    _run_pdb2gmx(gmx, cfg, run_dir)
    log.info("GROMACS: preparing ligand topology for %s", row.get("id"))
    lig_gro = _prepare_ligand_topology(cfg, row, run_dir)
    _merge_gro(
        run_dir / "protein.gro",
        lig_gro,
        run_dir / "complex.gro",
    )
    _patch_topology(run_dir / "topol.top", run_dir)
    log.info("GROMACS: solvating and adding ions for %s", row.get("id"))
    _solvate_and_ions(gmx, cfg, run_dir)
    _write_index(run_dir)
    log.info("GROMACS: running em/nvt/npt/production for %s", row.get("id"))
    _run_stages(gmx, cfg, run_dir)
    log.info("GROMACS: analyzing RMSD/RMSF for %s", row.get("id"))
    return _analyze_gromacs_output(cfg, run_dir)


def _run_pdb2gmx(gmx: str, cfg: ResolvedConfig, run_dir: Path) -> None:
    protein = run_dir / "protein.pdb"
    if not protein.exists():
        raise DockingError(f"protein.pdb missing in {run_dir}")
    cmd = [
        gmx,
        "pdb2gmx",
        "-f",
        "protein.pdb",
        "-o",
        "protein.gro",
        "-p",
        "topol.top",
        "-i",
        "posre.itp",
        "-ff",
        str(cfg.get("md_simulation", "protein_forcefield", "amber99sb-ildn")),
        "-water",
        str(cfg.get("md_simulation", "water", "tip3p")),
        "-ignh",
    ]
    _check_gmx(
        run_command(
            cmd,
            timeout=_timeout(cfg),
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
        ),
        "pdb2gmx",
    )


def _prepare_ligand_topology(
    cfg: ResolvedConfig,
    row,
    run_dir: Path,
) -> Path:
    lig_id = str(row.get("id", "ligand"))
    supplied = _copy_supplied_topology(cfg, lig_id, run_dir)
    if supplied is not None:
        _normalize_ligand_itp(run_dir / "ligand.itp")
        return supplied

    acpype = find_tool("acpype")
    if not acpype:
        raise ToolNotFoundError(
            "ACPYPE/acpype not found and md_simulation.topology_dir is empty; "
            "install acpype + AmberTools or place <id>.itp/<id>.gro in topology_dir"
        )
    obabel = find_tool("obabel")
    if not obabel:
        raise ToolNotFoundError("Open Babel is required to create ligand.mol2")
    lig_pdb = run_dir / "ligand.pdb"
    mol2 = run_dir / "ligand.mol2"
    result = run_command(
        [obabel, str(lig_pdb), "-O", str(mol2)],
        timeout=180,
        cwd=run_dir,
    )
    _check_gmx(result, "obabel ligand -> mol2")

    charge = int(cfg.get("md_simulation", "ligand_charge", 0) or 0)
    ff = str(cfg.get("md_simulation", "ligand_forcefield", "gaff2") or "gaff2")
    cmd = [
        acpype,
        "-i",
        "ligand.mol2",
        "-f",
        "mol2",
        "-a",
        ff,
        "-c",
        "bcc",
        "-n",
        str(charge),
        "-o",
        "gmx",
        "-d",
    ]
    _check_gmx(
        run_command(cmd, timeout=int(cfg.get("md_simulation", "timeout_seconds", 3600)), cwd=run_dir),
        "acpype",
    )

    itp = _find_acpype_file(run_dir, "*_GMX.itp")
    gro = _find_acpype_file(run_dir, "*_GMX.gro")
    if itp is None or gro is None:
        raise DockingError("acpype did not produce _GMX.itp/_GMX.gro files")
    shutil.copyfile(itp, run_dir / "ligand.itp")
    shutil.copyfile(gro, run_dir / "ligand.gro")
    _normalize_ligand_itp(run_dir / "ligand.itp")
    return run_dir / "ligand.gro"


def _copy_supplied_topology(
    cfg: ResolvedConfig,
    lig_id: str,
    run_dir: Path,
) -> Path | None:
    top_dir_value = cfg.get("md_simulation", "topology_dir")
    if not top_dir_value:
        return None
    top_dir = cfg._resolve(str(top_dir_value), cfg.workdir)
    candidates = [
        top_dir / f"{lig_id}.itp",
        top_dir / lig_id / f"{lig_id}.itp",
    ]
    itp = next((path for path in candidates if path.exists()), None)
    gro_candidates = [
        top_dir / f"{lig_id}.gro",
        top_dir / lig_id / f"{lig_id}.gro",
    ]
    gro = next((path for path in gro_candidates if path.exists()), None)
    if itp is None or gro is None:
        raise DockingError(
            f"topology_dir does not contain {lig_id}.itp and {lig_id}.gro"
        )
    shutil.copyfile(itp, run_dir / "ligand.itp")
    shutil.copyfile(gro, run_dir / "ligand.gro")
    return run_dir / "ligand.gro"


def _find_acpype_file(run_dir: Path, pattern: str) -> Path | None:
    matches = sorted(
        path
        for path in run_dir.rglob(pattern)
        if "protein" not in path.stem.lower() and "topol" not in path.stem.lower()
    )
    return matches[0] if matches else None


def _normalize_ligand_itp(itp_path: Path) -> None:
    text = itp_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    marker = None
    for idx, line in enumerate(lines):
        if line.strip() == "[ moleculetype ]":
            marker = idx + 1
            break
    if marker is None or marker >= len(lines):
        return
    data_idx = marker
    while data_idx < len(lines) and (
        not lines[data_idx].strip()
        or lines[data_idx].lstrip().startswith(";")
    ):
        data_idx += 1
    if data_idx >= len(lines):
        return
    parts = lines[data_idx].split()
    if parts and parts[0] != "LIG":
        parts[0] = "LIG"
        lines[data_idx] = " ".join(parts)
        itp_path.write_text("\n".join(lines), encoding="utf-8")


def _merge_gro(protein: Path, ligand: Path, out: Path) -> None:
    prot_lines = protein.read_text(encoding="utf-8", errors="replace").splitlines()
    lig_lines = ligand.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(prot_lines) < 3 or len(lig_lines) < 3:
        raise DockingError("invalid GRO file while merging protein and ligand")
    n_prot = int(prot_lines[1].strip())
    n_lig = int(lig_lines[1].strip())
    prot_atoms = prot_lines[2 : 2 + n_prot]
    lig_atoms = lig_lines[2 : 2 + n_lig]
    if not prot_atoms or not lig_atoms:
        raise DockingError("empty protein or ligand atom section in GRO")
    total = n_prot + n_lig
    merged = [
        line[:15] + f"{i:5d}" + line[20:]
        for i, line in enumerate(prot_atoms + lig_atoms, 1)
    ]
    box = prot_lines[2 + n_prot] if 2 + n_prot < len(prot_lines) else ""
    text = "\n".join(
        ["Protein-Ligand complex", f"{total:5d}", *merged, box]
    ).rstrip() + "\n"
    out.write_text(text, encoding="utf-8")


def _patch_topology(topol: Path, run_dir: Path) -> None:
    text = topol.read_text(encoding="utf-8", errors="replace")
    if "ligand.itp" in text and "\nLIG" in text.split("[ molecules ]")[-1]:
        return
    lines = text.splitlines()
    system_idx = next(
        (
            idx
            for idx, line in enumerate(lines)
            if line.strip().startswith("[ system ]")
        ),
        len(lines),
    )
    include = '\n; ligand parameters\n#include "ligand.itp"\n'
    lines.insert(system_idx, include)
    molecules_idx = next(
        (
            idx
            for idx, line in enumerate(lines)
            if line.strip().startswith("[ molecules ]")
        ),
        None,
    )
    if molecules_idx is not None and "LIG" not in " ".join(lines[molecules_idx + 1 :]):
        lines.append("LIG             1")
    text = "\n".join(lines).rstrip() + "\n"
    topol.write_text(text, encoding="utf-8")


def _solvate_and_ions(gmx: str, cfg: ResolvedConfig, run_dir: Path) -> None:
    timeout = _timeout(cfg)
    padding = float(cfg.get("md_simulation", "box_padding_nm", 1.2) or 1.2)
    box_type = str(cfg.get("md_simulation", "box_type", "cubic") or "cubic")
    _check_gmx(
        run_command(
            [
                gmx,
                "editconf",
                "-f",
                "complex.gro",
                "-o",
                "box.gro",
                "-c",
                "-d",
                str(padding),
                "-bt",
                box_type,
            ],
            timeout=timeout,
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
        ),
        "editconf",
    )
    _check_gmx(
        run_command(
            [
                gmx,
                "solvate",
                "-cp",
                "box.gro",
                "-cs",
                "spc216.gro",
                "-p",
                "topol.top",
                "-o",
                "solvated.gro",
            ],
            timeout=timeout,
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
        ),
        "solvate",
    )
    _check_gmx(
        run_command(
            [
                gmx,
                "grompp",
                "-f",
                "em.mdp",
                "-c",
                "solvated.gro",
                "-p",
                "topol.top",
                "-o",
                "ions.tpr",
                "-maxwarn",
                str(int(cfg.get("md_simulation", "maxwarn", 20) or 20)),
            ],
            timeout=timeout,
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
        ),
        "grompp ions",
    )
    result = run_command(
        [
            gmx,
            "genion",
            "-s",
            "ions.tpr",
            "-o",
            "solvated_ions.gro",
            "-p",
            "topol.top",
            "-pname",
            "NA",
            "-nname",
            "CL",
            "-neutral",
            "-conc",
            str(float(cfg.get("md_simulation", "ion_concentration", 0.15) or 0.15)),
        ],
        timeout=timeout,
        cwd=run_dir,
        env=_gmx_env(gmx, cfg),
        stdin_text="SOL\n",
    )
    _check_gmx(result, "genion")


def _write_index(run_dir: Path) -> None:
    gro_path = run_dir / "solvated_ions.gro"
    lines = gro_path.read_text(encoding="utf-8", errors="replace").splitlines()
    total = int(lines[1].strip())
    complex_gro = run_dir / "complex.gro"
    complex_lines = complex_gro.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    n_complex = int(complex_lines[1].strip())
    protein_gro = run_dir / "protein.gro"
    n_protein = int(
        protein_gro.read_text(encoding="utf-8", errors="replace").splitlines()[1].strip()
    )
    if n_complex > total:
        raise DockingError("solvated GRO has fewer atoms than complex GRO")
    groups = {
        "Protein": list(range(1, n_protein + 1)),
        "LIG": list(range(n_protein + 1, n_complex + 1)),
        "Protein_LIG": list(range(1, n_complex + 1)),
        "Water_and_ions": list(range(n_complex + 1, total + 1)),
        "System": list(range(1, total + 1)),
    }
    blocks: list[str] = []
    for name, indices in groups.items():
        blocks.append(f"[ {name} ]")
        for start in range(0, len(indices), 16):
            blocks.append(" ".join(str(i) for i in indices[start : start + 16]))
        blocks.append("")
    (run_dir / "index.ndx").write_text("\n".join(blocks), encoding="utf-8")


def _run_stages(gmx: str, cfg: ResolvedConfig, run_dir: Path) -> None:
    timeout = _timeout(cfg)
    stages = [
        ("em", "em.mdp", "solvated_ions.gro", None),
        ("nvt", "nvt.mdp", "em.gro", "em.gro"),
        ("npt", "npt.mdp", "nvt.gro", "nvt.gro"),
        ("md", "md.mdp", "npt.gro", "npt.gro"),
    ]
    for name, mdp, coord, ref in stages:
        grompp = [
            gmx,
            "grompp",
            "-f",
            mdp,
            "-c",
            coord,
            "-p",
            "topol.top",
            "-o",
            f"{name}.tpr",
            "-maxwarn",
            str(int(cfg.get("md_simulation", "maxwarn", 20) or 20)),
        ]
        if ref:
            grompp += ["-r", ref]
        if name in ("npt", "md"):
            grompp += ["-t", f"{'nvt' if name == 'npt' else 'npt'}.cpt"]
        if (run_dir / "index.ndx").exists():
            grompp += ["-n", "index.ndx"]
        _check_gmx(
            run_command(
                grompp,
                timeout=timeout,
                cwd=run_dir,
                env=_gmx_env(gmx, cfg),
            ),
            f"grompp {name}",
        )
        mdrun = [
            gmx,
            "mdrun",
            "-deffnm",
            name,
            "-ntmpi",
            "1",
            "-ntomp",
            str(int(cfg.get("md_simulation", "cpu", 4) or 4)),
        ]
        if cfg.get("md_simulation", "gpu", False):
            mdrun += ["-nb", "gpu", "-pme", "gpu"]
        _check_gmx(
            run_command(
                mdrun,
                timeout=timeout,
                cwd=run_dir,
                env=_gmx_env(gmx, cfg),
            ),
            f"mdrun {name}",
        )
    _make_whole_trajectory(gmx, cfg, run_dir)


def _make_whole_trajectory(gmx: str, cfg: ResolvedConfig, run_dir: Path) -> None:
    src = run_dir / "md.xtc"
    if not src.exists():
        return
    if (run_dir / "index.ndx").exists():
        result = run_command(
            [
                gmx,
                "trjconv",
                "-s",
                "md.tpr",
                "-f",
                "md.xtc",
                "-o",
                "md_nojump.xtc",
                "-pbc",
                "mol",
                "-center",
                "-n",
                "index.ndx",
            ],
            timeout=_timeout(cfg),
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
            stdin_text="Protein\nSystem\n",
        )
        if result.returncode == 0:
            return
    # Fallback: copy the raw trajectory when centering is not available.
    shutil.copyfile(src, run_dir / "md_nojump.xtc")


def _analyze_gromacs_output(cfg: ResolvedConfig, run_dir: Path) -> dict:
    metrics = {"time_ns": "", "rmsd_protein_mean_nm": "", "rmsd_ligand_mean_nm": "", "rmsf_ligand_mean_nm": ""}
    tpr = run_dir / "md.tpr"
    xtc = run_dir / "md_nojump.xtc"
    if not tpr.exists() or not xtc.exists():
        return metrics
    gmx = _find_gmx(cfg)
    if not gmx:
        return metrics
    timeout = _timeout(cfg)
    ndx = run_dir / "index.ndx"
    protein_last_time: float | None = None
    if ndx.exists():
        for label in ("protein", "ligand"):
            out_xvg = run_dir / f"rmsd_{label}.xvg"
            result = run_command(
                [
                    gmx,
                    "rms",
                    "-s",
                    str(tpr),
                    "-f",
                    str(xtc),
                    "-n",
                    str(ndx),
                    "-o",
                    str(out_xvg),
                    "-tu",
                    "ns",
                ],
                timeout=timeout,
                cwd=run_dir,
                env=_gmx_env(gmx, cfg),
                stdin_text=("Protein\n" if label == "protein" else "LIG\n")
                + ("Protein\n" if label == "protein" else "LIG\n"),
            )
            if result.returncode == 0:
                data = parse_xvg(out_xvg)
                if data is not None and len(data):
                    metrics[f"rmsd_{label}_mean_nm"] = float(np.mean(data[:, 1]))
                    if label == "protein":
                        protein_last_time = float(data[-1, 0])
                    if label == "ligand" and metrics.get("time_ns") == "":
                        metrics["time_ns"] = float(data[-1, 0]) if len(data) else ""
        out_rmsf = run_dir / "rmsf_ligand.xvg"
        result = run_command(
            [
                gmx,
                "rmsf",
                "-s",
                str(tpr),
                "-f",
                str(xtc),
                "-n",
                str(ndx),
                "-o",
                str(out_rmsf),
            ],
            timeout=timeout,
            cwd=run_dir,
            env=_gmx_env(gmx, cfg),
            stdin_text="LIG\n",
        )
        if result.returncode == 0:
            data = parse_xvg(out_rmsf)
            if data is not None and len(data):
                metrics["rmsf_ligand_mean_nm"] = float(np.mean(data[:, 1]))
    if metrics.get("time_ns") == "" and protein_last_time is not None:
        metrics["time_ns"] = protein_last_time
    if cfg.get("md_simulation", "figures", True):
        try:
            make_figures(run_dir)
        except Exception:
            pass
    return metrics


def parse_xvg(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        try:
            rows.append([float(value) for value in parts])
        except ValueError:
            continue
    if not rows:
        return None
    return np.asarray(rows, dtype=float)


def make_figures(run_dir: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=False)
    for ax, label, title, color in [
        (axes[0], "rmsd_protein.xvg", "Protein backbone RMSD", "#1665c0"),
        (axes[1], "rmsd_ligand.xvg", "Ligand RMSD", "#c0392b"),
        (axes[2], "rmsf_ligand.xvg", "Ligand atom RMSF", "#2e7d32"),
    ]:
        data = parse_xvg(run_dir / label)
        if data is None or not len(data):
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            continue
        ax.plot(data[:, 0], data[:, 1] * 10.0, color=color, linewidth=1.2)
        ax.set_title(title)
        ax.set_ylabel("Angstrom")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("Time (ns) or atom index")
    fig.tight_layout()
    out = run_dir / "md_rmsd_rmsf.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)


def _write_mdp_files(run_dir: Path, cfg: ResolvedConfig) -> None:
    em_steps = int(cfg.get("md_simulation", "em_steps", 5000) or 5000)
    equil_steps = int(cfg.get("md_simulation", "equil_steps", 5000) or 5000)
    prod_steps = int(cfg.get("md_simulation", "prod_steps", 250000) or 250000)
    dt = float(cfg.get("md_simulation", "dt_ps", 0.002) or 0.002)
    temp = float(cfg.get("md_simulation", "temperature", 300) or 300)
    pressure = float(cfg.get("md_simulation", "pressure", 1.0) or 1.0)
    (run_dir / "em.mdp").write_text(
        _em_mdp(em_steps),
        encoding="utf-8",
    )
    (run_dir / "nvt.mdp").write_text(
        _thermo_mdp("nvt", equil_steps, dt, temp, pressure),
        encoding="utf-8",
    )
    (run_dir / "npt.mdp").write_text(
        _thermo_mdp("npt", equil_steps, dt, temp, pressure),
        encoding="utf-8",
    )
    (run_dir / "md.mdp").write_text(
        _thermo_mdp("md", prod_steps, dt, temp, pressure),
        encoding="utf-8",
    )


def _em_mdp(steps: int) -> str:
    return f"""; energy minimization
integrator               = steep
nsteps                   = {steps}
emtol                    = 1000.0
emstep                   = 0.01
nstlist                  = 10
cutoff-scheme            = Verlet
ns_type                  = grid
rlist                    = 1.2
coulombtype              = PME
rcoulomb                 = 1.2
vdwtype                  = cutoff
vdw-modifier             = force-switch
rvdw-switch              = 1.0
rvdw                     = 1.2
pbc                      = xyz
DispCorr                 = no
"""


def _thermo_mdp(kind: str, steps: int, dt: float, temp: float, pressure: float) -> str:
    define = ""
    pcoupl = "no"
    continuation = "no"
    gen_vel = "yes"
    if kind == "nvt":
        define = "define                   = -DPOSRES\n"
    elif kind == "npt":
        pcoupl = "Parrinello-Rahman"
        continuation = "yes"
        gen_vel = "no"
    else:
        pcoupl = "Parrinello-Rahman"
        continuation = "yes"
        gen_vel = "no"
    title = {
        "nvt": "NVT equilibration",
        "npt": "NPT equilibration",
        "md": "Production MD",
    }[kind]
    if kind == "md":
        output_controls = "nstxout = 0\nnstvout = 0\n"
    else:
        output_controls = "nstxout = 500\nnstvout = 500\n"
    return f"""; {title}
integrator               = md
dt                       = {dt}
nsteps                   = {steps}
{define}{output_controls}
nstenergy                = 500
nstlog                   = 500
nstxout-compressed       = 5000

cutoff-scheme            = Verlet
ns_type                  = grid
nstlist                  = 10
rlist                    = 1.2
rcoulomb                 = 1.2
rvdw                     = 1.2
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
DispCorr                 = no

tcoupl                   = V-rescale
tc-grps                  = Protein_LIG Water_and_ions
tau_t                    = 0.1 0.1
ref_t                    = {temp} {temp}

pcoupl                   = {pcoupl}
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = {pressure}
compressibility          = 4.5e-5

pbc                      = xyz
constraints              = h-bonds
constraint_algorithm     = LINCS
lincs_iter               = 1
lincs_order              = 4

continuation             = {continuation}
gen_vel                  = {gen_vel}
gen_temp                 = {temp}
gen_seed                 = -1
"""


def _write_prepare_readme(run_dir: Path, cfg: ResolvedConfig) -> None:
    text = f"""# GROMACS MD inputs

This folder contains one top docking pose converted into a protein-ligand
complex. Run mode used: {cfg.get("md_simulation", "mode", "prepare")}.

To complete the simulation:
1. Install GROMACS and ACPYPE (acpype + AmberTools), or place topology files
   in md_simulation.topology_dir as <id>.itp and <id>.gro.
2. Run `python scripts\\run_docking.py md-simulation --md-mode auto`.
3. Check md_simulation_results.csv, md_simulation_summary.json and figures.
"""
    (run_dir / "README.md").write_text(text, encoding="utf-8")


def _write_report(out_dir: Path, rows: list[dict], summary: dict) -> None:
    table = "".join(
        "<tr>"
        f"<td>{html.escape(str(row.get('id', '')))}</td>"
        f"<td>{html.escape(str(row.get('status', '')))}</td>"
        f"<td>{html.escape(str(row.get('time_ns', '')))}</td>"
        f"<td>{html.escape(str(row.get('rmsd_ligand_mean_nm', '')))}</td>"
        f"<td>{html.escape(str(row.get('error', '')))}</td>"
        "</tr>"
        for row in rows
    )
    figures = sorted(
        path.relative_to(out_dir).as_posix()
        for path in out_dir.rglob("*.png")
    )
    images = "".join(
        f'<figure><img src="{html.escape(name)}"><figcaption>{html.escape(name)}</figcaption></figure>'
        for name in figures
    )
    (out_dir / "md_simulation_report.html").write_text(
        f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>MD Simulation Report</title>
<style>
body {{ font-family: "Segoe UI", Arial, sans-serif; margin: 24px; color: #1f2933; background: #f5f7fa; }}
.card {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 18px; margin-bottom: 16px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th, td {{ border: 1px solid #e5e7eb; padding: 7px 8px; text-align: left; }}
th {{ background: #eef2f7; }}
.gallery {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }}
img {{ width: 100%; border: 1px solid #e5e7eb; border-radius: 6px; }}
</style>
</head>
<body>
<h1>MD Simulation Report</h1>
<div class="card"><p>{html.escape(str(summary))}</p></div>
<div class="card"><h2>Results</h2><table><thead><tr><th>ID</th><th>Status</th><th>Time (ns)</th><th>Ligand RMSD mean (nm)</th><th>Error</th></tr></thead><tbody>{table}</tbody></table></div>
<div class="card"><h2>Figures</h2><div class="gallery">{images or '<p>No figures</p>'}</div></div>
</body>
</html>
""",
        encoding="utf-8",
    )


def _gmx_env(gmx: str, cfg: ResolvedConfig) -> dict:
    env = os.environ.copy()
    data_dir = cfg.get("md_simulation", "gmx_data_dir")
    if not data_dir:
        candidate = Path(gmx).resolve().parent.parent
        if (candidate / "share" / "gromacs" / "top").is_dir():
            data_dir = str(candidate)
    if data_dir:
        data_path = Path(str(data_dir)).expanduser().resolve()
        env["GMXDATA"] = str(data_path)
    return env


def _find_gmx(cfg: ResolvedConfig) -> str | None:
    configured = cfg.get("md_simulation", "executable")
    if configured:
        return str(configured)
    # Windows GROMACS installs sometimes expose a gmx.EXE wrapper with a
    # hard-coded data prefix; gmx.exe is the real binary.
    return find_tool("gmx.exe") or find_tool("gmx")


def _timeout(cfg: ResolvedConfig) -> int:
    return int(cfg.get("md_simulation", "timeout_seconds", 86400) or 86400)


def _check_gmx(result: subprocess.CompletedProcess, label: str) -> None:
    if result.returncode != 0:
        detail = _short_gmx_error(
            (result.stderr or "") + "\n" + (result.stdout or "")
        )
        if "No force fields found" in detail:
            detail += (
                " | set md_simulation.gmx_data_dir or GMXDATA to the "
                "GROMACS data root (contains share/gromacs/top)"
            )
        raise DockingError(
            f"{label} failed (code {result.returncode}): {detail}"
        )


def _short_gmx_error(text: str) -> str:
    lines = [line.strip() for line in (text or "").splitlines()]
    starts = [
        index
        for index, line in enumerate(lines)
        if line in ("Fatal error:", "Error in user input:")
        or line.startswith(("ERROR", "There were"))
    ]
    if starts:
        begin = starts[-1]
        end = len(lines)
        for index in range(begin + 1, len(lines)):
            if not lines[index]:
                end = index
                break
            if lines[index].startswith("For more information"):
                end = index
                break
        detail = "\n".join(lines[begin:end]).strip()
        if detail:
            return tail(detail, 1500)
    return tail(text)


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
