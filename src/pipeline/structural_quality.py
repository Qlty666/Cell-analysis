"""Structural-quality gates for docking and MD outputs."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from docking.utils import safe_name, write_json
from .pose_qc import run_pose_qc


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def _number(value: Any) -> float | None:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return None if pd.isna(parsed) else float(parsed)


def build_structural_quality(
    workdir: Path,
    out_dir: Path,
    *,
    protein_rmsd_std_cutoff_nm: float = 0.15,
    ligand_rmsd_std_cutoff_nm: float = 0.15,
) -> dict:
    """Summarize docking/MD validation gates for every successful target."""
    integration = workdir / "outputs" / "integration"
    docking_summary = _read_json(integration / "docking_summary.json")
    docking_targets = (
        pd.read_csv(integration / "docking_targets.csv")
        if (integration / "docking_targets.csv").exists()
        else pd.DataFrame()
    )
    rows: list[dict] = []
    for _, docking in docking_targets.iterrows():
        gene = str(docking.get("gene") or "").strip()
        if not gene:
            continue
        target = workdir / "work" / safe_name(gene, gene)
        analysis_summary = _read_json(
            target
            / "outputs"
            / "run_001"
            / "results"
            / "01_analysis"
            / "summary.json"
        )
        md_summary = _read_json(
            target
            / "outputs"
            / "run_001"
            / "results"
            / "06_md"
            / "md_simulation_summary.json"
        )
        md_results_path = (
            target
            / "outputs"
            / "run_001"
            / "results"
            / "06_md"
            / "md_simulation_results.csv"
        )
        md_results = (
            pd.read_csv(md_results_path)
            if md_results_path.exists()
            else pd.DataFrame()
        )
        protein_std = None
        ligand_std = None
        mmpbsa_available = False
        if not md_results.empty:
            if "rmsd_protein_tail_std_nm" in md_results.columns:
                protein_std = _number(
                    pd.to_numeric(
                        md_results["rmsd_protein_tail_std_nm"],
                        errors="coerce",
                    ).max()
                )
            if "rmsd_ligand_tail_std_nm" in md_results.columns:
                ligand_std = _number(
                    pd.to_numeric(
                        md_results["rmsd_ligand_tail_std_nm"],
                        errors="coerce",
                    ).max()
                )
            for column in md_results.columns:
                if not column.lower().startswith("mmpbsa"):
                    continue
                numeric = pd.to_numeric(md_results[column], errors="coerce")
                if numeric.notna().any() and (numeric.abs() > 1e-12).any():
                    mmpbsa_available = True
                    break
        protein_ok = bool(
            protein_std is not None
            and protein_std <= protein_rmsd_std_cutoff_nm
        )
        ligand_ok = bool(
            ligand_std is not None
            and ligand_std <= ligand_rmsd_std_cutoff_nm
        )
        md_rmsd_stable = bool(protein_ok and ligand_ok)
        rows.append(
            {
                "gene": gene,
                "docking_status": str(docking.get("status") or ""),
                "best_affinity": docking.get("best_affinity"),
                "docking_hits": docking.get("hits"),
                "replicate_consensus": bool(
                    analysis_summary.get("replicate_consensus", False)
                ),
                "md_mode": md_summary.get("mode", ""),
                "md_completed": int(md_summary.get("completed", 0) or 0),
                "md_prepared": int(md_summary.get("prepared", 0) or 0),
                "md_total_ns": md_summary.get("total_ns"),
                "protein_rmsd_tail_std_nm": protein_std,
                "ligand_rmsd_tail_std_nm": ligand_std,
                "mmpbsa_available": mmpbsa_available,
                "md_rmsd_stable": md_rmsd_stable,
            }
        )
    frame = pd.DataFrame(rows)
    positive_control = docking_summary.get("positive_control") or {}
    gates = {
        "positive_control": bool(positive_control.get("passed", False)),
        "replicate_consensus": bool(
            not frame.empty and frame["replicate_consensus"].any()
        ),
        "md_completed": bool(
            not frame.empty and (frame["md_completed"] > 0).any()
        ),
        "md_rmsd_stability": bool(
            not frame.empty and frame["md_rmsd_stable"].any()
        ),
        "mmpbsa_available": bool(
            not frame.empty and frame["mmpbsa_available"].any()
        ),
    }
    summary = {
        "targets": int(len(frame)),
        "docking_ok": int(
            (frame.get("docking_status", pd.Series(dtype=str)) == "ok").sum()
        ),
        "gates": gates,
        "all_structural_gates": all(gates.values()),
        "cutoffs": {
            "protein_rmsd_tail_std_nm": protein_rmsd_std_cutoff_nm,
            "ligand_rmsd_tail_std_nm": ligand_rmsd_std_cutoff_nm,
        },
    }
    pose_qc = run_pose_qc(workdir, out_dir)
    gates["posebusters_valid"] = bool(pose_qc.get("gate_passed", False))
    summary["pose_qc"] = pose_qc
    summary["all_publication_structural_gates"] = all(gates.values())
    out_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_dir / "structural_quality_targets.csv", index=False)
    write_json(out_dir / "structural_quality_summary.json", summary)
    return summary
