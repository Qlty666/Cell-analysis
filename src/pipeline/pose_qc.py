"""Optional PoseBusters validation for SDF docking poses."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

from docking.utils import write_json


def run_pose_qc(workdir: Path, out_dir: Path) -> dict:
    """Validate SDF poses when PoseBusters and receptor structures are available."""
    out_dir.mkdir(parents=True, exist_ok=True)
    sdf_files = sorted((workdir / "work").glob("*/**/*.sdf"))
    results_path = out_dir / "pose_qc_results.csv"
    if not sdf_files:
        summary = {
            "status": "skipped",
            "reason": "no SDF pose files found",
            "sdf_files": 0,
            "poses": 0,
            "pb_valid": 0,
            "pb_valid_rate": None,
            "gate_passed": False,
        }
        pd.DataFrame(
            columns=[
                "gene",
                "pose_file",
                "pose_index",
                "pb_valid",
                "source",
            ]
        ).to_csv(results_path, index=False)
        write_json(out_dir / "pose_qc_summary.json", summary)
        return summary
    if importlib.util.find_spec("posebusters") is None:
        summary = {
            "status": "unavailable",
            "reason": "PoseBusters is not installed",
            "sdf_files": len(sdf_files),
            "poses": 0,
            "pb_valid": 0,
            "pb_valid_rate": None,
            "gate_passed": False,
        }
        pd.DataFrame(
            columns=[
                "gene",
                "pose_file",
                "pose_index",
                "pb_valid",
                "source",
            ]
        ).to_csv(results_path, index=False)
        write_json(out_dir / "pose_qc_summary.json", summary)
        return summary

    from posebusters import PoseBusters

    bust = PoseBusters(config="dock")
    rows: list[dict] = []
    for sdf in sdf_files:
        gene = sdf.parents[1].name if len(sdf.parents) > 1 else ""
        receptors = sorted(sdf.parents[1].glob("**/data/receptors/*.pdb"))
        if not receptors:
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": -1,
                    "pb_valid": False,
                    "source": "no_receptor",
                }
            )
            continue
        try:
            frame = bust.bust(
                mol_pred=str(sdf),
                mol_cond=str(receptors[0]),
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": -1,
                    "pb_valid": False,
                    "source": f"error:{exc}",
                }
            )
            continue
        bool_columns = [
            column
            for column in frame.select_dtypes(include="bool").columns
            if not column.lower().startswith("rmsd")
        ]
        for index, row in frame.iterrows():
            valid = bool(row[bool_columns].all()) if bool_columns else False
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": int(index),
                    "pb_valid": valid,
                    "source": "posebusters",
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(results_path, index=False)
    valid_count = int(result["pb_valid"].sum()) if "pb_valid" in result else 0
    total = int(len(result))
    rate = valid_count / total if total else 0.0
    summary = {
        "status": "completed",
        "reason": "",
        "sdf_files": len(sdf_files),
        "poses": total,
        "pb_valid": valid_count,
        "pb_valid_rate": rate,
        "gate_passed": bool(total > 0 and rate >= 0.50),
    }
    write_json(out_dir / "pose_qc_summary.json", summary)
    return summary
