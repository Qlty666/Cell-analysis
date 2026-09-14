"""Five-axis target validation scoring with transparent GO/NO-GO tiers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from docking.utils import write_json


def _numeric(frame: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def _tier(score: float, has_evidence: bool) -> str:
    if score >= 75:
        return "GO"
    if score >= 50:
        return "CONDITIONAL_GO"
    if score >= 25 or not has_evidence:
        return "REVIEW"
    return "NO_GO"


def build_target_validation(
    workdir: Path,
    out_dir: Path,
) -> tuple[pd.DataFrame, dict]:
    """Build five evidence-axis validation scores for prioritized targets."""
    priority_path = out_dir / "integrated_target_priority.csv"
    if not priority_path.exists():
        priority_path = out_dir / "target_priority.csv"
    if not priority_path.exists():
        frame = pd.DataFrame()
    else:
        frame = pd.read_csv(priority_path)
    if frame.empty:
        summary = {
            "targets": 0,
            "decisions": {},
            "top_targets": [],
        }
        write_json(out_dir / "target_validation_summary.json", summary)
        return frame, summary

    gene_column = next(
        (
            column
            for column in ("gene", "target_symbol")
            if column in frame.columns
        ),
        None,
    )
    if gene_column is None:
        raise ValueError("target priority table has no gene column")
    frame = frame.rename(columns={gene_column: "gene"})
    frame["gene"] = frame["gene"].astype(str).str.strip().str.upper()
    frame = frame.drop_duplicates("gene", keep="first")

    legacy = pd.DataFrame()
    legacy_path = out_dir / "gene_evidence.csv"
    if legacy_path.exists():
        legacy = pd.read_csv(legacy_path)
        if "gene" in legacy.columns:
            legacy["gene"] = legacy["gene"].astype(str).str.strip().str.upper()
            legacy = legacy.drop_duplicates("gene", keep="first")
    if not legacy.empty:
        frame = frame.merge(legacy, on="gene", how="left", suffixes=("", "_legacy"))

    genetic = _numeric(frame, "genetic_association")
    disease = _numeric(frame, "disease_association")
    direct = _numeric(frame, "direct_experimental")
    curated = _numeric(frame, "curated_target")
    clinical = _numeric(frame, "clinical_precedent")
    structure = _numeric(frame, "structure")
    safety = _numeric(frame, "safety_risk")
    bioactivity = np.clip(
        np.log1p(_numeric(frame, "chembl_bioactivities"))
        / np.log1p(50.0),
        0.0,
        1.0,
    )
    pdb = _numeric(frame, "pdb_structures")
    alphafold = _numeric(frame, "alphafold_structures")

    disease_association = np.maximum(genetic, disease)
    druggability = np.maximum.reduce([direct, curated, bioactivity])
    chemical_matter = np.maximum(direct, bioactivity)
    structural_data = np.maximum.reduce(
        [
            structure,
            np.where(pdb > 0, np.clip(0.65 + 0.1 * pdb, 0.0, 1.0), 0.0),
            np.where(alphafold > 0, 0.5, 0.0),
        ]
    )
    axes = pd.DataFrame(
        {
            "gene": frame["gene"],
            "disease_association": disease_association * 20.0,
            "druggability": druggability * 20.0,
            "chemical_matter": chemical_matter * 20.0,
            "clinical_precedent": clinical * 20.0,
            "structural_data": structural_data * 20.0,
            "safety_risk": safety,
            "safety_penalty": safety * 15.0,
        }
    )
    raw_total = axes[
        [
            "disease_association",
            "druggability",
            "chemical_matter",
            "clinical_precedent",
            "structural_data",
        ]
    ].sum(axis=1)
    axes["raw_score"] = raw_total
    axes["adjusted_score"] = (raw_total - axes["safety_penalty"]).clip(
        0.0,
        100.0,
    )
    axes["decision"] = [
        _tier(float(score), float(raw) > 0)
        for score, raw in zip(
            axes["adjusted_score"],
            axes["raw_score"],
        )
    ]
    axes = axes.sort_values(
        ["adjusted_score", "raw_score", "gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    axes.insert(0, "validation_rank", np.arange(1, len(axes) + 1))
    axes.to_csv(out_dir / "target_validation_scores.csv", index=False)
    summary = {
        "targets": int(len(axes)),
        "decisions": {
            str(key): int(value)
            for key, value in axes["decision"].value_counts().items()
        },
        "median_adjusted_score": float(axes["adjusted_score"].median()),
        "top_targets": axes.head(20)[
            [
                "validation_rank",
                "gene",
                "decision",
                "adjusted_score",
                "disease_association",
                "druggability",
                "chemical_matter",
                "clinical_precedent",
                "structural_data",
                "safety_penalty",
            ]
        ].to_dict(orient="records"),
    }
    write_json(out_dir / "target_validation_summary.json", summary)
    return axes, summary
