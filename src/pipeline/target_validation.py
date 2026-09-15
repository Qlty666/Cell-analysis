"""Five-axis target validation scoring with transparent GO/NO-GO tiers."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from docking.utils import write_json


AXES = (
    "disease_association",
    "druggability",
    "chemical_matter",
    "clinical_precedent",
    "structural_data",
)

DRUGGABLE_FAMILY_TERMS = (
    "kinase",
    "gpcr",
    "g protein-coupled",
    "protease",
    "ion channel",
    "nuclear receptor",
    "phosphatase",
    "epigenetic",
    "transporter",
)


def _numeric(
    frame: pd.DataFrame,
    column: str,
    default: float = 0.0,
) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def _numeric_optional(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _first_optional(frame: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    for column in columns:
        if column in frame.columns:
            values = _numeric_optional(frame, column)
            if values.notna().any():
                return values
    return pd.Series(np.nan, index=frame.index, dtype=float)


def _split_sources(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [
        item.strip()
        for item in re.split(r"[;,|]", text)
        if item.strip()
    ]


def _row_sources(
    frame: pd.DataFrame,
    row: pd.Series,
    columns: Iterable[str],
) -> list[str]:
    sources: list[str] = []
    for column in columns:
        if column not in frame.columns:
            continue
        sources.extend(_split_sources(row.get(column)))
    return list(dict.fromkeys(sources))


def _fallback_sources(
    sources: list[str],
    values: pd.Series,
    row_index: object,
) -> list[str]:
    if sources:
        return sources
    value = values.loc[row_index]
    if pd.notna(value) and float(value) > 0:
        return ["legacy_aggregate"]
    return []


def _tier(score: float, has_evidence: bool) -> str:
    if score >= 75:
        return "GO"
    if score >= 50:
        return "CONDITIONAL_GO"
    if score >= 25 or not has_evidence:
        return "REVIEW"
    return "NO_GO"


def _axis_status(
    value: float,
    sources: list[str],
    target_evidence_status: str,
) -> str:
    if (pd.notna(value) and float(value) > 0) or sources:
        return "available"
    if target_evidence_status == "not_queried":
        return "not_queried"
    if target_evidence_status == "covered":
        return "not_found"
    return "not_available"


def build_target_validation(
    workdir: Path,
    out_dir: Path,
) -> tuple[pd.DataFrame, dict]:
    """Build non-overlapping five-axis validation scores for priority targets.

    The axes intentionally do not reuse the same underlying evidence merely to
    increase the total score. Disease genetics and disease associations are
    collapsed with ``max``; druggability uses curated tractability/family
    information; chemical matter uses experimental or bioactivity evidence;
    clinical and structural evidence remain independent axes.
    """
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
            "axis_availability": {},
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
    frame = frame.drop_duplicates("gene", keep="first").reset_index(drop=True)

    legacy = pd.DataFrame()
    legacy_path = out_dir / "gene_evidence.csv"
    if legacy_path.exists():
        legacy = pd.read_csv(legacy_path)
        if "gene" in legacy.columns:
            legacy["gene"] = legacy["gene"].astype(str).str.strip().str.upper()
            legacy = legacy.drop_duplicates("gene", keep="first")
    if not legacy.empty:
        frame = frame.merge(legacy, on="gene", how="left", suffixes=("", "_legacy"))

    genetic = _numeric_optional(frame, "genetic_association")
    disease = _numeric_optional(frame, "disease_association")
    direct = _numeric_optional(frame, "direct_experimental")
    curated = _numeric_optional(frame, "curated_target")
    clinical = _numeric_optional(frame, "clinical_precedent")
    structure = _numeric_optional(frame, "structure")
    safety = _numeric(frame, "safety_risk")
    bioactivity = np.clip(
        np.log1p(_numeric(frame, "chembl_bioactivities"))
        / np.log1p(50.0),
        0.0,
        1.0,
    )
    pdb = _numeric(frame, "pdb_structures")
    alphafold = _numeric(frame, "alphafold_structures")

    family_column = next(
        (
            column
            for column in (
                "protein_family",
                "target_family",
                "uniprot_family",
                "family",
            )
            if column in frame.columns
        ),
        None,
    )
    if family_column:
        family_text = frame[family_column].fillna("").astype(str).str.lower()
        family_tractability = family_text.map(
            lambda text: 0.7
            if any(term in text for term in DRUGGABLE_FAMILY_TERMS)
            else np.nan
        )
    else:
        family_tractability = pd.Series(np.nan, index=frame.index, dtype=float)
    explicit_tractability = _first_optional(
        frame,
        ("target_tractability", "tractability_score", "druggability_score"),
    )

    disease_association = np.fmax(
        genetic.to_numpy(dtype=float),
        disease.to_numpy(dtype=float),
    )
    druggability_signal = np.fmax.reduce(
        [
            curated.to_numpy(dtype=float),
            explicit_tractability.to_numpy(dtype=float),
            family_tractability.to_numpy(dtype=float),
        ]
    )
    chemical_signal = np.fmax(
        direct.to_numpy(dtype=float),
        bioactivity.to_numpy(dtype=float),
    )
    pdb_signal = np.where(
        pdb.to_numpy(dtype=float) > 0,
        np.clip(0.65 + 0.1 * pdb.to_numpy(dtype=float), 0.0, 1.0),
        np.nan,
    )
    alphafold_signal = np.where(
        alphafold.to_numpy(dtype=float) > 0,
        0.5,
        np.nan,
    )
    structural_signal = np.fmax.reduce(
        [
            structure.to_numpy(dtype=float),
            pdb_signal,
            alphafold_signal,
        ]
    )
    clinical_signal = clinical.to_numpy(dtype=float)

    axis_values = {
        "disease_association": disease_association,
        "druggability": druggability_signal,
        "chemical_matter": chemical_signal,
        "clinical_precedent": clinical_signal,
        "structural_data": structural_signal,
    }
    source_columns = {
        "disease_association": (
            "genetic_association_sources",
            "disease_association_sources",
        ),
        "druggability": (
            "curated_target_sources",
            "target_tractability_sources",
            "druggability_sources",
        ),
        "chemical_matter": (
            "direct_experimental_sources",
            "chembl_bioactivities_sources",
        ),
        "clinical_precedent": ("clinical_precedent_sources",),
        "structural_data": (
            "structure_sources",
            "pdb_structures_sources",
            "alphafold_structures_sources",
        ),
    }
    source_values: dict[str, list[list[str]]] = {}
    axis_statuses: dict[str, list[str]] = {}
    evidence_status = (
        frame["evidence_status"].fillna("").astype(str)
        if "evidence_status" in frame.columns
        else pd.Series("unknown", index=frame.index, dtype=str)
    )
    for axis in AXES:
        values = pd.Series(axis_values[axis], index=frame.index)
        axis_sources: list[list[str]] = []
        statuses: list[str] = []
        for index, row in frame.iterrows():
            sources = _fallback_sources(
                _row_sources(frame, row, source_columns[axis]),
                values,
                index,
            )
            axis_sources.append(sources)
            statuses.append(
                _axis_status(
                    float(values.loc[index])
                    if pd.notna(values.loc[index])
                    else np.nan,
                    sources,
                    str(evidence_status.loc[index]),
                )
            )
        source_values[axis] = axis_sources
        axis_statuses[axis] = statuses

    axes = pd.DataFrame({"gene": frame["gene"]})
    for axis in AXES:
        values = pd.Series(axis_values[axis], index=frame.index).fillna(0.0)
        axes[axis] = values.to_numpy(dtype=float) * 20.0
        axes[f"{axis}_sources"] = [
            ";".join(sources) for sources in source_values[axis]
        ]
        axes[f"{axis}_evidence_status"] = axis_statuses[axis]
    axes["safety_risk"] = safety.to_numpy(dtype=float)
    axes["safety_penalty"] = axes["safety_risk"] * 15.0
    axes["evidence_completeness"] = (
        sum(
            pd.Series(axis_statuses[axis]) == "available"
            for axis in AXES
        )
        / len(AXES)
    )
    axes["axis_evidence_status"] = [
        json.dumps(
            {
                axis: axis_statuses[axis][index]
                for axis in AXES
            },
            ensure_ascii=True,
            sort_keys=True,
        )
        for index in range(len(frame))
    ]
    clinical_sources = [
        ";".join(sources) for sources in source_values["clinical_precedent"]
    ]
    axes["clinical_evidence_quality"] = [
        (
            "text_search_only"
            if sources and "clinicaltrials" in sources.lower()
            else "source_reported"
            if sources
            else "missing"
        )
        for sources in clinical_sources
    ]
    text_search_only = axes["clinical_evidence_quality"] == "text_search_only"
    if text_search_only.any():
        axes.loc[text_search_only, "clinical_precedent"] = np.minimum(
            axes.loc[text_search_only, "clinical_precedent"],
            10.0,
        )
    axes["clinical_evidence_cap_applied"] = text_search_only
    raw_total = axes[list(AXES)].sum(axis=1)
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
        "axis_availability": {
            axis: int(
                (axes[f"{axis}_evidence_status"] == "available").sum()
            )
            for axis in AXES
        },
        "median_adjusted_score": float(axes["adjusted_score"].median()),
        "median_evidence_completeness": float(
            axes["evidence_completeness"].median()
        ),
        "top_targets": axes.head(20)[
            [
                "validation_rank",
                "gene",
                "decision",
                "adjusted_score",
                "evidence_completeness",
                *AXES,
                "safety_penalty",
            ]
        ].to_dict(orient="records"),
    }
    write_json(out_dir / "target_validation_summary.json", summary)
    return axes, summary
