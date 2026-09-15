"""Evidence-aware target prioritisation for the full pipeline.

This module keeps the DEG ranking and multi-source evidence ranking separate.
Observed components are combined with a neutral prior for unavailable
components, so missing evidence is neither treated as a negative result nor
allowed to inflate a target score.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd


DEFAULT_WEIGHTS = {
    "expression": 0.25,
    "evidence": 0.45,
    "knockout": 0.20,
    "advanced": 0.10,
}


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _percentile(series: pd.Series) -> pd.Series:
    values = _numeric(series)
    if values.notna().sum() == 0:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return values.rank(method="average", pct=True)


def _first_present(frame: pd.DataFrame, names: list[str]) -> pd.Series:
    for name in names:
        if name in frame.columns:
            return frame[name]
    return pd.Series(np.nan, index=frame.index, dtype=float)


def _normalise_weights(weights: Mapping[str, float] | None) -> dict[str, float]:
    raw = dict(DEFAULT_WEIGHTS)
    raw.update(
        {
            str(key): float(value)
            for key, value in (weights or {}).items()
            if key in DEFAULT_WEIGHTS
        }
    )
    total = sum(max(0.0, value) for value in raw.values())
    if total <= 0:
        return dict(DEFAULT_WEIGHTS)
    return {
        key: max(0.0, value) / total
        for key, value in raw.items()
    }


def _component_score(
    row: pd.Series,
    components: Mapping[str, str],
    weights: Mapping[str, float],
) -> tuple[float, str]:
    observed: list[tuple[str, float, float]] = []
    for component, column in components.items():
        value = pd.to_numeric(
            pd.Series([row.get(column)]),
            errors="coerce",
        ).iloc[0]
        if pd.notna(value):
            observed.append(
                (component, float(value), float(weights.get(component, 0.0)))
            )
    observed = [item for item in observed if item[2] > 0]
    weight_sum = sum(
        max(0.0, float(weight))
        for weight in weights.values()
    )
    if weight_sum <= 0:
        return 0.0, "no configured evidence components"
    observed_names = {name for name, _, _ in observed}
    score = 0.0
    for name, weight in weights.items():
        weight = max(0.0, float(weight))
        value = next(
            (value for item_name, value, _ in observed if item_name == name),
            0.5,
        )
        score += weight * value
    score = score / weight_sum
    used = "observed=" + (
        "+".join(name for name, _, _ in observed) or "none"
    )
    missing = [name for name in weights if name not in observed_names]
    if missing:
        used += ";neutral_prior=" + "+".join(missing)
    return float(np.clip(score, 0.0, 1.0)), used


def _decision(
    score: float,
    coverage: float,
    category_count: float,
    *,
    go_min_score: float,
    go_min_coverage: float,
    conditional_min_score: float,
) -> str:
    if (
        score >= go_min_score
        and coverage >= go_min_coverage
        and category_count >= 3
    ):
        return "GO"
    if score >= conditional_min_score and (
        coverage >= 0.25 or category_count >= 2
    ):
        return "CONDITIONAL_GO"
    return "REVIEW"


def build_target_priority(
    candidate_universe_csv: Path,
    evidence_priority_csv: Path | None,
    output_csv: Path,
    *,
    knockout_csv: Path | None = None,
    weights: Mapping[str, float] | None = None,
    evidence_queried: bool | None = None,
    go_min_score: float = 0.75,
    go_min_coverage: float = 0.50,
    conditional_min_score: float = 0.50,
) -> pd.DataFrame:
    """Merge expression, evidence and perturbation scores into one ranking."""
    candidates = pd.read_csv(candidate_universe_csv)
    if candidates.empty:
        raise ValueError("candidate universe is empty")
    if "gene" not in candidates.columns:
        raise ValueError("candidate universe requires a gene column")

    candidates = candidates.copy()
    candidates["gene"] = candidates["gene"].astype(str).str.strip().str.upper()
    candidates = candidates[candidates["gene"] != ""]
    candidates = candidates.drop_duplicates("gene", keep="first")
    if "candidate_origin" not in candidates.columns:
        candidates["candidate_origin"] = "DEG"
    candidates["candidate_origin"] = (
        candidates["candidate_origin"]
        .fillna("DEG")
        .astype(str)
        .replace({"": "DEG"})
    )

    if "deg_rank" in candidates.columns:
        rank = _numeric(candidates["deg_rank"])
        max_rank = max(1, int(rank.max(skipna=True)) if rank.notna().any() else 1)
        rank_score = 1.0 - (rank - 1.0) / max(1.0, max_rank - 1.0)
    else:
        rank_score = _percentile(
            -pd.Series(np.arange(len(candidates)), index=candidates.index)
        )
    effect_score = _percentile(
        _first_present(
            candidates,
            ["abs_log2fc", "avg_log2fc", "log2fc"],
        ).abs()
    )
    candidates["expression_score"] = (
        0.6 * rank_score.fillna(0.5) + 0.4 * effect_score.fillna(0.5)
    ).clip(0.0, 1.0)

    evidence = pd.DataFrame()
    if evidence_priority_csv is not None and Path(evidence_priority_csv).exists():
        try:
            evidence = pd.read_csv(evidence_priority_csv)
        except pd.errors.EmptyDataError:
            evidence = pd.DataFrame()
        if "target_symbol" in evidence.columns:
            evidence = evidence.rename(columns={"target_symbol": "gene"})
        if "gene" in evidence.columns:
            evidence["gene"] = (
                evidence["gene"].astype(str).str.strip().str.upper()
            )
            evidence = evidence.drop_duplicates("gene", keep="first")
        else:
            evidence = pd.DataFrame()

    if not evidence.empty:
        evidence = evidence.rename(
            columns={"priority_score": "evidence_priority_score"}
        )
        candidates = candidates.merge(evidence, on="gene", how="left")
    else:
        candidates["evidence_priority_score"] = np.nan

    if "evidence_score" not in candidates.columns:
        candidates["evidence_score"] = np.nan
    if "coverage_ratio" not in candidates.columns:
        candidates["coverage_ratio"] = np.nan
    if "category_count" not in candidates.columns:
        candidates["category_count"] = np.nan
    if "missing_categories" not in candidates.columns:
        candidates["missing_categories"] = ""
    if "rank_stability" not in candidates.columns:
        candidates["rank_stability"] = np.nan

    if evidence_queried is None:
        evidence_queried = not evidence.empty
    candidates["evidence_status"] = np.where(
        not evidence_queried,
        "not_queried",
        np.where(
            candidates["evidence_priority_score"].notna(),
            "covered",
            "not_found",
        ),
    )

    candidates["knockout_score"] = np.nan
    if knockout_csv is not None and Path(knockout_csv).exists():
        knockout = pd.read_csv(knockout_csv)
        gene_column = next(
            (
                column
                for column in ("gene", "target", "symbol")
                if column in knockout.columns
            ),
            None,
        )
        score_column = next(
            (
                column
                for column in (
                    "target_score",
                    "knockout_score",
                    "integrated_score",
                )
                if column in knockout.columns
            ),
            None,
        )
        if gene_column and score_column:
            knockout = knockout[[gene_column, score_column]].rename(
                columns={gene_column: "gene", score_column: "knockout_score"}
            )
            knockout["gene"] = (
                knockout["gene"].astype(str).str.strip().str.upper()
            )
            knockout["knockout_score"] = _numeric(knockout["knockout_score"])
            knockout = knockout.drop_duplicates("gene", keep="first")
            candidates = candidates.drop(columns=["knockout_score"]).merge(
                knockout,
                on="gene",
                how="left",
            )

    candidates["advanced_priority_score"] = _numeric(
        candidates.get("advanced_priority_score", np.nan)
    )
    weights_used = _normalise_weights(weights)
    components = {
        "expression": "expression_score",
        "evidence": "evidence_priority_score",
        "knockout": "knockout_score",
        "advanced": "advanced_priority_score",
    }
    integrated: list[float] = []
    reasons: list[str] = []
    for _, row in candidates.iterrows():
        score, used = _component_score(row, components, weights_used)
        integrated.append(score)
        reasons.append(used)
    candidates["integrated_score"] = integrated
    candidates["selection_reason"] = reasons
    candidates["decision"] = [
        _decision(
            float(row["integrated_score"]),
            float(row["coverage_ratio"])
            if pd.notna(row["coverage_ratio"])
            else 0.0,
            float(row["category_count"])
            if pd.notna(row["category_count"])
            else 0.0,
            go_min_score=go_min_score,
            go_min_coverage=go_min_coverage,
            conditional_min_score=conditional_min_score,
        )
        for _, row in candidates.iterrows()
    ]
    candidates = candidates.sort_values(
        ["integrated_score", "coverage_ratio", "expression_score", "gene"],
        ascending=[False, False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    candidates.insert(0, "priority_rank", np.arange(1, len(candidates) + 1))
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output_csv, index=False)
    return candidates


def write_target_priority_summary(
    frame: pd.DataFrame,
    output_json: Path,
    *,
    evidence_hub_used: bool,
    evidence_source_count: int = 0,
    evidence_targets_queried: int = 0,
    benchmark: dict | None = None,
) -> dict:
    """Write concise ranking and evidence-coverage diagnostics."""
    coverage = _numeric(frame.get("coverage_ratio", pd.Series(dtype=float)))
    summary = {
        "targets": int(len(frame)),
        "evidence_hub_used": bool(evidence_hub_used),
        "evidence_source_count": int(evidence_source_count),
        "evidence_targets_queried": int(evidence_targets_queried),
        "evidence_covered": int(
            (frame.get("evidence_status", pd.Series(dtype=str)) == "covered").sum()
        ),
        "evidence_not_found": int(
            (frame.get("evidence_status", pd.Series(dtype=str)) == "not_found").sum()
        ),
        "evidence_not_queried": int(
            (frame.get("evidence_status", pd.Series(dtype=str)) == "not_queried").sum()
        ),
        "median_coverage_ratio": (
            float(coverage.median()) if coverage.notna().any() else None
        ),
        "decisions": {
            str(key): int(value)
            for key, value in frame.get("decision", pd.Series(dtype=str))
            .value_counts()
            .items()
        },
        "top_targets": frame.head(20)[
            [
                column
                for column in (
                    "priority_rank",
                    "gene",
                    "candidate_origin",
                    "integrated_score",
                    "decision",
                    "evidence_status",
                    "coverage_ratio",
                )
                if column in frame.columns
            ]
        ].to_dict(orient="records"),
        "candidate_origins": {
            str(key): int(value)
            for key, value in frame.get(
                "candidate_origin",
                pd.Series(dtype=str),
            )
            .value_counts()
            .items()
        },
        "benchmark": benchmark or {},
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return summary
