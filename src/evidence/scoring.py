"""Source-aware target prioritisation with missing-evidence semantics."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Mapping

import numpy as np
import pandas as pd

from .models import EvidenceTier


DEFAULT_CATEGORY_WEIGHTS: dict[str, float] = {
    "direct_experimental": 0.24,
    "curated_target": 0.14,
    "predicted_target": 0.06,
    "genetic_association": 0.16,
    "disease_association": 0.12,
    "liver_context": 0.10,
    "dependency": 0.10,
    "pathway": 0.05,
    "structure": 0.03,
    "clinical_precedent": 0.08,
    "genetic_causal": 0.12,
}

CATEGORY_LABELS = {
    "direct_experimental": "Direct experimental compound-target evidence",
    "curated_target": "Curated compound-target evidence",
    "predicted_target": "Computational compound-target prediction",
    "genetic_association": "Human genetic disease association",
    "disease_association": "Curated or text-mined disease association",
    "liver_context": "Liver expression or cell-context evidence",
    "dependency": "Functional dependency evidence",
    "pathway": "Pathway membership or pathway evidence",
    "structure": "Experimental or predicted structure evidence",
    "clinical_precedent": "Clinical trial or translational precedent",
    "safety_risk": "Known toxicity or safety liability",
    "genetic_causal": "Human genetic causal or colocalization evidence",
}


@dataclass(frozen=True)
class ScoringConfig:
    weights: Mapping[str, float]
    evidence_weight: float = 0.8
    coverage_weight: float = 0.2
    safety_penalty_weight: float = 0.15

    @classmethod
    def from_dict(cls, value: Mapping | None) -> "ScoringConfig":
        raw = dict(value or {})
        weights = dict(DEFAULT_CATEGORY_WEIGHTS)
        weights.update(
            {
                key: float(weight)
                for key, weight in (raw.get("weights") or {}).items()
                if key in DEFAULT_CATEGORY_WEIGHTS
            }
        )
        total = sum(max(0.0, weight) for weight in weights.values())
        if total <= 0:
            raise ValueError("scoring weights must contain a positive value")
        weights = {
            key: max(0.0, weight) / total
            for key, weight in weights.items()
        }
        evidence_weight = float(raw.get("evidence_weight", 0.8))
        coverage_weight = float(raw.get("coverage_weight", 0.2))
        if evidence_weight < 0 or coverage_weight < 0:
            raise ValueError("scoring weights must be non-negative")
        total_mix = evidence_weight + coverage_weight
        if total_mix <= 0:
            raise ValueError("evidence_weight + coverage_weight must be positive")
        return cls(
            weights=weights,
            evidence_weight=evidence_weight / total_mix,
            coverage_weight=coverage_weight / total_mix,
            safety_penalty_weight=max(
                0.0,
                float(raw.get("safety_penalty_weight", 0.15)),
            ),
        )


def _evidence_category(row: Mapping) -> str:
    relation = str(row.get("relation") or "").strip().lower()
    subject_type = str(row.get("subject_type") or "").strip().lower()
    evidence_type = str(row.get("evidence_type") or "").strip().lower()
    tier_raw = row.get("tier")
    try:
        tier = EvidenceTier.parse(tier_raw)
    except (ValueError, KeyError):
        tier = EvidenceTier.CURATED

    if relation in {"targets", "binds", "perturbs"} and tier == EvidenceTier.EXPERIMENTAL:
        return "direct_experimental"
    if relation in {"targets", "binds", "perturbs"} and tier == EvidenceTier.CURATED:
        return "curated_target"
    if relation in {"targets", "binds", "perturbs"} and tier == EvidenceTier.PREDICTED:
        return "predicted_target"
    if (
        "safety" in evidence_type
        or "toxicity" in evidence_type
        or "adverse" in evidence_type
        or relation in {"toxic_to", "causes"}
    ):
        return "safety_risk"
    if (
        "clinical" in evidence_type
        or relation in {"studied_in", "clinical_precedent"}
    ):
        return "clinical_precedent"
    if (
        "coloc" in evidence_type
        or "eqtl" in evidence_type
        or "pqtl" in evidence_type
        or "sqlt" in evidence_type
        or "twass" in evidence_type
        or relation in {"colocalizes_with", "causally_associated_with"}
    ):
        return "genetic_causal"
    if tier == EvidenceTier.GENETIC or "genetic" in evidence_type:
        return "genetic_association"
    if subject_type == "disease" and relation in {
        "associated_with",
        "implicated_in",
    }:
        return "disease_association"
    if relation in {"expressed_in", "enriched_in"}:
        return "liver_context"
    if relation in {"dependent_on", "essential_in"}:
        return "dependency"
    if relation in {"participates_in", "member_of"}:
        return "pathway"
    if relation in {"has_structure", "structured_as"}:
        return "structure"
    return ""


def _safe_score(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return float(np.clip(number, 0.0, 1.0))


def _category_frame(evidence: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty:
        return pd.DataFrame(
            columns=[
                "target_symbol",
                "category",
                "category_score",
                "n_records",
                "n_source_groups",
                "source_groups",
            ]
        )
    frame = evidence.copy()
    frame["target_symbol"] = frame["target_symbol"].astype(str).str.upper().str.strip()
    frame = frame[frame["target_symbol"] != ""]
    frame["category"] = frame.apply(_evidence_category, axis=1)
    frame = frame[frame["category"] != ""]
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "target_symbol",
                "category",
                "category_score",
                "n_records",
                "n_source_groups",
                "source_groups",
            ]
        )
    frame["score"] = frame["score"].map(_safe_score)
    # Presence of an auditable relation is retained even when the source does
    # not expose a numeric score. It receives a neutral value rather than zero.
    frame["score"] = frame["score"].fillna(0.5)
    if "source_group" in frame.columns:
        source_group = frame["source_group"].fillna("").astype(str).str.strip()
    else:
        source_group = pd.Series("", index=frame.index, dtype=str)
    frame["source_group"] = source_group.replace("", np.nan)
    frame["source_group"] = frame["source_group"].fillna(
        frame["source"].astype(str)
    )
    source_level = (
        frame.groupby(
            ["target_symbol", "category", "source_group"],
            as_index=False,
        )
        .agg(source_score=("score", "max"), source_records=("score", "size"))
    )
    grouped = (
        source_level.groupby(
            ["target_symbol", "category"],
            as_index=False,
        )
        .agg(
            category_score=("source_score", "mean"),
            n_records=("source_records", "sum"),
            n_source_groups=("source_group", "nunique"),
            source_groups=(
                "source_group",
                lambda values: ";".join(sorted(set(map(str, values)))),
            ),
        )
    )
    return grouped


def _score_from_category_frame(
    category_frame: pd.DataFrame,
    config: ScoringConfig,
) -> pd.DataFrame:
    targets = sorted(category_frame["target_symbol"].dropna().unique())
    rows: list[dict] = []
    for target in targets:
        target_rows = category_frame[category_frame["target_symbol"] == target]
        scores = {
            str(row["category"]): float(row["category_score"])
            for row in target_rows.to_dict("records")
        }
        available_weight = sum(
            config.weights.get(category, 0.0) for category in scores
        )
        total_weight = sum(config.weights.values())
        evidence_score = (
            sum(
                config.weights.get(category, 0.0) * score
                for category, score in scores.items()
            )
            / available_weight
            if available_weight > 0
            else 0.0
        )
        coverage_ratio = available_weight / total_weight if total_weight else 0.0
        safety_risk = float(scores.get("safety_risk", 0.0) or 0.0)
        safety_penalty = (
            config.safety_penalty_weight * safety_risk
        )
        priority_score = max(
            0.0,
            config.evidence_weight * evidence_score
            + config.coverage_weight * coverage_ratio
            - safety_penalty,
        )
        row = {
            "target_symbol": target,
            "priority_score": priority_score,
            "evidence_score": evidence_score,
            "coverage_ratio": coverage_ratio,
            "category_count": int(len(scores)),
            "source_group_count": int(
                pd.to_numeric(
                    target_rows["n_source_groups"],
                    errors="coerce",
                ).fillna(0).max()
            ),
            "evidence_record_count": int(
                pd.to_numeric(
                    target_rows["n_records"],
                    errors="coerce",
                ).fillna(0).sum()
            ),
            "missing_categories": ";".join(
                category
                for category in config.weights
                if category not in scores
            ),
            "safety_risk": safety_risk,
            "safety_penalty": safety_penalty,
        }
        for category in config.weights:
            row[category] = scores.get(category, np.nan)
            row[f"{category}_sources"] = (
                target_rows.loc[
                    target_rows["category"] == category,
                    "source_groups",
                ].iloc[0]
                if category in scores
                else ""
            )
        rows.append(row)
    output = pd.DataFrame(rows)
    if output.empty:
        return output
    return output.sort_values(
        ["priority_score", "coverage_ratio", "target_symbol"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def _rank_stability(
    evidence: pd.DataFrame,
    config: ScoringConfig,
) -> pd.DataFrame:
    full = _score_from_category_frame(_category_frame(evidence), config)
    if full.empty:
        return pd.DataFrame(
            columns=[
                "target_symbol",
                "ablation_rank_sd",
                "max_rank_delta",
                "rank_stability",
                "ablation_runs",
            ]
        )
    full_ranks = {
        row.target_symbol: rank
        for rank, row in enumerate(full.itertuples(index=False), start=1)
    }
    source_groups = sorted(
        {
            str(value)
            for value in evidence.get("source_group", pd.Series(dtype=str))
            .fillna("")
            .astype(str)
            if str(value).strip()
        }
    )
    ranks_by_target: dict[str, list[int]] = {
        target: [rank] for target, rank in full_ranks.items()
    }
    for group in source_groups:
        subset = evidence[
            evidence["source_group"].fillna("").astype(str) != group
        ]
        ablated = _score_from_category_frame(_category_frame(subset), config)
        ablated_ranks = {
            row.target_symbol: rank
            for rank, row in enumerate(ablated.itertuples(index=False), start=1)
        }
        for target in full_ranks:
            ranks_by_target[target].append(
                ablated_ranks.get(target, len(ablated) + 1)
            )
    rows = []
    for target, ranks in ranks_by_target.items():
        max_rank = max(len(full_ranks), max(ranks))
        sd = float(np.std(ranks, ddof=1)) if len(ranks) > 1 else 0.0
        rows.append(
            {
                "target_symbol": target,
                "ablation_rank_sd": sd,
                "max_rank_delta": int(
                    max(abs(rank - full_ranks[target]) for rank in ranks)
                ),
                "rank_stability": float(
                    max(0.0, 1.0 - sd / max(1.0, max_rank / 2.0))
                ),
                "ablation_runs": len(ranks),
            }
        )
    return pd.DataFrame(rows)


def score_targets(
    evidence: pd.DataFrame,
    *,
    config: Mapping | ScoringConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return priority table, target-category matrix and source ablation table."""
    scoring_config = (
        config
        if isinstance(config, ScoringConfig)
        else ScoringConfig.from_dict(config or {})
    )
    category_frame = _category_frame(evidence)
    priority = _score_from_category_frame(category_frame, scoring_config)
    stability = _rank_stability(evidence, scoring_config)
    if not priority.empty:
        priority = priority.merge(stability, on="target_symbol", how="left")
        priority["rank"] = np.arange(1, len(priority) + 1)
        columns = ["rank", "target_symbol"] + [
            column
            for column in priority.columns
            if column not in {"rank", "target_symbol"}
        ]
        priority = priority[columns]
    matrix = category_frame.pivot_table(
        index="target_symbol",
        columns="category",
        values="category_score",
        aggfunc="mean",
    ).reset_index()
    return priority, matrix, stability


def benchmark_ranking(
    priority: pd.DataFrame,
    *,
    positive_targets: Iterable[str],
    negative_targets: Iterable[str] | None = None,
    top_n: int = 20,
    benchmark_source: str = "",
    benchmark_version: str = "",
    benchmark_independent: bool = False,
    exclude_sources: Iterable[str] | None = None,
    source_groups_by_target: Mapping[str, Iterable[str]] | None = None,
) -> dict:
    """Evaluate recovery of a reference target set without inventing labels."""
    positives = {
        str(value).upper().strip()
        for value in positive_targets
        if str(value).strip()
    }
    negatives = {
        str(value).upper().strip()
        for value in (negative_targets or [])
        if str(value).strip()
    }
    overlapping_reference_labels = sorted(positives & negatives)
    negatives = negatives - positives
    excluded_sources = {
        str(value).strip().lower()
        for value in (exclude_sources or [])
        if str(value).strip()
    }
    source_groups = {
        str(target).upper().strip(): {
            str(value).strip().lower()
            for value in values
            if str(value).strip()
        }
        for target, values in (source_groups_by_target or {}).items()
        if str(target).strip()
    }
    for values in source_groups.values():
        values.difference_update(excluded_sources)
    ranked = [
        str(value).upper()
        for value in priority.get("target_symbol", pd.Series(dtype=str)).tolist()
    ]
    rank_map = {target: rank for rank, target in enumerate(ranked, start=1)}
    score_column = next(
        (
            column
            for column in ("priority_score", "integrated_score", "score")
            if column in priority.columns
        ),
        None,
    )
    if score_column:
        score_values = pd.to_numeric(
            priority[score_column],
            errors="coerce",
        )
        score_map = {
            str(target).upper(): (
                float(score)
                if pd.notna(score)
                else 1.0 / math.log1p(rank_map.get(
                    str(target).upper(),
                    len(rank_map) + 1,
                ))
            )
            for target, score in zip(
                priority.get("target_symbol", pd.Series(dtype=str)),
                score_values,
            )
        }
        score_metric = str(score_column)
    else:
        score_map = {
            target: 1.0 / math.log1p(rank)
            for target, rank in rank_map.items()
        }
        score_metric = "reciprocal_rank"
    top_n = max(1, int(top_n))
    recovered_top_n = [target for target in ranked[:top_n] if target in positives]
    positive_ranks = [
        rank_map[target] for target in positives if target in rank_map
    ]
    result = {
        "n_positive_reference": len(positives),
        "n_positive_with_evidence": len(positive_ranks),
        "n_ranked_targets": len(ranked),
        "recall_at_n": (
            len(recovered_top_n) / len(positives) if positives else None
        ),
        "top_n": int(top_n),
        "recovered_top_n": recovered_top_n,
        "median_positive_rank": (
            float(np.median(positive_ranks)) if positive_ranks else None
        ),
        "precision_at_n": (
            len(recovered_top_n) / min(top_n, len(ranked))
            if ranked
            else None
        ),
        "f1_at_n": None,
        "enrichment_factor": None,
        "benchmark_source": str(benchmark_source or ""),
        "benchmark_version": str(benchmark_version or ""),
        "benchmark_independent": bool(benchmark_independent),
        "benchmark_independent_verified": False,
        "benchmark_evidence_mode": (
            "source_holdout" if excluded_sources else "not_established"
        ),
        "excluded_sources": sorted(excluded_sources),
        "source_leakage_checked": bool(source_groups or excluded_sources),
        "source_leakage": None,
        "source_overlap": [],
        "positive_source_groups": {
            target: sorted(source_groups.get(target, set()))
            for target in sorted(positives & set(rank_map))
        },
        "negative_source_groups": {
            target: sorted(source_groups.get(target, set()))
            for target in sorted(negatives & set(rank_map))
        },
        "overlapping_reference_labels": overlapping_reference_labels,
        "score_metric": score_metric,
    }
    labeled = [
        (target, 1)
        for target in sorted(positives)
        if target in rank_map
    ]
    labeled.extend(
        (target, 0)
        for target in sorted(negatives)
        if target in rank_map
    )
    result["n_positive"] = int(sum(label for _, label in labeled))
    result["n_negative"] = int(
        sum(1 for _, label in labeled if label == 0)
    )
    result["n_missing_positive"] = int(
        len([target for target in positives if target not in rank_map])
    )
    result["n_missing_negative"] = int(
        len([target for target in negatives if target not in rank_map])
    )
    reference_source_groups = {
        source
        for target, _ in labeled
        for source in source_groups.get(target, set())
    }
    ranking_source_groups = {
        source
        for target in ranked
        for source in source_groups.get(target, set())
    }
    source_overlap = sorted(
        reference_source_groups & ranking_source_groups
    )
    result["source_overlap"] = source_overlap
    if source_groups:
        result["source_leakage"] = bool(source_overlap)
    result["benchmark_independent_verified"] = bool(
        benchmark_independent and result["source_leakage"] is not True
    )
    if result["precision_at_n"] is not None:
        recall_at_n = result["recall_at_n"]
        if recall_at_n is not None and (
            result["precision_at_n"] + recall_at_n
        ) > 0:
            result["f1_at_n"] = float(
                2
                * result["precision_at_n"]
                * recall_at_n
                / (result["precision_at_n"] + recall_at_n)
            )
    if len({label for _, label in labeled}) < 2:
        result["auroc"] = None
        result["auprc"] = None
        result["auroc_ci_low"] = None
        result["auroc_ci_high"] = None
        result["auprc_ci_low"] = None
        result["auprc_ci_high"] = None
        result["permutation_p_value"] = None
        result["permutation_n"] = 0
        result["n_labeled"] = len(labeled)
        return result
    labels = np.asarray([label for _, label in labeled], dtype=int)
    scores = np.asarray(
        [score_map.get(target, 0.0) for target, _ in labeled],
        dtype=float,
    )
    auroc = _binary_auroc(labels, scores)
    result["auroc"] = auroc
    result["n_labeled"] = len(labeled)
    prevalence = float(np.mean(labels))
    if prevalence > 0 and result["precision_at_n"] is not None:
        result["enrichment_factor"] = float(
            result["precision_at_n"] / prevalence
        )

    order = np.argsort(-scores)
    sorted_labels = labels[order]
    cumulative_positives = np.cumsum(sorted_labels)
    ranks = np.arange(1, len(sorted_labels) + 1)
    precision = cumulative_positives / ranks
    recall = cumulative_positives / max(1, int(cumulative_positives[-1]))
    average_precision = 0.0
    previous_recall = 0.0
    for value, recall_value in zip(precision, recall):
        average_precision += (recall_value - previous_recall) * value
        previous_recall = recall_value
    result["auprc"] = float(average_precision)

    rng = np.random.default_rng(42)
    boot = []
    boot_auprc = []
    permutation_scores = []
    for _ in range(1000):
        indices = rng.integers(0, len(labels), size=len(labels))
        sampled_labels = labels[indices]
        sampled_scores = scores[indices]
        if len(set(sampled_labels)) < 2:
            continue
        boot.append(_binary_auroc(sampled_labels, sampled_scores))
        order = np.argsort(-sampled_scores)
        cumulative = np.cumsum(sampled_labels[order])
        precision_values = cumulative / np.arange(
            1,
            len(cumulative) + 1,
        )
        recall_values = cumulative / max(1, int(cumulative[-1]))
        ap = 0.0
        previous = 0.0
        for precision_value, recall_value in zip(
            precision_values,
            recall_values,
        ):
            ap += (recall_value - previous) * precision_value
            previous = recall_value
        boot_auprc.append(ap)
        permuted = rng.permutation(labels)
        permutation_scores.append(_binary_auroc(permuted, scores))
    if boot:
        result["auroc_ci_low"] = float(np.percentile(boot, 2.5))
        result["auroc_ci_high"] = float(np.percentile(boot, 97.5))
    else:
        result["auroc_ci_low"] = None
        result["auroc_ci_high"] = None
    if boot_auprc:
        result["auprc_ci_low"] = float(np.percentile(boot_auprc, 2.5))
        result["auprc_ci_high"] = float(np.percentile(boot_auprc, 97.5))
    else:
        result["auprc_ci_low"] = None
        result["auprc_ci_high"] = None
    if permutation_scores:
        result["permutation_p_value"] = float(
            (
                1
                + sum(
                    value >= auroc
                    for value in permutation_scores
                )
            )
            / (len(permutation_scores) + 1)
        )
        result["permutation_n"] = len(permutation_scores)
    else:
        result["permutation_p_value"] = None
        result["permutation_n"] = 0
    return result


def _binary_auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    order = np.argsort(scores)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    positive_ranks = ranks[labels == 1]
    negative_count = int(np.sum(labels == 0))
    if len(positive_ranks) == 0 or negative_count == 0:
        return float("nan")
    return float(
        (
            positive_ranks.sum()
            - len(positive_ranks)
            * (len(positive_ranks) + 1)
            / 2.0
        )
        / (len(positive_ranks) * negative_count)
    )
