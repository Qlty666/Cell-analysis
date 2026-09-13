"""Sample-level differential abundance statistics."""

from __future__ import annotations

import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd

from docking.utils import write_json

from .errors import IntegrationError

log = logging.getLogger("full_pipeline")


def _chi2_contingency(table) -> tuple[float, float]:
    """Pearson chi-square (Yates-corrected 2x2) with the erfc formula."""
    matrix = np.asarray(table, dtype=float)
    if matrix.shape != (2, 2):
        return float("nan"), float("nan")
    total = float(matrix.sum())
    if total == 0 or np.any(matrix < 0):
        return float("nan"), float("nan")
    row_totals = matrix.sum(axis=1)
    col_totals = matrix.sum(axis=0)
    expected = np.outer(row_totals, col_totals) / total
    with np.errstate(divide="ignore", invalid="ignore"):
        corrected = float(
            np.nansum(
                (np.abs(matrix - expected) - 0.5) ** 2 / expected
            )
        )
    chi2 = max(0.0, corrected)
    if chi2 == 0:
        return 0.0, 1.0
    p_value = math.erfc(math.sqrt(chi2 / 2.0))
    return chi2, min(max(p_value, 0.0), 1.0)


def _bh_adjust(p_values) -> list[float]:
    """Benjamini-Hochberg FDR adjustment without extra dependencies."""
    values = np.asarray(p_values, dtype=float)
    count = len(values)
    if count == 0:
        return []
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * count / np.arange(1, count + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return [float(value) for value in out]


def _betacf(
    a: float,
    b: float,
    x: float,
    max_iter: int = 200,
) -> float:
    """Continued fraction for the regularized incomplete beta function."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3e-12:
            break
    return h


def _regularized_incomplete_beta(
    a: float,
    b: float,
    x: float,
) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - _regularized_incomplete_beta(b, a, 1.0 - x)


def _welch_ttest(sample_a, sample_b) -> tuple[float, float]:
    """Two-sided Welch t-test without scipy; returns (t, p) or nan values."""
    a = np.asarray(sample_a, dtype=float)
    b = np.asarray(sample_b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return float("nan"), float("nan")
    var_a = float(a.var(ddof=1))
    var_b = float(b.var(ddof=1))
    standard_error = var_a / n_a + var_b / n_b
    if not math.isfinite(standard_error) or standard_error <= 0.0:
        return float("nan"), float("nan")
    t_value = float(
        (a.mean() - b.mean()) / math.sqrt(standard_error)
    )
    if t_value == 0.0:
        return 0.0, 1.0
    df = standard_error**2 / (
        (var_a / n_a) ** 2 / (n_a - 1)
        + (var_b / n_b) ** 2 / (n_b - 1)
    )
    if not math.isfinite(df) or df <= 0.0:
        return t_value, float("nan")
    x = df / (df + t_value * t_value)
    p_value = _regularized_incomplete_beta(df / 2.0, 0.5, x)
    return t_value, float(min(max(p_value, 0.0), 1.0))


def run_differential_abundance(
    single_cell_root: Path,
    out_dir: Path,
    config: dict | None = None,
) -> dict:
    """Pair DE with a sample-level cell-type composition shift test."""
    cfg = config or {}
    out_dir.mkdir(parents=True, exist_ok=True)
    ann_path = (
        single_cell_root
        / "results"
        / "data"
        / "04_annotation"
        / "fig_05_16_17_cell_annotations.csv"
    )
    if not ann_path.exists():
        return {
            "status": "skipped",
            "reason": "cell annotation CSV not found",
            "output_csv": "",
        }
    try:
        annotations = pd.read_csv(ann_path)
    except (OSError, ValueError) as exc:
        log.warning(
            "could not read cell annotation CSV %s: %s",
            ann_path,
            exc,
        )
        return {
            "status": "skipped",
            "reason": "cell annotation CSV unreadable",
            "output_csv": "",
        }
    required = {"celltype_annot", "condition", "sample"}
    if not required.issubset(annotations.columns):
        return {
            "status": "skipped",
            "reason": (
                "annotation CSV lacks celltype_annot/condition/sample columns"
            ),
            "output_csv": "",
        }
    annotations = annotations.dropna(
        subset=["celltype_annot", "condition", "sample"]
    ).copy()
    annotations["celltype_annot"] = annotations[
        "celltype_annot"
    ].astype(str)
    annotations["condition"] = annotations["condition"].astype(str)
    annotations["sample"] = annotations["sample"].astype(str)
    conditions = sorted(annotations["condition"].unique().tolist())
    if len(conditions) != 2:
        raise IntegrationError(
            "differential abundance requires exactly two conditions; found "
            f"{len(conditions)}: {conditions}"
        )
    cond0, cond1 = conditions

    multi_condition = annotations.groupby("sample")[
        "condition"
    ].nunique()
    ambiguous = multi_condition[multi_condition > 1].index.tolist()
    if ambiguous:
        log.warning(
            "excluding %s sample(s) that map to more than one condition: %s",
            len(ambiguous),
            ", ".join(str(name) for name in ambiguous[:5]),
        )
        annotations = annotations[
            ~annotations["sample"].isin(ambiguous)
        ]
    sample_condition = annotations.groupby("sample")["condition"].first()
    sample_counts = pd.crosstab(
        annotations["sample"],
        annotations["celltype_annot"],
    )
    sample_totals = sample_counts.sum(axis=1).replace(0, np.nan)
    sample_props = sample_counts.div(sample_totals, axis=0)

    counts = pd.crosstab(
        annotations["celltype_annot"],
        annotations["condition"],
    )
    if cond0 not in counts.columns or cond1 not in counts.columns:
        return {
            "status": "skipped",
            "reason": "condition columns missing from crosstab",
            "output_csv": "",
        }

    total0 = int(counts[cond0].sum())
    total1 = int(counts[cond1].sum())
    min_cells = int(cfg.get("min_cells", 5))
    mask0 = (sample_condition == cond0).reindex(
        sample_props.index,
        fill_value=False,
    )
    mask1 = (sample_condition == cond1).reindex(
        sample_props.index,
        fill_value=False,
    )
    rows: list[dict] = []
    for celltype, row in counts.iterrows():
        a = int(row.get(cond0, 0))
        b = int(row.get(cond1, 0))
        if a + b < min_cells:
            continue
        if celltype in sample_props.columns:
            props0 = sample_props.loc[mask0, celltype].dropna().tolist()
            props1 = sample_props.loc[mask1, celltype].dropna().tolist()
        else:
            props0, props1 = [], []
        mean0 = (
            float(np.mean(props0))
            if props0
            else (a / total0 if total0 else 0.0)
        )
        mean1 = (
            float(np.mean(props1))
            if props1
            else (b / total1 if total1 else 0.0)
        )
        t_stat, p_value = _welch_ttest(props0, props1)
        fraction0 = a / total0 if total0 else 0.0
        fraction1 = b / total1 if total1 else 0.0
        rows.append(
            {
                "celltype": str(celltype),
                f"{cond0}_cells": a,
                f"{cond1}_cells": b,
                f"{cond0}_fraction": round(fraction0, 6),
                f"{cond1}_fraction": round(fraction1, 6),
                "n_cells": a + b,
                f"n_samples_{cond0}": len(props0),
                f"n_samples_{cond1}": len(props1),
                "t_stat": (
                    None if math.isnan(t_stat) else round(t_stat, 6)
                ),
                "p_value": (
                    None if math.isnan(p_value) else p_value
                ),
                "direction": (
                    "enriched_in_" + cond0
                    if mean0 > mean1
                    else "enriched_in_" + cond1
                ),
            }
        )

    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["p_adjust"] = _bh_adjust(frame["p_value"].fillna(1.0))
        frame["significant"] = frame["p_adjust"] < float(
            cfg.get("fdr", 0.05)
        )
        frame = frame.sort_values(
            ["p_adjust", "p_value"],
            na_position="last",
        ).reset_index(drop=True)

    csv_path = out_dir / "differential_abundance.csv"
    frame.to_csv(csv_path, index=False)
    samples0 = int(mask0.sum())
    samples1 = int(mask1.sum())
    summary = {
        "status": "completed",
        "test": "welch_t_test_on_per_sample_proportions",
        "conditions": [cond0, cond1],
        "samples_per_condition": {
            cond0: samples0,
            cond1: samples1,
        },
        "celltypes_tested": int(len(frame)),
        "significant_celltypes": int(
            frame["significant"].sum() if not frame.empty else 0
        ),
        "output_csv": str(csv_path),
    }
    if samples0 < 2 or samples1 < 2:
        log.warning(
            "differential abundance has fewer than two samples in a condition "
            "(samples per condition: %s); p-values are not estimable",
            summary["samples_per_condition"],
        )
    write_json(
        out_dir / "differential_abundance_summary.json",
        summary,
    )
    log.info(
        "differential abundance: %s cell types tested across %s/%s samples, "
        "%s significant",
        summary["celltypes_tested"],
        samples0,
        samples1,
        summary["significant_celltypes"],
    )
    return summary
