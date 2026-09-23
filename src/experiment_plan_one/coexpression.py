"""Reference-informed co-expression module screening for bulk transcriptomes.

This is deliberately labelled as a WGCNA-style signed co-expression analysis,
not as the original R WGCNA implementation. It keeps the reproducible parts
used by the reference papers: variance filtering, soft-threshold selection,
module-trait association, module membership and hub export.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA

from .common import (
    bh_fdr,
    ensure_dir,
    read_expression,
    save_figure,
    write_json,
)

LOG = logging.getLogger("experiment_plan_one.coexpression")

CASE_TOKENS = ("nash", "nafld", "masld", "case", "disease", "ss")
CONTROL_TOKENS = ("hc", "healthy", "normal", "control")


def _condition_value(value: Any) -> int | None:
    text = str(value or "").strip().lower()
    if any(token in text for token in CONTROL_TOKENS):
        return 0
    if any(token in text for token in CASE_TOKENS):
        return 1
    return None


def _scale_free_fit(connectivity: np.ndarray, bins: int = 10) -> float:
    values = np.asarray(connectivity, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < bins * 2:
        return float("nan")
    counts, edges = np.histogram(values, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2.0
    valid = (counts > 0) & (centers > 0)
    if valid.sum() < 3:
        return float("nan")
    x = np.log10(centers[valid])
    y = np.log10(counts[valid])
    slope, intercept = np.polyfit(x, y, 1)
    residual = y - (slope * x + intercept)
    total = y - y.mean()
    if np.sum(total**2) <= 0:
        return float("nan")
    return float(1.0 - np.sum(residual**2) / np.sum(total**2))


def _select_soft_power(
    correlation: np.ndarray,
    candidate_powers: list[int],
) -> tuple[int, pd.DataFrame, bool]:
    rows: list[dict[str, float]] = []
    signed_correlation = np.clip((1.0 + correlation) / 2.0, 0.0, 1.0)
    for power in candidate_powers:
        adjacency = np.power(signed_correlation, power)
        connectivity = np.maximum(adjacency.sum(axis=1) - 1.0, 0.0)
        rows.append(
            {
                "power": int(power),
                "scale_free_r2": _scale_free_fit(connectivity),
                "mean_connectivity": float(np.mean(connectivity)),
            }
        )
    frame = pd.DataFrame(rows)
    eligible = frame[
        frame["scale_free_r2"].ge(0.8)
        & frame["mean_connectivity"].between(10.0, 1000.0)
    ]
    if not eligible.empty:
        selected = int(eligible.iloc[0]["power"])
        target_met = True
    else:
        low_connectivity = frame[frame["mean_connectivity"] <= 100.0]
        if not low_connectivity.empty:
            selected = int(
                low_connectivity.sort_values(
                    ["scale_free_r2", "power"],
                    ascending=[False, True],
                    na_position="last",
                ).iloc[0]["power"]
            )
        else:
            selected = int(
                frame.sort_values(
                    ["mean_connectivity", "power"],
                    ascending=[True, True],
                ).iloc[0]["power"]
            )
        target_met = False
    return selected, frame, target_met


def _module_eigengenes(
    matrix: np.ndarray,
    labels: np.ndarray,
) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    eigengenes: dict[int, np.ndarray] = {}
    scores: dict[int, np.ndarray] = {}
    for module in sorted(set(labels)):
        indices = np.flatnonzero(labels == module)
        if len(indices) < 3:
            continue
        module_matrix = matrix[indices, :]
        component = PCA(n_components=1, random_state=42).fit_transform(
            module_matrix.T
        )[:, 0]
        average = module_matrix.mean(axis=0)
        if np.corrcoef(component, average)[0, 1] < 0:
            component = -component
        eigengenes[int(module)] = component
        scores[int(module)] = component
    return eigengenes, scores


def run_coexpression_analysis(
    expression_path: Path,
    metadata_path: Path,
    output_dir: Path,
    *,
    condition_column: str = "condition",
    max_genes: int = 2000,
    min_module_size: int = 20,
    distance_threshold: float = 0.5,
    kme_threshold: float = 0.5,
    top_hubs_per_module: int = 25,
    candidate_powers: list[int] | None = None,
) -> dict[str, Any]:
    """Run a transparent signed co-expression module screen.

    The implementation is intentionally lightweight enough to run inside the
    local pipeline. It does not claim to reproduce WGCNA's topological
    overlap matrix or dynamic tree cut algorithm.
    """
    expression = read_expression(expression_path)
    metadata = pd.read_csv(metadata_path, index_col=0)
    samples = [
        sample
        for sample in expression.columns
        if sample in metadata.index
    ]
    if len(samples) < 6:
        raise ValueError("co-expression screening requires at least six matched samples")
    labels = metadata.loc[samples, condition_column].map(_condition_value)
    valid = labels.notna()
    samples = labels[valid].index.tolist()
    labels = labels.loc[samples].astype(int)
    if labels.nunique() != 2:
        raise ValueError(
            "co-expression screening requires both disease and control samples"
        )
    matrix = expression[samples].apply(pd.to_numeric, errors="coerce")
    variance = matrix.var(axis=1, skipna=True)
    matrix = matrix.loc[variance.notna() & variance.gt(0)]
    selected_count = min(max(50, int(max_genes)), len(matrix))
    if selected_count < 50:
        raise ValueError("fewer than 50 variable genes are available")
    mad = matrix.sub(matrix.median(axis=1), axis=0).abs().median(axis=1)
    matrix = matrix.loc[mad.sort_values(ascending=False).head(selected_count).index]
    values = matrix.to_numpy(dtype=float)
    values = (
        values - np.nanmean(values, axis=1, keepdims=True)
    ) / np.nanstd(values, axis=1, keepdims=True)
    values = np.nan_to_num(values, nan=0.0)
    correlation = np.corrcoef(values)
    correlation = np.nan_to_num(correlation, nan=0.0)
    np.fill_diagonal(correlation, 1.0)

    powers = candidate_powers or list(range(1, 21))
    selected_power, threshold_frame, scale_free_target_met = _select_soft_power(
        correlation,
        powers,
    )
    distance = np.clip(1.0 - correlation, 0.0, 2.0)
    np.fill_diagonal(distance, 0.0)
    hierarchy = linkage(squareform(distance, checks=False), method="average")
    raw_labels = fcluster(
        hierarchy,
        t=float(distance_threshold),
        criterion="distance",
    )
    counts = pd.Series(raw_labels).value_counts()
    keep = counts[counts >= int(min_module_size)].index
    remap = {int(module): index for index, module in enumerate(keep, start=1)}
    module_labels = np.array(
        [remap.get(int(module), 0) for module in raw_labels],
        dtype=int,
    )
    if module_labels.max(initial=0) < 1:
        raise RuntimeError("no co-expression module reached the minimum size")
    eigengenes, _ = _module_eigengenes(values, module_labels)
    if not eigengenes:
        raise RuntimeError("module eigengenes could not be calculated")

    trait = labels.to_numpy(dtype=float)
    trait_rows: list[dict[str, Any]] = []
    for module, eigengene in eigengenes.items():
        if np.std(eigengene) == 0:
            correlation_value = 0.0
            p_value = 1.0
        else:
            correlation_value, p_value = stats.pearsonr(eigengene, trait)
        trait_rows.append(
            {
                "module": f"M{module}",
                "n_genes": int(np.sum(module_labels == module)),
                "module_trait_r": float(correlation_value),
                "p_value": float(p_value),
            }
        )
    trait_frame = pd.DataFrame(trait_rows)
    trait_frame["fdr"] = bh_fdr(trait_frame["p_value"])
    trait_frame = trait_frame.sort_values(
        ["fdr", "module_trait_r"],
        ascending=[True, False],
    )
    selected_modules = set(
        trait_frame.loc[
            trait_frame["fdr"] < 0.05,
            "module",
        ].astype(str)
    )
    trait_lookup = trait_frame.set_index("module")

    gene_rows: list[dict[str, Any]] = []
    for index, gene in enumerate(matrix.index):
        module_number = int(module_labels[index])
        if module_number == 0:
            continue
        module = f"M{module_number}"
        eigengene = eigengenes[module_number]
        if np.std(values[index, :]) == 0 or np.std(eigengene) == 0:
            kme = 0.0
        else:
            kme = float(np.corrcoef(values[index, :], eigengene)[0, 1])
        gene_rows.append(
            {
                "gene": str(gene),
                "module": module,
                "kME": kme,
                "abs_kME": abs(kme),
                "module_trait_r": float(
                    trait_lookup.loc[module, "module_trait_r"]
                ),
                "module_fdr": float(trait_lookup.loc[module, "fdr"]),
                "disease_associated_module": module in selected_modules,
            }
        )
    gene_frame = pd.DataFrame(gene_rows)
    if gene_frame.empty:
        raise RuntimeError("no genes could be assigned to a retained module")
    hub_frame = (
        gene_frame[gene_frame["disease_associated_module"]]
        .sort_values(["module", "abs_kME"], ascending=[True, False])
        .groupby("module", as_index=False, group_keys=False)
        .head(max(1, int(top_hubs_per_module)))
        .reset_index(drop=True)
    )
    strong = hub_frame[hub_frame["abs_kME"] >= float(kme_threshold)]
    hub_frame = strong.reset_index(drop=True)
    module_gene_frame = (
        gene_frame[gene_frame["disease_associated_module"]]
        .sort_values(["module", "abs_kME"], ascending=[True, False])
        .reset_index(drop=True)
    )

    out_dir = ensure_dir(output_dir)
    threshold_frame.to_csv(
        out_dir / "coexpression_soft_threshold.csv",
        index=False,
    )
    trait_frame.to_csv(
        out_dir / "coexpression_module_trait.csv",
        index=False,
    )
    gene_frame.to_csv(
        out_dir / "coexpression_gene_modules.csv",
        index=False,
    )
    hub_frame.to_csv(
        out_dir / "coexpression_hubs.csv",
        index=False,
    )
    module_gene_frame.to_csv(
        out_dir / "coexpression_disease_module_genes.csv",
        index=False,
    )
    _plot_module_summary(
        threshold_frame,
        trait_frame,
        selected_power,
        hub_frame,
        out_dir / "fig_supp_coexpression_modules.png",
    )
    summary = {
        "status": "completed",
        "scientific_outcome": "significant_modules" if selected_modules else "valid_negative",
        "inference_role": "exploratory; not a preselected ML feature set",
        "method": "signed correlation adjacency, WGCNA-style",
        "implementation_note": (
            "This is a reproducible local co-expression screen, not the "
            "original R WGCNA implementation with topological overlap and "
            "dynamic tree cut."
        ),
        "n_samples": int(len(samples)),
        "n_genes_analyzed": int(len(matrix)),
        "selected_soft_power": int(selected_power),
        "scalefree_r2": float(
            threshold_frame.loc[
                threshold_frame["power"] == selected_power,
                "scale_free_r2",
            ].iloc[0]
        ),
        "scalefree_target_met": bool(scale_free_target_met),
        "n_modules": int(gene_frame["module"].nunique()),
        "disease_associated_modules": sorted(selected_modules),
        "n_disease_module_genes": int(len(module_gene_frame)),
        "n_hubs": int(len(hub_frame)),
        "outputs": {
            "soft_threshold": str(out_dir / "coexpression_soft_threshold.csv"),
            "module_trait": str(out_dir / "coexpression_module_trait.csv"),
            "gene_modules": str(out_dir / "coexpression_gene_modules.csv"),
            "hubs": str(out_dir / "coexpression_hubs.csv"),
            "disease_module_genes": str(
                out_dir / "coexpression_disease_module_genes.csv"
            ),
        },
    }
    write_json(out_dir / "coexpression_summary.json", summary)
    LOG.info(
        "co-expression screening retained %s modules and %s hubs",
        summary["n_modules"],
        summary["n_hubs"],
    )
    return summary


def _plot_module_summary(
    threshold: pd.DataFrame,
    trait: pd.DataFrame,
    selected_power: int,
    hubs: pd.DataFrame,
    output: Path,
) -> None:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.8),
        gridspec_kw={"width_ratios": [0.9, 1.1]},
    )
    axes[0].plot(
        threshold["power"],
        threshold["scale_free_r2"],
        marker="o",
        color="#2f6bb3",
        lw=1.5,
    )
    axes[0].axhline(0.8, color="#c05b4d", linestyle="--", lw=1)
    axes[0].axvline(selected_power, color="#c05b4d", linestyle=":", lw=1)
    axes[0].set_xlabel("Soft-threshold power")
    axes[0].set_ylabel("Scale-free fit R2")
    axes[0].set_title("Soft-threshold selection", fontweight="bold")

    ordered = trait.sort_values("module_trait_r")
    colors = ["#4f7fa8" if value < 0 else "#c05b4d" for value in ordered["module_trait_r"]]
    axes[1].barh(
        ordered["module"],
        ordered["module_trait_r"],
        color=colors,
    )
    axes[1].axvline(0, color="#8a949e", lw=0.8)
    axes[1].set_xlabel("Module-trait correlation")
    axes[1].set_ylabel("")
    axes[1].set_title("Disease-associated modules", fontweight="bold")
    if not hubs.empty:
        top = hubs.sort_values("abs_kME", ascending=False).head(8)
        axes[1].text(
            0.02,
            -0.22,
            "Top hubs: " + ", ".join(top["gene"].astype(str).tolist()),
            transform=axes[1].transAxes,
            ha="left",
            va="top",
            fontsize=6.5,
        )
    fig.tight_layout()
    save_figure(fig, output)
