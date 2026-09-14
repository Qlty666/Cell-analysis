"""Single-cell mouse discovery and human validation analyses."""

from __future__ import annotations

import gzip
import json
import logging
import os
import re
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from scipy.io import mmread

from .common import bh_fdr, ensure_dir, save_figure, slug, write_json

LOG = logging.getLogger("experiment_plan_one.single_cell")

MOUSE_MARKERS = {
    "Hepatocytes": ["Alb", "Apoa1", "Apoa2", "Apoe", "Ttr", "Hnf4a", "Asgr1"],
    "Cholangiocytes": ["Krt19", "Krt7", "Epcam", "Sox9"],
    "Kupffer cells": ["Clec4f", "Cd68", "Adgre1", "Csf1r"],
    "Monocytes": ["Lyz2", "Ccr2", "Itgam", "Fcgr3"],
    "Endothelial cells": ["Pecam1", "Kdr", "Cdh5", "Vwf"],
    "Hepatic stellate cells": ["Col1a1", "Col1a2", "Acta2", "Des", "Pdgfrb", "Rgs5"],
    "T/NK cells": ["Cd3d", "Cd3e", "Nkg7", "Klrb1c", "Trbc2"],
    "B cells": ["Cd79a", "Ms4a1", "Cd74"],
    "Neutrophils": ["S100a8", "S100a9", "Ly6g"],
}

HUMAN_MARKERS = {
    "Hepatocytes": ["ALB", "APOA1", "APOA2", "APOE", "TTR", "HNF4A", "ASGR1"],
    "Cholangiocytes": ["KRT19", "KRT7", "EPCAM", "SOX9"],
    "Macrophages": ["CD68", "CLEC4F", "ADGRE1", "CSF1R"],
    "Monocytes": ["LYZ", "CCR2", "ITGAM", "FCGR3A"],
    "Endothelial cells": ["PECAM1", "KDR", "CDH5", "VWF"],
    "Hepatic stellate cells": ["COL1A1", "COL1A2", "ACTA2", "DES", "PDGFRB", "RGS5"],
    "T/NK cells": ["CD3D", "CD3E", "NKG7", "KLRB1", "TRBC2"],
    "B/Plasma cells": ["CD79A", "MS4A1", "MZB1", "IGHG1"],
    "Neutrophils": ["S100A8", "S100A9", "LY6G"],
}

LIGAND_RECEPTOR_DB = [
    ("TNF", "TNF", "TNFRSF1A"),
    ("TNF", "TNF", "TNFRSF1B"),
    ("IL6", "IL6", "IL6R"),
    ("IL6", "IL6", "IL6ST"),
    ("TGFB", "TGFB1", "TGFBR1"),
    ("TGFB", "TGFB1", "TGFBR2"),
    ("CCL", "CCL2", "CCR2"),
    ("CCL", "CCL5", "CCR5"),
    ("CXCL", "CXCL12", "CXCR4"),
    ("PDGF", "PDGFB", "PDGFRB"),
    ("NOTCH", "JAG1", "NOTCH2"),
    ("WNT", "WNT2", "FZD1"),
    ("BMP", "BMP2", "BMPR1A"),
    ("Lipid transport", "APOE", "LRP1"),
    ("Lipid transport", "APOC3", "LRP1"),
]


def parse_geo_soft_samples(path: Path) -> dict[str, dict[str, Any]]:
    """Parse sample-level fields from a GEO family SOFT archive."""
    records: dict[str, dict[str, Any]] = {}
    current: dict[str, Any] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE ="):
                if current and current.get("gsm"):
                    records[str(current["gsm"])] = current
                current = {"gsm": line.split("=", 1)[1].strip()}
                continue
            if current is None or not line.startswith("!"):
                continue
            key, _, value = line[1:].partition("=")
            key = key.strip()
            value = value.strip()
            if key == "Sample_title":
                current["title"] = value
            elif key == "Sample_characteristics_ch1":
                if ":" in value:
                    feature, feature_value = value.split(":", 1)
                    current[slug(feature)] = feature_value.strip()
            elif key == "Sample_supplementary_file":
                current.setdefault("supplementary", []).append(value)
    if current and current.get("gsm"):
        records[str(current["gsm"])] = current
    return records


def load_mouse_dataset(
    extracted_dir: Path,
    soft_path: Path,
    *,
    min_genes: int = 200,
    max_mito_pct: float = 20.0,
) -> ad.AnnData:
    """Load the four GSE270583 10x libraries into one AnnData object."""
    metadata = parse_geo_soft_samples(soft_path)
    matrix_files = sorted(extracted_dir.rglob("*_filtered_feature_bc_matrix_matrix.mtx.gz"))
    if not matrix_files:
        matrix_files = sorted(extracted_dir.rglob("*filtered*matrix.mtx.gz"))
    if not matrix_files:
        raise FileNotFoundError(
            f"no filtered 10x matrix files under {extracted_dir}"
        )
    objects: list[ad.AnnData] = []
    for matrix_file in matrix_files:
        prefix = matrix_file.name[: -len("matrix.mtx.gz")]
        barcode_file = _find_sibling(matrix_file, prefix + "barcodes.tsv.gz")
        feature_file = _find_sibling(matrix_file, prefix + "features.tsv.gz")
        if barcode_file is None or feature_file is None:
            LOG.warning("skipping incomplete 10x triplet: %s", matrix_file)
            continue
        matrix = mmread(matrix_file).tocsr().T
        barcodes = pd.read_csv(barcode_file, header=None, sep="\t").iloc[:, 0].astype(str)
        features = pd.read_csv(feature_file, header=None, sep="\t", dtype=str)
        gene_ids = features.iloc[:, 0].astype(str)
        symbols = features.iloc[:, 1].astype(str).str.upper()
        obsm = pd.DataFrame(index=[f"{prefix}{barcode}" for barcode in barcodes])
        varm = pd.DataFrame({"gene_id": gene_ids.to_numpy()}, index=symbols.to_numpy())
        sample = ad.AnnData(X=matrix, obs=obsm, var=varm)
        sample.obs["sample_id"] = prefix
        gsm = prefix.split("_", 1)[0]
        record = metadata.get(gsm, {})
        treatment = str(record.get("treatment") or record.get("title") or "unknown")
        if "HFD" in treatment.upper():
            condition = "HFD"
        elif "NCD" in treatment.upper():
            condition = "NCD"
        elif "JQF" in treatment.upper():
            condition = "JQF"
        else:
            condition = treatment
        sample.obs["condition"] = condition
        sample.obs["batch"] = prefix
        sample.var_names_make_unique()
        sample.obs_names_make_unique()
        objects.append(sample)
    if not objects:
        raise RuntimeError("no complete 10x libraries were loaded")
    data = ad.concat(objects, join="outer", fill_value=0, index_unique=None)
    data.var_names_make_unique()

    data.var["mt"] = data.var_names.astype(str).str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], inplace=True, percent_top=None)
    sc.pp.filter_cells(data, min_genes=min_genes)
    sc.pp.filter_genes(data, min_cells=3)
    data = data[data.obs["pct_counts_mt"].astype(float) <= max_mito_pct].copy()
    data.layers["counts"] = data.X.copy()
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    data.raw = data
    sc.pp.highly_variable_genes(data, n_top_genes=min(2500, data.n_vars), batch_key="batch")
    hvg = data[:, data.var["highly_variable"]].copy()
    _integrate_batches(hvg)
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=min(40, hvg.n_vars - 1, hvg.n_obs - 1), svd_solver="arpack")
    sc.pp.neighbors(hvg, n_neighbors=15, n_pcs=min(30, hvg.obsm["X_pca"].shape[1]))
    sc.tl.umap(hvg, random_state=42)
    sc.tl.leiden(hvg, resolution=0.7, key_added="cluster", flavor="igraph", n_iterations=2)
    data.obsm["X_pca"] = hvg.obsm["X_pca"]
    data.obsm["X_umap"] = hvg.obsm["X_umap"]
    data.obs["cluster"] = hvg.obs["cluster"].astype(str)
    data.obs["cell_type"] = annotate_clusters(data, MOUSE_MARKERS)
    data.obs["cell_type"] = (
        data.obs["cell_type"]
        .astype(str)
        .str.replace(r"\s+\d+$", "", regex=True)
    )
    return data


def _find_sibling(source: Path, filename: str) -> Path | None:
    candidate = source.with_name(filename)
    if candidate.exists():
        return candidate
    matches = list(source.parent.glob(filename))
    return matches[0] if matches else None


def _integrate_batches(data: ad.AnnData) -> None:
    if data.obs["batch"].nunique() < 2:
        return
    try:
        sc.external.pp.harmony_integrate(data, "batch", max_iter_harmony=30)
        data.obsm["X_pca"] = data.obsm["X_pca_harmony"]
        LOG.info("Harmony integration completed")
    except Exception as exc:  # noqa: BLE001
        LOG.warning("Harmony unavailable/failed; using ComBat fallback: %s", exc)
        try:
            sc.pp.combat(data, key="batch")
        except Exception as combat_exc:  # noqa: BLE001
            LOG.warning("ComBat failed; continuing without batch correction: %s", combat_exc)


def annotate_clusters(
    data: ad.AnnData,
    marker_sets: dict[str, list[str]],
    *,
    cluster_column: str = "cluster",
) -> pd.Series:
    """Assign cluster labels from mean marker expression after z-scoring."""
    normalized = data.raw.to_adata() if data.raw is not None else data
    gene_lookup = {str(gene).upper(): str(gene) for gene in normalized.var_names}
    score_rows: list[dict[str, float]] = []
    clusters = normalized.obs[cluster_column].astype(str)
    for cluster in sorted(clusters.unique(), key=_natural_key):
        mask = clusters == cluster
        row: dict[str, float] = {"cluster": cluster, "n_cells": int(mask.sum())}
        for cell_type, markers in marker_sets.items():
            present = [
                gene_lookup[gene.upper()]
                for gene in markers
                if gene.upper() in gene_lookup
            ]
            if not present:
                row[cell_type] = np.nan
                continue
            matrix = normalized[mask, present].X
            if sparse.issparse(matrix):
                values = np.asarray(matrix.mean(axis=0)).ravel()
            else:
                values = np.asarray(matrix.mean(axis=0)).ravel()
            row[cell_type] = float(np.nanmean(values))
        score_rows.append(row)
    score_frame = pd.DataFrame(score_rows).set_index("cluster")
    numeric = score_frame.drop(columns=["n_cells"], errors="ignore").copy()
    for column in numeric.columns:
        values = pd.to_numeric(numeric[column], errors="coerce")
        mean = values.mean(skipna=True)
        std = values.std(skipna=True, ddof=0)
        numeric[column] = (
            (values - mean) / std if pd.notna(std) and std > 0 else 0.0
        )
    numeric = numeric.fillna(0.0)
    cluster_labels = numeric.idxmax(axis=1).to_dict()
    unique_labels = _make_unique_labels(cluster_labels)
    return clusters.map(unique_labels)


def _natural_key(value: str) -> tuple[Any, ...]:
    return tuple(int(token) if token.isdigit() else token for token in re.split(r"(\d+)", value))


def _make_unique_labels(mapping: dict[str, str]) -> dict[str, str]:
    counts: dict[str, int] = {}
    output: dict[str, str] = {}
    for cluster, label in mapping.items():
        counts[label] = counts.get(label, 0) + 1
        output[cluster] = (
            label if counts[label] == 1 else f"{label} {counts[label]}"
        )
    return output


def run_mouse_single_cell(
    extracted_dir: Path,
    soft_path: Path,
    core_genes: list[str],
    output_dir: Path,
    *,
    sample_key: str | None = None,
) -> dict[str, Any]:
    ensure_dir(output_dir)
    data_path = output_dir / "mouse_liver_processed.h5ad"
    if data_path.exists() and not sample_key:
        data = ad.read_h5ad(data_path)
    else:
        data = load_mouse_dataset(extracted_dir, soft_path)
        data.write_h5ad(data_path, compression="gzip")
    data.obs["cell_type"] = annotate_clusters(data, MOUSE_MARKERS)
    data.obs["cell_type"] = (
        data.obs["cell_type"]
        .astype(str)
        .str.replace(r"\s+\d+$", "", regex=True)
    )
    data.obs["condition"] = pd.Categorical(
        data.obs["condition"].astype(str),
        categories=["NCD", "HFD", "JQF"],
    )
    _plot_umap(
        data,
        color="cell_type",
        output=output_dir / "fig4a_umap_cell_types.png",
        title="Mouse liver single-cell atlas",
        legend_title="Cell type",
    )
    _plot_marker_dotplot(data, MOUSE_MARKERS, output_dir / "fig4b_cell_type_markers.png")
    core = _resolve_genes(data.var_names, core_genes)
    gene_frames: list[pd.DataFrame] = []
    for gene in core:
        values = _expression_vector(data, gene)
        frame = pd.DataFrame(
            {
                "cell": data.obs_names,
                "gene": gene,
                "expression": values,
                "cell_type": data.obs["cell_type"].astype(str).to_numpy(),
                "condition": data.obs["condition"].astype(str).to_numpy(),
                "sample_id": data.obs["sample_id"].astype(str).to_numpy(),
            }
        )
        gene_frames.append(frame)
    gene_expression = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    gene_expression.to_csv(output_dir / "core_gene_expression_long.csv.gz", index=False, compression="gzip")
    mouse_stats = _mouse_core_statistics(gene_expression)
    mouse_stats.to_csv(output_dir / "core_gene_condition_stats.csv", index=False)
    if core:
        _plot_core_violin(
            gene_expression,
            output_dir / "fig4c_core_gene_violin.png",
        )
        _plot_feature_grid(data, core, output_dir / "fig4d_core_gene_umap.png")
    composition = (
        data.obs.groupby(["sample_id", "condition", "cell_type"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    composition["fraction"] = composition["n_cells"] / composition.groupby(
        "sample_id",
        observed=True,
    )["n_cells"].transform("sum")
    composition.to_csv(output_dir / "cell_composition.csv", index=False)
    _plot_composition(composition, output_dir / "fig4e_cell_composition.png")

    rank = sc.tl.rank_genes_groups(
        data,
        groupby="condition",
        groups=["HFD"],
        reference="NCD",
        method="wilcoxon",
        use_raw=True,
        pts=True,
    )
    differential = sc.get.rank_genes_groups_df(data, group="HFD")
    differential.to_csv(output_dir / "hfd_vs_ncd_cell_level_deg.csv", index=False)

    interaction_path = output_dir / "cellchat_like_interactions.csv"
    pathway_path = output_dir / "cellchat_like_pathways.csv"
    if interaction_path.exists() and pathway_path.exists():
        lr = {
            "interactions": pd.read_csv(interaction_path),
            "pathways": pd.read_csv(pathway_path),
        }
    else:
        lr = _cellchat_like_analysis(data)
        lr["interactions"].to_csv(interaction_path, index=False)
        lr["pathways"].to_csv(pathway_path, index=False)
    _plot_cell_communication(
        lr["interactions"],
        lr["pathways"],
        output_dir / "fig4f_cellchat_network.png",
    )
    write_json(
        output_dir / "mouse_single_cell_summary.json",
        {
            "n_cells": int(data.n_obs),
            "n_genes": int(data.n_vars),
            "samples": data.obs["sample_id"].value_counts().to_dict(),
            "conditions": data.obs["condition"].value_counts().to_dict(),
            "cell_types": data.obs["cell_type"].value_counts().to_dict(),
            "core_genes": core,
            "cellchat_note": (
                "CellChat-like ligand-receptor scoring was performed with an "
                "explicit local receptor-pair table; it is a transparent "
                "screening approximation, not the full CellChat permutation test."
            ),
        },
    )
    return {
        "h5ad": data_path,
        "cellchat": lr,
        "core_genes": core,
        "composition": output_dir / "cell_composition.csv",
    }


def _resolve_genes(var_names: pd.Index, genes: list[str]) -> list[str]:
    lookup = {str(gene).upper(): str(gene) for gene in var_names}
    return [lookup[str(gene).upper()] for gene in genes if str(gene).upper() in lookup]


def _expression_vector(data: ad.AnnData, gene: str) -> np.ndarray:
    source = data
    if data.raw is not None and gene in data.raw.var_names:
        source = data.raw.to_adata()[:, gene]
    else:
        source = data[:, gene]
    values = source.X
    if sparse.issparse(values):
        values = values.toarray().ravel()
    return np.asarray(values).ravel()


def _plot_umap(
    data: ad.AnnData,
    *,
    color: str,
    output: Path,
    title: str,
    legend_title: str,
) -> None:
    frame = pd.DataFrame(data.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=data.obs_names)
    frame[color] = data.obs[color].astype(str)
    categories = sorted(frame[color].unique())
    palette = plt.get_cmap("tab20", max(len(categories), 1))
    color_map = {category: palette(index) for index, category in enumerate(categories)}
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for category in categories:
        subset = frame[frame[color] == category]
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            s=5,
            alpha=0.72,
            linewidths=0,
            color=color_map[category],
            label=category,
        )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(title, fontweight="bold")
    ax.legend(
        title=legend_title,
        markerscale=2.5,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=True,
        framealpha=0.9,
        edgecolor="#c7d0d8",
        borderpad=0.35,
        labelspacing=0.3,
    )
    save_figure(fig, output)


def _plot_marker_dotplot(
    data: ad.AnnData,
    marker_sets: dict[str, list[str]],
    output: Path,
) -> None:
    genes = [gene for values in marker_sets.values() for gene in values]
    genes = _resolve_genes(data.var_names, genes)
    if not genes:
        return
    matrix = data[:, genes]
    if sparse.issparse(matrix.X):
        values = matrix.X.toarray()
    else:
        values = np.asarray(matrix.X)
    cell_types = data.obs["cell_type"].astype(str)
    rows = []
    for marker_group, markers in marker_sets.items():
        for gene in _resolve_genes(data.var_names, markers):
            for cell_type in sorted(cell_types.unique()):
                subset = values[(cell_types == cell_type).to_numpy(), genes.index(gene)]
                rows.append(
                    {
                        "marker_group": marker_group,
                        "gene": gene,
                        "cell_type": cell_type,
                        "mean_expression": float(np.mean(subset)),
                        "pct_expression": float(np.mean(subset > 0) * 100),
                    }
                )
    frame = pd.DataFrame(rows)
    fig, ax = plt.subplots(
        figsize=(7.2, max(4.2, frame["cell_type"].nunique() * 0.42))
    )
    cell_order = sorted(frame["cell_type"].unique())
    gene_order = list(dict.fromkeys(frame["gene"]))
    x = {gene: index for index, gene in enumerate(gene_order)}
    y = {cell: index for index, cell in enumerate(cell_order)}
    ax.scatter(
        [x[gene] for gene in frame["gene"]],
        [y[cell] for cell in frame["cell_type"]],
        s=np.clip(frame["pct_expression"] * 1.5, 4, 150),
        c=frame["mean_expression"],
        cmap="YlOrRd",
        edgecolors="white",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=90, fontsize=6.2)
    ax.set_yticks(range(len(cell_order)))
    ax.set_yticklabels(cell_order, fontsize=6)
    ax.set_xlabel("Canonical marker gene")
    ax.set_ylabel("Cell type")
    ax.set_title("Canonical marker expression", fontweight="bold")
    save_figure(fig, output)


def _plot_core_violin(gene_expression: pd.DataFrame, output: Path) -> None:
    frame = gene_expression[gene_expression["condition"].isin(["NCD", "HFD"])].copy()
    if frame.empty:
        return
    genes = list(dict.fromkeys(frame["gene"]))
    fig, axes = plt.subplots(
        1,
        len(genes),
        figsize=(7.2, 3.6),
        squeeze=False,
    )
    for index, gene in enumerate(genes):
        ax = axes.flat[index]
        subset = frame[frame["gene"] == gene]
        groups = ["NCD", "HFD"]
        values = [subset.loc[subset["condition"] == group, "expression"].to_numpy() for group in groups]
        parts = ax.violinplot(values, showmeans=False, showextrema=False)
        for body, color in zip(parts["bodies"], ["#6b9ac4", "#c56b4a"]):
            body.set_facecolor(color)
            body.set_alpha(0.7)
        ax.boxplot(values, widths=0.18, showfliers=False)
        p_value = "NA"
        if min(map(len, values)) >= 2:
            p_value = f"{stats.mannwhitneyu(values[0], values[1]).pvalue:.2g}"
        ax.set_xticks([1, 2])
        ax.set_xticklabels(groups)
        ax.set_title(gene)
        ax.text(0.98, 0.98, f"p={p_value}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    fig.suptitle("Core-gene expression in mouse NAFLD", fontweight="bold")
    save_figure(fig, output)


def _plot_feature_grid(data: ad.AnnData, genes: list[str], output: Path) -> None:
    if not genes:
        return
    coords = data.obsm["X_umap"]
    fig, axes = plt.subplots(
        1,
        len(genes),
        figsize=(7.2, 3.6),
        squeeze=False,
    )
    for ax, gene in zip(axes.flat, genes):
        values = _expression_vector(data, gene)
        points = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=4, cmap="viridis", linewidths=0)
        ax.set_title(gene)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.colorbar(points, ax=ax, shrink=0.7)
    save_figure(fig, output)


def _plot_composition(composition: pd.DataFrame, output: Path) -> None:
    frame = (
        composition[composition["condition"].isin(["NCD", "HFD"])]
        .groupby(["condition", "cell_type"], observed=True)["n_cells"]
        .sum()
        .unstack(fill_value=0)
    )
    if frame.empty:
        return
    proportions = frame.div(frame.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    bottom = np.zeros(len(proportions))
    colors = plt.get_cmap("tab20", frame.shape[1])
    for index, cell_type in enumerate(frame.columns):
        values = proportions[cell_type].to_numpy()
        ax.bar(proportions.index, values, bottom=bottom, label=cell_type, color=colors(index))
        bottom += values
    ax.set_ylabel("Cell fraction")
    ax.set_title("Cell composition in mouse liver", fontweight="bold")
    ax.legend(
        bbox_to_anchor=(0.5, -0.22),
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=6,
    )
    save_figure(fig, output)


def _cellchat_like_analysis(data: ad.AnnData) -> dict[str, pd.DataFrame]:
    """Score explicit ligand-receptor pairs in NCD and HFD conditions."""
    normalized = data.raw.to_adata() if data.raw is not None else data
    expression = normalized.X
    cell_types = sorted(data.obs["cell_type"].astype(str).unique())
    conditions = ["NCD", "HFD"]
    records: list[dict[str, Any]] = []
    for pathway, ligand, receptor in LIGAND_RECEPTOR_DB:
        ligand_gene = _resolve_genes(normalized.var_names, [ligand])
        receptor_gene = _resolve_genes(normalized.var_names, [receptor])
        if not ligand_gene or not receptor_gene:
            continue
        ligand_values = _gene_matrix_column(expression, normalized.var_names, ligand_gene[0])
        receptor_values = _gene_matrix_column(expression, normalized.var_names, receptor_gene[0])
        for condition in conditions:
            condition_mask = (data.obs["condition"].astype(str) == condition).to_numpy()
            for source in cell_types:
                source_mask = condition_mask & (data.obs["cell_type"].astype(str) == source).to_numpy()
                if source_mask.sum() < 5:
                    continue
                ligand_mean = float(np.mean(ligand_values[source_mask]))
                for target in cell_types:
                    target_mask = condition_mask & (data.obs["cell_type"].astype(str) == target).to_numpy()
                    if target_mask.sum() < 5:
                        continue
                    receptor_mean = float(np.mean(receptor_values[target_mask]))
                    score = min(ligand_mean, receptor_mean)
                    records.append(
                        {
                            "pathway": pathway,
                            "ligand": ligand_gene[0],
                            "receptor": receptor_gene[0],
                            "condition": condition,
                            "source": source,
                            "target": target,
                            "communication_score": score,
                        }
                    )
    interactions = pd.DataFrame(records)
    if interactions.empty:
        return {"interactions": interactions, "pathways": pd.DataFrame()}
    wide = interactions.pivot_table(
        index=["pathway", "ligand", "receptor", "source", "target"],
        columns="condition",
        values="communication_score",
        fill_value=0.0,
    ).reset_index()
    for condition in conditions:
        if condition not in wide:
            wide[condition] = 0.0
    wide["delta_HFD_NCD"] = wide["HFD"] - wide["NCD"]
    pathways = (
        wide.groupby("pathway", as_index=False)
        .agg(
            NCD=("NCD", "sum"),
            HFD=("HFD", "sum"),
            delta_HFD_NCD=("delta_HFD_NCD", "sum"),
            n_pairs=("pathway", "size"),
        )
        .sort_values("delta_HFD_NCD", ascending=False)
    )
    return {"interactions": wide, "pathways": pathways}


def _gene_matrix_column(matrix: Any, var_names: pd.Index, gene: str) -> np.ndarray:
    index = int(var_names.get_loc(gene))
    values = matrix[:, index]
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def _plot_cell_communication(
    interactions: pd.DataFrame,
    pathways: pd.DataFrame,
    output: Path,
) -> None:
    if pathways.empty:
        return
    if interactions.empty or "delta_HFD_NCD" not in interactions.columns:
        return
    import networkx as nx

    frame = interactions.copy()
    frame["absolute_delta"] = pd.to_numeric(
        frame["delta_HFD_NCD"],
        errors="coerce",
    ).abs()
    frame = frame.dropna(subset=["absolute_delta"]).sort_values(
        "absolute_delta",
        ascending=False,
    )
    selected = frame.drop_duplicates(
        subset=["pathway", "source", "target"],
        keep="first",
    ).head(16)
    graph = nx.DiGraph()
    for row in selected.itertuples(index=False):
        source = str(row.source)
        target = str(row.target)
        pathway = str(row.pathway)
        delta = float(row.delta_HFD_NCD)
        pathway_node = f"P:{pathway}"
        graph.add_node(source, node_kind="cell")
        graph.add_node(target, node_kind="cell")
        graph.add_node(pathway_node, node_kind="pathway")
        width = 0.6 + 3.0 * min(abs(delta), 1.0)
        color = "#c05b4d" if delta >= 0 else "#4f7fa8"
        graph.add_edge(source, pathway_node, width=width, color=color)
        graph.add_edge(pathway_node, target, width=width, color=color)
    if graph.number_of_nodes() == 0:
        return
    position = nx.spring_layout(graph, seed=42, k=1.5 / (graph.number_of_nodes() ** 0.5))
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    cells = [
        node
        for node, data in graph.nodes(data=True)
        if data.get("node_kind") == "cell"
    ]
    pathway_nodes = [node for node in graph.nodes if node.startswith("P:")]
    cell_colors = {
        node: plt.get_cmap("tab10")(index % 10)
        for index, node in enumerate(sorted(cells))
    }
    nx.draw_networkx_nodes(
        graph,
        position,
        nodelist=cells,
        node_size=[360 + 80 * graph.degree(node) for node in cells],
        node_color=[cell_colors[node] for node in cells],
        edgecolors="white",
        linewidths=1.0,
        ax=ax,
    )
    nx.draw_networkx_nodes(
        graph,
        position,
        nodelist=pathway_nodes,
        node_size=[220 + 45 * graph.degree(node) for node in pathway_nodes],
        node_shape="D",
        node_color="#e8c66a",
        edgecolors="#6c5a2c",
        linewidths=0.8,
        ax=ax,
    )
    for source, target, data in graph.edges(data=True):
        ax.annotate(
            "",
            xy=position[target],
            xytext=position[source],
            arrowprops={
                "arrowstyle": "-|>",
                "color": data["color"],
                "linewidth": data["width"],
                "alpha": 0.55,
                "connectionstyle": "arc3,rad=0.08",
            },
        )
    label_texts = nx.draw_networkx_labels(
        graph,
        position,
        labels={
            node: node.replace("P:", "")[:28]
            for node in graph.nodes
        },
        font_size=6.2,
        font_family="Arial",
        font_color="#1f2933",
        ax=ax,
    )
    for text in label_texts.values():
        text.set_path_effects(
            [path_effects.withStroke(linewidth=1.8, foreground="white")]
        )
    handles = [
        plt.Line2D([0], [0], color="#c05b4d", lw=3, label="HFD > NCD"),
        plt.Line2D([0], [0], color="#4f7fa8", lw=3, label="HFD < NCD"),
    ]
    ax.legend(handles=handles, loc="upper right", frameon=False)
    ax.set_title(
        "Ligand-receptor communication network (CellChat-like scoring)",
        fontweight="bold",
    )
    ax.axis("off")
    save_figure(fig, output)


def run_human_single_cell(
    extracted_dir: Path,
    soft_path: Path,
    core_genes: list[str],
    output_dir: Path,
    *,
    max_cells_per_sample: int = 1200,
    seed: int = 42,
) -> dict[str, Any]:
    """Process GSE202379 raw count CSVs into a disease-spectrum atlas."""
    ensure_dir(output_dir)
    data_path = output_dir / "human_liver_processed.h5ad"
    if data_path.exists():
        data = ad.read_h5ad(data_path)
    else:
        metadata = parse_geo_soft_samples(soft_path)
        csv_files = sorted(extracted_dir.rglob("*raw_counts.csv.gz"))
        if not csv_files:
            csv_files = sorted(extracted_dir.rglob("*raw*counts.csv.gz"))
        if not csv_files:
            raise FileNotFoundError(f"no raw count CSVs under {extracted_dir}")
        objects: list[ad.AnnData] = []
        rng = np.random.default_rng(seed)
        for index, csv_path in enumerate(csv_files, start=1):
            match = re.match(r"(GSM\d+)", csv_path.name)
            gsm = match.group(1) if match else csv_path.name.split("_", 1)[0]
            record = metadata.get(gsm, {})
            frame = pd.read_csv(csv_path, index_col=0)
            if frame.empty:
                continue
            if frame.shape[1] > max_cells_per_sample:
                selected_columns = np.sort(
                    rng.choice(frame.shape[1], size=max_cells_per_sample, replace=False)
                )
                frame = frame.iloc[:, selected_columns]
            matrix = sparse.csr_matrix(frame.to_numpy(dtype=np.float32).T)
            obs = pd.DataFrame(
                index=[f"{gsm}_{cell}" for cell in frame.columns.astype(str)],
                data={
                    "sample_id": gsm,
                    "condition": _human_condition(record),
                    "patient_id": str(record.get("patient_id") or gsm),
                    "batch": gsm,
                },
            )
            var = pd.DataFrame(index=frame.index.astype(str))
            sample = ad.AnnData(X=matrix, obs=obs, var=var)
            sample.var_names_make_unique()
            sample.obs_names_make_unique()
            objects.append(sample)
            LOG.info("loaded human sample %s/%s: %s", index, len(csv_files), gsm)
        if not objects:
            raise RuntimeError("no human sample CSV files could be loaded")
        data = ad.concat(objects, join="outer", fill_value=0, index_unique=None)
        data.var_names_make_unique()
        data.layers["counts"] = data.X.copy()
        data.var["mt"] = data.var_names.astype(str).str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], inplace=True, percent_top=None)
        sc.pp.filter_cells(data, min_genes=200)
        sc.pp.filter_genes(data, min_cells=3)
        data = data[data.obs["pct_counts_mt"].astype(float) <= 20.0].copy()
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        data.raw = data
        sc.pp.highly_variable_genes(data, n_top_genes=min(2500, data.n_vars), batch_key="batch")
        hvg = data[:, data.var["highly_variable"]].copy()
        _integrate_batches(hvg)
        sc.pp.scale(hvg, max_value=10)
        sc.tl.pca(hvg, n_comps=min(40, hvg.n_vars - 1, hvg.n_obs - 1), svd_solver="arpack")
        sc.pp.neighbors(hvg, n_neighbors=15, n_pcs=min(30, hvg.obsm["X_pca"].shape[1]))
        sc.tl.umap(hvg, random_state=seed)
        sc.tl.leiden(hvg, resolution=0.8, key_added="cluster", flavor="igraph", n_iterations=2)
        data.obsm["X_pca"] = hvg.obsm["X_pca"]
        data.obsm["X_umap"] = hvg.obsm["X_umap"]
        data.obs["cluster"] = hvg.obs["cluster"].astype(str)
        data.obs["cell_type"] = annotate_clusters(data, HUMAN_MARKERS)
        data.write_h5ad(data_path, compression="gzip")
    data.obs["cell_type"] = (
        data.obs["cell_type"]
        .astype(str)
        .str.replace(r"\s+\d+$", "", regex=True)
    )

    _plot_umap(
        data,
        color="cell_type",
        output=output_dir / "fig4i_umap_human_cell_types.png",
        title="Human MASLD/MASH single-nucleus atlas",
        legend_title="Cell type",
    )
    core = _resolve_genes(data.var_names, core_genes)
    gene_frames = []
    for gene in core:
        values = _expression_vector(data, gene)
        gene_frames.append(
            pd.DataFrame(
                {
                    "cell": data.obs_names,
                    "gene": gene,
                    "expression": values,
                    "cell_type": data.obs["cell_type"].astype(str).to_numpy(),
                    "condition": data.obs["condition"].astype(str).to_numpy(),
                    "patient_id": data.obs["patient_id"].astype(str).to_numpy(),
                }
            )
        )
    expression_long = pd.concat(gene_frames, ignore_index=True) if gene_frames else pd.DataFrame()
    expression_long.to_csv(output_dir / "human_core_gene_expression_long.csv.gz", index=False, compression="gzip")
    human_stats = _human_core_statistics(expression_long)
    human_stats.to_csv(output_dir / "human_core_gene_cell_type_stats.csv", index=False)
    if not expression_long.empty:
        _plot_human_core_genes(
            expression_long,
            output_dir / "fig4j_human_core_gene_by_cell_type.png",
        )
    write_json(
        output_dir / "human_single_cell_summary.json",
        {
            "n_cells": int(data.n_obs),
            "n_genes": int(data.n_vars),
            "samples": int(data.obs["sample_id"].nunique()),
            "conditions": data.obs["condition"].value_counts().to_dict(),
            "cell_types": data.obs["cell_type"].value_counts().to_dict(),
            "core_genes": core,
            "sampling_note": (
                f"At most {max_cells_per_sample} cells were sampled per GEO sample "
                "to keep memory use reproducible on the local workstation."
            ),
        },
    )
    return {
        "h5ad": data_path,
        "core_genes": core,
        "expression": output_dir / "human_core_gene_expression_long.csv.gz",
    }


def _human_condition(record: dict[str, Any]) -> str:
    status = str(record.get("disease_status") or record.get("title") or "").lower()
    if "healthy" in status:
        return "Healthy"
    if "end stage" in status:
        return "End-stage"
    if "with cirrhosis" in status:
        return "MASH cirrhosis"
    if "nash" in status or "mash" in status:
        return "MASH"
    if "nafld" in status or "masld" in status:
        return "MASLD"
    return "Unknown"


def _plot_human_core_genes(expression: pd.DataFrame, output: Path) -> None:
    data = expression[
        ~expression["condition"].isin(["Unknown", "End-stage"])
    ].copy()
    if data.empty:
        return
    data["cell_type"] = data["cell_type"].astype(str).str.replace(
        r"\s+\d+$",
        "",
        regex=True,
    )
    genes = list(dict.fromkeys(data["gene"]))
    cell_types = (
        data.groupby("cell_type")["cell"]
        .nunique()
        .sort_values(ascending=False)
        .head(10)
        .index.tolist()
    )
    condition_order = ["Healthy", "MASLD", "MASH", "MASH cirrhosis"]
    disease_conditions = condition_order[1:]
    rows: list[dict[str, Any]] = []
    for gene in genes:
        for cell_type in cell_types:
            subset = data[(data["gene"] == gene) & (data["cell_type"] == cell_type)]
            values = {
                condition: pd.to_numeric(
                    subset.loc[subset["condition"] == condition, "expression"],
                    errors="coerce",
                ).dropna()
                for condition in condition_order
            }
            healthy = values["Healthy"]
            if healthy.empty:
                continue
            healthy_median = float(healthy.median())
            nonempty = [
                values_group
                for values_group in values.values()
                if len(values_group) >= 2
            ]
            p_value = (
                stats.kruskal(*nonempty).pvalue if len(nonempty) >= 2 else np.nan
            )
            for condition in disease_conditions:
                values_group = values[condition]
                rows.append(
                    {
                        "gene": gene,
                        "cell_type": cell_type,
                        "condition": condition,
                        "n_healthy": int(len(healthy)),
                        "n_condition": int(len(values_group)),
                        "healthy_median": healthy_median,
                        "condition_median": (
                            float(values_group.median())
                            if not values_group.empty
                            else np.nan
                        ),
                        "median_difference_vs_healthy": (
                            float(values_group.median()) - healthy_median
                            if not values_group.empty
                            else np.nan
                        ),
                        "kruskal_p_value": p_value,
                    }
                )
    statistics = pd.DataFrame(rows)
    if statistics.empty:
        return
    unique_tests = (
        statistics[["gene", "cell_type", "kruskal_p_value"]]
        .drop_duplicates(subset=["gene", "cell_type"])
        .copy()
    )
    unique_tests["kruskal_fdr"] = bh_fdr(unique_tests["kruskal_p_value"])
    statistics = statistics.drop(columns=["kruskal_fdr"], errors="ignore").merge(
        unique_tests[["gene", "cell_type", "kruskal_fdr"]],
        on=["gene", "cell_type"],
        how="left",
    )
    statistics.to_csv(output.with_suffix(".csv"), index=False)

    matrix = (
        statistics.pivot(
            index=["gene", "condition"],
            columns="cell_type",
            values="median_difference_vs_healthy",
        )
        .reindex(
            pd.MultiIndex.from_product(
                [genes, disease_conditions],
                names=["gene", "condition"],
            )
        )
        .reindex(columns=cell_types)
    )
    fdr = (
        statistics.pivot(
            index=["gene", "condition"],
            columns="cell_type",
            values="kruskal_fdr",
        )
        .reindex(matrix.index)
    )
    finite = np.abs(matrix.to_numpy(dtype=float))
    vmax = max(0.5, float(np.nanpercentile(finite, 95))) if np.isfinite(finite).any() else 1.0
    fig, ax = plt.subplots(figsize=(7.2, max(4.6, len(matrix) * 0.42)))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=5.5)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(
        [f"{gene} | {condition}" for gene, condition in matrix.index],
        fontsize=5.8,
    )
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = float(matrix.iloc[row, column])
            if not np.isfinite(value):
                continue
            fdr_value = float(fdr.iloc[row, column])
            marker = (
                "***"
                if np.isfinite(fdr_value) and fdr_value <= 0.001
                else "**"
                if np.isfinite(fdr_value) and fdr_value <= 0.01
                else "*"
                if np.isfinite(fdr_value) and fdr_value <= 0.05
                else ""
            )
            ax.text(
                column,
                row,
                marker,
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                color="white" if abs(value) > 0.55 * vmax else "#20262d",
            )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Core gene | disease stage")
    ax.set_title(
        "Core-gene expression shifts across human MASLD stages",
        fontweight="bold",
    )
    colorbar = fig.colorbar(image, ax=ax, shrink=0.72, pad=0.02)
    colorbar.set_label("Median expression difference vs Healthy")
    ax.text(
        1.0,
        -0.34,
        "*FDR<0.05  **FDR<0.01  ***FDR<0.001 (Kruskal-Wallis)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.2,
        color="#4a5560",
    )
    save_figure(fig, output)


def _mouse_core_statistics(expression: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if expression.empty:
        return pd.DataFrame()
    for (gene, cell_type), subset in expression.groupby(
        ["gene", "cell_type"], observed=True
    ):
        control = subset.loc[subset["condition"] == "NCD", "expression"]
        disease = subset.loc[subset["condition"] == "HFD", "expression"]
        if len(control) < 2 or len(disease) < 2:
            p_value = np.nan
        else:
            p_value = stats.mannwhitneyu(
                control,
                disease,
                alternative="two-sided",
            ).pvalue
        rows.append(
            {
                "gene": gene,
                "cell_type": cell_type,
                "n_NCD": int(len(control)),
                "n_HFD": int(len(disease)),
                "median_NCD": float(control.median()) if len(control) else np.nan,
                "median_HFD": float(disease.median()) if len(disease) else np.nan,
                "median_delta_HFD_NCD": (
                    float(disease.median() - control.median())
                    if len(control) and len(disease)
                    else np.nan
                ),
                "p_value": p_value,
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["fdr"] = bh_fdr(frame["p_value"])
    return frame.sort_values(["fdr", "gene", "cell_type"], na_position="last")


def _human_core_statistics(expression: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if expression.empty:
        return pd.DataFrame()
    frame = expression[
        ~expression["condition"].isin(["Unknown", "End-stage"])
    ].copy()
    frame["major_cell_type"] = frame["cell_type"].astype(str).str.replace(
        r"\s+\d+$",
        "",
        regex=True,
    )
    for (gene, cell_type), subset in frame.groupby(
        ["gene", "major_cell_type"],
        observed=True,
    ):
        groups = [
            group["expression"].to_numpy()
            for _, group in subset.groupby("condition", observed=True)
            if len(group) >= 2
        ]
        if len(groups) >= 2:
            p_value = stats.kruskal(*groups).pvalue
        else:
            p_value = np.nan
        rows.append(
            {
                "gene": gene,
                "cell_type": cell_type,
                "n_cells": int(len(subset)),
                "n_conditions": int(subset["condition"].nunique()),
                "p_value": p_value,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_fdr(result["p_value"])
    return result.sort_values(["fdr", "gene", "cell_type"], na_position="last")
