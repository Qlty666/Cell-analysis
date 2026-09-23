#!/usr/bin/env python3
"""Re-render experiment-plan-one figures from existing source tables.

The optimizer never changes raw result tables. It copies unavailable panels
unchanged and re-renders panels whose source data are available with the
publication-oriented plotting functions in ``src/experiment_plan_one``.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from experiment_plan_one.bulk import candidate_heatmap, validation_boxplots  # noqa: E402
from experiment_plan_one.classify import PANEL_ALIASES  # noqa: E402
from experiment_plan_one.common import (  # noqa: E402
    ensure_dir,
    read_json,
    save_figure,
    write_json,
)
from experiment_plan_one.docking_md import (  # noqa: E402
    _plot_docking_heatmap,
    _write_docking_pose_figure,
    _write_interaction_figure,
)
from experiment_plan_one.enrichment import _plot_enrichment  # noqa: E402
from experiment_plan_one.ml import (  # noqa: E402
    _plot_auc_heatmap,
    _plot_calibration,
    _plot_roc,
    _shap_bar,
)
from experiment_plan_one.pipeline import (  # noqa: E402
    _export_knockout_outputs,
    _plot_not_run_panel,
    _plot_source_coverage,
    _plot_target_source_counts,
)
from experiment_plan_one.ppi import (  # noqa: E402
    _bar_rank,
    _draw_network,
    ppi_hub_metrics,
)
from experiment_plan_one.single_cell import (  # noqa: E402
    _plot_cell_communication,
    _plot_composition,
    _plot_core_violin,
)
from experiment_plan_one.targets import (  # noqa: E402
    make_venn_figure,
    make_workflow_figure,
    write_compound_figures,
)


def _copy_existing_panels(source_root: Path, output_root: Path) -> None:
    for alias in PANEL_ALIASES:
        source = source_root / alias.source
        if not source.exists():
            continue
        for suffix in (".png", ".pdf", ".svg"):
            candidate = source.with_suffix(suffix)
            if candidate.exists():
                destination = output_root / alias.source
                destination = destination.with_suffix(suffix)
                ensure_dir(destination.parent)
                shutil.copy2(candidate, destination)
                if suffix == ".png":
                    _flatten_png(destination)


def _flatten_png(path: Path) -> None:
    with Image.open(path) as image:
        if image.mode != "RGBA":
            return
        alpha = image.getchannel("A")
        flattened = Image.new("RGB", image.size, "white")
        flattened.paste(image.convert("RGB"), mask=alpha)
        flattened.save(
            path,
            format="PNG",
            dpi=image.info.get("dpi", (600, 600)),
            optimize=True,
        )


def _render_compound_and_enrichment(source_root: Path, output_root: Path) -> None:
    source_dir = source_root / "01_compound_characterization"
    summary = read_json(source_dir / "target_collection_summary.json", {})
    properties = summary.get("compound") or {}
    descriptors = summary.get("compound_descriptors") or {}
    smiles = str(properties.get("canonical_smiles") or "")
    target_dir = ensure_dir(output_root / "01_compound_characterization")
    if smiles and descriptors:
        write_compound_figures(smiles, properties, descriptors, target_dir)
    make_workflow_figure(target_dir / "fig1a_workflow.png")

    compound_sources: dict[str, set[str]] = {}
    for path in sorted((source_dir / "sources").glob("*.csv")):
        frame = pd.read_csv(path)
        if "gene" not in frame.columns or frame.empty:
            continue
        compound_sources[path.stem] = set(
            frame["gene"].dropna().astype(str).str.upper()
        )
    _plot_target_source_counts(
        compound_sources,
        target_dir / "fig1d_compound_target_source_counts.png",
    )

    disease_dir = source_root / "02_disease_targets"
    disease_summary = read_json(disease_dir / "disease_target_summary.json", {})
    _plot_source_coverage(
        disease_summary.get("sources") or {},
        output_root / "02_disease_targets" / "fig1e_disease_target_source_status.png",
    )

    compound_path = source_dir / "compound_targets.csv"
    disease_path = disease_dir / "disease_targets.csv"
    if compound_path.exists() and disease_path.exists():
        compound = set(
            pd.read_csv(compound_path)["gene"].dropna().astype(str).str.upper()
        )
        disease = set(
            pd.read_csv(disease_path)["gene"].dropna().astype(str).str.upper()
        )
        make_venn_figure(
            {"Compound targets": compound, "Disease targets": disease},
            output_root
            / "03_intersection_ppi"
            / "fig1f_compound_disease_venn.png",
            title="6PPD-Q targets and NAFLD targets",
        )

    enrichment_dir = source_root / "03_intersection_ppi" / "enrichment"
    summary = read_json(enrichment_dir / "enrichment_summary.json", {})
    kegg_value = str(summary.get("kegg_csv") or "").strip()
    go_value = str(summary.get("go_csv") or "").strip()
    kegg_path = Path(kegg_value) if kegg_value else None
    go_path = Path(go_value) if go_value else None
    output_enrichment = ensure_dir(
        output_root / "03_intersection_ppi" / "enrichment"
    )
    if go_path is not None and go_path.exists():
        _plot_enrichment(
            pd.read_csv(go_path),
            output_enrichment / "fig1h_go_enrichment_bubble.png",
            "GO enrichment",
        )
    if kegg_path is not None and kegg_path.exists():
        _plot_enrichment(
            pd.read_csv(kegg_path),
            output_enrichment / "fig1g_kegg_enrichment_bar.png",
            "KEGG enrichment",
        )


def _render_ppi(source_root: Path, output_root: Path) -> None:
    source_dir = source_root / "03_intersection_ppi"
    output_dir = ensure_dir(output_root / "03_intersection_ppi")
    edges_path = source_dir / "string_edges.tsv"
    if not edges_path.exists():
        return
    edges = pd.read_csv(source_dir / "string_edges.tsv", sep="\t")
    edges.to_csv(output_dir / "string_edges.tsv", sep="\t", index=False)
    overlap_path = source_dir / "compound_disease_overlap.csv"
    genes = (
        pd.read_csv(overlap_path)["gene"].astype(str).tolist()
        if overlap_path.exists()
        else edges["node1"].astype(str).tolist()
    )
    metrics = ppi_hub_metrics(edges, genes=genes)
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        graph.add_edge(
            str(row.node1).upper(),
            str(row.node2).upper(),
            weight=float(row.score) / 1000.0,
        )
    summary = read_json(source_dir / "ppi_summary.json", {})
    modules = summary.get("modules") or {}
    node_modules = {
        str(gene): str(module)
        for module, genes_in_module in modules.items()
        for gene in genes_in_module
    }
    focus = set(metrics["gene"].astype(str))
    _draw_network(
        graph,
        metrics,
        output_dir / "fig2a_string_network.png",
        title=f"STRING PPI network | score >= {summary.get('required_score', 700) / 1000:.2f}",
        node_color="#3f7f93",
        focus_genes=focus,
    )
    _draw_network(
        graph,
        metrics,
        output_dir / "fig2b_cytoscape_module_network.png",
        title="PPI modules (Louvain communities)",
        node_color=None,
        focus_genes=focus,
        node_modules=node_modules,
    )
    _bar_rank(
        metrics,
        "degree",
        output_dir / "fig2c_degree_top20.png",
        "Degree",
        20,
    )
    _bar_rank(
        metrics,
        "betweenness",
        output_dir / "fig2d_betweenness_top20.png",
        "Betweenness centrality",
        20,
    )
    make_venn_figure(
        {
            "MCC Top 10": set(metrics.nlargest(10, "mcc")["gene"]),
            "Degree Top 10": set(metrics.nlargest(10, "degree")["gene"]),
        },
        output_dir / "fig2e_mcc_degree_venn.png",
        title="Consensus hub genes",
    )
    metrics.to_csv(output_dir / "ppi_hub_metrics.csv", index=False)


def _render_bulk(source_root: Path, output_root: Path) -> None:
    source_dir = source_root / "04_bulk_training"
    output_dir = ensure_dir(output_root / "04_bulk_training")
    processed = source_root / "00_data" / "processed"
    heatmap_path = source_dir / "fig2f_candidate_gene_heatmap.csv"
    expression_path = processed / "GSE89632" / "expression.csv.gz"
    metadata_path = processed / "GSE89632" / "metadata.csv"
    if heatmap_path.exists() and expression_path.exists() and metadata_path.exists():
        genes = pd.read_csv(heatmap_path, index_col=0).index.astype(str).tolist()
        candidate_heatmap(
            expression_path,
            metadata_path,
            genes,
            output_dir / "fig2f_candidate_gene_heatmap.png",
            group_order=["HC", "SS", "NASH"],
        )

    for prefix, accession, comparison in (
        (
            "fig2g_gse49541",
            "GSE49541",
            ("mild_fibrosis", "advanced_fibrosis"),
        ),
        (
            "fig2h_gse164441",
            "GSE164441",
            ("adjacent_normal", "tumor"),
        ),
    ):
        stats_path = source_dir / f"{prefix}_candidate_boxplot_stats.csv"
        expression_path = processed / accession / "expression.csv.gz"
        metadata_path = processed / accession / "metadata.csv"
        if not stats_path.exists() or not expression_path.exists():
            continue
        genes = pd.read_csv(stats_path)["gene"].astype(str).tolist()
        validation_boxplots(
            expression_path,
            metadata_path,
            genes,
            output_dir,
            comparisons=comparison,
            prefix=prefix,
        )


def _render_ml(source_root: Path, output_root: Path) -> None:
    source_dir = source_root / "05_machine_learning"
    output_dir = ensure_dir(output_root / "05_machine_learning")
    performance_path = source_dir / "model_performance.csv"
    if performance_path.exists():
        _plot_auc_heatmap(
            pd.read_csv(performance_path),
            output_dir / "fig3a_multimodel_auc_heatmap.png",
        )
    roc_specs = (
        (
            "best_model_cv_predictions.csv",
            "fig3b_best_model_roc_training.png",
            "Best model: out-of-fold training ROC",
        ),
        (
            "GSE49541_fibrosis_predictions.csv",
            "fig3c_gse49541_external_roc.png",
            "GSE49541: advanced vs mild fibrosis",
        ),
        (
            "GSE164441_tumor_predictions.csv",
            "fig3d_gse164441_external_roc.png",
            "GSE164441: tumor vs adjacent normal",
        ),
        (
            "GSE135251_NAFLD_predictions.csv",
            "fig3d_supplementary_gse135251_roc.png",
            "GSE135251: NAFLD vs control",
        ),
    )
    for filename, output_name, title in roc_specs:
        path = source_dir / filename
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        if not {"observed", "probability"}.issubset(frame.columns):
            continue
        _plot_roc(
            frame["observed"].to_numpy(dtype=int),
            frame["probability"].to_numpy(dtype=float),
            output_dir / output_name,
            title=title,
        )
    calibration_path = source_dir / "best_model_cv_predictions.csv"
    if calibration_path.exists():
        frame = pd.read_csv(calibration_path)
        _plot_calibration(
            frame["observed"].to_numpy(dtype=int),
            frame["probability"].to_numpy(dtype=float),
            output_dir / "fig3e_calibration_curve.png",
        )
    importance_path = source_dir / "fig3f_shap_importance.csv"
    if importance_path.exists():
        _shap_bar(
            pd.read_csv(importance_path),
            output_dir / "fig3f_shap_importance.png",
        )


def _render_single_cell(source_root: Path, output_root: Path) -> None:
    mouse_source = source_root / "06_single_cell_mouse"
    mouse_output = ensure_dir(output_root / "06_single_cell_mouse")
    expression_path = mouse_source / "core_gene_expression_long.csv.gz"
    if expression_path.exists():
        _plot_core_violin(
            pd.read_csv(expression_path),
            mouse_output / "fig4c_core_gene_violin.png",
        )
    composition_path = mouse_source / "cell_composition.csv"
    if composition_path.exists():
        _plot_composition(
            pd.read_csv(composition_path),
            mouse_output / "fig4e_cell_composition.png",
        )
    interaction_path = mouse_source / "cellchat_interactions.csv"
    pathway_path = mouse_source / "cellchat_pathways.csv"
    if not interaction_path.exists():
        interaction_path = mouse_source / "cellchat_like_interactions.csv"
    if not pathway_path.exists():
        pathway_path = mouse_source / "cellchat_like_pathways.csv"
    if interaction_path.exists() and pathway_path.exists():
        _plot_cell_communication(
            pd.read_csv(interaction_path),
            pd.read_csv(pathway_path),
            mouse_output / "fig4f_cellchat_network.png",
        )
    knockout_root = mouse_source / "virtual_knockout"
    if knockout_root.exists():
        _export_knockout_outputs(knockout_root, mouse_output)

    human_source = source_root / "07_single_cell_human"
    human_output = ensure_dir(output_root / "07_single_cell_human")
    stats_path = human_source / "fig4j_human_core_gene_by_cell_type.csv"
    if not stats_path.exists():
        return
    stats = pd.read_csv(stats_path)
    genes = list(dict.fromkeys(stats["gene"].astype(str)))
    conditions = ["MASLD", "MASH", "MASH cirrhosis"]
    cell_types = (
        stats.groupby("cell_type")["n_condition"]
        .sum()
        .sort_values(ascending=False)
        .head(10)
        .index.tolist()
    )
    matrix = (
        stats.pivot(
            index=["gene", "condition"],
            columns="cell_type",
            values="median_difference_vs_healthy",
        )
        .reindex(
            pd.MultiIndex.from_product(
                [genes, conditions],
                names=["gene", "condition"],
            )
        )
        .reindex(columns=cell_types)
    )
    fdr = (
        stats.pivot(
            index=["gene", "condition"],
            columns="cell_type",
            values="kruskal_fdr",
        )
        .reindex(matrix.index)
    )
    finite = np.abs(matrix.to_numpy(dtype=float))
    vmax = max(0.5, float(np.nanpercentile(finite, 95))) if np.isfinite(finite).any() else 1.0
    fig, ax = plt.subplots(figsize=(7.2, max(4.8, len(matrix) * 0.48)))
    image = ax.imshow(
        matrix.to_numpy(dtype=float),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=6.8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(
        [f"{gene} | {condition}" for gene, condition in matrix.index],
        fontsize=6.8,
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
                fontsize=7.2,
                fontweight="bold",
                color="white" if abs(value) > 0.55 * vmax else "#20262d",
            )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Core gene | disease stage")
    ax.set_title(
        "Core-gene expression shifts across human MASLD stages",
        fontweight="bold",
    )
    colorbar = fig.colorbar(
        image,
        ax=ax,
        location="bottom",
        shrink=0.82,
        pad=0.14,
    )
    colorbar.set_label("Median expression difference vs Healthy")
    ax.text(
        1.0,
        -0.52,
        "*FDR<0.05  **FDR<0.01  ***FDR<0.001 (Kruskal-Wallis)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
        color="#4a5560",
    )
    save_figure(fig, human_output / "fig4j_human_core_gene_by_cell_type.png")


def _render_docking(source_root: Path, output_root: Path) -> None:
    source_dir = source_root / "08_docking"
    output_dir = ensure_dir(output_root / "08_docking")
    scores_path = source_dir / "docking_scores.csv"
    if not scores_path.exists():
        return
    scores = pd.read_csv(scores_path)
    if scores.empty:
        return
    _plot_docking_heatmap(scores, output_dir / "fig5c_docking_affinity_heatmap.png")
    best = scores.sort_values("best_affinity_kcal_mol").iloc[0]
    target_dir = Path(str(best["workdir"]))
    if not target_dir.is_absolute():
        target_dir = source_root / target_dir
    _write_docking_pose_figure(
        target_dir,
        output_dir / "fig5a_docking_pose_3d.png",
        str(best["gene"]),
    )
    _write_interaction_figure(
        target_dir,
        output_dir / "fig5b_interaction_schematic.png",
        str(best["gene"]),
    )


def _render_md_placeholders(output_root: Path) -> None:
    output_dir = ensure_dir(output_root / "09_md_mmpbsa")
    for panel, filename in (
        ("d_rmsd", "fig5d_rmsd.png"),
        ("e_ligand_rmsd", "fig5e_ligand_rmsd.png"),
        ("f_rmsf", "fig5f_rmsf.png"),
        ("g_rg", "fig5g_rg.png"),
        ("h_mmpbsa", "fig5h_mmpbsa.png"),
    ):
        _plot_not_run_panel(panel, output_dir / filename)


def optimize(source_root: Path, output_root: Path) -> dict[str, Any]:
    source_root = source_root.resolve()
    output_root = output_root.resolve()
    if not source_root.exists():
        raise FileNotFoundError(f"result root not found: {source_root}")
    ensure_dir(output_root)
    _copy_existing_panels(source_root, output_root)
    _render_compound_and_enrichment(source_root, output_root)
    _render_ppi(source_root, output_root)
    _render_bulk(source_root, output_root)
    _render_ml(source_root, output_root)
    _render_single_cell(source_root, output_root)
    _render_docking(source_root, output_root)
    _render_md_placeholders(output_root)
    for path in output_root.rglob("*.png"):
        _flatten_png(path)
    expected_paths = {
        Path(alias.source).with_suffix(suffix).as_posix()
        for alias in PANEL_ALIASES
        for suffix in (".png", ".pdf", ".svg")
    }
    expected_paths.update(
        {
            f"09_md_mmpbsa/{stem}{suffix}"
            for stem in (
                "fig5d_rmsd",
                "fig5e_ligand_rmsd",
                "fig5f_rmsf",
                "fig5g_rg",
                "fig5h_mmpbsa",
            )
            for suffix in (".png", ".pdf", ".svg")
        }
    )
    files = sorted(
        relative
        for relative in expected_paths
        if (output_root / relative).exists()
    )
    manifest = {
        "source_root": str(source_root),
        "output_root": str(output_root),
        "files": files,
        "note": (
            "Existing panels are copied unchanged when they cannot be "
            "reconstructed from source tables; copied panels retain their "
            "original resolution and are not upsampled. The optimizer never "
            "changes raw result data."
        ),
    }
    write_json(output_root / "figure_optimization_manifest.json", manifest)
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Re-render experiment-plan-one figures from existing results."
    )
    parser.add_argument("--result-root", required=True)
    parser.add_argument(
        "--output-dir",
        default=None,
        help="default: <result-root>/10_reports/optimized_figures",
    )
    args = parser.parse_args(argv)
    source = Path(args.result_root).expanduser()
    output = (
        Path(args.output_dir).expanduser()
        if args.output_dir
        else source / "10_reports" / "optimized_figures"
    )
    manifest = optimize(source, output)
    print(f"optimized figures: {manifest['output_root']}")
    print(f"files: {len(manifest['files'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
