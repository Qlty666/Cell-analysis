"""GO/KEGG enrichment wrappers and publication-oriented summary plots."""

from __future__ import annotations

import shutil
import subprocess
import textwrap
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import cache
from .common import ensure_dir, read_json, save_figure, write_json


def run_go_kegg(
    genes: list[str],
    output_dir: Path,
    *,
    species: str = "hs",
    force: bool = False,
    timeout: int = 1200,
    universe: list[str] | None = None,
) -> dict[str, Any]:
    ensure_dir(output_dir)
    genes = sorted(set(str(gene).upper() for gene in genes if str(gene).strip()))
    script = Path(__file__).resolve().parents[1] / "docking" / "insilico_enrichment.R"
    universe = sorted(set(str(g).strip() for g in universe or [] if str(g).strip()))
    key = cache.signature(files=[script, Path(__file__)],
                          parameters={"genes": genes, "universe": universe, "species": species, "fdr": 0.05})
    state = output_dir / "enrichment.cache.json"
    summary_path = output_dir / "enrichment_summary.json"
    outputs = [
        summary_path,
        output_dir / "insilico_go_enrichment.csv",
        output_dir / "insilico_kegg_enrichment.csv",
        output_dir / "query_id_mapping.csv",
        output_dir / "background_id_mapping.csv",
        output_dir / "enrichment_database_versions.csv",
    ]
    cached = read_json(summary_path, None)
    if (
        not force
        and isinstance(cached, dict)
        and cached.get("status") == "completed"
        and cache.valid(state, key, outputs)
    ):
        return cached
    # Do not leave stale significant figures behind when a rerun is negative/failed.
    for stem in ("fig1h_go_enrichment_bubble", "fig1g_kegg_enrichment_bar"):
        for suffix in (".png", ".pdf", ".svg"):
            (output_dir / (stem + suffix)).unlink(missing_ok=True)
    if not universe:
        status = {"status": "blocked", "reason": "explicit testable background genes are required"}
        write_json(summary_path, status)
        return status
    if len(genes) < 3:
        status = {"status": "skipped", "reason": "fewer than three genes"}
        write_json(output_dir / "enrichment_summary.json", status)
        return status
    input_path = output_dir / "enrichment_input_genes.csv"
    pd.DataFrame({"gene": genes}).to_csv(input_path, index=False)
    background_path = output_dir / "enrichment_background_genes.csv"
    pd.DataFrame({"gene": universe}).to_csv(background_path, index=False)
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if not rscript:
        status = {"status": "failed", "reason": "Rscript unavailable"}
        write_json(output_dir / "enrichment_summary.json", status)
        return status
    script = Path(__file__).resolve().parents[1] / "docking" / "insilico_enrichment.R"
    proc = subprocess.run(
        [rscript, str(script), str(input_path), str(output_dir), species, str(background_path)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )
    if proc.returncode != 0:
        status = {
            "status": "failed",
            "reason": (proc.stderr or proc.stdout)[-3000:],
        }
        write_json(output_dir / "enrichment_summary.json", status)
        return status
    go_path = output_dir / "insilico_go_enrichment.csv"
    kegg_path = output_dir / "insilico_kegg_enrichment.csv"
    go = pd.read_csv(go_path) if go_path.exists() else pd.DataFrame()
    kegg = pd.read_csv(kegg_path) if kegg_path.exists() else pd.DataFrame()
    go_count = len(_significant(go))
    kegg_count = len(_significant(kegg))
    _plot_enrichment(go, output_dir / "fig1h_go_enrichment_bubble.png", "GO enrichment")
    _plot_enrichment(kegg, output_dir / "fig1g_kegg_enrichment_bar.png", "KEGG enrichment")
    status = {
        "status": "completed",
        "genes": len(genes),
        "background_genes": len(universe),
        "scientific_outcome": "significant_terms" if go_count or kegg_count else "valid_negative",
        "fdr_method": "Benjamini-Hochberg",
        "primary_fdr_cutoff": 0.05,
        "go_terms": go_count,
        "kegg_terms": kegg_count,
        "go_csv": str(go_path),
        "kegg_csv": str(kegg_path),
    }
    write_json(output_dir / "enrichment_summary.json", status)
    cache.save(state, key, outputs)
    return status


def _significant(frame: pd.DataFrame) -> pd.DataFrame:
    if "p.adjust" not in frame:
        return frame.iloc[:0]
    return frame.loc[pd.to_numeric(frame["p.adjust"], errors="coerce") < 0.05].copy()


def _plot_enrichment(frame: pd.DataFrame, output: Path, title: str) -> None:
    if frame.empty:
        return
    required = {"Description", "p.adjust"}
    if not required.issubset(frame.columns):
        return
    values = _significant(frame)
    values["p.adjust"] = pd.to_numeric(values["p.adjust"], errors="coerce")
    values = values.dropna(subset=["p.adjust"])
    if values.empty:
        return
    ontology_values = (
        [value for value in ("BP", "CC", "MF") if value in set(values["ONTOLOGY"])]
        if "ONTOLOGY" in values
        else []
    )
    if ontology_values:
        fig, axes = plt.subplots(
            1,
            len(ontology_values),
            figsize=(7.2, 4.8),
            squeeze=False,
        )
        scatter = None
        for ax, ontology in zip(axes.flat, ontology_values):
            subset = values[values["ONTOLOGY"] == ontology]
            subset = subset.sort_values("p.adjust").head(10)
            scatter = _plot_enrichment_axis(
                ax,
                subset,
                f"GO {ontology}",
                include_legend=ontology == ontology_values[0],
            )
        fig.suptitle(title, fontweight="bold")
        if scatter is not None:
            fig.colorbar(scatter, ax=axes.ravel().tolist(), label="-log10(FDR)", shrink=0.7)
        save_figure(fig, output)
        return
    values = values.sort_values("p.adjust").head(10)
    fig, ax = plt.subplots(figsize=(7.2, max(3.8, len(values) * 0.38)))
    scatter = _plot_enrichment_axis(ax, values, title)
    if scatter is not None:
        fig.colorbar(scatter, ax=ax, label="-log10(FDR)", shrink=0.7)
    save_figure(fig, output)


def _plot_enrichment_axis(
    ax,
    values: pd.DataFrame,
    title: str,
    *,
    include_legend: bool = True,
):
    if values.empty:
        ax.axis("off")
        ax.set_title(title, fontweight="bold")
        return None
    values = values.copy()
    if "Count" in values:
        count = pd.to_numeric(values["Count"], errors="coerce").fillna(1)
    else:
        count = pd.Series(np.ones(len(values)), index=values.index)
    size = 30 + 14 * count.to_numpy()
    negative_log10 = -np.log10(values["p.adjust"].clip(lower=1e-300))
    scatter = ax.scatter(
        negative_log10,
        np.arange(len(values)),
        s=size,
        c=negative_log10,
        cmap="cividis",
        edgecolors="#25313d",
        linewidths=0.35,
    )
    ax.set_yticks(np.arange(len(values)))
    ax.set_yticklabels(
        [
            textwrap.fill(str(value), width=48)
            for value in values["Description"].astype(str)
        ],
        fontsize=6.8,
    )
    ax.invert_yaxis()
    ax.set_xlabel("-log10(adjusted P)")
    ax.grid(axis="x", color="#dfe5ea", linewidth=0.6, alpha=0.8)
    unique_counts = sorted({int(value) for value in count if pd.notna(value)})
    legend_counts = unique_counts[:4]
    if len(unique_counts) > 4:
        legend_counts.append(unique_counts[-1])
    handles = [
        ax.scatter(
            [],
            [],
            s=30 + 14 * value,
            facecolors="none",
            edgecolors="#53606c",
            linewidths=0.8,
            label=str(value),
        )
        for value in legend_counts
    ]
    if handles and include_legend:
        ax.legend(
            handles=handles,
            title="Gene count",
            loc="lower right",
            frameon=False,
            fontsize=6.5,
            title_fontsize=7,
            labelspacing=0.8,
            borderpad=0.2,
        )
    ax.set_title(title, fontweight="bold")
    return scatter
