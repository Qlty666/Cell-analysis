"""Bulk transcriptome preparation, differential expression and visualization."""

from __future__ import annotations

import logging
import math
import re
import shutil
import subprocess
import gzip
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from .common import (
    LOG,
    benjamini_hochberg,
    collapse_probes,
    download_file,
    ensure_dir,
    extract_tar,
    parse_series_matrix,
    read_expression,
    save_expression,
    save_figure,
    split_gene_symbol,
    standardize_expression,
    write_json,
    zscore_rows,
)

GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
SERIES_MATRIX_URLS = {
    "GSE89632": f"{GEO}/series/GSE89nnn/GSE89632/matrix/GSE89632_series_matrix.txt.gz",
    "GSE49541": f"{GEO}/series/GSE49nnn/GSE49541/matrix/GSE49541_series_matrix.txt.gz",
}
GSE164441_COUNT_URL = (
    f"{GEO}/series/GSE164nnn/GSE164441/suppl/"
    "GSE164441_RNAseq_10Tvs10NT_cuffdiff_count.txt.gz"
)


def prepare_bulk_data(
    raw_dir: Path,
    processed_dir: Path,
    *,
    force: bool = False,
) -> dict[str, Path]:
    """Download, annotate and persist the three bulk datasets."""
    ensure_dir(raw_dir)
    ensure_dir(processed_dir)
    outputs: dict[str, Path] = {}
    platform_cache: dict[str, dict[str, str]] = {}

    for accession, url in SERIES_MATRIX_URLS.items():
        matrix_path = raw_dir / f"{accession}_series_matrix.txt.gz"
        download_file(url, matrix_path, timeout=180)
        dataset_dir = ensure_dir(processed_dir / accession)
        expression_path = dataset_dir / "expression.csv.gz"
        metadata_path = dataset_dir / "metadata.csv"
        if expression_path.exists() and metadata_path.exists() and not force:
            outputs[f"{accession}_expression"] = expression_path
            outputs[f"{accession}_metadata"] = metadata_path
            continue

        matrix, metadata, _ = parse_series_matrix(matrix_path)
        if accession == "GSE89632":
            mapping = _platform_mapping(
                matrix.index.astype(str).tolist(),
                processed_dir,
                "GPL14951",
                platform_cache,
            )
            mapped = collapse_probes(matrix, mapping, log=LOG)
            metadata["condition"] = [
                _normalize_gse89632_condition(row)
                for _, row in metadata.iterrows()
            ]
            metadata["fibrosis_stage"] = pd.to_numeric(
                metadata.get("fibrosis_stage", pd.Series(index=metadata.index, dtype=float)),
                errors="coerce",
            )
        elif accession == "GSE49541":
            mapping = _platform_mapping(
                matrix.index.astype(str).tolist(),
                processed_dir,
                "GPL570",
                platform_cache,
            )
            mapped = collapse_probes(matrix, mapping, log=LOG)
            metadata["condition"] = [
                _normalize_gse49541_condition(row)
                for _, row in metadata.iterrows()
            ]
            metadata["fibrosis_stage"] = metadata["condition"].map(
                {"mild_fibrosis": 1, "advanced_fibrosis": 3}
            )
        else:
            raise AssertionError(accession)
        mapped = mapped.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        save_expression(mapped, expression_path)
        metadata.to_csv(metadata_path, encoding="utf-8")
        outputs[f"{accession}_expression"] = expression_path
        outputs[f"{accession}_metadata"] = metadata_path

    count_path = raw_dir / "GSE164441_count.txt.gz"
    download_file(GSE164441_COUNT_URL, count_path, timeout=180)
    dataset_dir = ensure_dir(processed_dir / "GSE164441")
    expression_path = dataset_dir / "expression.csv.gz"
    metadata_path = dataset_dir / "metadata.csv"
    if force or not expression_path.exists() or not metadata_path.exists():
        frame = pd.read_csv(count_path, sep="\t")
        gene_column = _first_matching(frame.columns, ["tracking_id", "gene_id", "gene"])
        if gene_column is None:
            raise ValueError(f"GSE164441 count table has no gene ID column: {list(frame.columns)}")
        frame = frame.set_index(gene_column)
        sample_columns = [
            column for column in frame.columns if str(column).lower().endswith("_count")
        ]
        frame = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
        symbol_map = map_ensembl_symbols(
            frame.index.astype(str).tolist(),
            dataset_dir / "ensembl_to_symbol.csv",
        )
        mapped_rows: list[tuple[str, float, pd.Series]] = []
        for gene_id, row in frame.iterrows():
            clean_id = str(gene_id).split(".", 1)[0].upper()
            symbol = symbol_map.get(clean_id) or symbol_map.get(str(gene_id)) or ""
            if symbol:
                mapped_rows.append((symbol, float(row.var(ddof=0) or 0.0), row))
        if not mapped_rows:
            raise ValueError("GSE164441 could not map any Ensembl IDs to symbols")
        selected: dict[str, tuple[float, pd.Series]] = {}
        for symbol, variance, row in mapped_rows:
            if symbol not in selected or variance > selected[symbol][0]:
                selected[symbol] = (variance, row)
        mapped = pd.DataFrame(
            {symbol: value for symbol, (_, value) in selected.items()}
        ).T
        mapped.index.name = "gene"
        mapped = standardize_expression(mapped, already_log=False)
        metadata = pd.DataFrame(
            {
                "sample_id": sample_columns,
                "condition": [
                    "tumor" if str(sample).upper().endswith("T_COUNT") else "adjacent_normal"
                    for sample in sample_columns
                ],
                "patient": [
                    str(sample).rsplit("_", 2)[0] for sample in sample_columns
                ],
            }
        ).set_index("sample_id")
        save_expression(mapped, expression_path)
        metadata.to_csv(metadata_path, encoding="utf-8")
    outputs["GSE164441_expression"] = expression_path
    outputs["GSE164441_metadata"] = metadata_path

    gse135_archive = raw_dir / "GSE135251_RAW.tar"
    if gse135_archive.exists():
        extracted = raw_dir / "extracted" / "GSE135251"
        extract_tar(gse135_archive, extracted)
        dataset_dir = ensure_dir(processed_dir / "GSE135251")
        expression_path = dataset_dir / "expression.csv.gz"
        metadata_path = dataset_dir / "metadata.csv"
        if force or not expression_path.exists() or not metadata_path.exists():
            count_files = sorted(extracted.glob("*.counts.txt.gz"))
            samples: dict[str, pd.Series] = {}
            for count_path in count_files:
                gsm = count_path.name.split("_", 1)[0]
                frame = pd.read_csv(
                    count_path,
                    sep="\t",
                    header=None,
                    names=["gene_id", gsm],
                )
                series = frame.set_index("gene_id")[gsm].astype(float)
                samples[gsm] = series
            if not samples:
                raise FileNotFoundError(f"no GSE135251 count files under {extracted}")
            frame = pd.DataFrame(samples)
            frame.index = frame.index.astype(str)
            frame = frame.groupby(frame.index.str.split(".").str[0], sort=False).sum()
            frame.index.name = "gene_id"
            symbol_map = map_ensembl_symbols(
                frame.index.astype(str).tolist(),
                dataset_dir / "ensembl_to_symbol.csv",
            )
            mapped_rows: list[tuple[str, float, pd.Series]] = []
            for gene_id, row in frame.iterrows():
                symbol = symbol_map.get(str(gene_id).split(".", 1)[0])
                if symbol:
                    mapped_rows.append((symbol, float(row.var(ddof=0) or 0.0), row))
            if not mapped_rows:
                raise ValueError("GSE135251 could not map any Ensembl IDs to symbols")
            selected: dict[str, tuple[float, pd.Series]] = {}
            for symbol, variance, row in mapped_rows:
                if symbol not in selected or variance > selected[symbol][0]:
                    selected[symbol] = (variance, row)
            mapped = pd.DataFrame(
                {symbol: value for symbol, (_, value) in selected.items()}
            ).T
            mapped.index.name = "gene"
            mapped = standardize_expression(mapped, already_log=False)
            metadata = _parse_gse135251_metadata(raw_dir / "GSE135251_family.soft.gz")
            metadata = metadata.reindex(mapped.columns)
            save_expression(mapped, expression_path)
            metadata.to_csv(metadata_path, encoding="utf-8")
        outputs["GSE135251_expression"] = expression_path
        outputs["GSE135251_metadata"] = metadata_path
    return outputs


def _parse_gse135251_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"GSE135251 family SOFT not found: {path}")
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE ="):
                if current:
                    rows.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if current is None or not line.startswith("!Sample_characteristics_ch1"):
                continue
            value = line.split("=", 1)[1].strip()
            if ":" not in value:
                continue
            key, item = value.split(":", 1)
            key = re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")
            current[key] = item.strip()
    if current:
        rows.append(current)
    frame = pd.DataFrame(rows).set_index("sample_id")
    if "disease" not in frame.columns:
        raise ValueError("GSE135251 metadata has no disease field")
    disease = frame["disease"].astype(str).str.lower()
    frame["condition"] = np.where(disease.str.contains("control|healthy"), "control", "NAFLD")
    if "fibrosis_stage" in frame.columns:
        frame["fibrosis_stage"] = pd.to_numeric(frame["fibrosis_stage"], errors="coerce")
    return frame


def _platform_mapping(
    probe_ids: list[str],
    processed_dir: Path,
    platform: str,
    cache: dict[str, dict[str, str]],
) -> dict[str, str]:
    if platform in cache:
        return cache[platform]
    cache_path = processed_dir / f"{platform}_probe_to_symbol.csv"
    if cache_path.exists():
        frame = pd.read_csv(cache_path, dtype=str).fillna("")
        mapping = dict(zip(frame["probe_id"], frame["symbol"]))
        cache[platform] = mapping
        return mapping
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if not rscript:
        raise RuntimeError(
            f"Rscript is required to map {platform} probes with its Bioconductor annotation package"
        )
    input_path = cache_path.with_suffix(".input.txt")
    input_path.write_text("\n".join(sorted(set(probe_ids))), encoding="utf-8")
    script = Path(__file__).resolve().parent / "R" / "map_probes.R"
    proc = subprocess.run(
        [rscript, str(script), platform, str(input_path), str(cache_path)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=1800,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"probe mapping failed for {platform}: "
            f"{(proc.stderr or proc.stdout)[-2000:]}"
        )
    frame = pd.read_csv(cache_path, dtype=str).fillna("")
    mapping = dict(zip(frame["probe_id"], frame["symbol"]))
    cache[platform] = mapping
    return mapping


def _normalize_gse89632_condition(row: pd.Series) -> str:
    value = str(row.get("diagnosis") or row.get("Sample_title") or "").upper()
    if "NASH" in value or "STEATOHEPATITIS" in value:
        return "NASH"
    if re.search(r"\bSS\b|STEATOSIS", value):
        return "SS"
    if "HC" in value or "HEALTHY" in value or "CONTROL" in value:
        return "HC"
    title = str(row.get("Sample_title") or "").upper()
    if "_NASH_" in title:
        return "NASH"
    if "_SS_" in title:
        return "SS"
    if "_HC_" in title:
        return "HC"
    raise ValueError(f"cannot infer GSE89632 condition from {value!r}")


def _normalize_gse49541_condition(row: pd.Series) -> str:
    value = str(row.get("stage") or row.get("Sample_characteristics_ch1") or "").lower()
    if "advanced" in value or "3-4" in value:
        return "advanced_fibrosis"
    if "mild" in value or "0-1" in value:
        return "mild_fibrosis"
    raise ValueError(f"cannot infer GSE49541 condition from {value!r}")


def _first_matching(columns: list[Any], candidates: list[str]) -> str | None:
    lookup = {str(column).lower(): str(column) for column in columns}
    return next((lookup[item] for item in candidates if item in lookup), None)


def map_ensembl_symbols(
    ensembl_ids: list[str],
    output_path: Path,
    *,
    timeout: int = 1800,
) -> dict[str, str]:
    """Map Ensembl IDs to HGNC symbols through the local org.Hs.eg.db package."""
    if output_path.exists():
        frame = pd.read_csv(output_path, dtype=str).fillna("")
        return dict(zip(frame["ensembl_id"], frame["symbol"]))
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if not rscript:
        LOG.warning("Rscript unavailable; keeping Ensembl IDs for GSE164441")
        return {}
    input_path = output_path.with_suffix(".input.txt")
    input_path.write_text("\n".join(sorted(set(map(str, ensembl_ids)))), encoding="utf-8")
    script = Path(__file__).resolve().parent / "R" / "map_ensembl.R"
    proc = subprocess.run(
        [rscript, str(script), str(input_path), str(output_path)],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )
    if proc.returncode != 0:
        LOG.warning("Ensembl mapping failed: %s", (proc.stderr or proc.stdout)[-1500:])
        return {}
    frame = pd.read_csv(output_path, dtype=str).fillna("")
    return dict(zip(frame["ensembl_id"], frame["symbol"]))


def differential_expression_limma(
    expression_path: Path,
    metadata_path: Path,
    output_path: Path,
    *,
    condition_column: str,
    group1: list[str],
    group2: list[str],
    comparison: str,
    timeout: int = 1800,
) -> pd.DataFrame:
    """Run limma on a matrix where rows are genes and columns are samples."""
    ensure_dir(output_path.parent)
    if output_path.exists():
        return pd.read_csv(output_path)
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    if not rscript:
        raise RuntimeError("Rscript is required for limma differential expression")
    script = Path(__file__).resolve().parent / "R" / "limma_contrasts.R"
    proc = subprocess.run(
        [
            rscript,
            str(script),
            str(expression_path),
            str(metadata_path),
            str(output_path),
            condition_column,
            ",".join(group1),
            ",".join(group2),
            comparison,
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"limma failed for {comparison}: {(proc.stderr or proc.stdout)[-3000:]}"
        )
    return pd.read_csv(output_path)


def candidate_heatmap(
    expression_path: Path,
    metadata_path: Path,
    genes: list[str],
    output_path: Path,
    *,
    condition_column: str = "condition",
    group_order: list[str] | None = None,
    top_n: int = 30,
) -> dict[str, Any]:
    """Plot candidate genes across GSE89632 groups."""
    expression = read_expression(expression_path)
    metadata = pd.read_csv(metadata_path, index_col=0)
    wanted = [gene for gene in genes if gene in expression.index][:top_n]
    if not wanted:
        raise ValueError("no requested genes are present in the expression matrix")
    selected = zscore_rows(expression.loc[wanted])
    metadata = metadata.loc[[sample for sample in selected.columns if sample in metadata.index]]
    selected = selected[metadata.index]
    order = group_order or sorted(metadata[condition_column].astype(str).unique())
    sample_order: list[str] = []
    for group in order:
        group_samples = metadata.index[metadata[condition_column].astype(str) == group].tolist()
        if len(group_samples) > 2:
            distances = pdist(selected.loc[:, group_samples].T, metric="euclidean")
            groups = _single_linkage_order(len(group_samples), distances)
            group_samples = [group_samples[index] for index in groups]
        sample_order.extend(group_samples)
    selected = selected[sample_order]
    metadata = metadata.loc[sample_order]
    if selected.shape[0] > 2:
        row_order = leaves_list(linkage(selected.to_numpy(), method="average"))
        selected = selected.iloc[row_order]
    fig, (group_ax, ax) = plt.subplots(
        2,
        1,
        figsize=(7.2, max(3.8, selected.shape[0] * 0.18)),
        gridspec_kw={"height_ratios": [0.08, 1.0]},
        sharex=True,
    )
    group_colors = [plt.get_cmap("tab10")(idx) for idx, _ in enumerate(order)]
    group_lookup = {group: color for group, color in zip(order, group_colors)}
    group_strip = np.array(
        [
            [
                order.index(group) if group in order else -1
                for group in metadata[condition_column].astype(str)
            ]
        ],
        dtype=float,
    )
    group_ax.imshow(
        group_strip,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(group_colors),
        vmin=0,
        vmax=max(len(order) - 1, 1),
    )
    group_ax.set_yticks([])
    group_ax.set_ylabel("Group", rotation=0, ha="right", va="center")
    image = ax.imshow(selected.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(np.arange(selected.shape[1]))
    if selected.shape[1] <= 32:
        ax.set_xticklabels(selected.columns, rotation=90, fontsize=5.5)
    else:
        ax.set_xticklabels([])
        ax.set_xlabel("Samples ordered within condition")
    ax.set_yticks(np.arange(selected.shape[0]))
    ax.set_yticklabels(selected.index, fontsize=6)
    boundaries = []
    last = None
    for index, group in enumerate(metadata[condition_column].astype(str)):
        if last is not None and group != last:
            boundaries.append(index - 0.5)
        last = group
    for boundary in boundaries:
        ax.axvline(boundary, color="black", linewidth=0.8)
    handles = [
        plt.Line2D([0], [0], color=group_lookup[group], lw=5, label=group)
        for group in order
    ]
    group_sizes = metadata[condition_column].astype(str).value_counts()
    for handle, group in zip(handles, order):
        handle.set_label(f"{group} (n={int(group_sizes.get(group, 0))})")
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False,
        fontsize=5.5,
    )
    ax.set_title("Candidate-gene expression in GSE89632")
    fig.colorbar(image, ax=ax, shrink=0.4, label="Row z-score")
    save_figure(fig, output_path)
    data_path = output_path.with_suffix(".csv")
    selected.to_csv(data_path)
    return {"figure": output_path, "data": data_path, "genes": wanted}


def _single_linkage_order(count: int, distances: np.ndarray) -> list[int]:
    if count < 3:
        return list(range(count))
    linkage_matrix = linkage(distances, method="single")
    return [int(value) for value in leaves_list(linkage_matrix)]


def validation_boxplots(
    expression_path: Path,
    metadata_path: Path,
    genes: list[str],
    output_dir: Path,
    *,
    condition_column: str = "condition",
    comparisons: tuple[str, str] | None = None,
    prefix: str = "validation",
) -> dict[str, Any]:
    """Draw per-gene boxplots and report two-sided Mann-Whitney tests."""
    expression = read_expression(expression_path)
    metadata = pd.read_csv(metadata_path, index_col=0)
    wanted = [gene for gene in genes if gene in expression.index]
    if not wanted:
        raise ValueError("none of the requested genes are present for validation")
    conditions = sorted(metadata[condition_column].astype(str).unique())
    if comparisons is None:
        if len(conditions) != 2:
            raise ValueError(f"expected two conditions, got {conditions}")
        comparisons = (conditions[0], conditions[1])
    labels = [comparisons[0], comparisons[1]]
    n_cols = min(3, len(wanted))
    n_rows = math.ceil(len(wanted) / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(7.2, 2.9 * n_rows),
        squeeze=False,
    )
    rows: list[dict[str, Any]] = []
    for index, gene in enumerate(wanted):
        ax = axes.flat[index]
        values_by_group = []
        for group in labels:
            samples = metadata.index[metadata[condition_column].astype(str) == group]
            samples = [sample for sample in samples if sample in expression.columns]
            values = pd.to_numeric(expression.loc[gene, samples], errors="coerce").dropna()
            values_by_group.append(values)
            if not values.empty:
                ax.scatter(
                    np.full(len(values), group),
                    values,
                    color="#4a5568" if group == labels[0] else "#c05640",
                    s=12,
                    alpha=0.7,
                    zorder=3,
                )
        ax.boxplot(
            [values.to_numpy() for values in values_by_group],
            labels=labels,
            showfliers=False,
            widths=0.55,
            patch_artist=True,
            boxprops={"facecolor": "#eef2f7", "edgecolor": "#52606d"},
            medianprops={"color": "#1f2933"},
        )
        if len(values_by_group[0]) >= 2 and len(values_by_group[1]) >= 2:
            stat, p_value = stats.mannwhitneyu(
                values_by_group[0],
                values_by_group[1],
                alternative="two-sided",
            )
        else:
            stat, p_value = np.nan, np.nan
        rows.append(
            {
                "gene": gene,
                f"{labels[0]}_n": int(len(values_by_group[0])),
                f"{labels[1]}_n": int(len(values_by_group[1])),
                f"{labels[0]}_median": float(values_by_group[0].median()) if len(values_by_group[0]) else np.nan,
                f"{labels[1]}_median": float(values_by_group[1].median()) if len(values_by_group[1]) else np.nan,
                "mann_whitney_u": stat,
                "p_value": p_value,
            }
        )
        ax.set_title(gene, fontsize=7, fontweight="bold")
        ax.set_ylabel("Expression", fontsize=6)
        ax.tick_params(axis="x", rotation=25, labelsize=5.5)
        ax.tick_params(axis="y", labelsize=5.5)
        p_text = "NA" if not np.isfinite(p_value) else f"p={p_value:.3g}"
        ax.text(
            0.98,
            0.98,
            p_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=5.5,
        )
    for index in range(len(wanted), n_rows * n_cols):
        axes.flat[index].axis("off")
    save_figure(fig, output_dir / f"{prefix}_candidate_boxplots.png")
    results = benjamini_hochberg(pd.DataFrame(rows))
    data_path = output_dir / f"{prefix}_candidate_boxplot_stats.csv"
    results.to_csv(data_path, index=False)
    return {
        "figure": output_dir / f"{prefix}_candidate_boxplots.png",
        "statistics": data_path,
        "genes": wanted,
    }


def differential_summary_table(
    result: pd.DataFrame,
    *,
    fdr_cutoff: float = 0.05,
    logfc_cutoff: float = 0.5,
) -> pd.DataFrame:
    frame = benjamini_hochberg(result, p_column="P.Value")
    frame["significant"] = (frame["adj.P.Val"] <= fdr_cutoff) & (
        frame["logFC"].abs() >= logfc_cutoff
    )
    frame["direction"] = np.where(frame["logFC"] >= 0, "up", "down")
    return frame.sort_values(["adj.P.Val", "P.Value"], na_position="last")
