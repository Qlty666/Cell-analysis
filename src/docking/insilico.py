#!/usr/bin/env python3
"""CellOracle-inspired single-cell in-silico knockout analysis.

The implementation is deliberately lightweight: it uses the same building
blocks as CellOracle (KNN-imputed expression, a GRN coefficient matrix, and
iterative signal propagation), but builds candidate regulator edges from
expression-correlation unless a regulator CSV is supplied. It also combines
scTenifoldKnk-style differential-regulation summaries and optional GO/KEGG
enrichment into one report.

This is a screening/prediction heuristic, not a real knockout experiment.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .config import ResolvedConfig
from .utils import DockingError, write_json

APP_ROOT = Path(__file__).resolve().parent.parent.parent

_CELL_ID_COLS = {
    "cell",
    "barcode",
    "cell_id",
    "cell_barcode",
    "sample",
    "sample_id",
    "index",
}
_CELL_TYPE_COLS = {
    "cell_type",
    "celltype",
    "cell_types",
    "annotation",
    "annotations",
    "cell_annotation",
    "cluster_label",
    "cell_type_annotation",
    "louvain_annot",
}
_TARGET_EXCLUDE = re.compile(
    r"^(MT-|MTRNR|RPL|RPS|MRPL|MRPS|SNORD|SNORA|SCGB|IGH|IGK|IGL|TRA|TRB|TRG|"
    r"HLA-D|LINC|RP[0-9]|AC[0-9]|AL[0-9])",
    re.I,
)

# A compact, species-neutral set of commonly used transcription factors. The
# KO gene is always included even when it is not in this list.
DEFAULT_TRANSCRIPTION_FACTORS = [
    "ARNT", "ATF3", "ATF4", "BATF", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CTCF",
    "E2F1", "E2F2", "EGR1", "EGR2", "EP300", "ERG", "ETS1", "ETS2", "FLI1",
    "FOS", "FOXO1", "FOXO3", "GATA1", "GATA2", "GATA3", "GATA4", "GATA6",
    "GFI1", "GFI1B", "HIF1A", "IRF1", "IRF2", "IRF7", "IRF8", "JUN", "JUNB",
    "JUND", "KLF1", "KLF2", "KLF3", "KLF4", "KLF6", "LMO2", "MAFB", "MEIS1",
    "MYB", "MYBL1", "MYC", "NANOG", "NFE2", "NFE2L2", "NFKB1", "NR4A1",
    "NR4A2", "PBX1", "POU5F1", "RELA", "RUNX1", "RUNX2", "SOX2", "SOX4",
    "SOX9", "SPI1", "SPIB", "STAT1", "STAT3", "STAT5A", "STAT5B", "TAL1",
    "TP53", "TWIST1", "ZEB1", "ZEB2", "ZFP36",
]

# Keyed by upper-case symbol. Entries are intentionally conservative and are
# used only to provide a short biological context in the report.
_KNOWN_GENE_NOTES = {
    "GATA1": (
        "红细胞系和巨核细胞系分化的关键主效转录因子；其缺失通常抑制红系成熟 "
        "基因程序，并促进髓系/粒细胞方向的状态偏转。"
    )
}


def run_insilico_knockout(
    cfg: ResolvedConfig,
    log,
    ko_gene: str | None = None,
) -> dict:
    """Run a single-cell in-silico knockout analysis and write report outputs."""
    isko = cfg.data.get("insilico_knockout", {})
    knockout = cfg.data.get("knockout", {})
    enabled = bool(isko.get("enabled", False)) or bool(ko_gene)
    if not enabled:
        return {"status": "skipped", "reason": "in-silico knockout disabled"}

    gene = _first_gene(
        ko_gene,
        isko.get("ko_gene"),
        knockout.get("ko_gene"),
    )
    if not gene:
        raise DockingError(
            "in-silico knockout requires --insilico-gene or "
            "insilico_knockout.ko_gene"
        )

    expression_csv = _resolve_path(cfg, knockout.get("expression_csv"))
    metadata_csv = _resolve_path(cfg, knockout.get("metadata_csv"))
    if expression_csv is None or not expression_csv.exists():
        raise DockingError(
            "in-silico knockout expression CSV not found: "
            f"{expression_csv}"
        )

    matrix = _load_wide_expression(expression_csv, knockout)
    if matrix.shape[1] < int(isko.get("min_cells", 100)):
        raise DockingError(
            "in-silico knockout needs a cell-level expression matrix; "
            f"found {matrix.shape[1]} columns, minimum "
            f"{isko.get('min_cells', 100)} cells"
        )
    if gene not in matrix.index:
        raise DockingError(
            f"knockout gene {gene} is not in the expression matrix"
        )

    metadata = _load_cell_metadata(metadata_csv, matrix.columns, knockout)
    cell_types = _cell_type_series(metadata, matrix.columns)
    embedding = _load_embedding(cfg, isko, matrix.columns, metadata)
    species = _resolve_species(cfg, isko, metadata)

    max_cells = int(isko.get("max_cells", 5000))
    if len(matrix.columns) > max_cells:
        rng = np.random.default_rng(int(isko.get("seed", 123)))
        keep = sorted(
            rng.choice(len(matrix.columns), size=max_cells, replace=False)
        )
        matrix = matrix.iloc[:, keep]
        if cell_types is not None:
            cell_types = cell_types.iloc[keep]
        if embedding is not None:
            embedding = embedding.iloc[keep]

    log.info(
        "in-silico knockout: %s across %s cells x %s genes",
        gene,
        matrix.shape[1],
        matrix.shape[0],
    )

    log_mat = _to_log_normalized(matrix)
    selected = _select_genes(log_mat, gene, int(isko.get("max_genes", 1800)))
    log_mat = log_mat.loc[selected]
    if cell_types is not None:
        cell_types = cell_types.reindex(log_mat.columns)
    log_mat = log_mat.T

    embedding_df = _ensure_embedding(
        cfg,
        isko,
        log_mat,
        embedding,
        cell_names=log_mat.index,
    )
    if not {"cell", "umap_1", "umap_2"}.issubset(embedding_df.columns):
        raise DockingError("embedding must contain two numeric coordinates")
    embedding_coords = (
        embedding_df.set_index("cell")[["umap_1", "umap_2"]]
        .reindex(log_mat.index.astype(str))
        .dropna()
    )
    if len(embedding_coords) != len(log_mat):
        raise DockingError("embedding rows do not cover all cells")

    pca_model, pca_space = _pca_space(log_mat, isko)
    gem = _knn_impute(log_mat, pca_space, isko)
    regulators = _resolve_regulators(cfg, isko, log_mat, gene)

    log.info(
        "building GRN with %s regulator rows, %s genes",
        len(regulators),
        gem.shape[1],
    )
    coef = _coefficient_matrix(
        gem,
        regulators,
        max_edges=int(isko.get("network_edges_per_regulator", 300)),
        seed=int(isko.get("seed", 123)),
    )
    ko_sim, delta = _simulate_shift(
        gem,
        coef,
        gene,
        n_propagation=int(isko.get("n_propagation", 3)),
    )
    shift = _embedding_shift(
        gem,
        ko_sim,
        pca_model,
        embedding_coords.to_numpy(),
        n_neighbors=int(isko.get("n_neighbors", 100)),
        seed=int(isko.get("seed", 123)),
    )

    wt_mean = gem.mean(axis=0)
    ko_mean = ko_sim.mean(axis=0)
    changes = pd.DataFrame(
        {
            "gene": gem.columns,
            "wt_mean": wt_mean.to_numpy(),
            "ko_mean": ko_mean.to_numpy(),
            "delta": (ko_mean - wt_mean).to_numpy(),
        }
    )
    changes["delta"] = changes["delta"].fillna(0.0)
    changes["direction"] = np.where(changes["delta"] < 0, "down", "up")
    changes["abs_delta"] = changes["delta"].abs()
    changes = changes.sort_values(
        ["abs_delta", "gene"], ascending=[False, True]
    ).reset_index(drop=True)

    out_dir = cfg.knockout_dir() / "in_silico"
    data_dir = out_dir / "data"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    figures: list[str] = []
    if isko.get("figures", True):
        figures = _write_figures(
            fig_dir,
            data_dir,
            gem,
            ko_sim,
            changes,
            coef,
            gene,
            cell_types,
            embedding_coords,
            shift,
            isko,
        )

    data_files = _write_data_files(
        data_dir,
        gem,
        ko_sim,
        changes,
        coef,
        gene,
        cell_types,
        embedding_coords,
        shift,
    )

    enrichment = _run_enrichment(
        cfg,
        isko,
        data_dir,
        changes,
        species,
        log,
    )
    if enrichment and enrichment.get("figures"):
        figures.extend(enrichment["figures"])
    if enrichment and enrichment.get("status") == "completed":
        for name in ("insilico_go_enrichment.csv", "insilico_kegg_enrichment.csv"):
            path = data_dir / name
            if path.exists():
                data_files[f"{Path(name).stem}_csv"] = path

    report_path = _write_html_report(
        out_dir,
        fig_dir,
        data_dir,
        gene,
        species,
        matrix,
        cell_types,
        changes,
        enrichment,
        isko,
    )
    summary = {
        "status": "completed",
        "ko_gene": gene,
        "cells": int(len(gem)),
        "genes_modeled": int(gem.shape[1]),
        "regulators": [str(g) for g in regulators],
        "n_propagation": int(isko.get("n_propagation", 3)),
        "embedding": str(embedding_coords.columns.tolist()),
        "species": species,
        "top_targets": changes.head(
            int(isko.get("target_top_n", 15))
        )["gene"].astype(str).tolist(),
        "data": {name: str(path) for name, path in data_files.items()},
        "figures": figures,
        "report": str(report_path),
        "enrichment": enrichment or {},
        "warning": (
            "Prediction is a network heuristic; wet-lab validation is required."
        ),
    }
    write_json(out_dir / "insilico_summary.json", summary)
    _copy_photo_outputs(cfg, isko, summary)
    log.info("in-silico knockout complete: output %s", out_dir)
    return summary


def _first_gene(*values) -> str | None:
    for value in values:
        if value:
            text = str(value).strip()
            if text and str(text).lower() not in ("nan", "none", "auto"):
                return text
    return None


def _resolve_path(cfg: ResolvedConfig, value) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = cfg.workdir / path
    return path.resolve()


def _load_wide_expression(path: Path, ko: dict) -> pd.DataFrame:
    df = pd.read_csv(path)
    if df.empty:
        raise DockingError(f"expression CSV is empty: {path}")
    gene_col = ko.get("gene_column") or "gene"
    sample_col = ko.get("sample_column") or "sample"
    value_col = ko.get("value_column") or "value"
    if value_col in df.columns and gene_col in df.columns:
        matrix = df.pivot_table(
            index=gene_col,
            columns=sample_col,
            values=value_col,
            aggfunc="mean",
        )
    else:
        first_col = df.columns[0]
        matrix = df.set_index(first_col).apply(pd.to_numeric, errors="coerce")
    matrix = matrix.dropna(how="all").fillna(0.0)
    matrix.columns = matrix.columns.astype(str)
    matrix.index = matrix.index.astype(str)
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    return matrix


def _load_cell_metadata(
    path: Path | None,
    cell_names: pd.Index,
    ko: dict,
) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    df = pd.read_csv(path)
    if df.empty:
        return None
    cell_col = next(
        (c for c in df.columns if c.lower() in _CELL_ID_COLS),
        None,
    )
    if cell_col is None:
        return None
    df = df[df[cell_col].astype(str).isin(cell_names.astype(str))]
    return df


def _cell_type_series(
    metadata: pd.DataFrame | None,
    cell_names: pd.Index,
) -> pd.Series | None:
    if metadata is None or metadata.empty:
        return None
    cell_col = next(
        (c for c in metadata.columns if c.lower() in _CELL_ID_COLS),
        None,
    )
    type_col = next(
        (c for c in metadata.columns if c.lower() in _CELL_TYPE_COLS),
        None,
    )
    if cell_col is None or type_col is None:
        return None
    pairs = metadata[[cell_col, type_col]].dropna().drop_duplicates(cell_col)
    mapping = dict(zip(pairs[cell_col].astype(str), pairs[type_col].astype(str)))
    labels = pd.Series(
        [mapping.get(name, "Unannotated") for name in cell_names],
        index=cell_names,
    )
    if labels.nunique() < 1:
        return None
    return labels


def _load_embedding(
    cfg: ResolvedConfig,
    isko: dict,
    cell_names: pd.Index,
    metadata: pd.DataFrame | None,
) -> pd.DataFrame | None:
    path = _resolve_path(cfg, isko.get("embedding_csv"))
    df = None
    if path is not None and path.exists():
        df = pd.read_csv(path)
    elif metadata is not None:
        xcol = next((c for c in metadata.columns if c.lower() in ("umap_1", "umap1", "x")), None)
        ycol = next((c for c in metadata.columns if c.lower() in ("umap_2", "umap2", "y")), None)
        cell_col = next(
            (c for c in metadata.columns if c.lower() in _CELL_ID_COLS),
            None,
        )
        if cell_col is not None and xcol is not None and ycol is not None:
            df = metadata[[cell_col, xcol, ycol]].copy()
            df.columns = ["cell", "umap_1", "umap_2"]
    if df is None or df.empty:
        return None
    cell_col = next(
        (c for c in df.columns if c.lower() in _CELL_ID_COLS),
        None,
    )
    if cell_col is None:
        df = df.reset_index().rename(columns={"index": "cell"})
        cell_col = "cell"
    xcol = next(
        (c for c in df.columns if c.lower() in ("umap_1", "umap1", "x", "dim1")),
        None,
    )
    ycol = next(
        (c for c in df.columns if c.lower() in ("umap_2", "umap2", "y", "dim2")),
        None,
    )
    if xcol is None or ycol is None:
        return None
    df = df[[cell_col, xcol, ycol]].copy()
    df.columns = ["cell", "umap_1", "umap_2"]
    df["umap_1"] = pd.to_numeric(df["umap_1"], errors="coerce")
    df["umap_2"] = pd.to_numeric(df["umap_2"], errors="coerce")
    df = df.dropna().drop_duplicates("cell")
    df = df.set_index("cell").reindex(cell_names.astype(str)).dropna()
    df = df.reset_index()
    if "index" in df.columns:
        df = df.rename(columns={"index": "cell"})
    return df


def _resolve_species(cfg: ResolvedConfig, isko: dict, metadata) -> str:
    value = str(isko.get("species") or cfg.get("knockout", "species") or "").lower()
    if value in ("hs", "mm", "human", "mouse", "mus_musculus", "homo_sapiens"):
        return "mm" if value.startswith(("mm", "mouse", "mus")) else "hs"
    if metadata is not None:
        text = " ".join(
            metadata[[c for c in metadata.columns if c.lower() in ("organism", "species")]]
            .astype(str)
            .to_numpy()
            .ravel()
        ).lower()
        if "mouse" in text or "mus musculus" in text:
            return "mm"
    return "hs"


def _to_log_normalized(matrix: pd.DataFrame) -> pd.DataFrame:
    """Return a log-scale expression matrix, handling raw count input."""
    df = matrix.astype(float)
    flat = df.to_numpy()
    finite = np.isfinite(flat)
    integer_like = finite & (flat == np.round(flat))
    fraction_integer = float(integer_like.sum() / max(1, finite.sum()))
    if fraction_integer > 0.6 and np.nanmax(flat) > 20:
        totals = df.sum(axis=0).replace(0, np.nan)
        scale = 10_000.0 / totals
        return np.log1p(df.multiply(scale, axis=1)).fillna(0.0)
    if np.nanmax(flat) > 20:
        return np.log1p(df).fillna(0.0)
    return df.fillna(0.0)


def _select_genes(
    log_mat: pd.DataFrame,
    ko_gene: str,
    max_genes: int,
) -> pd.Index:
    detection = (log_mat > 0).mean(axis=1)
    variance = log_mat.var(axis=1)
    score = (detection * (log_mat.mean(axis=1) + 0.01) * variance).fillna(0.0)
    score = score.sort_values(ascending=False)
    selected = list(score.head(max_genes).index)
    if ko_gene not in selected:
        selected.append(ko_gene)
    return log_mat.index[log_mat.index.isin(selected)]


def _pca_space(log_mat: pd.DataFrame, isko: dict):
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    n_comp = min(
        30,
        log_mat.shape[0],
        log_mat.shape[1] - 1,
        max(2, len(log_mat.columns) - 1),
    )
    scaler = StandardScaler()
    scaled = scaler.fit_transform(log_mat.to_numpy())
    pca = PCA(n_components=max(2, int(n_comp)), random_state=int(isko.get("seed", 123)))
    pcs = pca.fit_transform(scaled)
    return pca, pcs


def _knn_impute(
    log_mat: pd.DataFrame,
    pca_space: np.ndarray,
    isko: dict,
) -> pd.DataFrame:
    from sklearn.neighbors import NearestNeighbors

    k = int(isko.get("knn_impute_neighbors", 66))
    k = min(k, len(log_mat.index), max(2, len(log_mat.index) - 1))
    if k < 2:
        return log_mat.copy()
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(pca_space)
    _, idx = nn.kneighbors(pca_space)
    imputed = log_mat.to_numpy()[idx].mean(axis=1)
    return pd.DataFrame(
        imputed,
        index=log_mat.index,
        columns=log_mat.columns,
    )


def _resolve_regulators(
    cfg: ResolvedConfig,
    isko: dict,
    log_mat: pd.DataFrame,
    ko_gene: str,
) -> list[str]:
    genes = list(log_mat.columns)
    upper_to_gene = {str(g).upper(): str(g) for g in genes}
    regulators: list[str] = []
    path = _resolve_path(cfg, isko.get("regulators_csv"))
    if path is not None and path.exists():
        df = pd.read_csv(path)
        col = next(
            (
                c
                for c in df.columns
                if c.lower() in ("regulator", "symbol", "gene", "tf", "hgnc")
            ),
            df.columns[0],
        )
        regulators = df[col].astype(str).tolist()
    else:
        regulators = [
            upper_to_gene[sym.upper()]
            for sym in DEFAULT_TRANSCRIPTION_FACTORS
            if sym.upper() in upper_to_gene
        ]
    regulators = list(dict.fromkeys(regulators))
    if ko_gene not in regulators:
        regulators.append(ko_gene)
    return [g for g in regulators if g in genes]


def _ensure_embedding(
    cfg: ResolvedConfig,
    isko: dict,
    log_mat: pd.DataFrame,
    provided: pd.DataFrame | None,
    cell_names: pd.Index | None = None,
) -> pd.DataFrame:
    """Return cell/umap_1/umap_2 coordinates aligned to log_mat columns."""
    names = log_mat.index if cell_names is not None else log_mat.columns
    names = pd.Index(names).astype(str)
    if provided is not None and not provided.empty:
        df = provided.copy()
        if "cell" not in df.columns:
            cell_col = next(
                (c for c in df.columns if c.lower() in _CELL_ID_COLS),
                "index" if "index" in df.columns else df.columns[0],
            )
            if cell_col != "cell":
                df = df.rename(columns={cell_col: "cell"})
            if "index" in df.columns and df["index"].equals(df.index):
                df = df.drop(columns=["index"])
        df = (
            df.set_index("cell")
            .reindex(names)
            .dropna()
        )
        if len(df) == len(names):
            df = df.reset_index()
            if "index" in df.columns:
                df = df.rename(columns={"index": "cell"})
            return df

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    n_comp = min(
        30,
        len(log_mat.index),
        len(log_mat.columns),
        max(2, len(log_mat.columns) - 1),
    )
    scaled = StandardScaler().fit_transform(log_mat.to_numpy())
    pca = PCA(n_components=max(2, int(n_comp)), random_state=42)
    pcs = pca.fit_transform(scaled)
    try:
        import umap

        reducer = umap.UMAP(
            n_neighbors=min(30, max(5, len(log_mat.index) - 1)),
            min_dist=0.3,
            random_state=int(isko.get("seed", 123)),
        )
        coords = reducer.fit_transform(pcs)
    except Exception:
        coords = pcs[:, :2]
    return pd.DataFrame(
        {
            "cell": names.to_numpy(),
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
        }
    )


def _coefficient_matrix(
    gem: pd.DataFrame,
    regulators: list[str],
    max_edges: int,
    seed: int,
) -> pd.DataFrame:
    """Build a sparse directed GRN from expression correlations.

    Each row is a regulator. Pairwise slopes preserve the sign of the
    correlation and are then sparse-filtered and spectrally damped to keep
    iterative propagation stable.
    """
    genes = list(gem.columns)
    regulators = [r for r in regulators if r in genes]
    if not regulators:
        return pd.DataFrame(0.0, index=genes, columns=genes)
    reg_matrix = gem[regulators].to_numpy()
    target_matrix = gem.to_numpy()
    reg_center = reg_matrix - reg_matrix.mean(axis=0)
    target_center = target_matrix - target_matrix.mean(axis=0)
    reg_sd = np.linalg.norm(reg_center, axis=0)
    target_sd = np.linalg.norm(target_center, axis=0)
    denom = reg_sd[:, None] * target_sd[None, :]
    corr = reg_center.T @ target_center / np.where(denom == 0, np.nan, denom)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    beta = corr * (
        np.nan_to_num(gem.std(axis=0).to_numpy()[None, :])
        / np.nan_to_num(gem[regulators].std(axis=0).to_numpy()[:, None])
    )
    beta = np.nan_to_num(beta, nan=0.0, posinf=0.0, neginf=0.0)
    for i, reg in enumerate(regulators):
        if reg in genes:
            beta[i, genes.index(reg)] = 0.0
        row = np.abs(beta[i])
        k = min(max_edges, len(row))
        if k < len(row):
            threshold = np.partition(row, len(row) - k)[len(row) - k]
            beta[i, row < threshold] = 0.0

    coef = pd.DataFrame(0.0, index=genes, columns=genes)
    coef.loc[regulators] = beta
    np.fill_diagonal(coef.values, 0.0)
    singular = float(np.linalg.svd(coef.to_numpy(), compute_uv=False)[0])
    if singular > 1.0:
        coef = coef * (0.95 / singular)
    return coef


def _simulate_shift(
    gem: pd.DataFrame,
    coef: pd.DataFrame,
    ko_gene: str,
    n_propagation: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Iteratively propagate a KO perturbation through the GRN."""
    from scipy import sparse

    simulation_input = gem.copy()
    simulation_input[ko_gene] = 0.0
    delta_input = simulation_input - gem
    delta = delta_input.copy()
    coef_sparse = sparse.csr_matrix(coef.to_numpy())
    for _ in range(max(1, n_propagation)):
        values = delta.to_numpy() @ coef_sparse
        delta = pd.DataFrame(values, index=gem.index, columns=gem.columns)
        delta.loc[:, ko_gene] = delta_input.loc[:, ko_gene]
        simulated = gem + delta
        simulated[simulated < 0] = 0.0
        delta = simulated - gem
    return gem + delta, delta


def _embedding_shift(
    gem: pd.DataFrame,
    ko_sim: pd.DataFrame,
    pca_model,
    embedding: np.ndarray,
    n_neighbors: int,
    seed: int,
) -> pd.DataFrame:
    """Approximate cell state movement by transition to nearest WT cells."""
    from sklearn.neighbors import NearestNeighbors

    mean = gem.mean(axis=0).to_numpy()
    sd = np.where(gem.std(axis=0).to_numpy() == 0, 1.0, gem.std(axis=0).to_numpy())
    wt_std = (gem.to_numpy() - mean) / sd
    future_std = (ko_sim.to_numpy() - mean) / sd
    future_pcs = pca_model.transform(future_std)
    k = min(n_neighbors + 1, len(gem), max(3, len(gem) - 1))
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean")
    nn.fit(pca_model.transform(wt_std))
    distances, indices = nn.kneighbors(future_pcs)
    sigma = np.median(distances[:, -1]) if distances.size else 1.0
    if sigma <= 0:
        sigma = 1.0
    weights = np.exp(-(distances**2) / (2.0 * sigma**2))
    weights /= weights.sum(axis=1, keepdims=True)
    shifts = np.zeros((len(gem), 2))
    for i in range(len(gem)):
        neighbors = indices[i]
        # Exclude the cell itself if present, then fall back to all neighbors.
        mask = neighbors != i
        if mask.any():
            neighbors = neighbors[mask]
            w = weights[i][mask]
        else:
            w = weights[i]
        w = w / w.sum()
        shifts[i] = (w[:, None] * (embedding[neighbors] - embedding[i])).sum(axis=0)
    return pd.DataFrame(
        shifts,
        index=gem.index,
        columns=["shift_1", "shift_2"],
    )


def _display_changes(changes: pd.DataFrame, blacklist: bool = True) -> pd.DataFrame:
    df = changes.copy()
    if blacklist:
        df = df[~df["gene"].astype(str).str.match(_TARGET_EXCLUDE)]
    return df.reset_index(drop=True)


def _write_figures(
    fig_dir: Path,
    data_dir: Path,
    gem: pd.DataFrame,
    ko_sim: pd.DataFrame,
    changes: pd.DataFrame,
    coef: pd.DataFrame,
    ko_gene: str,
    cell_types: pd.Series | None,
    embedding: pd.DataFrame,
    shift: pd.DataFrame,
    isko: dict,
) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures: list[str] = []
    n_targets = int(isko.get("target_top_n", 15))
    display = _display_changes(changes).head(n_targets)
    table_path = fig_dir / "fig_63_ko_target_expression_table.png"
    _plot_target_table(display, ko_gene, table_path)
    figures.append(table_path.name)
    plt.close("all")

    bar_top = int(isko.get("bar_top_n", 8))
    bar = _display_changes(changes).head(max(2, bar_top))
    bar_path = fig_dir / "fig_64_ko_target_expression_bar.png"
    _plot_target_bar(gem, ko_sim, bar, ko_gene, bar_path)
    figures.append(bar_path.name)
    plt.close("all")

    network_path = fig_dir / "fig_65_ko_regulatory_network.png"
    if _plot_network(coef, ko_gene, network_path):
        figures.append(network_path.name)
    plt.close("all")

    umap_path = fig_dir / "fig_66_ko_shift_umap.png"
    _plot_umap_shift(embedding, shift, cell_types, ko_gene, umap_path)
    figures.append(umap_path.name)
    plt.close("all")
    return figures


def _plot_target_table(
    frame: pd.DataFrame,
    ko_gene: str,
    path: Path,
) -> None:
    import matplotlib.pyplot as plt

    frame = frame.head(15).copy()
    rows = [
        [
            str(row.gene),
            f"{row.wt_mean:.4f}",
            f"{row.ko_mean:.4f}",
            f"{row.delta:.4f}",
            "up" if row.delta > 0 else "down",
        ]
        for row in frame.itertuples(index=False)
    ]
    fig, ax = plt.subplots(figsize=(10, max(3.0, 0.42 * len(frame) + 2)))
    ax.axis("off")
    cell_colours = [["#f7f7f7"] * 5 for _ in range(len(frame))]
    for i, row in enumerate(frame.itertuples(index=False)):
        colour = "#d9534f" if row.delta < 0 else "#3d9a50"
        cell_colours[i][2] = colour
        cell_colours[i][3] = colour
    table = ax.table(
        cellText=rows,
        colLabels=["Gene", "WT mean", "KO mean", "Delta", "Trend"],
        cellLoc="center",
        loc="center",
        cellColours=cell_colours,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.35)
    ax.set_title(f"{ko_gene} KO Target Expression Changes (Top 15)", pad=18)
    fig.tight_layout()
    fig.savefig(path, dpi=160)


def _plot_target_bar(
    gem: pd.DataFrame,
    ko_sim: pd.DataFrame,
    bar: pd.DataFrame,
    ko_gene: str,
    path: Path,
) -> None:
    import matplotlib.pyplot as plt

    genes = bar["gene"].astype(str).tolist()
    wt = bar["wt_mean"].to_numpy()
    ko = bar["ko_mean"].to_numpy()
    x = np.arange(len(genes))
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.bar(x - 0.2, wt, width=0.4, color="#8c8c8c", label="Wild Type (WT)")
    ax.bar(x + 0.2, ko, width=0.4, color="#2f6db0", label="Knockout (KO)")
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Normalized Log Expression")
    ax.set_xlabel("Gene Symbol")
    ax.set_title(f"Target Gene Expression Changes Post-{ko_gene} Knockout")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)


def _plot_network(
    coef: pd.DataFrame,
    ko_gene: str,
    path: Path,
    top_n: int = 20,
) -> bool:
    import matplotlib.pyplot as plt

    if ko_gene not in coef.index:
        return False
    row = coef.loc[ko_gene].drop(index=[ko_gene])
    row = row[row != 0]
    row = row.reindex(row.abs().sort_values(ascending=False).index).head(top_n)
    if row.empty:
        return False
    fig, ax = plt.subplots(figsize=(9, 9))
    angles = np.linspace(0, 2 * math.pi, len(row), endpoint=False)
    radius = 1.0
    positions = {gene: (radius * math.cos(a), radius * math.sin(a)) for gene, a in zip(row.index, angles)}
    positions[ko_gene] = (0.0, 0.0)
    max_w = float(row.abs().max())
    for gene, weight in row.items():
        x0, y0 = positions[ko_gene]
        x1, y1 = positions[gene]
        colour = "#2e8b57" if weight > 0 else "#d1495b"
        ax.plot(
            [x0, x1],
            [y0, y1],
            color=colour,
            alpha=0.75,
            linewidth=0.8 + 5.0 * abs(weight) / max_w,
        )
    for gene, (x, y) in positions.items():
        if gene == ko_gene:
            ax.scatter(x, y, s=850, color="#2f6db0", edgecolor="white", zorder=4)
            ax.text(x, y, gene, ha="center", va="center", color="white", fontsize=10, fontweight="bold")
        else:
            weight = row.get(gene, 0.0)
            ax.scatter(x, y, s=380, color="#dce8f5", edgecolor="#2f6db0", zorder=3)
            ax.text(x, y, gene, ha="center", va="center", fontsize=7.5)
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"{ko_gene} Local Regulatory Network Weights")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    return True


def _plot_umap_shift(
    embedding: pd.DataFrame,
    shift: pd.DataFrame,
    cell_types: pd.Series | None,
    ko_gene: str,
    path: Path,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.5, 7))
    if cell_types is not None and cell_types.nunique() <= 30:
        from matplotlib import colormaps

        cmap = colormaps["tab20"]
        labels = cell_types.astype(str).to_numpy()
        for i, label in enumerate(sorted(pd.unique(labels))):
            mask = labels == label
            ax.scatter(
                embedding.loc[mask, "umap_1"],
                embedding.loc[mask, "umap_2"],
                s=6,
                color=cmap(i % 20),
                label=label,
                alpha=0.8,
            )
        ax.legend(frameon=False, fontsize=7, loc="lower left", markerscale=1.5)
    else:
        ax.scatter(embedding["umap_1"], embedding["umap_2"], s=6, color="#9aa7b0", alpha=0.8)
    max_len = np.percentile(np.linalg.norm(shift[["shift_1", "shift_2"]].to_numpy(), axis=1), 95)
    if max_len > 0:
        u = shift["shift_1"].to_numpy() / max_len
        v = shift["shift_2"].to_numpy() / max_len
        if len(u) > 700:
            n_grid = 26
            x = embedding["umap_1"].to_numpy()
            y = embedding["umap_2"].to_numpy()
            pad_x = max(1e-6, (x.max() - x.min()) * 0.04)
            pad_y = max(1e-6, (y.max() - y.min()) * 0.04)
            hist, x_edges, y_edges = np.histogram2d(
                x,
                y,
                bins=n_grid,
                range=[[x.min() - pad_x, x.max() + pad_x], [y.min() - pad_y, y.max() + pad_y]],
            )
            sum_u, _, _ = np.histogram2d(
                x,
                y,
                bins=n_grid,
                range=[[x.min() - pad_x, x.max() + pad_x], [y.min() - pad_y, y.max() + pad_y]],
                weights=u,
            )
            sum_v, _, _ = np.histogram2d(
                x,
                y,
                bins=n_grid,
                range=[[x.min() - pad_x, x.max() + pad_x], [y.min() - pad_y, y.max() + pad_y]],
                weights=v,
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                grid_u = np.where(hist > 0, sum_u / np.maximum(hist, 1), 0.0)
                grid_v = np.where(hist > 0, sum_v / np.maximum(hist, 1), 0.0)
            gx = (x_edges[:-1] + x_edges[1:]) / 2
            gy = (y_edges[:-1] + y_edges[1:]) / 2
            Xg, Yg = np.meshgrid(gx, gy)
            valid = (hist >= 2) & (np.hypot(grid_u, grid_v) > 1e-6)
            ax.quiver(
                Xg[valid],
                Yg[valid],
                grid_u[valid],
                grid_v[valid],
                color="#122a45",
                alpha=0.62,
                angles="xy",
                scale_units="xy",
                scale=1.0,
                width=0.004,
                headwidth=3.0,
                headlength=3.2,
            )
        else:
            step = max(1, len(u) // 700)
            ax.quiver(
                embedding["umap_1"].to_numpy()[::step],
                embedding["umap_2"].to_numpy()[::step],
                u[::step],
                v[::step],
                color="#122a45",
                alpha=0.45,
                angles="xy",
                scale_units="xy",
                scale=1.0,
                width=0.002,
            )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP Projection of {ko_gene} Knockout Shift Vectors")
    fig.tight_layout()
    fig.savefig(path, dpi=160)


def _write_data_files(
    data_dir: Path,
    gem: pd.DataFrame,
    ko_sim: pd.DataFrame,
    changes: pd.DataFrame,
    coef: pd.DataFrame,
    ko_gene: str,
    cell_types: pd.Series | None,
    embedding_coords: pd.DataFrame,
    shift: pd.DataFrame,
) -> dict[str, Path]:
    files: dict[str, Path] = {}
    cell_df = pd.DataFrame(
        {
            "cell": gem.index,
            "ko_expression": ko_sim[ko_gene].to_numpy(),
            "wt_expression": gem[ko_gene].to_numpy(),
            "delta_magnitude": (ko_sim.to_numpy() - gem.to_numpy()).std(axis=1),
        }
    )
    if cell_types is not None:
        cell_df["cell_type"] = cell_types.astype(str).to_numpy()
    cell_df = cell_df.set_index("cell")
    cell_df = cell_df.join(embedding_coords)
    cell_df = cell_df.join(shift)
    cell_df = cell_df.reset_index()
    cell_path = data_dir / "insilico_cell_shift.csv"
    cell_df.to_csv(cell_path, index=False)
    files["cell_shift_csv"] = cell_path

    change_path = data_dir / "insilico_target_changes.csv"
    changes.to_csv(change_path, index=False)
    files["target_changes_csv"] = change_path

    top = _display_changes(changes).head(15)
    top_path = data_dir / "fig_63_ko_target_top15.csv"
    top.to_csv(top_path, index=False)
    files["top15_csv"] = top_path

    edge_records = []
    if ko_gene in coef.index:
        row = coef.loc[ko_gene].drop(index=[ko_gene]).sort_values(key=abs, ascending=False)
        row = row[row != 0]
        for target, weight in row.head(200).items():
            edge_records.append(
                {
                    "regulator": ko_gene,
                    "target": target,
                    "weight": float(weight),
                    "regulation": "activation" if weight > 0 else "repression",
                }
            )
    edge_path = data_dir / "insilico_regulatory_edges.csv"
    pd.DataFrame(edge_records).to_csv(edge_path, index=False)
    files["regulatory_edges_csv"] = edge_path
    return files


def _run_enrichment(
    cfg: ResolvedConfig,
    isko: dict,
    data_dir: Path,
    changes: pd.DataFrame,
    species: str,
    log,
) -> dict:
    if not isko.get("run_enrichment", True):
        return {"status": "skipped", "reason": "enrichment disabled"}
    genes = (
        _display_changes(changes)
        .head(int(isko.get("enrichment_genes", 150)))["gene"]
        .astype(str)
        .tolist()
    )
    if len(genes) < 3:
        return {"status": "skipped", "reason": "too few target genes"}
    gene_csv = data_dir / "insilico_enrichment_input.csv"
    pd.DataFrame({"gene": genes}).to_csv(gene_csv, index=False)
    script = APP_ROOT / "src" / "docking" / "insilico_enrichment.R"
    if not script.exists():
        return {"status": "skipped", "reason": "enrichment R script missing"}
    rscript = _find_rscript()
    if rscript is None:
        return {"status": "skipped", "reason": "Rscript not found"}
    proc = subprocess.run(
        [rscript, str(script), str(gene_csv), str(data_dir), species],
        cwd=str(APP_ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=int(isko.get("enrichment_timeout", 900)),
    )
    if proc.returncode != 0:
        log.warning(
            "in-silico enrichment R failed: %s",
            (proc.stderr or proc.stdout)[-1500:],
        )
        return {"status": "skipped", "reason": proc.stderr[-500:]}
    go_csv = data_dir / "insilico_go_enrichment.csv"
    kegg_csv = data_dir / "insilico_kegg_enrichment.csv"
    fig_dir = data_dir.parent / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[str] = []
    if go_csv.exists():
        go_path = fig_dir / "fig_67_ko_go_enrichment.png"
        if _plot_enrichment_bubble(go_csv, go_path, "GO Enrichment", group_col="ONTOLOGY"):
            figures.append(go_path.name)
    if kegg_csv.exists():
        kegg_path = fig_dir / "fig_68_ko_kegg_enrichment.png"
        if _plot_enrichment_bubble(kegg_csv, kegg_path, "KEGG Enrichment"):
            figures.append(kegg_path.name)
    return {"status": "completed", "figures": figures}


def _find_rscript() -> str | None:
    import shutil

    for name in ("Rscript", "Rscript.exe"):
        found = shutil.which(name)
        if found:
            return found
    for base in (
        Path(r"C:\Program Files\R"),
        Path(r"C:\Program Files\Microsoft\R Open"),
        Path.home() / "AppData" / "Local" / "Programs" / "R",
    ):
        if base.exists():
            candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
            if candidates:
                return str(candidates[0])
    return None


def _plot_enrichment_bubble(
    csv_path: Path,
    out_path: Path,
    title: str,
    group_col: str | None = None,
) -> bool:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    df = pd.read_csv(csv_path)
    if df.empty or "Description" not in df.columns:
        return False
    if "pvalue" not in df.columns and "p.adjust" not in df.columns:
        return False
    p_col = "p.adjust" if "p.adjust" in df.columns else "pvalue"
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
    df = df.dropna(subset=[p_col]).sort_values(p_col).head(30)
    if df.empty:
        return False
    if "Count" not in df.columns and "count" in df.columns:
        df = df.rename(columns={"count": "Count"})
    if "Count" not in df.columns:
        df["Count"] = 5
    count = pd.to_numeric(df["Count"], errors="coerce").fillna(5)
    colour = -np.log10(pd.to_numeric(df[p_col], errors="coerce").clip(lower=1e-300))
    if "GeneRatio" in df.columns:
        ratios = df["GeneRatio"].astype(str).str.split("/", expand=True)
        if ratios.shape[1] == 2:
            x = pd.to_numeric(ratios[0], errors="coerce") / pd.to_numeric(
                ratios[1], errors="coerce"
            ) * 100
        else:
            x = colour
        xlabel = "Gene Ratio (%)"
    else:
        x = colour
        xlabel = "-log10 adjusted p"
    fig, ax = plt.subplots(figsize=(9, max(4.5, 0.34 * len(df) + 2)))
    y = np.arange(len(df))
    scatter = ax.scatter(x, y, s=count * 8, c=colour, cmap="RdYlBu_r", edgecolor="black", linewidth=0.3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["Description"].astype(str), fontsize=7.5)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    cb = fig.colorbar(scatter, ax=ax, pad=0.02)
    cb.set_label("-log10 p")
    legend_sizes = sorted(
        {
            int(np.ceil(np.percentile(count, p)))
            for p in (10, 50, 90)
        }
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=max(2.0, 2.0 * math.sqrt(size)),
            markerfacecolor="#d9d9d9",
            markeredgecolor="black",
            markeredgewidth=0.4,
        )
        for size in legend_sizes
    ]
    ax.legend(
        handles,
        [str(size) for size in legend_sizes],
        title="Count",
        loc="lower right",
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        labelspacing=0.7,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    return True


def _write_html_report(
    out_dir: Path,
    fig_dir: Path,
    data_dir: Path,
    ko_gene: str,
    species: str,
    matrix: pd.DataFrame,
    cell_types: pd.Series | None,
    changes: pd.DataFrame,
    enrichment: dict,
    isko: dict,
) -> Path:
    from datetime import datetime

    top = _display_changes(changes).head(15)
    rows = "".join(
        "<tr><td>{gene}</td><td>{wt:.4f}</td><td>{ko:.4f}</td><td>{d:.4f}</td>"
        "<td>{direction}</td></tr>".format(
            gene=row.gene,
            wt=row.wt_mean,
            ko=row.ko_mean,
            d=row.delta,
            direction="up" if row.delta > 0 else "down",
        )
        for row in top.itertuples(index=False)
    )
    note = str(isko.get("gene_note") or "").strip() or _KNOWN_GENE_NOTES.get(
        ko_gene.upper(), ""
    )
    def _fig(name: str, caption: str) -> str:
        if (fig_dir / name).exists():
            return (
                f'<figure><img src="figures/{name}">'
                f"<figcaption>{caption}</figcaption></figure>"
            )
        return f"<p>{caption}：未生成。</p>"

    go_fig = _fig(
        "fig_67_ko_go_enrichment.png",
        "图4：敲除后改变靶基因 GO 富集气泡图",
    )
    kegg_fig = _fig(
        "fig_68_ko_kegg_enrichment.png",
        "图5：敲除后改变靶基因 KEGG 通路富集气泡图",
    )
    html = f"""<!doctype html>
<html lang="zh-CN">
<head><meta charset="utf-8">
<title>{ko_gene} 虚拟敲除分析报告</title>
<style>
body{{font-family:'Microsoft YaHei',Arial,sans-serif;max-width:1100px;margin:32px auto;padding:0 24px;color:#222}}
h1{{border-bottom:3px solid #1f6f8b;padding-bottom:10px}}
h2{{margin-top:32px;color:#1f4e5f}}
.meta{{color:#555;font-size:14px}}
figure img{{max-width:100%;border:1px solid #ddd}}
figcaption{{font-size:12px;color:#666;margin:6px 0 22px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
td,th{{border:1px solid #ccc;padding:5px 7px}}
th{{background:#f1f5f6}}
.alert{{background:#fff7e6;border:1px solid #e7c270;padding:10px 14px;border-radius:4px}}
.note{{background:#eef7f3;padding:10px 14px;border-radius:4px}}
</style></head>
<body>
<h1>{ko_gene} 虚拟敲除分析报告 <small>(In Silico Knockout)</small></h1>
<div class="meta">分析时间：{datetime.now():%Y-%m-%d %H:%M}；物种：{species}；分析方式：单细胞调控网络扰动模拟</div>
<div class="alert">本报告为基于表达共调控网络的预测性结果，不等同于真实敲除实验；下游结论需经湿实验验证。</div>
<h2>1. 数据集来源与分析规模</h2>
<p>本次模拟使用输入表达矩阵中的 <b>{matrix.shape[1]}</b> 个细胞和 <b>{matrix.shape[0]}</b> 个基因。经高变基因筛选后，模型纳入 <b>{len(changes)}</b> 个基因；若提供细胞类型注释，UMAP 向量场按注释着色。</p>
<h2>2. 敲除基因生物学背景</h2>
<div class="note">{note or "未配置该基因的生物学背景说明。可在配置项 insilico_knockout.gene_note 中补充。"}</div>
<h2>3. 细胞命运偏转轨迹预测 (UMAP 向量场分析)</h2>
<p>下图箭头表示基因 {ko_gene} 功能缺失后预测的细胞状态位移。箭头方向指示状态转变方向，长度代表转变动力强弱。</p>
<figure><img src="figures/fig_66_ko_shift_umap.png"><figcaption>图1：{ko_gene} 虚拟敲除后细胞在 UMAP 空间中的命运偏转矢量场。</figcaption></figure>
<h2>4. 基因调控网络作用权重 (TF-Target Network)</h2>
<p>网络以 {ko_gene} 为中心，展示表达调控系数最强的下游靶基因。绿色边代表正调控/激活，红色边代表负调控/抑制。</p>
<figure><img src="figures/fig_65_ko_regulatory_network.png"><figcaption>图2：{ko_gene} 局部转录调控网络拓扑图。</figcaption></figure>
<h2>5. 下游靶基因定量变化 (Expression Changes)</h2>
<p>下图展示受敲除影响最明显的靶基因在野生型 (WT) 与敲除型 (KO) 条件下的平均表达量变化预测。</p>
<figure><img src="figures/fig_64_ko_target_expression_bar.png"><figcaption>图3：野生型与 {ko_gene} 虚拟敲除型细胞中关键靶基因表达水平的对比柱状图。</figcaption></figure>
<h3>靶基因表达定量变化数据表 (Top 15)</h3>
<table><tr><th>靶基因</th><th>野生型表达均值</th><th>敲除型表达均值</th><th>表达变化值</th><th>调控倾向</th></tr>{rows}</table>
<h2>6. GO / KEGG 富集分析</h2>
{go_fig}
{kegg_fig}
</body></html>"""
    path = out_dir / "in_silico_knockout_report.html"
    path.write_text(html, encoding="utf-8")
    return path


def _copy_photo_outputs(cfg: ResolvedConfig, isko: dict, summary: dict) -> None:
    value = isko.get("photo_output_dir") or cfg.get(
        "insilico_knockout", "photo_output_dir"
    )
    if not value:
        return
    target = Path(str(value)).expanduser()
    if not target.is_absolute():
        target = cfg.workdir / target
    target.mkdir(parents=True, exist_ok=True)
    data_target = target / "data"
    fig_target = target / "figures"
    data_target.mkdir(parents=True, exist_ok=True)
    fig_target.mkdir(parents=True, exist_ok=True)
    for name, path in summary.get("data", {}).items():
        source = Path(path)
        if source.exists():
            shutil.copy2(source, data_target / source.name)
    if summary.get("report"):
        source = Path(summary["report"])
        if source.exists():
            shutil.copy2(source, target / source.name)
    for figure in summary.get("figures", []):
        fig_source = cfg.knockout_dir() / "in_silico" / "figures" / figure
        if fig_source.exists():
            shutil.copy2(fig_source, fig_target / fig_source.name)
