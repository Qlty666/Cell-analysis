#!/usr/bin/env python3
"""Virtual knockout and multi-dimensional target prioritization.

The module combines tumor-vs-normal expression shift, proliferation
co-expression, co-expression network hubness and optional DepMap CRISPR
dependency into a single 0-1 knockout priority score. When additional
evidence is provided it also scores disease-signature reversal, pathway
control, cell-type specificity, clinical prognosis and druggability, then
merges everything into a target score and classifies candidates into core
drivers, microenvironment regulators or biomarkers.

It is a screening heuristic for prioritising genes before wet-lab knockout,
not a mechanistic simulation of the knockout phenotype.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import ResolvedConfig
from .network_toxicology import ppi_hub_scores
from .provenance import write_run_manifest
from .utils import DockingError, write_json

PROLIFERATION_SIGNATURE = [
    "MKI67",
    "PCNA",
    "TOP2A",
    "AURKA",
    "CDC20",
    "CDK1",
    "CCNB1",
    "BUB1",
    "BIRC5",
    "CENPA",
]

_META_SAMPLE_COLS = {
    "sample",
    "sample_id",
    "cell",
    "barcode",
    "cell_id",
    "cell_barcode",
}
_META_GROUP_COLS = {
    "condition",
    "group",
    "group_label",
    "label",
    "type",
    "status",
    "disease",
}
_DEPMAP_GENE_COLS = {"gene", "symbol", "hgnc", "hugo_symbol"}
_DEPMAP_EFFECT_COLS = {"gene_effect", "effect", "dependency", "score", "depmap"}
_DEPMAP_MODEL_COLS = {"modelid", "model_id", "model"}

_DEFAULT_WEIGHTS = {
    "expression": 0.35,
    "proliferation": 0.25,
    "network": 0.20,
    "depmap": 0.20,
}

DEFAULT_PATHWAY_GENES = {
    "cell_cycle": [
        "CDK1",
        "CDK2",
        "CDK4",
        "CCNA2",
        "CCNB1",
        "CCND1",
        "CCNE1",
        "CDC20",
        "CDC25A",
        "E2F1",
    ],
    "apoptosis": [
        "BCL2",
        "BAX",
        "BCL2L1",
        "MCL1",
        "CASP3",
        "CASP9",
        "BIRC5",
        "FAS",
        "TP53",
        "MDM2",
    ],
    "emt": [
        "VIM",
        "SNAI1",
        "SNAI2",
        "TWIST1",
        "ZEB1",
        "ZEB2",
        "CDH1",
        "CDH2",
        "FN1",
        "MMP2",
    ],
    "p53": ["TP53", "MDM2", "CDKN1A", "BAX", "BBC3", "PMAIP1", "GADD45A"],
    "pi3k_akt": [
        "PIK3CA",
        "PIK3R1",
        "AKT1",
        "AKT2",
        "MTOR",
        "PTEN",
        "RPS6KB1",
        "EIF4EBP1",
    ],
}

DEFAULT_TARGET_WEIGHTS = {
    "base": 0.35,
    "reversal": 0.20,
    "pathway": 0.15,
    "specificity": 0.10,
    "prognosis": 0.10,
    "druggability": 0.10,
    "ppi_hub": 0.10,
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
}


def run_knockout(cfg: ResolvedConfig, log) -> dict:
    """Score virtual knockout candidates and write ranking outputs."""
    ko = cfg.data.get("knockout", {})
    expression_csv = _resolve_path(cfg, ko.get("expression_csv"))
    if expression_csv is None or not expression_csv.exists():
        raise DockingError(
            "knockout expression CSV not found: "
            f"{expression_csv}; pass --expression-csv or set "
            "knockout.expression_csv"
        )

    matrix, long_meta = _load_expression(expression_csv, ko)
    if len(matrix) < 5:
        raise DockingError(
            f"expression matrix has only {len(matrix)} genes; need at least 5"
        )
    if matrix.shape[1] < 2:
        raise DockingError(
            f"expression matrix has only {matrix.shape[1]} sample column"
        )

    metadata_csv = _resolve_path(cfg, ko.get("metadata_csv"))
    metadata = _load_metadata(metadata_csv, ko.get("group_column") or "condition")
    if metadata is None and long_meta is not None:
        metadata = long_meta
    case_cols, normal_cols, group_names = _assign_groups(matrix, metadata, ko)
    if case_cols and normal_cols:
        log.info(
            "knockout groups: case=%s (%s samples), normal=%s (%s samples)",
            group_names["case"],
            len(case_cols),
            group_names["normal"],
            len(normal_cols),
        )
    else:
        log.info(
            "no case/normal groups detected; using expression level, "
            "proliferation and network scores only"
        )

    mat = np.log2(matrix + 1.0)
    expression = _expression_stat(mat, case_cols, normal_cols)
    expression_score = _percentile_rank(expression)
    prolif_corr, prolif_score = _proliferation_scores(mat)
    hub_degree, hub_score = _network_hub_scores(mat, ko, log)
    ppi_hub_frame, ppi_hub_score = _ppi_hub_scores(cfg, mat.index, ko, log)
    depmap_effect, depmap_score, depmap_lines = _depmap_scores(
        cfg, mat.index, ko, log
    )

    scores = {
        "expression": expression_score,
        "proliferation": prolif_score,
        "network": hub_score,
        "depmap": depmap_score,
    }
    weights = _normalize_weights(ko.get("weights"), scores)
    knockout_score = pd.Series(0.0, index=mat.index)
    for name, weight in weights.items():
        knockout_score = knockout_score + weight * scores[name]

    reversal_corr, reversal_score = _reversal_scores(
        mat, expression, case_cols, normal_cols, ko, log
    )
    pathway_control, pathway_score = _pathway_scores(mat, ko, log)
    cell_types = _load_cell_types(
        metadata if metadata is not None else metadata_csv,
        ko,
    )
    cell_type_specificity, specificity_score = _specificity_scores(
        mat, cell_types, log
    )
    prognosis_hr, prognosis_score = _prognosis_scores(
        cfg, mat.index, expression, ko, log
    )
    druggable_hits, druggability_score = _druggability_scores(
        cfg, mat.index, ko, log
    )
    off_target_paralogs, safety_concern = _off_target_flags(
        cfg, mat.index, ko, log
    )

    target_sources = {
        "base": knockout_score,
        "reversal": reversal_score,
        "pathway": pathway_score,
        "specificity": specificity_score,
        "prognosis": prognosis_score,
        "druggability": druggability_score,
        "ppi_hub": ppi_hub_score,
    }
    target_weights = _normalize_target_weights(
        ko.get("target_weights"), target_sources
    )
    target_score = pd.Series(0.0, index=mat.index)
    for name, weight in target_weights.items():
        target_score = target_score + weight * target_sources[name]
    if safety_concern is not None:
        penalty = float(ko.get("off_target_penalty", 0.05))
        target_score = target_score - penalty * safety_concern.fillna(0.0)
    target_score = target_score.clip(0.0, 1.0)

    frame = pd.DataFrame({"gene": mat.index}, index=mat.index)
    if case_cols and normal_cols:
        frame["case_mean"] = mat[case_cols].mean(axis=1)
        frame["normal_mean"] = mat[normal_cols].mean(axis=1)
        frame["log2fc"] = expression
    else:
        frame["case_mean"] = np.nan
        frame["normal_mean"] = np.nan
        frame["log2fc"] = np.nan
        frame["mean_expression"] = expression
    frame["expression_score"] = expression_score
    frame["proliferation_correlation"] = prolif_corr
    frame["proliferation_score"] = prolif_score
    frame["hub_degree"] = hub_degree
    frame["hub_score"] = hub_score
    if ppi_hub_frame is not None:
        frame["ppi_degree"] = _as_series(ppi_hub_frame.get("ppi_degree"))
        frame["ppi_betweenness"] = _as_series(
            ppi_hub_frame.get("ppi_betweenness")
        )
        frame["ppi_clustering"] = _as_series(
            ppi_hub_frame.get("ppi_clustering")
        )
        frame["ppi_hub_score"] = _as_series(ppi_hub_score)
    frame["depmap_effect"] = depmap_effect
    frame["depmap_score"] = depmap_score
    frame["knockout_score"] = knockout_score
    frame["reversal_correlation"] = _as_series(reversal_corr)
    frame["reversal_score"] = _as_series(reversal_score)
    frame["pathway_control"] = _as_series(pathway_control)
    frame["pathway_score"] = _as_series(pathway_score)
    frame["cell_type_specificity"] = _as_series(cell_type_specificity)
    frame["specificity_score"] = _as_series(specificity_score)
    frame["prognosis_hr"] = _as_series(prognosis_hr)
    frame["prognosis_score"] = _as_series(prognosis_score)
    frame["druggable_hits"] = _as_series(druggable_hits)
    frame["druggability_score"] = _as_series(druggability_score)
    frame["off_target_paralogs"] = _as_series(off_target_paralogs)
    frame["safety_concern"] = _as_series(safety_concern)
    frame["target_score"] = target_score
    frame["target_priority"] = np.where(
        target_score >= 0.75,
        "high",
        np.where(target_score >= 0.5, "medium", "low"),
    )
    frame["target_class"] = _classify_targets(frame)
    frame["priority"] = np.where(
        knockout_score >= 0.75,
        "high",
        np.where(knockout_score >= 0.5, "medium", "low"),
    )
    frame = frame.sort_values("target_score", ascending=False).reset_index(
        drop=True
    )
    frame.insert(0, "rank", np.arange(1, len(frame) + 1))

    top_n = int(ko.get("top_n", 50))
    ko_dir = cfg.knockout_dir()
    data_dir = ko_dir / "data"
    figures_dir = ko_dir / "figures"
    ko_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = data_dir / "fig_52_53_ranked_knockout.csv"
    frame.to_csv(ranked_path, index=False)

    candidate_cols = [
        "rank",
        "gene",
        "target_class",
        "target_score",
        "knockout_score",
        "reversal_score",
        "pathway_score",
        "specificity_score",
        "prognosis_score",
        "druggability_score",
        "ppi_hub_score",
        "off_target_paralogs",
        "safety_concern",
    ]
    candidate_cols = [c for c in candidate_cols if c in frame.columns]
    candidates = frame[candidate_cols].head(top_n)
    candidates_path = data_dir / "fig_52_target_candidates.csv"
    candidates.to_csv(candidates_path, index=False)

    figures: list[str] = []
    if ko.get("figures", True):
        try:
            figures = _make_figures(frame, figures_dir, top_n)
        except Exception as exc:
            log.warning("knockout figures failed: %s", exc)

    summary = {
        "genes_scored": int(len(frame)),
        "case_samples": len(case_cols),
        "normal_samples": len(normal_cols),
        "case_group": group_names["case"],
        "normal_group": group_names["normal"],
        "depmap_lines": depmap_lines,
        "depmap_included": depmap_score is not None,
        "ppi_hub_included": ppi_hub_score is not None,
        "top_n": min(top_n, len(frame)),
        "weights": weights,
        "top_genes": frame["gene"].astype(str).head(top_n).tolist(),
        "high_priority": int((frame["priority"] == "high").sum()),
        "output_csv": str(ranked_path),
        "target_candidates_csv": str(candidates_path),
        "target_weights": target_weights,
        "target_class_counts": frame["target_class"].value_counts().to_dict(),
        "multidimensional_scoring": any(
            score is not None
            for score in [
                reversal_score,
                pathway_score,
                specificity_score,
                prognosis_score,
                druggability_score,
            ]
        ),
        "figures": figures,
    }
    write_json(ko_dir / "summary.json", summary)
    manifest_path = write_run_manifest(
        ko_dir,
        cfg,
        "virtual-knockout",
        {
            "expression_csv": _resolve_path(cfg, ko.get("expression_csv")),
            "metadata_csv": _resolve_path(cfg, ko.get("metadata_csv")),
            "depmap_csv": _resolve_path(cfg, ko.get("depmap_csv")),
            "prognosis_csv": _resolve_path(cfg, ko.get("prognosis_csv")),
            "druggability_csv": _resolve_path(cfg, ko.get("druggability_csv")),
            "off_target_csv": _resolve_path(cfg, ko.get("off_target_csv")),
            "ppi_network_csv": _resolve_path(cfg, ko.get("ppi_network_csv")),
        },
        {
            "weights": weights,
            "target_weights": target_weights,
            "top_n": top_n,
            "off_target_penalty": ko.get("off_target_penalty", 0.05),
        },
    )
    summary["manifest"] = str(manifest_path)
    _write_markdown_report(data_dir, frame, summary, top_n)
    target_report = _write_target_report(data_dir, frame, summary, top_n)
    summary["target_report"] = str(target_report)
    log.info(
        "virtual knockout complete: %s genes scored, output %s",
        len(frame),
        ko_dir,
    )
    return summary


def _resolve_path(cfg: ResolvedConfig, value) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = cfg.workdir / path
    return path.resolve()


def _load_expression(path: Path, ko: dict) -> tuple[pd.DataFrame, pd.DataFrame | None]:
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
        long_meta = None
        group_col = ko.get("group_column") or "condition"
        if group_col in df.columns:
            meta_cols = [sample_col, group_col]
            extra_cols = [
                c
                for c in df.columns
                if c not in {gene_col, sample_col, value_col, group_col}
            ]
            meta_cols += extra_cols
            long_meta = (
                df[meta_cols]
                .dropna(subset=[sample_col, group_col])
                .drop_duplicates()
                .rename(columns={sample_col: "sample", group_col: "group"})
            )
    else:
        first_col = df.columns[0]
        matrix = df.set_index(first_col).apply(pd.to_numeric, errors="coerce")
        long_meta = None
    matrix = matrix.dropna(how="all").fillna(0.0)
    matrix.columns = matrix.columns.astype(str)
    matrix.index = matrix.index.astype(str)
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    return matrix, long_meta


def _load_metadata(path: Path | None, group_column: str) -> pd.DataFrame | None:
    if path is None:
        return None
    if not path.exists():
        raise DockingError(f"metadata CSV not found: {path}")
    meta = pd.read_csv(path)
    if meta.empty:
        return None
    sample_col = next(
        (c for c in meta.columns if c.lower() in _META_SAMPLE_COLS),
        meta.columns[0],
    )
    group_col = next(
        (
            c
            for c in meta.columns
            if c == group_column or c.lower() in _META_GROUP_COLS
        ),
        None,
    )
    if group_col is None:
        return None
    cell_col = next(
        (c for c in meta.columns if c.lower() in _CELL_TYPE_COLS),
        None,
    )
    keep = [sample_col, group_col]
    if cell_col is not None and cell_col not in keep:
        keep.append(cell_col)
    out = meta[keep].dropna(subset=[sample_col, group_col]).drop_duplicates()
    return out.rename(columns={sample_col: "sample", group_col: "group"})


def _load_cell_types(source, ko: dict) -> dict[str, str] | None:
    """Map sample names to cell-type labels from metadata or a CSV path."""
    if source is None:
        return None
    if isinstance(source, pd.DataFrame):
        df = source
    else:
        path = Path(str(source))
        if not path.exists():
            return None
        df = pd.read_csv(path)
    if df.empty:
        return None
    sample_col = next(
        (c for c in df.columns if c.lower() in _META_SAMPLE_COLS),
        None,
    )
    if sample_col is None:
        return None
    configured = ko.get("cell_type_column")
    cell_col = (
        configured
        if configured
        else next(
            (c for c in df.columns if c.lower() in _CELL_TYPE_COLS),
            None,
        )
    )
    if cell_col is None or cell_col not in df.columns:
        return None
    pairs = df[[sample_col, cell_col]].dropna()
    return {
        str(sample): str(label)
        for sample, label in zip(pairs[sample_col], pairs[cell_col])
    }


def _as_series(value):
    return np.nan if value is None else value


def _reversal_scores(
    mat: pd.DataFrame,
    expression: pd.Series,
    case_cols: list[str],
    normal_cols: list[str],
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Score how strongly a gene's shift aligns with the disease module."""
    if not case_cols or not normal_cols:
        return None, None
    up = [
        g
        for g in (ko.get("disease_up_genes") or PROLIFERATION_SIGNATURE)
        if g in mat.index
    ]
    down = [
        g for g in (ko.get("disease_down_genes") or []) if g in mat.index
    ]
    if not up:
        return None, None
    down_module = mat.loc[down].mean(axis=0) if down else None
    corr = pd.Series(np.nan, index=mat.index)
    for gene in mat.index:
        up_set = [g for g in up if g != gene]
        if not up_set:
            continue
        module = mat.loc[up_set].mean(axis=0)
        if down_module is not None:
            module = module - down_module
        corr.loc[gene] = mat.loc[gene].corr(module)
    corr = corr.fillna(0.0)
    signal = expression * corr
    return corr, _percentile_rank(signal)


def _pathway_scores(
    mat: pd.DataFrame,
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Mean absolute correlation with curated disease pathway gene sets."""
    pathways = ko.get("pathway_genes") or DEFAULT_PATHWAY_GENES
    if not isinstance(pathways, dict) or not pathways:
        return None, None
    genes: list[str] = []
    for items in pathways.values():
        if isinstance(items, (list, tuple)):
            genes.extend(str(g) for g in items)
    genes = [g for g in dict.fromkeys(genes) if g in mat.index]
    if not genes or mat.shape[1] < 4:
        return None, None
    max_genes = int(ko.get("max_genes", 2000))
    variances = mat.var(axis=1)
    top = mat.loc[
        variances.nlargest(min(max_genes, len(variances))).index
    ]
    pathway_top = [g for g in genes if g in top.index]
    if not pathway_top:
        return None, None
    corr = top.T.corr().reindex(
        index=top.index, columns=pathway_top
    ).fillna(0.0)
    control = corr.abs().mean(axis=1).reindex(mat.index).fillna(0.0)
    return control, _percentile_rank(control)


def _specificity_scores(
    mat: pd.DataFrame,
    cell_types: dict[str, str] | None,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Entropy-based expression specificity across annotated cell types."""
    if not cell_types:
        return None, None
    types = sorted(
        {
            label
            for sample in mat.columns
            if (label := cell_types.get(str(sample))) is not None
        }
    )
    if len(types) < 2:
        return None, None
    type_means = pd.DataFrame(
        {
            label: mat[
                [
                    c
                    for c in mat.columns
                    if cell_types.get(str(c)) == label
                ]
            ].mean(axis=1)
            for label in types
        }
    )
    total = type_means.sum(axis=1).replace(0, np.nan)
    proportions = type_means.div(total, axis=0).clip(lower=1e-12)
    log2_prop = np.log2(proportions)
    entropy = -(proportions * log2_prop).sum(axis=1)
    specificity = (1 - entropy / np.log2(len(types))).clip(0.0, 1.0)
    specificity = specificity.fillna(0.5)
    return specificity, _percentile_rank(specificity)


def _prognosis_scores(
    cfg: ResolvedConfig,
    genes: pd.Index,
    expression: pd.Series,
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Direction-aware hazard-ratio score from an optional prognosis CSV."""
    path = _resolve_path(cfg, ko.get("prognosis_csv"))
    if path is None or not path.exists():
        return None, None
    df = pd.read_csv(path)
    if df.empty:
        return None, None
    gene_col = next(
        (c for c in df.columns if c.lower() in _DEPMAP_GENE_COLS),
        df.columns[0],
    )
    hr_col = next(
        (
            c
            for c in df.columns
            if c.lower()
            in {
                "hr",
                "hazard_ratio",
                "hazardratio",
                "cox_hr",
                "log2hr",
            }
        ),
        None,
    )
    if hr_col is None:
        log.warning("prognosis CSV has no hazard ratio column: %s", path)
        return None, None
    hr = (
        df[[gene_col, hr_col]]
        .dropna()
        .set_index(gene_col)[hr_col]
        .astype(float)
    )
    hr = hr[hr > 0]
    if hr.empty:
        return None, None
    sign = np.sign(expression).reindex(genes).fillna(0.0)
    effect = sign * np.log(hr.reindex(genes).fillna(1.0))
    return hr.reindex(genes), _percentile_rank(effect)


def _druggability_scores(
    cfg: ResolvedConfig,
    genes: pd.Index,
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Combine known-ligand, structure and bioactivity evidence into a score."""
    path = _resolve_path(cfg, ko.get("druggability_csv"))
    if path is None or not path.exists():
        return None, None
    df = pd.read_csv(path)
    if df.empty:
        return None, None
    gene_col = next(
        (c for c in df.columns if c.lower() in _DEPMAP_GENE_COLS),
        df.columns[0],
    )
    count_cols = []
    aliases = {
        "known_ligands": {"known_ligands", "ligands", "ligand_count"},
        "pdb_structures": {
            "pdb_structures",
            "structures",
            "pdb_count",
            "structures_count",
        },
        "chembl_bioactivities": {
            "chembl_bioactivities",
            "bioactivities",
            "assays",
            "bioactivity_count",
        },
    }
    for canonical, names in aliases.items():
        col = next((c for c in df.columns if c.lower() in names), None)
        if col is not None:
            count_cols.append((canonical, col))
    if not count_cols:
        log.warning("druggability CSV has no count columns: %s", path)
        return None, None
    counts = pd.DataFrame(index=df[gene_col].astype(str))
    for canonical, col in count_cols:
        counts[canonical] = pd.to_numeric(
            df[col], errors="coerce"
        ).fillna(0.0).to_numpy()
    hits = counts.sum(axis=1).apply(lambda x: np.log1p(max(0.0, x)))
    hits = hits.groupby(level=0).max()
    hits = hits.reindex(genes.astype(str)).fillna(0.0)
    return hits, _percentile_rank(hits)


def _off_target_flags(
    cfg: ResolvedConfig,
    genes: pd.Index,
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None]:
    """Paralog count and safety-flag columns from an optional CSV."""
    path = _resolve_path(cfg, ko.get("off_target_csv"))
    if path is None or not path.exists():
        return None, None
    df = pd.read_csv(path)
    if df.empty:
        return None, None
    gene_col = next(
        (c for c in df.columns if c.lower() in _DEPMAP_GENE_COLS),
        df.columns[0],
    )
    paralog_col = next(
        (
            c
            for c in df.columns
            if c.lower()
            in {"off_target_paralogs", "paralogs", "paralog_count"}
        ),
        None,
    )
    concern_col = next(
        (
            c
            for c in df.columns
            if c.lower()
            in {"safety_concern", "safety", "off_target_concern", "warning"}
        ),
        None,
    )
    paralogs = (
        pd.to_numeric(df[paralog_col], errors="coerce").fillna(0.0)
        if paralog_col
        else pd.Series(0.0, index=df.index)
    )
    paralogs = (
        pd.Series(paralogs.values, index=df[gene_col])
        .groupby(level=0)
        .max()
        .reindex(genes)
        .fillna(0.0)
    )
    if concern_col:
        raw = df[concern_col].astype(str).str.strip().str.lower()
        concern = raw.map(
            {
                "1": 1,
                "true": 1,
                "yes": 1,
                "y": 1,
                "0": 0,
                "false": 0,
                "no": 0,
                "n": 0,
            }
        )
        concern = concern.fillna(
            pd.to_numeric(df[concern_col], errors="coerce").fillna(0)
        ).astype(float)
        concern = (
            pd.Series(concern.values, index=df[gene_col])
            .groupby(level=0)
            .max()
            .reindex(genes)
            .fillna(0.0)
        )
    else:
        concern = pd.Series(0.0, index=genes)
    return paralogs, concern


def _assign_groups(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame | None,
    ko: dict,
) -> tuple[list[str], list[str], dict[str, str | None]]:
    case_label = str(ko.get("case_label") or "").strip()
    normal_label = str(ko.get("normal_label") or "").strip()
    sample_names = [str(c) for c in matrix.columns]

    if metadata is not None and "group" in metadata.columns:
        mapping = dict(
            zip(metadata["sample"].astype(str), metadata["group"].astype(str))
        )
        case_group = case_label or None
        normal_group = normal_label or None
        if not case_group or not normal_group:
            values = sorted(
                {
                    v
                    for v in mapping.values()
                    if str(v).strip().lower() not in ("", "na", "nan")
                }
            )
            if len(values) >= 2:
                normal_group = normal_group or values[0]
                case_group = case_group or values[1]
        if case_group and normal_group:
            case_cols = [
                s for s in sample_names if mapping.get(s) == case_group
            ]
            normal_cols = [
                s for s in sample_names if mapping.get(s) == normal_group
            ]
            if case_cols and normal_cols:
                return case_cols, normal_cols, {
                    "case": case_group,
                    "normal": normal_group,
                }

    if case_label and normal_label:
        case_cols = [s for s in sample_names if case_label.lower() in s.lower()]
        normal_cols = [
            s for s in sample_names if normal_label.lower() in s.lower()
        ]
        if case_cols and normal_cols:
            return case_cols, normal_cols, {
                "case": case_label,
                "normal": normal_label,
            }

    return sample_names, [], {"case": None, "normal": None}


def _expression_stat(
    mat: pd.DataFrame,
    case_cols: list[str],
    normal_cols: list[str],
) -> pd.Series:
    if case_cols and normal_cols:
        return mat[case_cols].mean(axis=1) - mat[normal_cols].mean(axis=1)
    return mat.mean(axis=1)


def _percentile_rank(values: pd.Series | np.ndarray) -> pd.Series:
    series = pd.Series(values)
    if series.notna().sum() == 0:
        return series.fillna(0.5)
    return series.rank(pct=True).fillna(0.5)


def _proliferation_scores(mat: pd.DataFrame):
    signature = [g for g in PROLIFERATION_SIGNATURE if g in mat.index]
    if not signature or mat.shape[1] < 3:
        return pd.Series(np.nan, index=mat.index), None
    signature_expr = mat.loc[signature].mean(axis=0)
    varying = mat.index[mat.var(axis=1) > 0]
    corr = pd.Series(np.nan, index=mat.index)
    if len(varying):
        part = mat.loc[varying].T.corrwith(signature_expr)
        corr.loc[part.index] = part
    corr = corr.replace([np.inf, -np.inf], np.nan)
    return corr, _percentile_rank(corr.fillna(0.0))


def _network_hub_scores(
    mat: pd.DataFrame,
    ko: dict,
    log,
) -> tuple[pd.Series, pd.Series]:
    empty_degree = pd.Series(0.0, index=mat.index)
    empty_score = pd.Series(0.5, index=mat.index)
    if mat.shape[1] < 4:
        return empty_degree, empty_score

    max_samples = int(ko.get("max_samples", 5000))
    if mat.shape[1] > max_samples:
        step = mat.shape[1] // max_samples
        mat = mat.iloc[:, ::step].iloc[:, :max_samples]
    max_genes = int(ko.get("max_genes", 2000))
    variances = mat.var(axis=1)
    top = mat.loc[variances.nlargest(min(max_genes, len(variances))).index]
    corr = top.T.corr()
    corr = corr.replace([np.inf, -np.inf], 0.0).fillna(0.0)
    np.fill_diagonal(corr.values, 0.0)

    cutoff = float(ko.get("corr_cutoff", 0.7))
    degree = (corr.abs() >= cutoff).sum(axis=1)
    if degree.max() <= 0:
        degree = corr.abs().mean(axis=1)
        log.debug("no co-expression edges above %s; using mean correlation", cutoff)
    full_degree = pd.Series(0.0, index=mat.index)
    full_degree.loc[degree.index] = degree
    return full_degree, _percentile_rank(full_degree)


def _ppi_hub_scores(
    cfg: ResolvedConfig,
    genes: pd.Index,
    ko: dict,
    log,
) -> tuple[pd.DataFrame | None, pd.Series | None]:
    """Score PPI hubness from an optional STRING-style edge table."""
    path = _resolve_path(cfg, ko.get("ppi_network_csv"))
    if path is None or not path.exists():
        return None, None
    try:
        hub = ppi_hub_scores(path, genes=genes)
    except Exception as exc:
        log.warning("PPI hub scoring skipped: %s", exc)
        return None, None
    lookup = {str(g).upper(): str(g) for g in genes}
    hub["gene"] = hub["gene"].map(lookup).fillna(hub["gene"])
    hub = (
        hub.set_index("gene")
        .reindex([str(g) for g in genes])
        .reset_index()
    )
    if "index" in hub.columns and "gene" not in hub.columns:
        hub = hub.rename(columns={"index": "gene"})
    hub = hub.set_index("gene")
    score = hub["ppi_hub_score"]
    log.info(
        "PPI hub scoring included: %s/%s genes matched",
        int(score.notna().sum()),
        len(genes),
    )
    return hub, score


def _depmap_scores(
    cfg: ResolvedConfig,
    genes: pd.Index,
    ko: dict,
    log,
) -> tuple[pd.Series | None, pd.Series | None, int]:
    path = _resolve_path(cfg, ko.get("depmap_csv"))
    if path is None:
        return None, None, 0
    if not path.exists():
        log.warning("DepMap CSV not found, skipping dependency score: %s", path)
        return None, None, 0

    df = pd.read_csv(path)
    if df.empty:
        log.warning("DepMap CSV is empty: %s", path)
        return None, None, 0
    lower = {str(c).lower(): c for c in df.columns}
    lineage_filter = str(ko.get("liver_lineage", "liver")).strip().lower()

    if "modelid" in lower:
        model_col = lower["modelid"]
        lineage_col = lower.get("lineage")
        if lineage_col and lineage_filter:
            mask = (
                df[lineage_col]
                .astype(str)
                .str.lower()
                .str.contains(lineage_filter, na=False)
            )
            if mask.any():
                df = df[mask]
        gene_cols = [
            c
            for c in df.columns
            if c != model_col
            and c != lineage_col
            and pd.api.types.is_numeric_dtype(df[c])
        ]
        if not gene_cols:
            raise DockingError(
                "DepMap wide CSV has no numeric gene columns after ModelID"
            )
        effects = df[gene_cols].mean(axis=0)
        n_lines = len(df)
    else:
        gene_col = next(
            (lower[c] for c in _DEPMAP_GENE_COLS if c in lower), None
        )
        effect_col = next(
            (lower[c] for c in _DEPMAP_EFFECT_COLS if c in lower), None
        )
        if gene_col is None or effect_col is None:
            raise DockingError(
                "DepMap CSV must be wide (ModelID + gene columns) or long "
                "(gene/effect columns)"
            )
        lineage_col = lower.get("lineage")
        if lineage_col and lineage_filter:
            mask = (
                df[lineage_col]
                .astype(str)
                .str.lower()
                .str.contains(lineage_filter, na=False)
            )
            if mask.any():
                df = df[mask]
        effects = (
            df.groupby(gene_col)[effect_col].mean().astype(float)
        )
        n_lines = int(df[gene_col].nunique())

    common = [g for g in genes if g in effects.index]
    if not common:
        log.warning("no DepMap genes matched the expression matrix")
        return None, None, n_lines
    effect = effects.reindex(genes)
    score = _percentile_rank(-effect)
    return effect, score, n_lines


def _normalize_weights(
    weights,
    scores: dict[str, pd.Series | None],
) -> dict[str, float]:
    merged = dict(_DEFAULT_WEIGHTS)
    if isinstance(weights, dict):
        for key, value in weights.items():
            if value is not None:
                merged[key] = float(value)
    available = [
        key for key in merged if scores.get(key) is not None
    ]
    total = sum(merged[key] for key in available)
    if total <= 0:
        return {key: 1.0 / len(available) for key in available}
    return {key: merged[key] / total for key in available}


def _normalize_target_weights(
    weights,
    scores: dict[str, pd.Series | None],
) -> dict[str, float]:
    merged = dict(DEFAULT_TARGET_WEIGHTS)
    if isinstance(weights, dict):
        for key, value in weights.items():
            if value is not None:
                merged[key] = float(value)
    available = [key for key in merged if scores.get(key) is not None]
    total = sum(merged[key] for key in available)
    if total <= 0:
        return {key: 1.0 / len(available) for key in available}
    return {key: merged[key] / total for key in available}


def _classify_targets(frame: pd.DataFrame) -> pd.Series:
    """Classify candidates into driver, microenvironment or biomarker roles."""

    def classify(row) -> str:
        base = row["knockout_score"]
        dep = row.get("depmap_score")
        reversal = row.get("reversal_score")
        if base >= 0.75 and (
            (pd.notna(dep) and dep >= 0.6)
            or (pd.notna(reversal) and reversal >= 0.6)
        ):
            return "core_driver"
        if base >= 0.75:
            return "high_priority"
        pathway = row.get("pathway_score")
        network = row.get("hub_score")
        specificity = row.get("specificity_score")
        if (
            pd.notna(pathway)
            and pathway >= 0.6
        ) or (
            pd.notna(network)
            and pd.notna(specificity)
            and network >= 0.6
            and specificity >= 0.6
        ):
            return "microenvironment_regulator"
        if pd.notna(row.get("expression_score")) and row["expression_score"] >= 0.6:
            return "biomarker"
        return "low_priority"

    return frame.apply(classify, axis=1)


def _make_figures(frame: pd.DataFrame, fig_dir: Path, top_n: int) -> list[str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir.mkdir(parents=True, exist_ok=True)
    figures: list[str] = []
    top = frame.head(top_n).iloc[::-1]
    if not top.empty:
        fig, ax = plt.subplots(figsize=(8, max(3.0, len(top) * 0.28)))
        ax.barh(
            top["gene"].astype(str),
            top["knockout_score"],
            color="#1d6f42",
        )
        ax.set_xlim(0, 1)
        ax.set_xlabel("Knockout priority score")
        ax.set_title(f"Top {len(top)} virtual knockout candidates")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_52_knockout_top_candidates.png", dpi=150)
        plt.close(fig)
        figures.append("fig_52_knockout_top_candidates.png")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(
        frame["knockout_score"],
        bins=40,
        color="#4c7bb8",
        edgecolor="white",
    )
    ax.set_xlabel("Knockout priority score")
    ax.set_ylabel("Gene count")
    ax.set_title("Virtual knockout score distribution")
    fig.tight_layout()
    fig.savefig(fig_dir / "fig_53_knockout_score_distribution.png", dpi=150)
    plt.close(fig)
    figures.append("fig_53_knockout_score_distribution.png")
    return figures


def _write_target_report(
    ko_dir: Path,
    frame: pd.DataFrame,
    summary: dict,
    top_n: int,
) -> Path:
    """Write the multi-dimensional target prioritization report."""
    lines = [
        "# Target Prioritization Report",
        "",
        "The target score combines the base knockout priority with optional "
        "evidence dimensions when they are available. Missing dimensions are "
        "excluded and the remaining weights are renormalized to 0-1.",
        "",
        "| Dimension | Meaning |",
        "| --- | --- |",
        "| knockout_score | expression shift, proliferation, network hub, DepMap |",
        "| reversal_score | expression shift aligned with the disease module |",
        "| pathway_score | mean correlation with disease pathway gene sets |",
        "| specificity_score | entropy-based cell-type expression specificity |",
        "| prognosis_score | direction-aware hazard ratio from clinical cohorts |",
        "| druggability_score | known ligands, structures and bioactivity counts |",
        "| safety_concern | off-target/paralog flag that lowers the target score |",
        "",
        "## Candidate classes",
        "",
    ]
    counts = summary.get("target_class_counts") or {}
    for name in [
        "core_driver",
        "high_priority",
        "microenvironment_regulator",
        "biomarker",
        "low_priority",
    ]:
        lines.append(f"- {name}: {counts.get(name, 0)}")
    lines += [
        "",
        "## Top candidates",
        "",
        "| rank | gene | class | target | knockout | reversal | pathway | "
        "specificity | prognosis | druggability | safety |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in frame.head(top_n).iterrows():
        lines.append(
            "| {rank} | {gene} | {cls} | {target:.3f} | {ko:.3f} | {rev:.3f} | "
            "{path:.3f} | {spec:.3f} | {prog:.3f} | {drug:.3f} | {safe:.0f} |".format(
                rank=row["rank"],
                gene=row["gene"],
                cls=row["target_class"],
                target=row["target_score"],
                ko=row["knockout_score"],
                rev=_na_or(row, "reversal_score"),
                path=_na_or(row, "pathway_score"),
                spec=_na_or(row, "specificity_score"),
                prog=_na_or(row, "prognosis_score"),
                drug=_na_or(row, "druggability_score"),
                safe=_na_or(row, "safety_concern"),
            )
        )
    lines += [
        "",
        "## Handoff",
        "",
        "Run `python scripts\\run_docking.py export-validation` to generate the staged "
        "wet-lab validation plan for these candidates.",
        "",
    ]
    path = ko_dir / "target_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _na_or(row, key: str) -> float:
    value = row.get(key)
    return float(value) if pd.notna(value) else 0.0


def _write_markdown_report(
    ko_dir: Path,
    frame: pd.DataFrame,
    summary: dict,
    top_n: int,
) -> None:
    lines = [
        "# 虚拟敲除优先级分析报告",
        "",
        "> 本报告综合表达差异、增殖共表达、共表达网络 hub 程度和可选 DepMap "
        "CRISPR 依赖数据对候选基因排序，用于实验前筛选；分数是统计优先级，"
        "不等同于真实敲除表型。",
        "",
        "## 分析概况",
        "",
        f"- 评分基因数：{summary['genes_scored']}",
        f"- 肿瘤/病例样本：{summary['case_samples']}",
        f"- 正常样本：{summary['normal_samples']}",
        f"- DepMap 细胞系：{summary['depmap_lines']}",
        f"- 高优先级基因：{summary['high_priority']}",
        "",
        "## Top 候选基因",
        "",
        "| 排名 | 基因 | 敲除分 | 表达分 | 增殖分 | 网络分 | DepMap 分 | 优先级 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    top = frame.head(top_n)
    for _, row in top.iterrows():
        lines.append(
            "| {rank} | {gene} | {score:.3f} | {expr:.3f} | {prolif:.3f} | "
            "{net:.3f} | {dep:.3f} | {pri} |".format(
                rank=row["rank"],
                gene=row["gene"],
                score=row["knockout_score"],
                expr=row["expression_score"],
                prolif=(
                    row["proliferation_score"]
                    if pd.notna(row["proliferation_score"])
                    else 0.0
                ),
                net=row["hub_score"],
                dep=(
                    row["depmap_score"]
                    if pd.notna(row["depmap_score"])
                    else 0.0
                ),
                pri=row["priority"],
            )
        )
    lines += ["", "## 输出文件", "", "- `fig_52_53_ranked_knockout.csv`：全部基因评分表"]
    if summary["figures"]:
        for name in summary["figures"]:
            lines.append(f"- `{name}`：结果图")
    lines.append("- `summary.json`：分析汇总")
    (ko_dir / "knockout_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
