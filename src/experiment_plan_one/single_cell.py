"""Single-cell mouse discovery and human validation analyses."""

from __future__ import annotations

import gzip
import json
import logging
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any
from urllib.parse import quote

try:  # Keep bulk/target stages importable without the optional scRNA stack.
    import anndata as ad
    import scanpy as sc
except Exception:  # pragma: no cover - optional stack may be absent or incompatible
    ad = None
    sc = None
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread, mmwrite

from . import cache
from .common import bh_fdr, ensure_dir, get_json, save_figure, slug, write_json

LOG = logging.getLogger("experiment_plan_one.single_cell")
MYGENE_API = "https://mygene.info/v3"


def _require_scanpy() -> None:
    if ad is None or sc is None:
        raise RuntimeError(
            "single-cell stages require anndata and scanpy; install the "
            "single-cell optional dependencies before running mouse/human stages"
        )

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

UNKNOWN_IDS = {"", "na", "n/a", "nan", "none", "unknown", "unassigned"}


def _clean_identifier(value: Any) -> str:
    text = str(value or "").strip()
    return "" if text.lower() in UNKNOWN_IDS else text


def _qc_snapshot(data: ad.AnnData) -> dict[str, Any]:
    return {
        "cells": int(data.n_obs),
        "genes": int(data.n_vars),
        "median_genes_per_cell": (
            float(data.obs["n_genes_by_counts"].median())
            if "n_genes_by_counts" in data.obs
            else None
        ),
        "median_counts_per_cell": (
            float(data.obs["total_counts"].median())
            if "total_counts" in data.obs
            else None
        ),
        "median_mitochondrial_percent": (
            float(data.obs["pct_counts_mt"].median())
            if "pct_counts_mt" in data.obs
            else None
        ),
    }


def _write_qc_report(
    output_dir: Path,
    before: dict[str, Any],
    after: dict[str, Any],
    *,
    doublet_status: str,
) -> None:
    write_json(
        output_dir / "qc_and_replication_report.json",
        {
            "before_qc": before,
            "after_qc": after,
            "doublet_assessment": doublet_status,
            "ambient_rna_assessment": (
                "not_evaluated_without_soupx_or_cellbender; no ambient-RNA "
                "correction is claimed"
            ),
            "statistical_unit_rule": (
                "Cells are never treated as biological replicates for group "
                "inference; use donor/sample pseudobulk."
            ),
        },
    )


def _optional_doublet_assessment(data: ad.AnnData, output_dir: Path) -> str:
    """Attempt Scrublet through Scanpy, but keep failure explicit."""
    if "sample_id" not in data.obs or data.obs["sample_id"].nunique() < 1:
        return "not_evaluated_missing_sample_ids"
    try:
        counts = data.layers.get("counts")
        if counts is None:
            return "not_evaluated_missing_raw_count_layer"
        work = ad.AnnData(
            X=counts.copy(),
            obs=data.obs.copy(),
            var=data.var.copy(),
        )
        sc.pp.scrublet(work, batch_key="sample_id")
        data.obs["doublet_score"] = work.obs["doublet_score"].reindex(
            data.obs_names
        )
        data.obs["predicted_doublet"] = work.obs["predicted_doublet"].reindex(
            data.obs_names
        )
        scores = pd.to_numeric(
            data.obs.get("doublet_score", pd.Series(index=data.obs_names)),
            errors="coerce",
        )
        labels = data.obs.get(
            "predicted_doublet",
            pd.Series(False, index=data.obs_names),
        )
        pd.DataFrame(
            {
                "cell": data.obs_names.astype(str),
                "sample_id": data.obs["sample_id"].astype(str).to_numpy(),
                "doublet_score": scores.to_numpy(),
                "predicted_doublet": np.asarray(labels, dtype=bool),
            }
        ).to_csv(output_dir / "doublet_scores.csv.gz", index=False, compression="gzip")
        return (
            f"completed_scrublet; predicted_doublets="
            f"{int(np.asarray(labels, dtype=bool).sum())}"
        )
    except Exception as exc:  # noqa: BLE001
        LOG.warning("doublet assessment unavailable: %s", exc)
        return f"unavailable_without_verified_doublet_tool: {exc}"


def _choose_inference_unit(data: ad.AnnData) -> tuple[str, dict[str, str]]:
    """Return a donor column and a library-to-unit mapping."""
    for column in ("donor_id", "animal_id", "patient_id"):
        if column not in data.obs:
            continue
        values = data.obs[column].map(_clean_identifier)
        if values.nunique() >= 2 and values.ne("").all():
            return column, {
                str(library): str(unit)
                for library, unit in data.obs.groupby("sample_id", observed=True)[
                    column
                ]
                .first()
                .items()
            }
    return "sample_id", {
        str(sample): str(sample)
        for sample in data.obs["sample_id"].dropna().astype(str).unique()
    }


def map_human_to_mouse_homologs(
    genes: list[str],
    output_path: Path,
    *,
    allow_network: bool = True,
    timeout: int = 60,
) -> dict[str, str]:
    """Map human symbols to one-to-one mouse orthologs with an auditable cache."""
    if output_path.exists():
        cached = pd.read_csv(output_path, sep="\t", dtype=str).fillna("")
        cached_mapping = {
            str(row["human_gene"]): str(row["mouse_gene"])
            for _, row in cached.iterrows()
            if row.get("mapping_status") == "one_to_one"
            and str(row.get("mouse_gene") or "").strip()
        }
        if cached_mapping or not allow_network:
            return cached_mapping
    rows: list[dict[str, str]] = []
    mapping: dict[str, str] = {}
    for gene in sorted(set(str(value).upper() for value in genes if value)):
        row = {
            "human_gene": gene,
            "mouse_gene": "",
            "mouse_gene_id": "",
            "homology_type": "",
            "mapping_status": (
                "network_disabled" if not allow_network else "unmapped"
            ),
            "reason": "",
            "source": "NCBI HomoloGene via MyGene.info",
            "source_version": "",
        }
        if allow_network:
            try:
                payload = get_json(
                    (
                        f"{MYGENE_API}/query?q={quote('symbol:' + gene)}"
                        "&species=human&fields=symbol,homologene&size=10"
                    ),
                    timeout=timeout,
                )
                hits = (
                    payload.get("hits") or []
                    if isinstance(payload, dict)
                    else []
                )
                exact = next(
                    (
                        hit
                        for hit in hits
                        if str(hit.get("symbol") or "").upper() == gene
                    ),
                    None,
                )
                homologene = (
                    exact.get("homologene") if isinstance(exact, dict) else {}
                ) or {}
                gene_groups = homologene.get("genes") or []
                mouse_members = [
                    group
                    for group in gene_groups
                    if len(group) >= 2 and int(group[0]) == 10090
                ]
                human_members = [
                    group
                    for group in gene_groups
                    if len(group) >= 2 and int(group[0]) == 9606
                ]
                if len(mouse_members) == 1 and len(human_members) == 1:
                    mouse_id = str(mouse_members[0][1])
                    lookup = get_json(
                        (
                            f"{MYGENE_API}/gene/{quote(mouse_id)}"
                            "?fields=symbol"
                        ),
                        timeout=timeout,
                    )
                    mouse_symbol = (
                        str(lookup.get("symbol") or "").upper()
                        if isinstance(lookup, dict)
                        else ""
                    )
                    if mouse_symbol:
                        mapping[gene] = mouse_symbol
                        row.update(
                            {
                                "mouse_gene": mouse_symbol,
                                "mouse_gene_id": mouse_id,
                                "homology_type": "homologene_one_to_one_group",
                                "mapping_status": "one_to_one",
                                "source_version": (
                                    f"HomoloGene:{homologene.get('id')}"
                                ),
                            }
                        )
                    else:
                        row["reason"] = "mouse symbol lookup failed"
                else:
                    row["reason"] = (
                        "HomoloGene group is missing or not one-to-one"
                    )
            except Exception as exc:  # noqa: BLE001
                row["reason"] = str(exc)
        rows.append(row)
    ensure_dir(output_path.parent)
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
    return mapping


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
    _require_scanpy()
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
    data.uns["qc_before"] = _qc_snapshot(data)
    sc.pp.filter_cells(data, min_genes=min_genes)
    sc.pp.filter_genes(data, min_cells=3)
    data = data[data.obs["pct_counts_mt"].astype(float) <= max_mito_pct].copy()
    data.uns["qc_after"] = _qc_snapshot(data)
    data.layers["counts"] = data.X.copy()
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    data.raw = data
    sc.pp.highly_variable_genes(data, n_top_genes=min(2500, data.n_vars), batch_key="batch")
    hvg = data[:, data.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=min(40, hvg.n_vars - 1, hvg.n_obs - 1), svd_solver="arpack")
    representation = _integrate_batches(hvg)
    sc.pp.neighbors(hvg, n_neighbors=15, n_pcs=min(30, hvg.obsm[representation].shape[1]), use_rep=representation)
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


def _integrate_batches(data: ad.AnnData) -> str:
    """Integrate an existing PCA representation; never change the count layer."""
    if "X_pca" not in data.obsm:
        raise ValueError("PCA must precede Harmony")
    if data.obs["batch"].nunique() < 2:
        return "X_pca"
    try:
        sc.external.pp.harmony_integrate(data, "batch", max_iter_harmony=30)
        return "X_pca_harmony"
    except Exception as exc:
        LOG.warning("Harmony failed; explicitly retaining uncorrected PCA: %s", exc)
        data.uns["integration_limitation"] = str(exc)
        return "X_pca"


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


def _rscript_path() -> str | None:
    return shutil.which("Rscript") or shutil.which("Rscript.exe")


def _cellchat_r_available() -> bool:
    executable = _rscript_path()
    if not executable:
        return False
    process = subprocess.run(
        [
            executable,
            "--vanilla",
            "-e",
            'cat(requireNamespace("CellChat", quietly=TRUE))',
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    return process.returncode == 0 and "TRUE" in process.stdout.upper()


def _export_cellchat_inputs(
    data: ad.AnnData,
    input_dir: Path,
    *,
    group_column: str,
    cell_type_column: str,
) -> tuple[Path, Path, Path, Path]:
    if "counts" not in data.layers:
        raise RuntimeError(
            "the original R CellChat path requires a raw-count layer named "
            "'counts'; no normalized matrix was substituted"
        )
    observation = data.obs.copy()
    if group_column not in observation or cell_type_column not in observation:
        raise RuntimeError(
            "CellChat metadata must contain condition and annotated cell type"
        )
    keep = (
        observation[group_column].astype(str).str.strip().ne("")
        & observation[cell_type_column].astype(str).str.strip().ne("")
    )
    subset = data[keep.to_numpy()].copy()
    if subset.n_obs < 20:
        raise RuntimeError("too few annotated cells for R CellChat")
    ensure_dir(input_dir)
    matrix_path = input_dir / "counts.mtx"
    genes_path = input_dir / "genes.txt"
    cells_path = input_dir / "cells.txt"
    metadata_path = input_dir / "metadata.tsv"
    counts = subset.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.transpose().tocoo()
    else:
        counts = sparse.coo_matrix(
            np.asarray(counts, dtype=float).transpose()
        )
    mmwrite(matrix_path, counts)
    genes_path.write_text(
        "\n".join(subset.var_names.astype(str)),
        encoding="utf-8",
    )
    cells = subset.obs_names.astype(str)
    cells_path.write_text("\n".join(cells), encoding="utf-8")
    metadata = pd.DataFrame(
        {
            "cell": cells,
            "condition": subset.obs[group_column].astype(str).to_numpy(),
            "cell_type": subset.obs[cell_type_column].astype(str).to_numpy(),
        }
    )
    metadata.to_csv(metadata_path, sep="\t", index=False)
    return matrix_path, genes_path, cells_path, metadata_path


def _run_r_cellchat_analysis(
    data: ad.AnnData,
    output_dir: Path,
    *,
    species: str = "mm",
    group_column: str = "condition",
    cell_type_column: str = "cell_type",
    min_cells: int = 10,
    n_bootstrap: int = 100,
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    """Run the original R CellChat implementation on raw counts."""
    executable = _rscript_path()
    script = Path(__file__).resolve().parent / "R" / "cellchat_analysis.R"
    if not executable or not script.exists():
        raise RuntimeError(
            "R CellChat is configured but Rscript or the bundled R script is "
            "unavailable; Figure 4f is blocked rather than substituted"
        )
    input_dir = ensure_dir(output_dir / "cellchat_r_inputs")
    matrix_path, genes_path, cells_path, metadata_path = (
        _export_cellchat_inputs(
            data,
            input_dir,
            group_column=group_column,
            cell_type_column=cell_type_column,
        )
    )
    process = subprocess.run(
        [
            executable,
            "--vanilla",
            str(script),
            f"--matrix={matrix_path}",
            f"--genes={genes_path}",
            f"--cells={cells_path}",
            f"--metadata={metadata_path}",
            f"--output={output_dir}",
            f"--species={species}",
            f"--group-column={group_column}",
            f"--celltype-column={cell_type_column}",
            f"--min-cells={int(min_cells)}",
            f"--nboot={int(n_bootstrap)}",
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=3600,
    )
    interaction_path = output_dir / "cellchat_r_interactions.csv"
    summary_path = output_dir / "cellchat_r_summary.json"
    if process.returncode != 0 or not summary_path.exists():
        detail = (process.stderr or process.stdout or "").strip()
        raise RuntimeError(
            "R CellChat failed; Figure 4f remains blocked and no substitute "
            f"method was run. {detail[-1200:]}"
        )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    raw = (
        pd.read_csv(interaction_path)
        if interaction_path.exists() and interaction_path.stat().st_size
        else pd.DataFrame()
    )
    if raw.empty:
        return {
            "interactions": pd.DataFrame(),
            "pathways": pd.DataFrame(),
            "status": {
                **summary,
                "status": "valid_negative",
                "reason": "R CellChat returned no communication interactions",
            },
        }
    raw = raw.rename(
        columns={
            "pathway_name": "pathway",
            "prob": "score",
            "pval": "cell_level_pvalue",
        }
    )
    required = {"source", "target", "ligand", "receptor", "group", "score"}
    if not required.issubset(raw.columns):
        raise RuntimeError(
            "R CellChat output is missing required communication columns"
        )
    raw["source"] = raw["source"].astype(str)
    raw["target"] = raw["target"].astype(str)
    raw["ligand"] = raw["ligand"].astype(str)
    raw["receptor"] = raw["receptor"].astype(str)
    raw["group"] = raw["group"].astype(str)
    raw["score"] = pd.to_numeric(raw["score"], errors="coerce")
    keys = ["pathway", "ligand", "receptor", "source", "target"]
    grouped = (
        raw.groupby([*keys, "group"], dropna=False, as_index=False)["score"]
        .max()
    )
    wide = (
        grouped.pivot_table(
            index=keys,
            columns="group",
            values="score",
            aggfunc="max",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for column in ("NCD", "HFD"):
        if column not in wide:
            wide[column] = np.nan
    wide["delta_HFD_NCD"] = wide["HFD"] - wide["NCD"]
    wide["inference_status"] = "cell_level_probability_only"
    wide["fdr"] = np.nan
    wide["significant"] = False
    pathways = (
        wide.groupby("pathway", as_index=False)
        .agg(
            NCD=("NCD", "sum"),
            HFD=("HFD", "sum"),
            delta_HFD_NCD=("delta_HFD_NCD", "sum"),
            n_pairs=("pathway", "size"),
            n_significant=("significant", "sum"),
        )
        .sort_values(["n_significant", "delta_HFD_NCD"], ascending=False)
    )
    return {
        "interactions": wide,
        "pathways": pathways,
        "status": {
            **summary,
            "status": "descriptive_only",
            "reason": (
                "R CellChat probabilities were computed per group; "
                "cell-level probabilities are not treated as independent "
                "biological replicates"
            ),
            "n_interactions": int(len(wide)),
            "n_significant": 0,
            "significant_inference": False,
        },
    }


def run_mouse_single_cell(
    extracted_dir: Path,
    soft_path: Path,
    core_genes: list[str],
    output_dir: Path,
    *,
    sample_key: str | None = None,
    cellchat_backend: str = "explicit_lr",
    cellchat_permutations: int = 100,
    cellchat_seed: int = 123,
    allow_network: bool = True,
) -> dict[str, Any]:
    _require_scanpy()
    ensure_dir(output_dir)
    data_path = output_dir / "mouse_liver_processed.h5ad"
    state_path = output_dir / "single_cell.cache.json"
    input_key = cache.signature(files=[soft_path, Path(__file__), *sorted(extracted_dir.rglob("*.gz"))],
                                parameters={"sample_key": sample_key, "species": "mouse"})
    if cache.valid(state_path, input_key, [data_path]) and not sample_key:
        data = ad.read_h5ad(data_path)
    else:
        data = load_mouse_dataset(extracted_dir, soft_path)
        data.write_h5ad(data_path, compression="gzip")
        cache.save(state_path, input_key, [data_path])
    doublet_status = _optional_doublet_assessment(data, output_dir)
    _write_qc_report(
        output_dir,
        dict(data.uns.get("qc_before") or {}),
        dict(data.uns.get("qc_after") or _qc_snapshot(data)),
        doublet_status=doublet_status,
    )
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
    homolog_path = output_dir / "human_mouse_homolog_mapping.tsv"
    homolog_map = map_human_to_mouse_homologs(
        core_genes,
        homolog_path,
        allow_network=allow_network,
    )
    requested_core = list(core_genes)
    mouse_core = [homolog_map.get(gene.upper(), gene) for gene in core_genes]
    core = _resolve_genes(data.var_names, mouse_core)
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

    write_json(output_dir / "disease_inference_status.json", {
        "status": "descriptive_only",
        "reason": "GSE270583 has one NCD library; cell counts cannot replace biological replicates",
        "legacy_cell_level_DE": "not used for disease inference",
        "jqf_role": "separate intervention arm; never merged into normal control",
    })
    interaction_path = output_dir / "cellchat_like_interactions.csv"
    pathway_path = output_dir / "cellchat_like_pathways.csv"
    backend = str(cellchat_backend or "explicit_lr").strip().lower()
    if backend in {"r_cellchat", "cellchat", "original"}:
        lr = _run_r_cellchat_analysis(
            data,
            output_dir,
            species="mm",
            group_column="condition",
            cell_type_column="cell_type",
            n_bootstrap=cellchat_permutations,
        )
    elif backend in {"explicit_lr", "python", "legacy"}:
        lr = _cellchat_like_analysis(
            data,
            n_permutations=cellchat_permutations,
            seed=cellchat_seed,
        )
    else:
        raise ValueError(
            "cellchat_backend must be r_cellchat or explicit_lr"
        )
    lr["interactions"].to_csv(interaction_path, index=False)
    lr["pathways"].to_csv(pathway_path, index=False)
    communication_status = dict(lr.get("status") or {})
    write_json(
        output_dir / "cellchat_permutation_summary.json",
        {
            **communication_status,
            "n_bootstrap": (
                int(cellchat_permutations)
                if backend in {"r_cellchat", "cellchat", "original"}
                else 0
            ),
            "n_biounit_permutations": (
                0
                if backend in {"r_cellchat", "cellchat", "original"}
                else int(cellchat_permutations)
            ),
            "n_permutations": (
                None
                if backend in {"r_cellchat", "cellchat", "original"}
                else int(cellchat_permutations)
            ),
            "n_interactions": len(lr["interactions"]),
            "n_significant": (
                int(lr["interactions"]["significant"].sum())
                if "significant" in lr["interactions"].columns
                else 0
            ),
            "note": (
                "Original R CellChat was executed on raw counts."
                if backend in {"r_cellchat", "cellchat", "original"}
                else (
                    "This is an explicit ligand-receptor scoring procedure, "
                    "not the R CellChat implementation. Significant group "
                    "inference requires at least two independent biological "
                    "units per group."
                )
            ),
        },
    )
    _plot_cell_communication(
        lr["interactions"],
        lr["pathways"],
        output_dir / "fig4f_cellchat_network.png",
        title=(
            "R CellChat communication network"
            if backend in {"r_cellchat", "cellchat", "original"}
            else "Explicit ligand-receptor scoring network"
        ),
    )
    unit_column, unit_map = _choose_inference_unit(data)
    library_manifest = (
        data.obs.groupby(["sample_id", "condition"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    library_manifest["biological_unit"] = library_manifest["sample_id"].map(
        unit_map
    )
    library_manifest["verified_biological_replicate"] = unit_column in {
        "donor_id",
        "animal_id",
    }
    library_manifest["role"] = np.where(
        library_manifest["condition"].astype(str).eq("JQF"),
        "separate_intervention_arm",
        "descriptive_library",
    )
    library_manifest_path = output_dir / "mouse_library_manifest.tsv"
    library_manifest.to_csv(library_manifest_path, sep="\t", index=False)
    cohort_rows: list[dict[str, Any]] = []
    for (sample_id, condition), subset in data.obs.groupby(
        ["sample_id", "condition"], observed=True
    ):
        unit = unit_map.get(str(sample_id), str(sample_id))
        verified_donor = unit_column in {"donor_id", "animal_id", "patient_id"}
        cohort_rows.append(
            {
                "accession": "GSE270583",
                "sample_id": str(sample_id),
                "library_id": str(sample_id),
                "donor_id": str(unit) if verified_donor else "",
                "species": "Mus musculus",
                "platform": "10x scRNA-seq",
                "tissue": "liver",
                "condition": str(condition),
                "original_diagnosis": str(condition),
                "paired_patient": "no",
                "batch": str(sample_id),
                "exposure_status": "unknown",
                "included": "yes",
                "inclusion_reason": (
                    "descriptive localization only; insufficient animal-level "
                    "replication for disease significance"
                ),
                "data_version": "current GEO RAW",
                "source": str(soft_path),
            }
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
            "requested_human_core_genes": requested_core,
            "homolog_mapping": str(homolog_path),
            "biological_unit_column": unit_column,
            "biological_units": len(set(unit_map.values())),
            "library_manifest": str(library_manifest_path),
            "disease_inference_status": "descriptive_only",
            "cellchat_note": (
                "Ligand-receptor communication was tested with an explicit "
                "local receptor-pair table. Group labels are permuted across "
                "biological units only when at least two units per condition "
                "exist; otherwise the output is descriptive."
            ),
        },
    )
    return {
        "h5ad": data_path,
        "cellchat": lr,
        "core_genes": core,
        "composition": output_dir / "cell_composition.csv",
        "cohort_rows": cohort_rows,
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
    gene_order: list[str] = []
    for marker_group in marker_sets:
        marker_group_genes = [
            gene
            for gene in frame["gene"].drop_duplicates().astype(str)
            if gene in {
                str(value)
                for value in _resolve_genes(data.var_names, marker_sets[marker_group])
            }
        ]
        for gene in marker_group_genes[:3]:
            if gene not in gene_order:
                gene_order.append(gene)
    frame = frame[frame["gene"].isin(gene_order)].copy()
    if frame.empty:
        return
    x = {gene: index for index, gene in enumerate(gene_order)}
    y = {cell: index for index, cell in enumerate(cell_order)}
    ax.scatter(
        [x[gene] for gene in frame["gene"]],
        [y[cell] for cell in frame["cell_type"]],
        s=np.clip(frame["pct_expression"] * 1.5, 4, 150),
        c=frame["mean_expression"],
        cmap="cividis",
        edgecolors="#2f3a45",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=50, ha="right", fontsize=6.5)
    ax.set_yticks(range(len(cell_order)))
    ax.set_yticklabels(cell_order, fontsize=6.5)
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
        sharey=True,
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
        ax.set_xticks([1, 2])
        ax.set_xticklabels(groups)
        ax.set_title(gene)
        ax.text(0.98, 0.98, "Descriptive cells; no replicate inference", transform=ax.transAxes,
                ha="right", va="top", fontsize=6)
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
        ax.set_title(gene, fontsize=8.2, fontweight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        colorbar = fig.colorbar(points, ax=ax, shrink=0.72, pad=0.02)
        colorbar.ax.tick_params(labelsize=6)
    save_figure(fig, output)


def _plot_composition(composition: pd.DataFrame, output: Path) -> None:
    frame = (
        composition[composition["condition"].isin(["NCD", "HFD"])]
        .groupby(["sample_id", "cell_type"], observed=True)["n_cells"]
        .sum()
        .unstack(fill_value=0)
    )
    if frame.empty:
        return
    proportions = frame.div(frame.sum(axis=1), axis=0)
    condition_lookup = (
        composition.drop_duplicates("sample_id")
        .set_index("sample_id")["condition"]
        .astype(str)
    )
    labels = [
        f"{sample}\n{condition_lookup.get(str(sample), '')}"
        for sample in proportions.index.astype(str)
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    bottom = np.zeros(len(proportions))
    colors = plt.get_cmap("tab20", frame.shape[1])
    for index, cell_type in enumerate(frame.columns):
        values = proportions[cell_type].to_numpy()
        ax.bar(labels, values, bottom=bottom, label=cell_type, color=colors(index))
        bottom += values
    ax.set_ylabel("Cell fraction")
    ax.set_title("Cell composition by independent library", fontweight="bold")
    ax.legend(
        bbox_to_anchor=(0.5, -0.22),
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=7,
    )
    save_figure(fig, output)


def _cellchat_like_analysis(
    data: ad.AnnData,
    *,
    n_permutations: int = 100,
    seed: int = 123,
    max_cells: int = 5000,
) -> dict[str, pd.DataFrame]:
    """Score LR pairs and permute condition labels at biological-unit level.

    This is an explicit ligand-receptor scoring procedure, not an imitation of
    R CellChat. If fewer than two independent biological units are available in
    either condition, observed scores are returned as descriptive only and no
    p-values or FDR claims are made.
    """
    normalized = data.raw.to_adata() if data.raw is not None else data
    observation = data.obs.copy()
    conditions = ["NCD", "HFD"]
    donor_column = next(
        (
            column
            for column in ("donor_id", "animal_id")
            if column in observation
        ),
        "",
    )
    verified_biological_units = bool(donor_column)
    if not donor_column:
        donor_column = "sample_id" if "sample_id" in observation else ""
        if not donor_column:
            return {
                "interactions": pd.DataFrame(),
                "pathways": pd.DataFrame(),
                "status": {
                    "status": "not_run",
                    "reason": "no donor/sample identifier available for permutation",
                },
            }
        observation = observation.copy()
        observation["sample_id"] = observation["sample_id"].map(
            _clean_identifier
        )
    work = pd.DataFrame(
        {
            "__position": np.arange(len(observation), dtype=int),
            "cell_type": observation["cell_type"].astype(str).to_numpy(),
            "condition": observation["condition"].astype(str).to_numpy(),
            "donor": observation[donor_column].map(_clean_identifier).to_numpy(),
        },
        index=observation.index,
    )
    work = work[
        work["condition"].isin(conditions)
        & work["donor"].ne("")
    ]
    if work.empty:
        return {
            "interactions": pd.DataFrame(),
            "pathways": pd.DataFrame(),
            "status": {
                "status": "not_run",
                "reason": "no labelled donor/sample observations in NCD/HFD",
            },
        }
    donor_conditions = (
        work.groupby("donor", observed=True)["condition"]
        .agg(lambda values: sorted(set(values.astype(str))))
        .to_dict()
    )
    if any(len(values) != 1 for values in donor_conditions.values()):
        raise ValueError("a biological unit has conflicting NCD/HFD annotations")
    donor_to_condition = {
        str(donor): str(values[0])
        for donor, values in donor_conditions.items()
    }
    donors = sorted(donor_to_condition)
    donor_index = {donor: index for index, donor in enumerate(donors)}
    donor_condition = np.asarray(
        [donor_to_condition[donor] == "HFD" for donor in donors],
        dtype=int,
    )
    n_by_condition = {
        condition: int(np.sum(donor_condition == int(condition == "HFD")))
        for condition in conditions
    }
    sufficient_replication = bool(
        verified_biological_units
        and all(value >= 2 for value in n_by_condition.values())
    )

    rng = np.random.default_rng(seed)
    sampled_indices: list[int] = []
    for _, group in work.groupby(["cell_type", "condition"], observed=True):
        take = max(1, int(round(max_cells / max(1, work.groupby(["cell_type", "condition"]).ngroups))))
        sampled_indices.extend(
            rng.choice(
                group["__position"].to_numpy(),
                size=min(take, len(group)),
                replace=False,
            ).tolist()
        )
    sampled = work.loc[work["__position"].isin(sampled_indices)]
    if len(sampled) > max_cells:
        sampled = sampled.sample(n=max_cells, random_state=seed)
    cell_types = sorted(sampled["cell_type"].unique())
    cell_type_codes = pd.Categorical(
        sampled["cell_type"],
        categories=cell_types,
    ).codes
    donor_codes = sampled["donor"].map(donor_index).to_numpy(dtype=int)
    n_cell_types = len(cell_types)
    expression = normalized.X
    pair_data: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
    for pathway, ligand, receptor in LIGAND_RECEPTOR_DB:
        ligand_gene = _resolve_genes(normalized.var_names, [ligand])
        receptor_gene = _resolve_genes(normalized.var_names, [receptor])
        if not ligand_gene or not receptor_gene:
            continue
        ligand_values = _gene_matrix_column(
            expression,
            normalized.var_names,
            ligand_gene[0],
        )[sampled["__position"].to_numpy()]
        receptor_values = _gene_matrix_column(
            expression,
            normalized.var_names,
            receptor_gene[0],
        )[sampled["__position"].to_numpy()]
        pair_data[(pathway, ligand_gene[0], receptor_gene[0])] = {
            "ligand": np.asarray(ligand_values, dtype=float),
            "receptor": np.asarray(receptor_values, dtype=float),
        }

    def means(
        values: np.ndarray,
        condition_labels: np.ndarray,
    ) -> np.ndarray:
        """Average donor x cell-type means with each donor weighted equally."""
        donor_cell_sums = np.zeros((len(donors), n_cell_types), dtype=float)
        donor_cell_counts = np.zeros((len(donors), n_cell_types), dtype=float)
        np.add.at(
            donor_cell_sums,
            (donor_codes, cell_type_codes),
            values,
        )
        np.add.at(
            donor_cell_counts,
            (donor_codes, cell_type_codes),
            1.0,
        )
        donor_cell_means = np.divide(
            donor_cell_sums,
            donor_cell_counts,
            out=np.full_like(donor_cell_sums, np.nan),
            where=donor_cell_counts > 0,
        )
        output = np.full((2, n_cell_types), np.nan, dtype=float)
        for condition_index in (0, 1):
            mask = condition_labels == condition_index
            if not mask.any():
                continue
            condition_values = donor_cell_means[mask]
            valid = np.isfinite(condition_values)
            counts = valid.sum(axis=0)
            sums = np.nansum(condition_values, axis=0)
            output[condition_index, counts > 0] = (
                sums[counts > 0] / counts[counts > 0]
            )
        return output

    permutation_deltas: dict[tuple[str, str, str], np.ndarray] = {
        key: np.full(
            (max(1, n_permutations), n_cell_types, n_cell_types),
            np.nan,
            dtype=float,
        )
        for key in pair_data
    }
    if sufficient_replication:
        for permutation_index in range(max(1, n_permutations)):
            permuted = rng.permutation(donor_condition)
            for key, values in pair_data.items():
                ligand_means = means(values["ligand"], permuted)
                receptor_means = means(values["receptor"], permuted)
                scores = np.fmin(
                    ligand_means[:, :, None],
                    receptor_means[:, None, :],
                )
                permutation_deltas[key][permutation_index] = (
                    scores[1] - scores[0]
                )

    records: list[dict[str, Any]] = []
    for (pathway, ligand, receptor), values in pair_data.items():
        ligand_means = means(values["ligand"], donor_condition)
        receptor_means = means(values["receptor"], donor_condition)
        scores = np.fmin(
            ligand_means[:, :, None],
            receptor_means[:, None, :],
        )
        delta = scores[1] - scores[0]
        null = permutation_deltas[(pathway, ligand, receptor)]
        for source_index, source in enumerate(cell_types):
            for target_index, target in enumerate(cell_types):
                observed_delta = float(delta[source_index, target_index])
                if not np.isfinite(observed_delta):
                    continue
                null_values = null[:, source_index, target_index]
                null_values = null_values[np.isfinite(null_values)]
                p_value = (
                    float(
                        (
                            1
                            + np.sum(
                                np.abs(null_values)
                                >= abs(observed_delta)
                            )
                        )
                        / (len(null_values) + 1)
                    )
                    if sufficient_replication and len(null_values)
                    else np.nan
                )
                records.append(
                    {
                        "pathway": pathway,
                        "ligand": ligand,
                        "receptor": receptor,
                        "source": source,
                        "target": target,
                        "NCD": (
                            float(scores[0, source_index, target_index])
                            if np.isfinite(scores[0, source_index, target_index])
                            else np.nan
                        ),
                        "HFD": (
                            float(scores[1, source_index, target_index])
                            if np.isfinite(scores[1, source_index, target_index])
                            else np.nan
                        ),
                        "delta_HFD_NCD": observed_delta,
                        "p_value": p_value,
                        "inference_status": (
                            "biological_unit_permutation"
                            if sufficient_replication
                            else "insufficient_biological_replicates"
                        ),
                        "n_NCD_units": n_by_condition["NCD"],
                        "n_HFD_units": n_by_condition["HFD"],
                        "donor_column": donor_column,
                    }
                )
    interactions = pd.DataFrame(records)
    if interactions.empty:
        return {
            "interactions": interactions,
            "pathways": pd.DataFrame(),
            "status": {
                "status": "valid_negative",
                "reason": "no matched ligand-receptor pairs in the local table",
                "donor_column": donor_column,
            },
        }
    finite_p = interactions["p_value"].map(
        lambda value: float(value) if pd.notna(value) else 1.0
    )
    interactions["fdr"] = bh_fdr(finite_p.tolist())
    if not sufficient_replication:
        interactions["fdr"] = np.nan
        interactions["significant"] = False
    else:
        interactions["significant"] = (
            pd.to_numeric(interactions["fdr"], errors="coerce") < 0.05
        )
    interactions["n_permutations"] = int(max(1, n_permutations))
    pathways = (
        interactions.groupby("pathway", as_index=False)
        .agg(
            NCD=("NCD", "sum"),
            HFD=("HFD", "sum"),
            delta_HFD_NCD=("delta_HFD_NCD", "sum"),
            n_pairs=("pathway", "size"),
            min_fdr=("fdr", "min"),
            n_significant=("significant", "sum"),
        )
        .sort_values(["n_significant", "delta_HFD_NCD"], ascending=False)
    )
    status = {
        "status": (
            "completed" if sufficient_replication else "descriptive_only"
        ),
        "method": "explicit ligand-receptor scoring with biological-unit label permutation",
        "statistical_unit": donor_column,
        "verified_biological_units": verified_biological_units,
        "n_NCD_units": n_by_condition["NCD"],
        "n_HFD_units": n_by_condition["HFD"],
        "reason": (
            ""
            if sufficient_replication
            else "fewer than two independent biological units in at least one condition"
            if verified_biological_units
            else "library IDs are not treated as independent animals/donors"
        ),
    }
    return {
        "interactions": interactions,
        "pathways": pathways,
        "status": status,
    }


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
    *,
    title: str = "Cell-cell communication",
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
    ).head(10)
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
                "alpha": 0.38,
                "connectionstyle": "arc3,rad=0.08",
            },
        )
    label_texts = nx.draw_networkx_labels(
        graph,
        position,
        labels={
            node: node.replace("P:", "")[:24]
            for node in graph.nodes
            if not node.startswith("P:")
            or node.replace("P:", "")
            in set(
                str(value).replace("P:", "")
                for value in pathway_nodes[:5]
            )
        },
        font_size=6.8,
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
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=7)
    ax.set_title(title, fontweight="bold")
    ax.axis("off")
    save_figure(fig, output)


def run_human_single_cell(
    extracted_dir: Path,
    soft_path: Path,
    core_genes: list[str],
    output_dir: Path,
    *,
    max_cells_per_sample: int = 0,
    seed: int = 42,
) -> dict[str, Any]:
    """Process GSE202379 raw count CSVs into a disease-spectrum atlas."""
    _require_scanpy()
    ensure_dir(output_dir)
    data_path = output_dir / "human_liver_processed.h5ad"
    state_path = output_dir / "single_cell.cache.json"
    input_key = cache.signature(files=[soft_path, Path(__file__), *sorted(extracted_dir.rglob("*.gz"))],
                                parameters={"max_cells_per_sample": max_cells_per_sample, "seed": seed})
    metadata = parse_geo_soft_samples(soft_path)
    if cache.valid(state_path, input_key, [data_path]):
        data = ad.read_h5ad(data_path)
    else:
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
            if max_cells_per_sample > 0 and frame.shape[1] > max_cells_per_sample:
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
                    "patient_id": _patient_id(record),
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
        data.uns["qc_before"] = _qc_snapshot(data)
        sc.pp.filter_cells(data, min_genes=200)
        sc.pp.filter_genes(data, min_cells=3)
        data = data[data.obs["pct_counts_mt"].astype(float) <= 20.0].copy()
        data.uns["qc_after"] = _qc_snapshot(data)
        sc.pp.normalize_total(data, target_sum=1e4)
        sc.pp.log1p(data)
        data.raw = data
        sc.pp.highly_variable_genes(data, n_top_genes=min(2500, data.n_vars), batch_key="batch")
        hvg = data[:, data.var["highly_variable"]].copy()
        sc.pp.scale(hvg, max_value=10)
        sc.tl.pca(hvg, n_comps=min(40, hvg.n_vars - 1, hvg.n_obs - 1), svd_solver="arpack")
        representation = _integrate_batches(hvg)
        sc.pp.neighbors(hvg, n_neighbors=15, n_pcs=min(30, hvg.obsm[representation].shape[1]), use_rep=representation)
        sc.tl.umap(hvg, random_state=seed)
        sc.tl.leiden(hvg, resolution=0.8, key_added="cluster", flavor="igraph", n_iterations=2)
        data.obsm["X_pca"] = hvg.obsm["X_pca"]
        data.obsm["X_umap"] = hvg.obsm["X_umap"]
        data.obs["cluster"] = hvg.obs["cluster"].astype(str)
        data.obs["cell_type"] = annotate_clusters(data, HUMAN_MARKERS)
        data.write_h5ad(data_path, compression="gzip")
        cache.save(state_path, input_key, [data_path])
    doublet_status = _optional_doublet_assessment(data, output_dir)
    _write_qc_report(
        output_dir,
        dict(data.uns.get("qc_before") or {}),
        dict(data.uns.get("qc_after") or _qc_snapshot(data)),
        doublet_status=doublet_status,
    )
    data.obs["cell_type"] = (
        data.obs["cell_type"]
        .astype(str)
        .str.replace(r"\s+\d+$", "", regex=True)
    )
    manifest_path, cohort_rows = _write_human_manifest(
        data,
        output_dir,
        soft_path,
        metadata,
    )
    included_samples = set(
        pd.read_csv(manifest_path, sep="\t")
        .loc[lambda frame: frame["included"].eq("yes"), "sample_id"]
        .astype(str)
    )
    unknown_or_ambiguous_libraries = len(set(data.obs["sample_id"].astype(str)) - included_samples)

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
            "donors": int(
                data.obs["patient_id"]
                .map(_clean_identifier)
                .loc[lambda values: values.ne("")]
                .nunique()
            ),
            "unknown_or_ambiguous_libraries": unknown_or_ambiguous_libraries,
            "conditions": data.obs["condition"].value_counts().to_dict(),
            "cell_types": data.obs["cell_type"].value_counts().to_dict(),
            "core_genes": core,
            "sampling_note": (
                "No per-library cell cap is applied by default; all retained "
                "nuclei are used unless a positive max_cells_per_sample is configured."
            ),
            "statistical_unit": "donor",
            "manifest": str(manifest_path),
        },
    )
    return {
        "h5ad": data_path,
        "core_genes": core,
        "expression": output_dir / "human_core_gene_expression_long.csv.gz",
        "cohort_rows": cohort_rows,
    }


def _patient_id(record: dict[str, Any]) -> str:
    for field in ("patient_id", "patient", "donor_id", "donor"):
        value = str(record.get(field) or "").strip()
        if value and value.lower() not in {"nan", "unknown", "na"}:
            return value
    # Dataset-specific documented title prefix, not GSM-as-patient substitution.
    match = re.match(r"^(P(?:CL|HL)?\d+)-", str(record.get("title") or ""), re.IGNORECASE)
    return match.group(1).upper() if match else ""


def _donor_expression(frame: pd.DataFrame, donor: str) -> pd.DataFrame:
    if frame.empty or donor not in frame:
        return frame.iloc[:0].copy()
    ids = frame[donor].fillna("").astype(str).str.strip()
    clean = frame.assign(**{donor: ids})
    clean = clean.loc[~ids.str.lower().isin(UNKNOWN_IDS)].copy()
    if clean.empty:
        return clean
    if clean.groupby(donor)["condition"].nunique().gt(1).any():
        raise ValueError("a biological donor has conflicting conditions")
    return clean.groupby(
        [donor, "gene", "cell_type", "condition"],
        observed=True,
        as_index=False,
    )["expression"].mean()


def _write_human_manifest(
    data: ad.AnnData,
    output_dir: Path,
    soft_path: Path,
    metadata: dict[str, dict[str, Any]],
) -> tuple[Path, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    for sample_id, subset in data.obs.groupby("sample_id", observed=True):
        patient_ids = subset["patient_id"].map(_clean_identifier)
        patients = sorted(set(patient_ids.loc[patient_ids.ne("")].astype(str)))
        conditions = sorted(set(subset["condition"].astype(str)))
        included = len(patients) == 1 and len(conditions) == 1
        reason = (
            "verified donor-level metadata"
            if included
            else "unknown_conflicting_or_ambiguous_patient_id"
        )
        rows.append(
            {
                "accession": "GSE202379",
                "sample_id": str(sample_id),
                "library_id": str(sample_id),
                "donor_id": patients[0] if len(patients) == 1 else "",
                "species": "Homo sapiens",
                "platform": "snRNA-seq",
                "tissue": "liver",
                "condition": conditions[0] if len(conditions) == 1 else "conflicting",
                "original_diagnosis": "|".join(conditions),
                "paired_patient": "yes" if included else "unknown",
                "batch": str(sample_id),
                "exposure_status": "unknown",
                "included": "yes" if included else "no",
                "inclusion_reason": reason,
                "data_version": "current GEO SOFT/RAW",
                "source": str(soft_path),
                "n_cells": len(subset),
            }
        )
    manifest = pd.DataFrame(rows)
    manifest_path = output_dir / "human_donor_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    correction_rows = []
    for gsm in ("GSM6112262", "GSM6112263"):
        record = metadata.get(gsm, {})
        correction_rows.append(
            {
                "gsm": gsm,
                "title": str(record.get("title") or ""),
                "patient_id": _patient_id(record),
                "disease_status": str(record.get("disease_status") or ""),
                "official_revision_date": "2025-06-06",
                "status": "current_SOFT_parsed_without_manual_swap",
            }
        )
    pd.DataFrame(correction_rows).to_csv(
        output_dir / "gse202379_metadata_correction_audit.tsv",
        sep="\t",
        index=False,
    )
    return manifest_path, rows


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
    expression = _donor_expression(expression, "patient_id")
    if expression.empty:
        return
    expression["cell"] = expression["patient_id"]
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
                (stats.kruskal(*nonempty).pvalue if np.ptp(np.concatenate(nonempty)) > 0 else 1.0)
                if len(nonempty) >= 2 and all(len(v) >= 2 for v in values.values() if len(v)) else np.nan
            )
            for condition in disease_conditions:
                values_group = values[condition]
                rows.append(
                    {
                        "gene": gene,
                        "cell_type": cell_type,
                        "condition": condition,
                        "n_healthy": len(healthy),
                        "n_condition": len(values_group),
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
    colorbar = fig.colorbar(
        image,
        ax=ax,
        location="bottom",
        shrink=0.82,
        pad=0.13,
    )
    colorbar.set_label("Median donor-mean expression difference vs Healthy")
    ax.text(
        1.0,
        -0.34,
        "*FDR<0.05  **FDR<0.01  ***FDR<0.001 (donor-mean Kruskal-Wallis; exploratory)",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=5.2,
        color="#4a5560",
    )
    save_figure(fig, output)


def _mouse_core_statistics(expression: pd.DataFrame) -> pd.DataFrame:
    expression = _donor_expression(expression, "sample_id")
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
                "inference_level": "sample_mean_exploratory" if min(len(control), len(disease)) >= 2 else "insufficient_biological_replicates",
                "n_NCD": len(control),
                "n_HFD": len(disease),
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
    expression = _donor_expression(expression, "patient_id")
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
            p_value = stats.kruskal(*groups).pvalue if np.ptp(np.concatenate(groups)) > 0 else 1.0
        else:
            p_value = np.nan
        rows.append(
            {
                "gene": gene,
                "cell_type": cell_type,
                "n_donors": len(subset),
                "inference_level": "donor_mean_exploratory",
                "n_conditions": int(subset["condition"].nunique()),
                "p_value": p_value,
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_fdr(result["p_value"])
    return result.sort_values(["fdr", "gene", "cell_type"], na_position="last")
