#!/usr/bin/env python3
"""Compound-disease network toxicology analysis.

The module implements the shared front half of the network toxicology
workflow reviewed from the six reference articles:

1. load compound targets and disease genes;
2. compute the compound-disease overlap;
3. score PPI hubness from an optional STRING-style edge table;
4. export compound-target-pathway-disease style nodes/edges and a Venn plot.

It is intentionally lightweight and only requires pandas/matplotlib for the
base outputs; networkx is used when available for betweenness/clustering.
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .cytoscape_network import export_network_to_cytoscape, write_xgmml_network
from .html_utils import esc
from .utils import DockingError, write_json

_GENE_COLS = (
    "gene",
    "symbol",
    "hgnc",
    "gene_symbol",
    "hugo",
    "hugo_symbol",
    "gene_name",
    "geneid",
    "entrez",
    "entrezgene",
    "target",
    "node2",
    "gene2",
    "string_id2",
    "to",
    "target_node",
    "source",
    "node1",
    "gene1",
    "string_id1",
    "from",
    "source_node",
    "protein1",
    "protein2",
)
_SOURCE_COLS = {
    "source",
    "database",
    "data_source",
    "target_source",
    "origin",
}
_EDGE_COL1_ALIASES = (
    "protein1",
    "node1",
    "source",
    "gene1",
    "string_id1",
    "from",
    "source_node",
)
_EDGE_COL2_ALIASES = (
    "protein2",
    "node2",
    "target",
    "gene2",
    "string_id2",
    "to",
    "target_node",
)
_CONFIDENCE_ALIASES = (
    "combined_score",
    "combinedscore",
    "string_combined_score",
    "combined score",
    "score",
    "confidence",
    "confidence_score",
    "evidence_score",
    "experimental_score",
    "database_score",
)


def _confidence_column(edges: pd.DataFrame) -> str | None:
    lower = {str(column).strip().lower(): column for column in edges.columns}
    return next(
        (lower[alias] for alias in _CONFIDENCE_ALIASES if alias in lower),
        None,
    )


def _select_ppi_edges(
    edges: pd.DataFrame,
    overlap_genes: set[str],
    max_ppi_edges: int | None = None,
) -> tuple[pd.DataFrame, str | None]:
    """Keep deduplicated STRING edges connecting overlap genes."""
    if edges is None or edges.empty:
        return edges.copy() if edges is not None else pd.DataFrame(), None
    col1, col2 = _edge_columns(edges)
    if col1 is None or col2 is None:
        col1, col2 = edges.columns[0], edges.columns[1]
    a = edges[col1].astype(str).str.strip().str.upper()
    b = edges[col2].astype(str).str.strip().str.upper()
    mask = (
        a.isin(overlap_genes)
        & b.isin(overlap_genes)
        & (a != "")
        & (b != "")
        & (a != b)
    )
    selected = edges.loc[mask].copy()
    if selected.empty:
        return selected, _confidence_column(selected)
    a_selected = a.loc[selected.index].to_numpy()
    b_selected = b.loc[selected.index].to_numpy()
    left = np.where(a_selected <= b_selected, a_selected, b_selected)
    right = np.where(a_selected <= b_selected, b_selected, a_selected)
    selected["_edge_key"] = left + "|" + right
    score_col = _confidence_column(selected)
    if score_col is not None:
        selected["_score"] = pd.to_numeric(
            selected[score_col],
            errors="coerce",
        )
        selected = selected.sort_values(
            ["_score", "_edge_key"],
            ascending=[False, True],
        )
    else:
        selected = selected.sort_values("_edge_key")
    selected = selected.drop_duplicates(subset="_edge_key", keep="first")
    if max_ppi_edges is not None and len(selected) > int(max_ppi_edges):
        selected = selected.head(int(max_ppi_edges))
    return (
        selected.drop(columns=["_edge_key", "_score"], errors="ignore").reset_index(
            drop=True
        ),
        score_col,
    )


def _edge_columns(edges: pd.DataFrame) -> tuple[str | None, str | None]:
    """Return the two node columns in an edge table when recognizable."""
    lower = {str(c).lower(): c for c in edges.columns}
    col1 = next(
        (lower[c] for c in _EDGE_COL1_ALIASES if c in lower),
        None,
    )
    col2 = next(
        (lower[c] for c in _EDGE_COL2_ALIASES if c in lower),
        None,
    )
    return col1, col2


def _pick_gene_column(df: pd.DataFrame) -> str | None:
    lower = {str(c).strip().lower(): c for c in df.columns}
    return next((lower[c] for c in _GENE_COLS if c in lower), None)


def read_gene_list(
    path: str | Path,
    gene_column: str | None = None,
) -> set[str]:
    """Read a CSV/TSV gene list and return normalized uppercase symbols."""
    p = Path(path)
    if not p.exists():
        raise DockingError(f"gene list not found: {p}")
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(p, sep=sep, dtype=str)
    if df.empty:
        return set()
    column = gene_column
    if column is None:
        column = _pick_gene_column(df) or df.columns[0]
    if column not in df.columns:
        raise DockingError(
            f"gene column '{column}' not found in {p.name}: {list(df.columns)}"
        )
    values = df[column].astype(str).str.strip()
    return set(values[values != ""].str.upper())


def read_target_table(
    path: str | Path,
    gene_column: str | None = None,
    source_name: str | None = None,
) -> pd.DataFrame:
    """Read compound targets and keep gene/source columns."""
    p = Path(path)
    if not p.exists():
        raise DockingError(f"compound target file not found: {p}")
    suffix = p.suffix.lower()
    if suffix == ".json":
        payload = json.loads(p.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = (
                payload.get("targets")
                or payload.get("records")
                or payload.get("data")
                or []
            )
        if not isinstance(payload, list):
            raise DockingError(f"JSON target file must contain a list: {p}")
        df = pd.DataFrame(payload)
    elif suffix in (".jsonl", ".ndjson"):
        rows = [
            json.loads(line)
            for line in p.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        df = pd.DataFrame(rows)
    else:
        sep = "\t" if suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(p, sep=sep, dtype=str)
    if df.empty:
        raise DockingError(f"compound target file is empty: {p}")
    column = gene_column
    if column is None:
        column = _pick_gene_column(df) or df.columns[0]
    if column not in df.columns:
        raise DockingError(
            f"gene column '{column}' not found in {p.name}: {list(df.columns)}"
        )
    source_col = next(
        (c for c in df.columns if c.strip().lower() in _SOURCE_COLS),
        None,
    )
    source = source_name or p.stem
    out = pd.DataFrame(
        {
            "gene": df[column].astype(str).str.strip().str.upper(),
            "source": (
                df[source_col].astype(str).str.strip()
                if source_col is not None
                else source
            ),
        }
    )
    return out[out["gene"] != ""].drop_duplicates()


def load_target_sources(
    cfg,
) -> dict[str, pd.DataFrame]:
    """Load one or more compound-target files from the config."""
    section = cfg.data.get("network_toxicology", {}) or {}
    mapping = section.get("target_sources") or {}
    single = section.get("compound_targets_csv")
    if isinstance(single, (list, tuple)):
        mapping.update({str(Path(str(p)).stem): str(p) for p in single})
    elif single:
        mapping = {"compound": str(single)}
    if section.get("target_sources_dir"):
        mapping = dict(mapping)
        mapping["auto_discovered"] = str(section["target_sources_dir"])

    if not mapping:
        raise DockingError(
            "network_toxicology.target_sources or "
            "network_toxicology.compound_targets_csv is required"
        )
    loaded: dict[str, pd.DataFrame] = {}
    for source, path in mapping.items():
        if path:
            p = Path(str(path)).expanduser()
            if not p.is_absolute():
                p = cfg.workdir / p
            if p.is_dir():
                candidates = []
                for pattern in ("*.csv", "*.tsv", "*.txt", "*.json", "*.jsonl"):
                    candidates.extend(sorted(p.glob(pattern)))
                if not candidates:
                    raise DockingError(f"no target tables found in directory: {p}")
                frames = [
                    read_target_table(
                        candidate,
                        source_name=candidate.stem,
                    )
                    for candidate in candidates
                ]
                loaded[str(source)] = pd.concat(frames, ignore_index=True)
            else:
                loaded[str(source)] = read_target_table(
                    p,
                    source_name=str(source),
                )
    if not loaded:
        raise DockingError("no compound target files could be loaded")
    return loaded


def load_disease_genes(cfg) -> set[str]:
    """Load disease genes from config or the full-pipeline key-gene table."""
    section = cfg.data.get("network_toxicology", {}) or {}
    path = section.get("disease_genes_csv")
    if not path:
        candidate = cfg.workdir / "outputs" / "integration" / "key_genes.csv"
        if candidate.exists():
            path = str(candidate)
    if not path:
        raise DockingError(
            "network_toxicology.disease_genes_csv is required "
            "(or run the full pipeline key-targets stage first)"
        )
    p = Path(str(path)).expanduser()
    if not p.is_absolute():
        p = cfg.workdir / p
    return read_gene_list(p, section.get("disease_gene_column"))


def overlap_analysis(
    target_sources: dict[str, pd.DataFrame],
    disease_genes: set[str],
) -> pd.DataFrame:
    """Return genes shared by compound targets and the disease gene set."""
    per_source = {
        name: set(frame["gene"]) for name, frame in target_sources.items()
    }
    all_targets: set[str] = set()
    for genes in per_source.values():
        all_targets |= genes
    rows = []
    for gene in sorted(all_targets & disease_genes):
        sources = sorted(name for name, genes in per_source.items() if gene in genes)
        rows.append(
            {
                "gene": gene,
                "n_sources": len(sources),
                "sources": ";".join(sources),
            }
        )
    frame = pd.DataFrame(
        rows,
        columns=["gene", "n_sources", "sources"],
    )
    frame["in_compound"] = True
    frame["in_disease"] = True
    return frame.sort_values(
        ["n_sources", "gene"],
        ascending=[False, True],
    ).reset_index(drop=True)


def ppi_hub_scores(
    edge_path: str | Path,
    genes: list[str] | pd.Index | None = None,
) -> pd.DataFrame:
    """Compute degree/betweenness/clustering from a STRING-style edge table.

    The edge table should contain two node columns such as
    ``protein1/protein2``, ``node1/node2`` or ``source/target``. networkx is
    used for betweenness/clustering; a pandas degree fallback is used when
    networkx is not installed.
    """
    p = Path(edge_path)
    if not p.exists():
        raise DockingError(f"PPI edge table not found: {p}")
    sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
    edges = pd.read_csv(p, sep=sep, dtype=str)
    if edges.empty:
        raise DockingError(f"PPI edge table is empty: {p}")
    col1, col2 = _edge_columns(edges)
    if col1 is None or col2 is None:
        raise DockingError(
            "PPI edge table must contain two node columns, e.g. "
            "protein1/protein2, node1/node2 or source/target"
        )
    pair = edges[[col1, col2]].astype(str).apply(
        lambda s: s.str.strip().str.upper()
    )
    pair = pair[(pair[col1] != "") & (pair[col2] != "")]
    pair = pair[pair[col1] != pair[col2]].drop_duplicates()

    degree = (
        pd.concat([pair[col1], pair[col2]])
        .value_counts()
        .rename("ppi_degree")
    )
    result = degree.to_frame()
    result["gene"] = result.index

    try:
        import networkx as nx

        graph = nx.Graph()
        graph.add_edges_from(
            zip(pair[col1].astype(str), pair[col2].astype(str))
        )

        def add_metric(name: str, compute) -> None:
            try:
                values = compute()
                result[name] = result["gene"].map(
                    pd.Series(values, dtype=float)
                )
            except Exception:
                result[name] = np.nan

        add_metric(
            "ppi_betweenness",
            lambda: nx.betweenness_centrality(graph),
        )
        add_metric(
            "ppi_closeness",
            lambda: nx.closeness_centrality(graph),
        )
        add_metric(
            "ppi_eigenvector",
            lambda: nx.eigenvector_centrality_numpy(graph),
        )
        add_metric(
            "ppi_pagerank",
            lambda: nx.pagerank(graph),
        )
        add_metric(
            "ppi_clustering",
            lambda: nx.clustering(graph),
        )

        # Maximum-clique centrality is expensive for large networks. Compute
        # it only when the graph remains tractable; otherwise retain NaN and
        # let the consensus score use the available topology metrics.
        if graph.number_of_nodes() <= 300:
            mcc = {node: 0.0 for node in graph.nodes}
            try:
                for clique in nx.find_cliques(graph):
                    weight = math.factorial(max(0, len(clique) - 1))
                    for node in clique:
                        mcc[node] += weight
            except Exception:
                mcc = {node: float("nan") for node in graph.nodes}
            add_metric("ppi_mcc", lambda: mcc)
        else:
            result["ppi_mcc"] = np.nan
    except ImportError:
        for name in (
            "ppi_betweenness",
            "ppi_closeness",
            "ppi_eigenvector",
            "ppi_pagerank",
            "ppi_mcc",
            "ppi_clustering",
        ):
            result[name] = np.nan

    if genes is not None:
        wanted = pd.Index(genes).astype(str).str.upper()
        result = result.set_index("gene").reindex(wanted).reset_index()
        if "index" in result.columns and "gene" not in result.columns:
            result = result.rename(columns={"index": "gene"})
    result = result.fillna(
        {
            "ppi_degree": 0,
            "ppi_betweenness": 0.0,
            "ppi_clustering": 0.0,
        }
    )
    metric_columns = [
        column
        for column in (
            "ppi_degree",
            "ppi_betweenness",
            "ppi_closeness",
            "ppi_eigenvector",
            "ppi_pagerank",
            "ppi_mcc",
            "ppi_clustering",
        )
        if column in result.columns
    ]
    rank_frame = pd.DataFrame(index=result.index)
    for column in metric_columns:
        values = pd.to_numeric(result[column], errors="coerce")
        rank_frame[column] = values.rank(pct=True)
    result["ppi_hub_score"] = rank_frame.mean(axis=1).fillna(0.5).clip(0.0, 1.0)
    return result.reset_index(drop=True)


def make_venn_figure(
    sets: dict[str, set[str]],
    out_path: str | Path,
) -> str:
    """Draw a two- or three-set Venn diagram with matplotlib only."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = list(sets.keys())[:3]
    values = [sets[name] for name in labels]
    if len(values) < 2:
        raise DockingError("Venn diagram needs at least two gene sets")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect("equal")
    ax.axis("off")
    colors = ["#4C72B0", "#DD8452", "#55A868"]
    if len(values) == 2:
        radius = 1.1
        centers = [(-0.55, 0.0), (0.55, 0.0)]
        regions = _venn_regions(values)
        _label_region(ax, regions.get(("set0",), 0), -1.15, 0.0)
        _label_region(ax, regions.get(("set0", "set1"), 0), 0.0, 0.0)
        _label_region(ax, regions.get(("set1",), 0), 1.15, 0.0)
        for (x, y), color, name in zip(centers, colors, labels):
            circle = plt.Circle(
                (x, y),
                radius,
                fill=False,
                edgecolor=color,
                linewidth=2.5,
            )
            ax.add_patch(circle)
            ax.text(x, y + 1.35, name, ha="center", fontsize=12, fontweight="bold")
    else:
        radius = 1.25
        angle = 90.0
        centers = []
        for i in range(3):
            rad = math.radians(angle)
            centers.append((0.95 * math.cos(rad), 0.95 * math.sin(rad)))
            angle += 120.0
        regions = _venn_regions(values)
        for key, count in regions.items():
            if not count:
                continue
            names = set(key)
            if names == {"set0"}:
                x, y = -1.75, 0.9
            elif names == {"set1"}:
                x, y = 1.75, 0.9
            elif names == {"set2"}:
                x, y = 0.0, -1.55
            elif names == {"set0", "set1"}:
                x, y = 0.0, 1.45
            elif names == {"set0", "set2"}:
                x, y = -0.95, -0.15
            elif names == {"set1", "set2"}:
                x, y = 0.95, -0.15
            else:
                x, y = 0.0, 0.0
            _label_region(ax, count, x, y)
        for (x, y), color, name in zip(centers, colors, labels):
            circle = plt.Circle(
                (x, y),
                radius,
                fill=False,
                edgecolor=color,
                linewidth=2.5,
            )
            ax.add_patch(circle)
            ax.text(x * 1.55, y * 1.55, name, ha="center", fontsize=11, fontweight="bold")

    fig.tight_layout()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return str(path)


def _venn_regions(sets: list[set[str]]) -> dict[tuple[str, ...], int]:
    n = len(sets)
    regions: dict[tuple[str, ...], int] = {}
    for mask in range(1, 1 << n):
        members = {
            f"set{i}" for i in range(n) if mask & (1 << i)
        }
        present = set.intersection(
            *(sets[i] for i in range(n) if mask & (1 << i))
        )
        excluded: set[str] = set()
        for i in range(n):
            if not (mask & (1 << i)):
                excluded |= sets[i]
        regions[tuple(sorted(members))] = len(present - excluded)
    return regions


def _label_region(ax, value, x, y) -> None:
    ax.text(
        x,
        y,
        str(int(value)),
        ha="center",
        va="center",
        fontsize=13,
    )


def write_ctpd_network(
    overlap: pd.DataFrame,
    target_sources: dict[str, pd.DataFrame],
    ppi_edges: pd.DataFrame | None,
    compound_name: str,
    disease_name: str,
    out_dir: Path,
    max_ppi_edges: int | None = None,
) -> dict[str, Path]:
    """Write compound-target-disease nodes/edges as CSV, XGMML and HTML."""
    out_dir.mkdir(parents=True, exist_ok=True)
    metric_columns = [
        column
        for column in (
            "ppi_degree",
            "ppi_betweenness",
            "ppi_closeness",
            "ppi_eigenvector",
            "ppi_pagerank",
            "ppi_mcc",
            "ppi_clustering",
            "ppi_hub_score",
            "n_sources",
        )
        if column in overlap.columns
    ]
    zero_metrics = {column: 0 for column in metric_columns}
    nodes = [
        {
            "node_id": "compound",
            "label": compound_name or "Compound",
            "node_type": "compound",
            **zero_metrics,
        },
        {
            "node_id": "disease",
            "label": disease_name or "Disease",
            "node_type": "disease",
            **zero_metrics,
        },
    ]
    edges = []
    source_map: dict[str, set[str]] = {}
    for source, frame in target_sources.items():
        source_map[source] = set(frame["gene"])
    overlap_genes = set(overlap["gene"])
    for _, row in overlap.iterrows():
        gene = str(row["gene"])
        node = {
            "node_id": gene,
            "label": gene,
            "node_type": "target",
        }
        for column in metric_columns:
            value = row.get(column)
            if pd.isna(value):
                value = 0
            node[column] = value
        nodes.append(node)
        edges.append(
            {
                "source": "compound",
                "target": gene,
                "edge_type": "compound_target",
                "sources": ";".join(
                    s for s, genes in source_map.items() if gene in genes
                ),
            }
        )
        edges.append(
            {
                "source": gene,
                "target": "disease",
                "edge_type": "disease_target",
                "sources": "disease",
            }
        )
    if ppi_edges is not None and not ppi_edges.empty:
        col1, col2 = _edge_columns(ppi_edges)
        if col1 is None or col2 is None:
            col1, col2 = ppi_edges.columns[0], ppi_edges.columns[1]
        selected, score_col = _select_ppi_edges(
            ppi_edges,
            overlap_genes,
            max_ppi_edges=max_ppi_edges,
        )
        for _, row in selected.iterrows():
            a = str(row[col1]).upper()
            b = str(row[col2]).upper()
            edge: dict = {
                "source": a,
                "target": b,
                "edge_type": "ppi",
                "sources": "STRING",
            }
            if score_col is not None and not pd.isna(row.get(score_col)):
                edge["combined_score"] = row[score_col]
            for column in ppi_edges.columns:
                if (
                    column in (col1, col2, score_col)
                    or column in ("source", "target", "sources", "edge_type")
                    or str(column).startswith("_")
                ):
                    continue
                if not pd.isna(row.get(column)):
                    edge[str(column)] = row[column]
            edges.append(edge)
    node_path = out_dir / "ctpd_nodes.csv"
    edge_path = out_dir / "ctpd_edges.csv"
    xgmml_path = out_dir / "ctpd_network.xgmml"
    pd.DataFrame(nodes).to_csv(node_path, index=False)
    edge_frame = pd.DataFrame(edges)
    if edge_frame.empty:
        edge_frame = pd.DataFrame(
            columns=["source", "target", "edge_type", "sources"]
        )
    edge_frame.to_csv(edge_path, index=False)
    network_label = (
        f"{compound_name or 'Compound'} / {disease_name or 'Disease'} "
        "C-T-P-D Network"
    )
    write_xgmml_network(
        node_path,
        edge_path,
        xgmml_path,
        network_label=network_label,
    )

    html_path = out_dir / "ctpd_network.html"
    rows = "".join(
        f"<tr><td>{esc(row['source'])}</td><td>{esc(row['target'])}</td>"
        f"<td>{esc(row['edge_type'])}</td><td>{esc(row['sources'])}</td></tr>"
        for row in edges[:500]
    )
    html = f"""<!doctype html>
<html lang="zh">
<head><meta charset="utf-8"><title>C-T-P-D network</title></head>
<body>
<h1>Compound-Target-Disease Network</h1>
<p>{esc(compound_name or "Compound")} / {esc(disease_name or "Disease")} /
{len(overlap)} overlapping targets / {len(edges)} edges</p>
<table border="1" cellpadding="4" cellspacing="0">
<thead><tr><th>Source</th><th>Target</th><th>Type</th><th>Evidence</th></tr></thead>
<tbody>{rows}</tbody>
</table>
</body>
</html>
"""
    html_path.write_text(html, encoding="utf-8")
    return {
        "nodes": node_path,
        "edges": edge_path,
        "xgmml": xgmml_path,
        "html": html_path,
    }


def run_network_toxicology(cfg, log) -> dict:
    """Run the compound-disease overlap, PPI hub and C-T-P-D export."""
    section = cfg.data.get("network_toxicology", {}) or {}
    compound_name = section.get("compound_name") or "Compound"
    disease_name = section.get("disease_name") or "Disease"
    cytoscape_mode = str(section.get("cytoscape") or "auto").lower()
    if cytoscape_mode not in ("auto", "on", "off"):
        cytoscape_mode = "auto"
    cytoscape_url = section.get("cytoscape_url") or "http://127.0.0.1:1234"
    cytoscape_layout = section.get("cytoscape_layout") or "cose"
    cytoscape_session = bool(section.get("cytoscape_save_session", False))
    max_ppi_edges_value = section.get("max_ppi_edges")
    max_ppi_edges = (
        int(max_ppi_edges_value)
        if max_ppi_edges_value not in (None, "")
        else None
    )
    run_enrichment = bool(section.get("run_enrichment", False))
    out_dir = cfg._resolve(
        section.get("output_dir") or "outputs/run_001/network_toxicology",
        cfg.workdir,
    )
    (out_dir / "data").mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    target_sources = load_target_sources(cfg)
    disease_genes = load_disease_genes(cfg)
    log.info(
        "network toxicology: %s target files, %s disease genes",
        len(target_sources),
        len(disease_genes),
    )

    overlap = overlap_analysis(target_sources, disease_genes)
    overlap_path = out_dir / "data" / "compound_disease_overlap.csv"
    overlap.to_csv(overlap_path, index=False)

    hub_path = out_dir / "data" / "ppi_hub_scores.csv"
    hub_frame = None
    edge_path = section.get("ppi_network_csv")
    if edge_path:
        edge_resolved = cfg._resolve(edge_path, cfg.workdir)
        edge_sep = (
            "\t"
            if edge_resolved.suffix.lower() in (".tsv", ".txt")
            else ","
        )
        raw_edges = pd.read_csv(
            edge_resolved,
            sep=edge_sep,
            dtype=str,
        )
        hub_frame = ppi_hub_scores(
            edge_resolved,
            genes=overlap["gene"].tolist(),
        )
        hub_frame.to_csv(hub_path, index=False)
        overlap = overlap.merge(
            hub_frame[
                ["gene", "ppi_degree", "ppi_betweenness", "ppi_clustering", "ppi_hub_score"]
            ],
            on="gene",
            how="left",
        )
        overlap.to_csv(overlap_path, index=False)

    venn_path = None
    if section.get("venn", True):
        sets = {"disease": disease_genes}
        for source, frame in target_sources.items():
            sets[str(source)] = set(frame["gene"])
        venn_path = make_venn_figure(
            sets,
            out_dir / "figures" / "compound_disease_venn.png",
        )

    ctpd = write_ctpd_network(
        overlap,
        target_sources,
        raw_edges if edge_path else None,
        compound_name,
        disease_name,
        out_dir / "data",
        max_ppi_edges=max_ppi_edges,
    )
    enrichment_result = {"status": "disabled"}
    if run_enrichment and len(overlap) >= 3:
        enrichment_result = _run_network_enrichment(
            overlap["gene"].tolist(),
            out_dir,
            cfg,
            log,
        )
    cytoscape_result: dict = {"status": "off"}
    if cytoscape_mode != "off":
        cytoscape_result = export_network_to_cytoscape(
            ctpd["nodes"],
            ctpd["edges"],
            out_dir,
            network_title=(
                f"{compound_name or 'Compound'} / "
                f"{disease_name or 'Disease'}"
            ),
            base_url=str(cytoscape_url),
            layout=str(cytoscape_layout),
            save_session=cytoscape_session,
            mode=cytoscape_mode,
        )
        for key, value in list(cytoscape_result.items()):
            if isinstance(value, Path):
                cytoscape_result[key] = str(value)
    all_compound_targets: set[str] = set()
    for frame in target_sources.values():
        all_compound_targets.update(frame["gene"].astype(str))
    summary = {
        "compound_name": compound_name,
        "disease_name": disease_name,
        "compound_targets": int(len(all_compound_targets)),
        "disease_genes": int(len(disease_genes)),
        "overlap_genes": int(len(overlap)),
        "ppi_hub_scored": bool(hub_frame is not None),
        "cytoscape": cytoscape_result,
        "outputs": {
            "overlap_csv": str(overlap_path),
            "ppi_hub_csv": str(hub_path) if hub_frame is not None else "",
            "venn": venn_path or "",
            "ctpd": {k: str(v) for k, v in ctpd.items()},
            "enrichment": enrichment_result,
        },
    }
    write_json(out_dir / "network_toxicology_summary.json", summary)
    log.info(
        "network toxicology complete: %s overlap genes -> %s",
        len(overlap),
        out_dir,
    )
    return summary


def _run_network_enrichment(
    genes: list[str],
    out_dir: Path,
    cfg,
    log,
) -> dict:
    """Run the existing clusterProfiler helper on the overlap gene set."""
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    script = Path(__file__).resolve().parent / "insilico_enrichment.R"
    if not rscript or not script.exists():
        return {"status": "skipped", "reason": "Rscript or enrichment script unavailable"}
    data_dir = out_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    gene_csv = data_dir / "network_enrichment_input.csv"
    pd.DataFrame({"gene": genes}).to_csv(gene_csv, index=False)
    species = str((cfg.data.get("species") or "hs")).lower()
    species = "mm" if species.startswith("mm") else "hs"
    try:
        proc = subprocess.run(
            [rscript, str(script), str(gene_csv), str(data_dir), species],
            cwd=cfg.workdir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=int(
                (cfg.data.get("network_toxicology", {}) or {}).get(
                    "enrichment_timeout",
                    900,
                )
            ),
        )
    except Exception as exc:
        return {"status": "failed", "reason": str(exc)}
    if proc.returncode != 0:
        return {
            "status": "failed",
            "reason": (proc.stderr or proc.stdout)[-2000:],
        }
    go_path = data_dir / "insilico_go_enrichment.csv"
    kegg_path = data_dir / "insilico_kegg_enrichment.csv"
    return {
        "status": "completed",
        "go_csv": str(go_path) if go_path.exists() else "",
        "kegg_csv": str(kegg_path) if kegg_path.exists() else "",
    }
