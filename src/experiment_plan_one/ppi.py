"""STRING PPI retrieval, topology ranking and module detection."""

from __future__ import annotations

import io
import math
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import cm

from .common import LOG, ensure_dir, save_figure, write_json
from .targets import make_venn_figure

STRING_NETWORK_URL = "https://string-db.org/api/tsv/network"


def fetch_string_network(
    genes: list[str],
    output_path: Path,
    *,
    species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 50,
    force: bool = False,
    timeout: int = 180,
) -> pd.DataFrame:
    """Fetch a high-confidence STRING network and normalize its columns."""
    ensure_dir(output_path.parent)
    if output_path.exists() and not force:
        return pd.read_csv(output_path, sep="\t")
    identifiers = "\r".join(sorted(set(gene.upper() for gene in genes if gene)))
    params = urllib.parse.urlencode(
        {
            "identifiers": identifiers,
            "species": str(species),
            "required_score": str(required_score),
            "add_nodes": str(max(0, int(add_nodes))),
            "caller_identity": "liver-cancer-experiment-plan-one",
        }
    )
    request = urllib.request.Request(
        f"{STRING_NETWORK_URL}?{params}",
        headers={"User-Agent": "Mozilla/5.0"},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        text = response.read().decode("utf-8", "replace")
    frame = pd.read_csv(io.StringIO(text), sep="\t")
    if frame.empty:
        raise RuntimeError("STRING returned no network edges")
    column_lookup = {str(column).lower(): str(column) for column in frame.columns}
    left = column_lookup.get("preferredname_a") or column_lookup.get("stringid_a")
    right = column_lookup.get("preferredname_b") or column_lookup.get("stringid_b")
    score = column_lookup.get("score") or column_lookup.get("combined_score")
    if not left or not right or not score:
        raise RuntimeError(f"unexpected STRING columns: {list(frame.columns)}")
    output = pd.DataFrame(
        {
            "node1": frame[left].astype(str).str.upper(),
            "node2": frame[right].astype(str).str.upper(),
            "score": pd.to_numeric(frame[score], errors="coerce"),
            "source": "STRING",
        }
    ).dropna(subset=["node1", "node2", "score"])
    output = output[output["node1"] != output["node2"]].drop_duplicates()
    output.to_csv(output_path, sep="\t", index=False)
    return output


def ppi_hub_metrics(edges: pd.DataFrame, genes: list[str] | None = None) -> pd.DataFrame:
    """Calculate common centralities and an interpretable consensus rank."""
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        graph.add_edge(
            str(row.node1).upper(),
            str(row.node2).upper(),
            weight=float(row.score),
        )
    if graph.number_of_nodes() == 0:
        raise ValueError("PPI graph has no nodes")
    nodes = sorted(set(genes or graph.nodes()))
    nodes = [node for node in nodes if node in graph]
    if not nodes:
        raise ValueError("none of the requested genes are present in the PPI graph")
    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, weight=None)
    closeness = nx.closeness_centrality(graph)
    try:
        eigenvector = nx.eigenvector_centrality_numpy(graph)
    except Exception:
        eigenvector = {node: 0.0 for node in graph}
    pagerank = nx.pagerank(graph, weight="weight")
    clustering = nx.clustering(graph)
    mcc = _maximum_clique_centrality(graph)
    frame = pd.DataFrame(
        {
            "gene": nodes,
            "degree": [degree.get(node, 0) for node in nodes],
            "betweenness": [betweenness.get(node, 0.0) for node in nodes],
            "closeness": [closeness.get(node, 0.0) for node in nodes],
            "eigenvector": [eigenvector.get(node, 0.0) for node in nodes],
            "pagerank": [pagerank.get(node, 0.0) for node in nodes],
            "clustering": [clustering.get(node, 0.0) for node in nodes],
            "mcc": [mcc.get(node, 0.0) for node in nodes],
        }
    )
    rank_columns = ["degree", "betweenness", "closeness", "eigenvector", "pagerank", "mcc"]
    ranks = frame[rank_columns].rank(ascending=False, method="average", pct=True)
    frame["ppi_rank_score"] = 1.0 - ranks.mean(axis=1)
    frame["degree_rank"] = frame["degree"].rank(ascending=False, method="min").astype(int)
    frame["betweenness_rank"] = frame["betweenness"].rank(ascending=False, method="min").astype(int)
    return frame.sort_values(
        ["ppi_rank_score", "degree", "betweenness", "gene"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def _maximum_clique_centrality(graph: nx.Graph) -> dict[str, float]:
    scores = {node: 0.0 for node in graph}
    try:
        for clique in nx.find_cliques(graph):
            weight = math.factorial(max(len(clique) - 1, 0))
            for node in clique:
                scores[node] += weight
    except Exception:
        pass
    return scores


def louvain_clusters(graph: nx.Graph) -> list[set[str]]:
    """Return deterministic Louvain communities for a weighted PPI network."""
    if graph.number_of_nodes() == 0:
        return []
    communities = nx.community.louvain_communities(
        graph,
        weight="weight",
        seed=42,
    )
    return [
        set(str(node) for node in community)
        for community in sorted(communities, key=len, reverse=True)
    ]


def run_ppi_analysis(
    genes: list[str],
    output_dir: Path,
    *,
    required_score: int = 700,
    top_n: int = 20,
    add_nodes: int = 50,
    force: bool = False,
) -> dict[str, Any]:
    ensure_dir(output_dir)
    edges = fetch_string_network(
        genes,
        output_dir / "string_edges.tsv",
        required_score=required_score,
        add_nodes=add_nodes,
        force=force,
    )
    metrics = ppi_hub_metrics(edges, genes=genes)
    metrics.to_csv(output_dir / "ppi_hub_metrics.csv", index=False)
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        graph.add_edge(
            str(row.node1).upper(),
            str(row.node2).upper(),
            weight=float(row.score) / 1000.0,
        )
    clusters = louvain_clusters(graph)
    major_clusters = [cluster for cluster in clusters if len(cluster) >= 3]
    assigned_major = set().union(*major_clusters) if major_clusters else set()
    clusters_for_output = major_clusters + [
        {node} for node in graph.nodes if node not in assigned_major
    ]
    cluster_map = {
        node: f"M{index + 1}"
        for index, cluster in enumerate(major_clusters)
        for node in cluster
    }
    node_modules = {
        node: cluster_map.get(node, "Other")
        for node in graph.nodes()
    }
    metrics["module"] = metrics["gene"].map(cluster_map).fillna("Other")
    metrics.to_csv(output_dir / "ppi_hub_metrics.csv", index=False)

    _draw_network(
        graph,
        metrics,
        output_dir / "fig2a_string_network.png",
        title=f"STRING PPI network | score >= {required_score / 1000:.2f}",
        node_color="#3f7f93",
        focus_genes=set(metrics["gene"]),
    )
    _draw_network(
        graph,
        metrics,
        output_dir / "fig2b_cytoscape_module_network.png",
        title="PPI modules (Louvain communities)",
        node_color=None,
        focus_genes=set(metrics["gene"]),
        node_modules=node_modules,
    )
    _bar_rank(
        metrics,
        "degree",
        output_dir / "fig2c_degree_top20.png",
        "Degree",
        top_n,
    )
    _bar_rank(
        metrics,
        "betweenness",
        output_dir / "fig2d_betweenness_top20.png",
        "Betweenness centrality",
        top_n,
    )
    mcc_top = set(metrics.nlargest(10, "mcc")["gene"])
    degree_top = set(metrics.nlargest(10, "degree")["gene"])
    make_venn_figure(
        {"MCC Top 10": mcc_top, "Degree Top 10": degree_top},
        output_dir / "fig2e_mcc_degree_venn.png",
        title="Consensus hub genes",
    )
    metrics.to_csv(output_dir / "ppi_hub_metrics.csv", index=False)
    write_json(
        output_dir / "ppi_summary.json",
        {
            "input_genes": len(set(genes)),
            "network_nodes": int(graph.number_of_nodes()),
            "network_edges": int(graph.number_of_edges()),
            "required_score": required_score,
            "module_method": "Louvain",
            "modules": {
                f"M{i + 1}": sorted(cluster)
                for i, cluster in enumerate(major_clusters)
            },
            "top_consensus": metrics.head(10).to_dict(orient="records"),
        },
    )
    return {
        "edges": output_dir / "string_edges.tsv",
        "metrics": output_dir / "ppi_hub_metrics.csv",
        "top_genes": metrics["gene"].head(top_n).tolist(),
        "top3": metrics["gene"].head(3).tolist(),
        "modules": clusters_for_output,
    }


def _draw_network(
    graph: nx.Graph,
    metrics: pd.DataFrame,
    output: Path,
    *,
    title: str,
    node_color: str | None,
    focus_genes: set[str] | None = None,
    node_modules: dict[str, str] | None = None,
) -> None:
    if graph.number_of_nodes() == 0:
        raise ValueError("cannot draw an empty PPI graph")
    position = nx.spring_layout(graph, seed=42, weight="weight", k=1.2 / np.sqrt(max(graph.number_of_nodes(), 1)))
    degree = dict(graph.degree())
    focus_genes = set(focus_genes or ())
    sizes = [180 + 42 * degree.get(node, 0) for node in graph.nodes()]
    if node_color is None:
        module_column = (
            pd.Series(node_modules)
            if node_modules
            else metrics.set_index("gene")["module"]
        )
        module_values = [
            module_column.get(node, "Other") for node in graph.nodes()
        ]
        modules = [
            module
            for module, count in pd.Series(module_values).value_counts().items()
            if module == "Other" or count >= 3
        ][:10]
        palette = cm.get_cmap("tab20", max(len(modules), 1))
        color_map = {
            module: palette(index)
            for index, module in enumerate(modules)
        }
        colors = [
            color_map.get(module_column.get(node, "Other"), "#b8bec5")
            for node in graph.nodes()
        ]
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=color_map[module],
                markersize=7,
                label=module,
            )
            for module in modules
        ]
    else:
        colors = [node_color] * graph.number_of_nodes()
        legend_handles = []
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    widths = [
        0.8 + 4.0 * float(data.get("weight", 0.7))
        for _, _, data in graph.edges(data=True)
    ]
    nx.draw_networkx_edges(
        graph,
        position,
        ax=ax,
        width=widths,
        alpha=0.42,
        edge_color="#66737f",
    )
    node_order = list(graph.nodes())
    if node_color is not None and focus_genes:
        background = [node for node in node_order if node not in focus_genes]
        highlighted = [node for node in node_order if node in focus_genes]
        if background:
            nx.draw_networkx_nodes(
                graph,
                position,
                ax=ax,
                nodelist=background,
                node_size=[sizes[node_order.index(node)] for node in background],
                node_color="#aebac3",
                linewidths=0.5,
                edgecolors="white",
            )
        nx.draw_networkx_nodes(
            graph,
            position,
            ax=ax,
            nodelist=highlighted,
            node_size=[sizes[node_order.index(node)] for node in highlighted],
            node_color=node_color,
            linewidths=0.7,
            edgecolors="white",
        )
        legend_handles = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=node_color,
                markersize=7,
                label="Seeded overlap genes",
            ),
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor="#aebac3",
                markersize=7,
                label="STRING context nodes (unlabeled)",
            ),
        ]
    else:
        nx.draw_networkx_nodes(
            graph,
            position,
            ax=ax,
            node_size=sizes,
            node_color=colors,
            linewidths=0.6,
            edgecolors="white",
        )
    if focus_genes:
        focus_texts = nx.draw_networkx_labels(
            graph,
            position,
            ax=ax,
            labels={node: node for node in graph.nodes() if node in focus_genes},
            font_size=6.8,
            font_family="Arial",
            font_weight="bold",
            font_color="#102b3d",
        )
        for text in focus_texts.values():
            text.set_path_effects(
                [path_effects.withStroke(linewidth=1.8, foreground="white")]
            )
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc="lower left",
            bbox_to_anchor=(0.0, -0.03),
            frameon=False,
            fontsize=6,
        )
    ax.set_title(title, fontsize=8, fontweight="bold")
    ax.axis("off")
    save_figure(fig, output)


def _bar_rank(
    metrics: pd.DataFrame,
    column: str,
    output: Path,
    title: str,
    top_n: int,
) -> None:
    values = metrics.nlargest(top_n, column).sort_values(column, ascending=True)
    fig, ax = plt.subplots(figsize=(7, max(4, len(values) * 0.3)))
    colors = ["#c05640" if index < 5 else "#627d98" for index in range(len(values) - 1, -1, -1)]
    ax.barh(values["gene"], values[column], color=colors)
    ax.set_xlabel(title)
    ax.set_ylabel("")
    ax.set_title(f"{title}: Top {min(top_n, len(values))}", fontweight="bold")
    save_figure(fig, output)
