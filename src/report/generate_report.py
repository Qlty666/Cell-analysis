#!/usr/bin/env python3
"""Generate a self-contained HTML report from pipeline outputs."""

import csv
import html
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_ROOT = Path(
    os.environ.get("LIVER_OUTPUT_ROOT", str(ROOT.parent / "liver_cancer"))
).resolve()
RES = OUTPUT_ROOT / "results"
FIG = RES / "figures"
DATA = RES / "data"


def esc(value) -> str:
    return html.escape(str(value))


def read_table(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def render_table(rows: list[dict], max_rows: int = 25, cols=None) -> str:
    if not rows:
        return "<p class='muted'>No data.</p>"
    if cols is None:
        cols = list(rows[0].keys())
    cols = cols[:10]
    rows = rows[:max_rows]

    header = "".join(f"<th>{esc(c)}</th>" for c in cols)
    body = ""
    for row in rows:
        body += "<tr>" + "".join(f"<td>{esc(row.get(c, ''))}</td>" for c in cols) + "</tr>"
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def image_card(name: str, title: str) -> str:
    path = FIG / name
    if not path.exists():
        return ""
    return (
        f"<figure><img src='figures/{name}' alt='{esc(title)}'>"
        f"<figcaption>{esc(title)}</figcaption></figure>"
    )


def data_link(name: str, label: str) -> str:
    path = DATA / name
    if not path.exists():
        return ""
    return f"<li><a href='data/{name}'>{esc(label)}</a></li>"


def main() -> int:
    summary_path = RES / "summary.json"
    if not summary_path.exists():
        (RES / "result_report.html").write_text(
            "<html><body><h1>Report not available yet</h1></body></html>",
            encoding="utf-8",
        )
        return 1

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    dataset = summary.get("dataset", "GSE")
    condition_counts = summary.get("condition_counts", {})
    condition_html = "".join(
        f"<span class='chip'>{esc(k)}: {esc(v)}</span>"
        for k, v in condition_counts.items()
    )

    top_degs = summary.get("top_degs", [])
    if isinstance(top_degs, dict):
        top_degs = [top_degs]
    deg_columns = ["gene", "avg_log2FC", "p_val_adj"]

    go_up = read_table(DATA / "enrichment_up_go.csv")
    go_down = read_table(DATA / "enrichment_down_go.csv")
    kegg_up = read_table(DATA / "enrichment_up_kegg.csv")
    kegg_down = read_table(DATA / "enrichment_down_kegg.csv")

    enrichment_cols = ["ID", "Description", "GeneRatio", "pvalue", "p.adjust"]
    if go_up and "Description" in go_up[0]:
        go_html = (
            "<h3>上调基因 GO BP</h3>" + render_table(go_up, 10, enrichment_cols) +
            "<h3>下调基因 GO BP</h3>" + render_table(go_down, 10, enrichment_cols)
        )
    else:
        go_html = "<p class='muted'>未找到显著 GO 富集结果。</p>"

    if kegg_up and "Description" in kegg_up[0]:
        kegg_html = (
            "<h3>上调基因 KEGG</h3>" + render_table(kegg_up, 10, enrichment_cols) +
            "<h3>下调基因 KEGG</h3>" + render_table(kegg_down, 10, enrichment_cols)
        )
    else:
        kegg_html = "<p class='muted'>未找到显著 KEGG 富集结果。</p>"

    links = "".join([
        data_link("deg_all.csv", "全部差异表达基因"),
        data_link("deg_significant.csv", "显著差异表达基因"),
        data_link("qc_metrics.csv", "QC 指标"),
        data_link("doublet_results.csv", "双细胞结果"),
        data_link("cell_annotations.csv", "细胞注释"),
        data_link("annotation_confusion.csv", "注释混淆矩阵"),
        data_link("umap_coordinates.csv", "UMAP 坐标"),
        data_link("cluster_composition.csv", "聚类组成"),
        data_link("enrichment_up_go.csv", "上调 GO 富集"),
        data_link("enrichment_down_go.csv", "下调 GO 富集"),
        data_link("enrichment_up_kegg.csv", "上调 KEGG 富集"),
        data_link("enrichment_down_kegg.csv", "下调 KEGG 富集"),
        data_link("cell_cycle_scores.csv", "细胞周期打分"),
        data_link("doublet_rate_by_sample.csv", "样本双细胞率"),
        data_link("cluster_markers.csv", "聚类 marker 基因"),
        data_link("signature_scores.csv", "功能签名打分"),
        data_link("cnv_heatmap.csv", "推断 CNV 热图矩阵"),
        data_link("singleR_annotations.csv", "SingleR 注释"),
        data_link("singleR_confusion.csv", "SingleR 混淆矩阵"),
        data_link("trajectory_pseudotime.csv", "拟时序结果"),
        data_link("cellchat_communication.csv", "CellChat 通讯矩阵"),
        data_link("cellchat_pathways.csv", "CellChat 通路表"),
        data_link("ml_classification_report.csv", "ML 分类报告"),
        data_link("sample_annotations.csv", "样本分组"),
        data_link("liver_cancer_seurat.rds", "Seurat 对象"),
    ])

    html_doc = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>肝癌单细胞分析报告</title>
<style>
body {{ font-family: "Segoe UI", Arial, sans-serif; margin: 0; background: #f5f7fa; color: #1f2933; }}
.wrap {{ max-width: 1180px; margin: 0 auto; padding: 24px; }}
h1 {{ font-size: 28px; margin: 0 0 8px; }}
h2 {{ font-size: 20px; margin: 32px 0 12px; border-bottom: 2px solid #dbe2ea; padding-bottom: 6px; }}
h3 {{ font-size: 16px; margin: 18px 0 8px; }}
.sub {{ color: #52606d; margin-bottom: 16px; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin: 18px 0; }}
.card {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 14px; }}
.card .num {{ font-size: 24px; font-weight: 700; color: #1665c0; }}
.card .label {{ color: #52606d; font-size: 13px; }}
.chip {{ display: inline-block; background: #eef2f7; border-radius: 12px; padding: 4px 10px; margin: 3px; font-size: 13px; }}
figure {{ margin: 14px 0; }}
img {{ max-width: 100%; border: 1px solid #e4e7eb; border-radius: 8px; background: #fff; }}
figcaption {{ color: #52606d; font-size: 13px; margin-top: 6px; }}
table {{ width: 100%; border-collapse: collapse; background: #fff; font-size: 12px; }}
th, td {{ border: 1px solid #e4e7eb; padding: 6px 8px; text-align: left; }}
th {{ background: #f0f4f8; }}
tr:nth-child(even) {{ background: #fafbfc; }}
ul {{ line-height: 1.8; }}
a {{ color: #1665c0; text-decoration: none; }}
.muted {{ color: #7b8794; }}
footer {{ margin-top: 40px; color: #7b8794; font-size: 12px; }}
</style>
</head>
<body>
<div class="wrap">
<h1>单细胞分析报告</h1>
<p class="sub">{esc(summary.get('title', 'single-cell analysis'))} | {esc(dataset)}</p>

<div class="cards">
  <div class="card"><div class="num">{esc(summary.get('n_cells_raw', 'NA'))}</div><div class="label">原始细胞</div></div>
  <div class="card"><div class="num">{esc(summary.get('n_cells_after_qc', 'NA'))}</div><div class="label">QC 后细胞</div></div>
  <div class="card"><div class="num">{esc(summary.get('n_cells_after_doublet_removal', 'NA'))}</div><div class="label">去双细胞后</div></div>
  <div class="card"><div class="num">{esc(summary.get('n_genes', 'NA'))}</div><div class="label">基因数</div></div>
  <div class="card"><div class="num">{esc(summary.get('n_clusters', 'NA'))}</div><div class="label">聚类数</div></div>
  <div class="card"><div class="num">{esc(summary.get('n_celltypes', 'NA'))}</div><div class="label">细胞类型</div></div>
  <div class="card"><div class="num">{esc(summary.get('deg_up', 'NA'))}</div><div class="label">HCC 上调 DEG</div></div>
  <div class="card"><div class="num">{esc(summary.get('deg_down', 'NA'))}</div><div class="label">HCC 下调 DEG</div></div>
</div>

<p>样本分组：{condition_html}</p>

<h2>1. 质控</h2>
{image_card('fig_01_qc_raw_violin.png', '原始 QC 指标')}
{image_card('fig_01_qc_filtered_violin.png', '过滤后 QC 指标')}

<h2>2. 双细胞检测</h2>
{image_card('fig_02_doublet_scores.png', '双细胞得分')}

<h2>3. 聚类与 UMAP</h2>
{image_card('fig_03_umap_clusters.png', 'Seurat 聚类 UMAP')}
{image_card('fig_04_umap_condition.png', 'HCC vs iCCA UMAP')}

<h2>4. 细胞注释</h2>
{image_card('fig_05_umap_annotation.png', '注释与发表细胞类型 UMAP')}
{image_card('fig_06_dotplot_markers.png', 'Marker 基因 DotPlot')}
{image_card('fig_07_annotation_confusion_heatmap.png', '注释混淆矩阵')}
{image_card('fig_16_featureplot_markers.png', 'Marker 基因 FeaturePlot')}
{image_card('fig_17_marker_violin.png', 'Marker 基因小提琴图')}
{image_card('fig_18_celltype_proportion.png', '细胞类型比例堆叠图')}
{image_card('fig_19_condition_proportion.png', '分组构成比例图')}

<h2>5. 差异表达与火山图</h2>
{image_card('fig_08_volcano.png', '差异表达火山图')}
{image_card('fig_09_deg_heatmap.png', 'Top DEG 热图')}
{image_card('fig_14_pca.png', 'PCA 分组图')}
{image_card('fig_15_elbow.png', '主成分 Elbow 图')}
<h3>Top 差异表达基因</h3>
{render_table(top_degs, 25, deg_columns)}

<h2>6. 富集分析</h2>
{image_card('fig_10_go_up.png', '上调基因 GO BP')}
{image_card('fig_11_go_down.png', '下调基因 GO BP')}
{image_card('fig_12_kegg_up.png', '上调基因 KEGG')}
{image_card('fig_13_kegg_down.png', '下调基因 KEGG')}
{image_card('fig_20_gsea_go.png', 'GSEA GO BP')}
{image_card('fig_21_gsea_kegg.png', 'GSEA KEGG')}
{image_card('fig_22_go_network.png', 'GO BP 通路网络')}
{image_card('fig_23_kegg_network.png', 'KEGG 通路网络')}

<h2>7. 机器学习分析</h2>
{image_card('fig_24_ml_feature_importance.png', 'ML 特征重要性')}
{image_card('fig_25_ml_shap.png', 'SHAP 可解释性')}
{image_card('fig_43_ml_confusion_matrix.png', 'ML 混淆矩阵')}
{image_card('fig_44_ml_roc_pr.png', 'ML ROC 与 PR 曲线')}
{image_card('fig_45_ml_cv_scores.png', 'ML 交叉验证得分')}
{go_html}
{kegg_html}

<h2>8. 高级分析与发表图</h2>
{image_card('fig_26_cellcycle_umap.png', '细胞周期 UMAP')}
{image_card('fig_27_cellcycle_proportion.png', '细胞周期比例')}
{image_card('fig_28_umap_sample.png', 'UMAP 按样本')}
{image_card('fig_29_doublet_rate_sample.png', '样本双细胞率')}
{image_card('fig_30_sample_proportion.png', '样本细胞类型比例')}
{image_card('fig_31_cluster_marker_heatmap.png', '聚类 marker 热图')}
{image_card('fig_32_cluster_marker_dotplot.png', '聚类 marker DotPlot')}
{image_card('fig_33_signature_scores_umap.png', '功能签名 UMAP')}
{image_card('fig_34_signature_scores_boxplot.png', '功能签名箱线图')}
{image_card('fig_35_celltype_abundance_effect.png', '细胞类型丰度变化')}
{image_card('fig_36_cnv_heatmap.png', '推断 CNV 热图')}
{image_card('fig_37_singler_umap.png', 'SingleR 注释 UMAP')}
{image_card('fig_38_singler_confusion_heatmap.png', 'SingleR 混淆矩阵')}
{image_card('fig_39_trajectory_umap.png', '拟时序轨迹')}

<h2>9. CellChat 细胞通讯</h2>
{image_card('fig_40_cellchat_network.png', 'CellChat 通讯网络')}
{image_card('fig_41_cellchat_heatmap.png', 'CellChat 通讯热图')}
{image_card('fig_42_cellchat_bubble.png', 'CellChat 配体受体气泡图')}

<h2>10. 数据文件</h2>
<ul>{links}</ul>

<footer>Pipeline: R / Seurat / scDblFinder / clusterProfiler | 自动生成于 {esc(summary.get('finished_at', ''))}</footer>
</div>
</body>
</html>
"""

    (RES / "result_report.html").write_text(html_doc, encoding="utf-8")
    print("HTML report written:", RES / "result_report.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
