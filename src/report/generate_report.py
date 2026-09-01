#!/usr/bin/env python3
"""Generate a self-contained HTML report from pipeline outputs."""

import csv
import html
import json
import os
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_ROOT = Path(
    os.environ.get("LIVER_OUTPUT_ROOT", str(ROOT.parent / "liver_cancer"))
).resolve()
RES = OUTPUT_ROOT / "results"
FIG = RES / "figures"
DATA = RES / "data"

MAX_ANALYSIS_SAMPLE_ROWS = 5000
MAX_PREVIEW_ROWS = 6
MAX_PREVIEW_COLS = 8
MAX_NUMERIC_COLS = 20
MAX_STAT_COLS = 5
IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".svg", ".gif"}
DATA_SUFFIXES = {".csv", ".tsv", ".json", ".xlsx", ".xls", ".rds", ".txt"}

STAGE_LABELS = {
    "01_qc": "质控",
    "02_doublets": "双细胞检测",
    "03_cluster": "聚类",
    "04_annotation": "细胞注释",
    "05_deg": "差异表达",
    "06_enrichment": "富集分析",
    "07_ml": "机器学习",
    "08_publication": "高级分析与发表图",
    "09_cellchat": "细胞通讯",
}

FIGURE_GUIDE = {
    "fig_01_qc_raw_violin.png": ("QC 小提琴图（原始）", "过滤前 nFeature、nCount、percent.mt 的分布，用于识别低质量细胞和数据批次差异。"),
    "fig_01_qc_filtered_violin.png": ("QC 小提琴图（过滤后）", "过滤后的 QC 指标分布，确认阈值是否合理、主要细胞群是否保留。"),
    "fig_48_qc_pvalue_comparison.png": ("QC 质控差异度 P 值图", "以 Wilcoxon 秩和检验 P 值衡量原始/过滤后 QC 指标在不同条件间的差异程度。"),
    "fig_02_doublet_scores.png": ("双细胞得分图", "scDblFinder 双细胞得分的分布，用于判断双细胞分类边界。"),
    "fig_03_umap_clusters.png": ("UMAP 聚类图", "Seurat 聚类在 UMAP 上的结构，检查分群是否清晰、是否存在过度分割。"),
    "fig_04_umap_condition.png": ("UMAP 分组图", "不同条件下细胞在 UMAP 上的分布，检查分组偏移和批次效应。"),
    "fig_05_umap_annotation.png": ("UMAP 注释图", "细胞类型注释在 UMAP 上的分布，验证注释与聚类结构的一致性。"),
    "fig_06_dotplot_markers.png": ("Marker 基因 DotPlot", "marker 基因在各细胞类型中的表达比例和表达量，用于验证注释。"),
    "fig_07_annotation_confusion_heatmap.png": ("注释混淆矩阵热图", "自动注释与发布注释之间的一致性，识别容易混淆的细胞类型。"),
    "fig_08_volcano.png": ("差异表达火山图", "以 log2FC 和校正 P 值展示差异基因，观察上调和下调基因的显著性分布。"),
    "fig_09_deg_heatmap.png": ("Top DEG 热图", "Top 差异基因在样本/条件下的表达模式，用于查看差异基因的整体趋势。"),
    "fig_09_deg_horizontal_violin.png": ("Top DEG 横向小提琴图（P 值）", "按校正 P 值排序的差异最显著基因，横向对比两组表达分布并标注 P 值。"),
    "fig_10_go_up.png": ("GO BP 富集图（上调）", "上调基因的 GO 生物过程富集结果，查看显著富集的通路。"),
    "fig_11_go_down.png": ("GO BP 富集图（下调）", "下调基因的 GO 生物过程富集结果，查看显著富集的通路。"),
    "fig_12_kegg_up.png": ("KEGG 富集图（上调）", "上调基因的 KEGG 通路富集结果。"),
    "fig_13_kegg_down.png": ("KEGG 富集图（下调）", "下调基因的 KEGG 通路富集结果。"),
    "fig_14_pca.png": ("PCA 分组图", "主成分空间中的样本/细胞分布，观察主要变异来源。"),
    "fig_15_elbow.png": ("主成分 Elbow 图", "各主成分解释的方差比例，用于选择后续分析的主成分数。"),
    "fig_16_featureplot_markers.png": ("Marker 基因 FeaturePlot", "marker 基因在 UMAP 上的表达位置，验证细胞类型注释。"),
    "fig_17_marker_violin.png": ("Marker 基因小提琴图", "marker 基因在不同细胞类型中的表达分布。"),
    "fig_18_celltype_proportion.png": ("细胞类型比例堆叠图", "不同样本/条件下细胞类型构成的变化。"),
    "fig_19_condition_proportion.png": ("分组构成比例图", "不同分组中细胞类型的比例对比。"),
    "fig_20_gsea_go.png": ("GSEA GO BP 富集图", "GO 生物过程的 GSEA 富集曲线，查看通路的整体上调/下调方向。"),
    "fig_21_gsea_kegg.png": ("GSEA KEGG 富集图", "KEGG 通路的 GSEA 富集曲线。"),
    "fig_22_go_network.png": ("GO BP 通路网络图（Top5 核心 + 5 延伸）", "上调基因 GO 富集结果经 p.adjust 筛选后保留前 5 个核心通路，并追加最佳 5 个延伸通路的基因-通路关系网络。"),
    "fig_23_kegg_network.png": ("KEGG 通路网络图（Top5 核心 + 5 延伸）", "上调基因 KEGG 富集结果经 p.adjust 筛选后保留前 5 个核心通路，并追加最佳 5 个延伸通路的基因-通路关系网络。"),
    "fig_46_go_top5.png": ("GO BP 筛选后 Top5 富集图", "上调基因 GO 富集结果按 p.adjust 排序的前 5 条通路气泡图。"),
    "fig_47_kegg_top5.png": ("KEGG 筛选后 Top5 富集图", "上调基因 KEGG 富集结果按 p.adjust 排序的前 5 条通路气泡图。"),
    "fig_24_ml_feature_importance.png": ("ML 特征重要性图", "分类模型中的重要基因/特征，用于筛选关键变量。"),
    "fig_25_ml_shap.png": ("SHAP 可解释性图", "特征对模型预测的方向和贡献，解释分类结果。"),
    "fig_26_cellcycle_umap.png": ("细胞周期 UMAP", "细胞周期阶段在 UMAP 上的分布，检查细胞周期对聚类的干扰。"),
    "fig_27_cellcycle_proportion.png": ("细胞周期比例图", "不同细胞类型/分组中的细胞周期比例。"),
    "fig_28_umap_sample.png": ("UMAP 按样本", "样本来源在 UMAP 上的分布，评估样本混合和批次效应。"),
    "fig_29_doublet_rate_sample.png": ("样本双细胞率图", "每个样本的双细胞率，识别异常样本。"),
    "fig_30_sample_proportion.png": ("样本细胞类型比例图", "每个样本的细胞类型组成，检查样本间异质性。"),
    "fig_31_cluster_marker_heatmap.png": ("聚类 Marker 热图", "各聚类 marker 基因的表达模式，用于重新定义聚类标签。"),
    "fig_32_cluster_marker_dotplot.png": ("聚类 Marker DotPlot", "各聚类 marker 的表达比例和水平。"),
    "fig_33_signature_scores_umap.png": ("功能签名 UMAP", "增殖、EMT、缺氧等功能签名在 UMAP 上的分布。"),
    "fig_34_signature_scores_boxplot.png": ("功能签名箱线图", "不同条件下功能签名得分的差异。"),
    "fig_35_celltype_abundance_effect.png": ("细胞类型丰度变化图", "细胞类型丰度在组间的变化方向和显著性。"),
    "fig_36_cnv_heatmap.png": ("推断 CNV 热图", "基于表达推断的拷贝数变异模式，用于区分恶性/非恶性细胞。"),
    "fig_37_singler_umap.png": ("SingleR 注释 UMAP", "SingleR 参考注释在 UMAP 上的分布。"),
    "fig_38_singler_confusion_heatmap.png": ("SingleR 混淆矩阵热图", "SingleR 注释与现有注释的一致性。"),
    "fig_39_trajectory_umap.png": ("拟时序轨迹图", "拟时序分析结果，查看细胞分化轨迹。"),
    "fig_40_cellchat_network.png": ("CellChat 通讯网络图", "细胞类型之间的配体受体通讯网络。"),
    "fig_41_cellchat_heatmap.png": ("CellChat 通讯热图", "细胞类型对之间的通讯强度。"),
    "fig_42_cellchat_bubble.png": ("CellChat 配体受体气泡图", "关键配体受体对在细胞类型间的通讯关系。"),
    "fig_43_ml_confusion_matrix.png": ("ML 混淆矩阵", "分类模型的预测混淆情况，查看具体错分。"),
    "fig_44_ml_roc_pr.png": ("ML ROC 与 PR 曲线", "分类模型的 ROC 和 PR 性能，评估分类能力。"),
    "fig_45_ml_cv_scores.png": ("ML 交叉验证得分图", "交叉验证得分的分布，评估模型稳定性。"),
}

DATA_GUIDE = {
    "fig_01_qc_metrics.csv": ("QC 指标表", "每细胞的 nFeature、nCount、percent.mt 以及样本/分组信息，用于检查过滤前后细胞质量。"),
    "qc_metrics.csv": ("QC 指标表", "每细胞的 nFeature、nCount、percent.mt 以及样本/分组信息，用于检查过滤前后细胞质量。"),
    "fig_02_doublet_results.csv": ("双细胞结果表", "每细胞的双细胞得分和分类结果，用于评估双细胞去除。"),
    "doublet_results.csv": ("双细胞结果表", "每细胞的双细胞得分和分类结果，用于评估双细胞去除。"),
    "fig_05_16_17_cell_annotations.csv": ("细胞注释表", "每细胞的聚类、细胞类型注释和样本分组信息。"),
    "cell_annotations.csv": ("细胞注释表", "每细胞的聚类、细胞类型注释和样本分组信息。"),
    "fig_07_annotation_confusion.csv": ("注释混淆矩阵", "自动注释与发布注释的交叉计数。"),
    "annotation_confusion.csv": ("注释混淆矩阵", "自动注释与发布注释的交叉计数。"),
    "fig_08_deg_all.csv": ("全部差异表达基因", "所有基因的差异表达统计，包括 log2FC、P 值和校正 P 值。"),
    "deg_all.csv": ("全部差异表达基因", "所有基因的差异表达统计，包括 log2FC、P 值和校正 P 值。"),
    "fig_09_deg_significant.csv": ("显著差异表达基因", "达到显著性阈值的差异基因列表。"),
    "deg_significant.csv": ("显著差异表达基因", "达到显著性阈值的差异基因列表。"),
    "fig_09_deg_horizontal_violin.csv": ("Top DEG 横向小提琴图数据", "横向小提琴图对应的基因、差异倍数和校正 P 值。"),
    "fig_10_enrichment_up_go.csv": ("上调基因 GO 富集表", "上调基因的 GO 生物过程富集结果。"),
    "enrichment_up_go.csv": ("上调基因 GO 富集表", "上调基因的 GO 生物过程富集结果。"),
    "fig_11_enrichment_down_go.csv": ("下调基因 GO 富集表", "下调基因的 GO 生物过程富集结果。"),
    "enrichment_down_go.csv": ("下调基因 GO 富集表", "下调基因的 GO 生物过程富集结果。"),
    "fig_12_enrichment_up_kegg.csv": ("上调基因 KEGG 富集表", "上调基因的 KEGG 通路富集结果。"),
    "enrichment_up_kegg.csv": ("上调基因 KEGG 富集表", "上调基因的 KEGG 通路富集结果。"),
    "fig_13_enrichment_down_kegg.csv": ("下调基因 KEGG 富集表", "下调基因的 KEGG 通路富集结果。"),
    "enrichment_down_kegg.csv": ("下调基因 KEGG 富集表", "下调基因的 KEGG 通路富集结果。"),
    "fig_18_19_celltype_proportion_stats.csv": ("细胞类型比例统计", "细胞类型比例的统计和组间比较结果。"),
    "fig_26_27_cell_cycle_scores.csv": ("细胞周期打分", "每细胞的 S 和 G2M 期打分。"),
    "cell_cycle_scores.csv": ("细胞周期打分", "每细胞的 S 和 G2M 期打分。"),
    "fig_29_doublet_rate_by_sample.csv": ("样本双细胞率", "每个样本的双细胞数量、比例和去双细胞前后细胞数。"),
    "doublet_rate_by_sample.csv": ("样本双细胞率", "每个样本的双细胞数量、比例和去双细胞前后细胞数。"),
    "fig_16_17_31_32_cluster_markers.csv": ("聚类 marker 基因", "每个聚类的 marker 基因及统计信息。"),
    "cluster_markers.csv": ("聚类 marker 基因", "每个聚类的 marker 基因及统计信息。"),
    "fig_33_34_35_signature_scores.csv": ("功能签名打分", "增殖、EMT、缺氧等签名在每细胞中的得分。"),
    "signature_scores.csv": ("功能签名打分", "增殖、EMT、缺氧等签名在每细胞中的得分。"),
    "fig_36_cnv_heatmap.csv": ("推断 CNV 热图矩阵", "基于染色体窗口均值推断的 CNV 矩阵。"),
    "cnv_heatmap.csv": ("推断 CNV 热图矩阵", "基于染色体窗口均值推断的 CNV 矩阵。"),
    "fig_37_singleR_annotations.csv": ("SingleR 注释表", "SingleR 参考注释和打分。"),
    "singleR_annotations.csv": ("SingleR 注释表", "SingleR 参考注释和打分。"),
    "fig_38_singleR_confusion.csv": ("SingleR 混淆矩阵", "SingleR 注释与现有注释的交叉计数。"),
    "singleR_confusion.csv": ("SingleR 混淆矩阵", "SingleR 注释与现有注释的交叉计数。"),
    "fig_39_trajectory_pseudotime.csv": ("拟时序结果", "每细胞的拟时序值和轨迹信息。"),
    "fig_40_cellchat_communication.csv": ("CellChat 通讯矩阵", "细胞类型之间的通讯数量或强度。"),
    "fig_40_cellchat_communication_weight.csv": ("CellChat 通讯权重", "细胞类型之间的通讯权重。"),
    "fig_42_cellchat_pathways.csv": ("CellChat 通路表", "CellChat 识别的配体受体通路。"),
    "fig_24_ml_feature_importance.csv": ("ML 特征重要性", "分类模型的特征重要性排序。"),
    "ml_feature_importance.csv": ("ML 特征重要性", "分类模型的特征重要性排序。"),
    "fig_43_44_45_ml_classification_results.csv": ("ML 分类预测", "每个样本或细胞的预测结果。"),
    "ml_classification_results.csv": ("ML 分类预测", "每个样本或细胞的预测结果。"),
    "fig_43_44_45_ml_classification_report.csv": ("ML 分类报告", "分类模型的精确率、召回率、F1 等指标。"),
    "ml_classification_report.csv": ("ML 分类报告", "分类模型的精确率、召回率、F1 等指标。"),
    "sample_annotations.csv": ("样本分组表", "样本与疾病分组/条件的对应关系。"),
    "liver_cancer_seurat.rds": ("Seurat 对象", "完整单细胞 Seurat 对象，可用于下游重分析和复现。"),
}

FIGURE_DATA_MAP = {
    "fig_01_qc_raw_violin.png": ["fig_01_qc_metrics.csv", "qc_metrics.csv"],
    "fig_01_qc_filtered_violin.png": ["fig_01_qc_metrics.csv", "qc_metrics.csv"],
    "fig_48_qc_pvalue_comparison.png": ["fig_48_qc_pvalue_comparison.csv"],
    "fig_02_doublet_scores.png": [
        "fig_02_doublet_results.csv",
        "doublet_results.csv",
        "fig_29_doublet_rate_by_sample.csv",
        "doublet_rate_by_sample.csv",
    ],
    "fig_03_umap_clusters.png": [
        "fig_03_04_05_umap_coordinates.csv",
        "umap_coordinates.csv",
    ],
    "fig_04_umap_condition.png": [
        "fig_03_04_05_umap_coordinates.csv",
        "umap_coordinates.csv",
    ],
    "fig_05_umap_annotation.png": [
        "fig_03_04_05_umap_coordinates.csv",
        "umap_coordinates.csv",
        "fig_05_16_17_cell_annotations.csv",
        "cell_annotations.csv",
        "fig_07_annotation_confusion.csv",
        "annotation_confusion.csv",
    ],
    "fig_06_dotplot_markers.png": [
        "fig_05_16_17_cell_annotations.csv",
        "cell_annotations.csv",
        "fig_16_17_31_32_cluster_markers.csv",
        "cluster_markers.csv",
    ],
    "fig_07_annotation_confusion_heatmap.png": [
        "fig_07_annotation_confusion.csv",
        "annotation_confusion.csv",
        "fig_38_singleR_confusion.csv",
        "singleR_confusion.csv",
    ],
    "fig_08_volcano.png": [
        "fig_08_deg_all.csv",
        "deg_all.csv",
        "fig_09_deg_significant.csv",
        "deg_significant.csv",
    ],
    "fig_09_deg_heatmap.png": [
        "fig_09_deg_significant.csv",
        "deg_significant.csv",
        "fig_09_deg_horizontal_violin.csv",
    ],
    "fig_09_deg_horizontal_violin.png": [
        "fig_09_deg_horizontal_violin.csv",
        "fig_09_deg_significant.csv",
        "deg_significant.csv",
    ],
    "fig_10_go_up.png": ["fig_10_enrichment_up_go.csv", "enrichment_up_go.csv"],
    "fig_11_go_down.png": [
        "fig_11_enrichment_down_go.csv",
        "enrichment_down_go.csv",
    ],
    "fig_12_kegg_up.png": ["fig_12_enrichment_up_kegg.csv", "enrichment_up_kegg.csv"],
    "fig_13_kegg_down.png": [
        "fig_13_enrichment_down_kegg.csv",
        "enrichment_down_kegg.csv",
    ],
    "fig_14_pca.png": [
        "fig_03_04_05_umap_coordinates.csv",
        "umap_coordinates.csv",
        "sample_annotations.csv",
    ],
    "fig_15_elbow.png": ["sample_annotations.csv"],
    "fig_16_featureplot_markers.png": [
        "fig_05_16_17_cell_annotations.csv",
        "cell_annotations.csv",
        "fig_16_17_31_32_cluster_markers.csv",
        "cluster_markers.csv",
    ],
    "fig_17_marker_violin.png": [
        "fig_05_16_17_cell_annotations.csv",
        "cell_annotations.csv",
        "fig_16_17_31_32_cluster_markers.csv",
        "cluster_markers.csv",
    ],
    "fig_18_celltype_proportion.png": [
        "fig_18_19_celltype_proportion_stats.csv",
        "fig_18_19_30_cluster_composition.csv",
        "cluster_composition.csv",
        "fig_05_16_17_cell_annotations.csv",
        "cell_annotations.csv",
    ],
    "fig_19_condition_proportion.png": [
        "fig_18_19_celltype_proportion_stats.csv",
        "fig_18_19_30_cluster_composition.csv",
        "cluster_composition.csv",
    ],
    "fig_20_gsea_go.png": [
        "fig_20_gsea_go.csv",
        "gsea_go.csv",
        "fig_08_deg_all.csv",
        "deg_all.csv",
    ],
    "fig_21_gsea_kegg.png": [
        "fig_21_gsea_kegg.csv",
        "gsea_kegg.csv",
        "fig_08_deg_all.csv",
        "deg_all.csv",
    ],
    "fig_22_go_network.png": [
        "fig_10_enrichment_up_go.csv",
        "enrichment_up_go.csv",
    ],
    "fig_23_kegg_network.png": [
        "fig_12_enrichment_up_kegg.csv",
        "enrichment_up_kegg.csv",
    ],
    "fig_24_ml_feature_importance.png": [
        "fig_24_ml_feature_importance.csv",
        "ml_feature_importance.csv",
        "ml_model_summary.json",
    ],
    "fig_25_ml_shap.png": [
        "ml_model_summary.json",
        "fig_43_44_45_ml_classification_results.csv",
        "ml_classification_results.csv",
    ],
    "fig_26_cellcycle_umap.png": [
        "fig_26_27_cell_cycle_scores.csv",
        "cell_cycle_scores.csv",
    ],
    "fig_27_cellcycle_proportion.png": [
        "fig_26_27_cell_cycle_scores.csv",
        "cell_cycle_scores.csv",
    ],
    "fig_28_umap_sample.png": [
        "fig_03_04_05_umap_coordinates.csv",
        "umap_coordinates.csv",
        "sample_annotations.csv",
    ],
    "fig_29_doublet_rate_sample.png": [
        "fig_29_doublet_rate_by_sample.csv",
        "doublet_rate_by_sample.csv",
    ],
    "fig_30_sample_proportion.png": [
        "fig_18_19_30_cluster_composition.csv",
        "cluster_composition.csv",
        "sample_annotations.csv",
    ],
    "fig_31_cluster_marker_heatmap.png": [
        "fig_16_17_31_32_cluster_markers.csv",
        "cluster_markers.csv",
    ],
    "fig_32_cluster_marker_dotplot.png": [
        "fig_16_17_31_32_cluster_markers.csv",
        "cluster_markers.csv",
    ],
    "fig_33_signature_scores_umap.png": [
        "fig_33_34_35_signature_scores.csv",
        "signature_scores.csv",
    ],
    "fig_34_signature_scores_boxplot.png": [
        "fig_33_34_35_signature_scores.csv",
        "signature_scores.csv",
    ],
    "fig_35_celltype_abundance_effect.png": [
        "differential_abundance.csv",
        "fig_33_34_35_signature_scores.csv",
        "signature_scores.csv",
    ],
    "fig_36_cnv_heatmap.png": ["fig_36_cnv_heatmap.csv", "cnv_heatmap.csv"],
    "fig_37_singler_umap.png": [
        "fig_37_singleR_annotations.csv",
        "singleR_annotations.csv",
    ],
    "fig_38_singler_confusion_heatmap.png": [
        "fig_38_singleR_confusion.csv",
        "singleR_confusion.csv",
    ],
    "fig_39_trajectory_umap.png": ["fig_39_trajectory_pseudotime.csv"],
    "fig_40_cellchat_network.png": [
        "fig_40_cellchat_communication.csv",
        "fig_40_cellchat_communication_weight.csv",
    ],
    "fig_41_cellchat_heatmap.png": [
        "fig_40_cellchat_communication.csv",
        "fig_40_cellchat_communication_weight.csv",
    ],
    "fig_42_cellchat_bubble.png": ["fig_42_cellchat_pathways.csv"],
    "fig_43_ml_confusion_matrix.png": [
        "fig_43_44_45_ml_classification_report.csv",
        "ml_classification_report.csv",
        "fig_43_44_45_ml_classification_results.csv",
        "ml_classification_results.csv",
        "ml_model_summary.json",
    ],
    "fig_44_ml_roc_pr.png": [
        "fig_43_44_45_ml_classification_report.csv",
        "ml_classification_report.csv",
        "fig_43_44_45_ml_classification_results.csv",
        "ml_classification_results.csv",
        "ml_model_summary.json",
    ],
    "fig_45_ml_cv_scores.png": [
        "fig_43_44_45_ml_classification_report.csv",
        "ml_classification_report.csv",
        "fig_43_44_45_ml_classification_results.csv",
        "ml_classification_results.csv",
        "ml_model_summary.json",
    ],
    "fig_45_ml_calibration_curve.png": [
        "fig_43_44_45_ml_classification_report.csv",
        "ml_classification_report.csv",
    ],
    "fig_46_go_top5.png": ["fig_10_enrichment_up_go.csv", "enrichment_up_go.csv"],
    "fig_47_kegg_top5.png": ["fig_12_enrichment_up_kegg.csv", "enrichment_up_kegg.csv"],
}


def esc(value) -> str:
    return html.escape(str(value))


def read_table(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def find_result(base: Path, name: str) -> Path | None:
    direct = base / name
    if direct.is_file():
        return direct
    for path in base.rglob(name):
        if path.is_file():
            return path
    return None


def read_named_table(name: str) -> list[dict]:
    path = find_result(DATA, name)
    return read_table(path) if path is not None else []


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
    path = find_result(FIG, name)
    if path is None:
        return ""
    rel = path.relative_to(RES).as_posix()
    return (
        f"<figure><img src='{rel}' alt='{esc(title)}'>"
        f"<figcaption>{esc(title)}</figcaption></figure>"
    )


def data_link(name: str, label: str) -> str:
    path = find_result(DATA, name)
    if path is None:
        return ""
    rel = path.relative_to(RES).as_posix()
    return f"<li><a href='{rel}'>{esc(label)}</a></li>"


def fmt_size(num: float) -> str:
    value = float(num)
    for unit in ("B", "KB", "MB", "GB"):
        if value < 1024 or unit == "GB":
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GB"


def fmt_num(value: float, digits: int = 4) -> str:
    if value == 0:
        return "0"
    if abs(value) < 1e-4 or abs(value) >= 1e6:
        return f"{value:.{digits}e}"
    return f"{value:.{digits}g}"


def iter_result_files():
    for base, kind in ((FIG, "figure"), (DATA, "data")):
        if not base.exists():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file():
                yield path.relative_to(RES).as_posix(), kind, path


def stage_label(rel: str) -> str:
    parts = Path(rel).parts
    if len(parts) >= 2 and parts[0] in ("figures", "data"):
        stage = STAGE_LABELS.get(parts[1])
        if stage:
            return stage
    match = re.match(r"^fig_(\d+)", Path(rel).name)
    if match:
        num = int(match.group(1))
        if num <= 1 or num == 48:
            return STAGE_LABELS["01_qc"]
        if num == 2:
            return STAGE_LABELS["02_doublets"]
        if num in (3, 4, 14, 15):
            return STAGE_LABELS["03_cluster"]
        if num in (5, 6, 7, 16, 17, 18, 19):
            return STAGE_LABELS["04_annotation"]
        if num in (8, 9):
            return STAGE_LABELS["05_deg"]
        if num in (10, 11, 12, 13, 20, 21, 22, 23, 46, 47):
            return STAGE_LABELS["06_enrichment"]
        if num in (24, 25, 43, 44, 45):
            return STAGE_LABELS["07_ml"]
        if 26 <= num <= 39:
            return STAGE_LABELS["08_publication"]
        if 40 <= num <= 42:
            return STAGE_LABELS["09_cellchat"]
    return "其他"


def figure_title(name: str) -> str:
    stem = Path(name).stem
    stem = re.sub(r"^fig_\d+_", "", stem)
    stem = stem.replace("_", " ").replace("-", " ").strip()
    return stem.title()


def analyze_delimited(path: Path, delimiter: str = ",") -> dict:
    findings = []
    preview = []
    total_rows = 0
    fieldnames = []
    numeric_values = {}
    counts = Counter()
    sample_values = {}

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        fieldnames = [str(c or "") for c in (reader.fieldnames or [])]
        lookup = {c.lower(): c for c in fieldnames}
        numeric_fields = fieldnames[:MAX_NUMERIC_COLS]
        numeric_values = {c: [] for c in numeric_fields}

        def col(*aliases):
            for alias in aliases:
                if alias.lower() in lookup:
                    return lookup[alias.lower()]
            return None

        p_col = col("p_val_adj", "padj", "fdr", "p.adjust", "pvalue", "p_value", "p.value", "PValue")
        fc_col = col("avg_log2FC", "log2FoldChange", "log2FC", "logFC")
        sig_col = col("significant")
        direction_col = col("direction", "change")
        sample_col = col("sample", "sample_id", "Sample")
        condition_col = col("condition", "group", "Group", "Condition")
        cluster_col = col("cluster", "seurat_clusters", "Cluster")
        celltype_col = col(
            "celltype", "cell_type", "celltype_annot",
            "celltype_annot_cell", "CellType",
        )
        doublet_col = col("doublet_call", "DF.classifications")
        desc_col = col("Description", "description")

        if sample_col:
            sample_values["sample"] = set()
        if condition_col:
            sample_values["condition"] = set()
        if cluster_col:
            sample_values["cluster"] = set()
        if celltype_col:
            sample_values["celltype"] = set()

        p_values = []
        fc_values = []
        for row in reader:
            total_rows += 1
            if len(preview) < MAX_PREVIEW_ROWS:
                preview.append(
                    {c: row.get(c, "") for c in fieldnames[:MAX_PREVIEW_COLS]}
                )
            for c in numeric_fields:
                try:
                    numeric_values[c].append(float(row.get(c)))
                except (TypeError, ValueError):
                    pass
            if p_col:
                try:
                    p_values.append(float(row.get(p_col)))
                except (TypeError, ValueError):
                    pass
            if fc_col:
                try:
                    fc_values.append(float(row.get(fc_col)))
                except (TypeError, ValueError):
                    pass
            if sig_col:
                counts[("significant", str(row.get(sig_col)).strip().lower())] += 1
            if direction_col:
                counts[("direction", str(row.get(direction_col)).strip().lower())] += 1
            if sample_col:
                value = row.get(sample_col)
                if value not in (None, ""):
                    sample_values["sample"].add(value)
            if condition_col:
                value = row.get(condition_col)
                if value not in (None, ""):
                    sample_values["condition"].add(value)
            if cluster_col:
                value = row.get(cluster_col)
                if value not in (None, ""):
                    sample_values["cluster"].add(value)
            if celltype_col:
                value = row.get(celltype_col)
                if value not in (None, ""):
                    sample_values["celltype"].add(value)
            if doublet_col:
                value = row.get(doublet_col)
                if value not in (None, ""):
                    counts[("doublet", str(value).strip().lower())] += 1

    if not fieldnames:
        findings.append("文件为空或没有表头。")
    else:
        findings.append(f"表格共 {total_rows} 行、{len(fieldnames)} 列。")
        preview_cols = "、".join(fieldnames[:10])
        if len(fieldnames) > 10:
            preview_cols += " 等"
        findings.append(f"主要字段：{preview_cols}。")

    for c in numeric_fields[:MAX_STAT_COLS]:
        values = numeric_values[c]
        if len(values) >= 2:
            findings.append(
                f"{c}：有效值 {len(values)}，均值 {fmt_num(statistics.mean(values))}，"
                f"范围 {fmt_num(min(values))}–{fmt_num(max(values))}。"
            )

    if p_values:
        lt_005 = sum(v < 0.05 for v in p_values)
        lt_001 = sum(v < 0.01 for v in p_values)
        findings.append(
            f"P 值/校正 P 值：有效值 {len(p_values)}，<0.05 共 {lt_005} 个，"
            f"<0.01 共 {lt_001} 个。"
        )
    if fc_values:
        up = sum(v > 0 for v in fc_values)
        down = sum(v < 0 for v in fc_values)
        findings.append(
            f"log2FC：上调（>0）{up} 个，下调（<0）{down} 个，"
            f"范围 {fmt_num(min(fc_values))}–{fmt_num(max(fc_values))}。"
        )
    if counts["significant"]:
        sig_parts = "，".join(f"{k}={v}" for k, v in counts["significant"].items())
        findings.append(f"显著性标记：{sig_parts}。")
    if counts["direction"]:
        dir_parts = "，".join(f"{k}={v}" for k, v in counts["direction"].items())
        findings.append(f"方向标记：{dir_parts}。")
    if "sample" in sample_values:
        findings.append(f"样本数：{len(sample_values['sample'])}。")
    if "condition" in sample_values:
        findings.append(f"分组数：{len(sample_values['condition'])}。")
    if "cluster" in sample_values:
        findings.append(f"聚类数：{len(sample_values['cluster'])}。")
    if "celltype" in sample_values:
        findings.append(f"细胞类型数：{len(sample_values['celltype'])}。")
    if counts["doublet"]:
        dt_parts = "，".join(f"{k}={v}" for k, v in counts["doublet"].items())
        findings.append(f"双细胞分类：{dt_parts}。")
    if desc_col and total_rows:
        findings.append(f"富集条目数：{total_rows} 条。")

    return {
        "kind": "data",
        "title": "",
        "rows": total_rows,
        "columns": fieldnames,
        "findings": findings,
        "preview": preview,
        "size": fmt_size(path.stat().st_size),
    }


def analyze_json(path: Path) -> dict:
    findings = []
    try:
        data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        findings.append("JSON 文件无法解析，建议人工检查格式。")
        data = None
    if isinstance(data, dict):
        findings.append(f"JSON 对象包含 {len(data)} 个顶层键。")
        for key, value in list(data.items())[:20]:
            if isinstance(value, (str, int, float, bool)) or value is None:
                findings.append(f"{key} = {value}")
            elif isinstance(value, list):
                findings.append(f"{key}：列表，共 {len(value)} 项。")
            elif isinstance(value, dict):
                findings.append(f"{key}：嵌套对象，共 {len(value)} 项。")
    elif isinstance(data, list):
        findings.append(f"JSON 数组共 {len(data)} 项。")
        if data and isinstance(data[0], dict):
            findings.append("首项字段：" + "、".join(list(data[0])[:10]))
    return {
        "kind": "data",
        "title": "",
        "rows": None,
        "columns": [],
        "findings": findings,
        "preview": [],
        "size": fmt_size(path.stat().st_size),
    }


def analyze_data_file(path: Path, rel: str) -> dict:
    name = path.name
    guide = DATA_GUIDE.get(name)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        result = analyze_delimited(path, ",")
    elif suffix == ".tsv":
        result = analyze_delimited(path, "\t")
    elif suffix == ".json":
        result = analyze_json(path)
    else:
        result = {
            "kind": "data",
            "title": "",
            "rows": None,
            "columns": [],
            "findings": [f"结果数据文件，格式为 {suffix or '未知'}。"],
            "preview": [],
            "size": fmt_size(path.stat().st_size),
        }
        if suffix == ".rds":
            result["findings"].append(
                "该文件是 R 序列化对象，建议用 readRDS() 加载后检查。"
            )
    if guide:
        result["title"] = guide[0]
        result["findings"].insert(0, guide[1])
    else:
        result["title"] = figure_title(name)
    result["rel"] = rel
    result["stage"] = stage_label(rel)
    return result


def analyze_figure(path: Path, rel: str) -> dict:
    name = path.name
    guide = FIGURE_GUIDE.get(name)
    findings = []
    if guide:
        title, description = guide
        findings.append(description)
    else:
        title = figure_title(name)
    findings.append(f"文件大小：{fmt_size(path.stat().st_size)}。")
    try:
        from PIL import Image
        with Image.open(path) as img:
            findings.append(f"图像尺寸：{img.width} × {img.height} 像素。")
    except Exception:
        pass
    lower = name.lower()
    if "volcano" in lower:
        findings.append("重点观察显著上调和下调基因的数量、分布对称性和离群点。")
    elif "umap" in lower or "tsne" in lower or "pca" in lower:
        findings.append("重点检查分群边界、样本混合程度和条件之间的分布偏移。")
    elif "heatmap" in lower or "dotplot" in lower or "violin" in lower:
        findings.append("重点检查表达模式是否与细胞类型或差异方向一致。")
    elif "enrich" in lower or "gsea" in lower or "kegg" in lower or "go_" in lower:
        findings.append("重点查看显著通路的富集方向和基因数。")
    elif "ml_" in lower or "shap" in lower:
        findings.append("重点检查模型性能指标、特征重要性和分类稳定性。")
    findings.append(f"所属阶段：{stage_label(rel)}。")
    return {
        "kind": "figure",
        "title": title,
        "rows": None,
        "columns": [],
        "findings": findings,
        "preview": [],
        "size": fmt_size(path.stat().st_size),
        "rel": rel,
        "stage": stage_label(rel),
    }


def companion_data_analyses(fig_name: str, analyses: list[dict]) -> list[dict]:
    aliases = FIGURE_DATA_MAP.get(fig_name, [])
    by_basename = {}
    for analysis in analyses:
        if analysis["kind"] == "data":
            by_basename.setdefault(Path(analysis["rel"]).name, analysis)
    matched = []
    seen = set()
    for alias in aliases:
        analysis = by_basename.get(alias)
        if analysis is not None and analysis["rel"] not in seen:
            matched.append(analysis)
            seen.add(analysis["rel"])
    return matched


def data_findings_brief(analysis: dict, limit: int = 6) -> list[str]:
    findings = list(analysis.get("findings", []))
    guide = DATA_GUIDE.get(Path(analysis["rel"]).name)
    if guide and findings and findings[0] == guide[1]:
        findings = findings[1:]
    return findings[:limit]


def num_str(value) -> str:
    try:
        return fmt_num(float(value))
    except (TypeError, ValueError):
        return str(value or "")


def top_table_summary(name: str, n: int = 3) -> str:
    rows = read_named_table(name)
    if not rows:
        alt = re.sub(r"^fig_\d+_", "", name)
        if alt != name:
            rows = read_named_table(alt)
    if not rows:
        return ""
    first = rows[0]
    label_col = next(
        (c for c in ("Description", "ID", "pathway", "gene") if c in first),
        None,
    )
    p_col = next((c for c in ("p.adjust", "qvalue", "pvalue") if c in first), None)
    if not label_col:
        return ""
    ranked = rows
    if p_col:
        try:
            ranked = sorted(rows, key=lambda r: float(r.get(p_col, 1)))[:n]
        except (TypeError, ValueError):
            ranked = rows[:n]
    parts = []
    for row in ranked:
        label = row.get(label_col, "")
        p = num_str(row.get(p_col, "")) if p_col else ""
        if p_col and p:
            parts.append(f"{label}（{p_col}={p}）")
        else:
            parts.append(str(label))
    return "；".join(parts)


def ml_metrics_line() -> str:
    rows = read_named_table("fig_43_44_45_ml_classification_report.csv")
    if not rows:
        rows = read_named_table("ml_classification_report.csv")
    if not rows:
        return ""
    parts = []
    for row in rows[:8]:
        label = str(
            row.get("label") or row.get("class") or row.get("") or ""
        )
        metrics = []
        for key in ("precision", "recall", "f1-score", "f1", "accuracy", "support"):
            if key in row:
                try:
                    value = fmt_num(float(row.get(key)))
                except (TypeError, ValueError):
                    value = str(row.get(key))
                metrics.append(f"{key}={value}")
        if metrics:
            parts.append(f"{label or '指标'}: " + "，".join(metrics))
    return "；".join(parts[:4])


def top_deg_line(summary: dict) -> str:
    top = summary.get("top_degs", [])
    if isinstance(top, dict):
        top = [top]
    if not top:
        return ""
    parts = []
    for item in top[:3]:
        gene = item.get("gene", "")
        fc = num_str(item.get("avg_log2FC", ""))
        p = num_str(item.get("p_val_adj", ""))
        parts.append(f"{gene}（log2FC={fc}, padj={p}）")
    return "Top 差异基因：" + "；".join(parts)


def joint_conclusions(
    fig_name: str,
    companions: list[dict],
    summary: dict,
) -> list[str]:
    items = []
    if fig_name in ("fig_01_qc_raw_violin.png", "fig_01_qc_filtered_violin.png"):
        raw = summary.get("n_cells_raw")
        qc = summary.get("n_cells_after_qc")
        if isinstance(raw, (int, float)) and isinstance(qc, (int, float)) and raw:
            items.append(
                f"与 QC 汇总一致：原始 {raw} 细胞，QC 后保留 {qc} 细胞"
                f"（{qc / raw * 100:.1f}%）。"
            )
    if fig_name == "fig_02_doublet_scores.png":
        qc = summary.get("n_cells_after_qc")
        dbl = summary.get("n_cells_after_doublet_removal")
        if isinstance(qc, (int, float)) and isinstance(dbl, (int, float)) and qc:
            items.append(
                f"去双细胞后保留 {dbl} / {qc} 细胞（{dbl / qc * 100:.1f}%）。"
            )
    if fig_name == "fig_08_volcano.png" and summary.get("deg_total") is not None:
        items.append(
            f"差异表达汇总：共 {summary.get('deg_total')} 个，"
            f"上调 {summary.get('deg_up', 0)} 个，下调 {summary.get('deg_down', 0)} 个。"
        )
    if fig_name in (
        "fig_08_volcano.png",
        "fig_09_deg_heatmap.png",
        "fig_09_deg_horizontal_violin.png",
    ):
        line = top_deg_line(summary)
        if line:
            items.append(line)
    enrichment_names = {
        "fig_10_go_up.png": "fig_10_enrichment_up_go.csv",
        "fig_11_go_down.png": "fig_11_enrichment_down_go.csv",
        "fig_12_kegg_up.png": "fig_12_enrichment_up_kegg.csv",
        "fig_13_kegg_down.png": "fig_13_enrichment_down_kegg.csv",
        "fig_20_gsea_go.png": "fig_20_gsea_go.csv",
        "fig_21_gsea_kegg.png": "fig_21_gsea_kegg.csv",
        "fig_22_go_network.png": "fig_10_enrichment_up_go.csv",
        "fig_23_kegg_network.png": "fig_12_enrichment_up_kegg.csv",
        "fig_46_go_top5.png": "fig_10_enrichment_up_go.csv",
        "fig_47_kegg_top5.png": "fig_12_enrichment_up_kegg.csv",
    }
    if fig_name in enrichment_names:
        terms = top_table_summary(enrichment_names[fig_name], 3)
        if terms:
            items.append(f"对应富集数据 Top 通路：{terms}。")
    if fig_name.startswith(("fig_24", "fig_25", "fig_43", "fig_44", "fig_45")):
        metrics = ml_metrics_line()
        if metrics:
            items.append(f"分类指标：{metrics}。")
    if companions:
        names = "、".join(Path(a["rel"]).name for a in companions[:5])
        more = " 等" if len(companions) > 5 else ""
        items.append(
            f"已联合 {len(companions)} 个结果数据文件：{names}{more}。"
        )
    else:
        items.append("未在结果数据目录中找到该图的对应数据文件，当前仅能做图面检查。")
    return items


def render_joint_analysis_block(
    fig_analysis: dict,
    companions: list[dict],
    summary: dict,
) -> str:
    rel = fig_analysis["rel"]
    title = fig_analysis.get("title", "")
    figure_points = "".join(
        f"<li>{esc(item)}</li>" for item in fig_analysis.get("findings", [])[:4]
    )
    data_blocks = ""
    for companion in companions:
        brief = data_findings_brief(companion, 6)
        findings = "".join(f"<li>{esc(item)}</li>" for item in brief)
        data_blocks += (
            "<div class='companion-data'>"
            f"<p class='companion-title'><a href='{esc(companion['rel'])}'>"
            f"{esc(Path(companion['rel']).name)}</a>"
            f"<span class='file-size'>{esc(companion.get('size', ''))}</span></p>"
            f"<ul class='findings'>{findings}</ul></div>"
        )
    conclusion_html = "".join(
        f"<li>{esc(item)}</li>"
        for item in joint_conclusions(Path(rel).name, companions, summary)
    )
    data_html = data_blocks or "<p class='muted'>未找到对应数据文件。</p>"
    return (
        "<div class='joint-analysis'>"
        f"<h3>{esc(title)}</h3>"
        f"<p class='file-head'><span class='badge-fig'>结果图</span>"
        f"<code>{esc(rel)}</code>"
        f"<span class='file-size'>{esc(fig_analysis.get('size', ''))}</span></p>"
        f"<figure class='analysis-thumb'><img src='{esc(rel)}' alt='{esc(title)}'>"
        f"<figcaption>{esc(title)}</figcaption></figure>"
        f"<h4>图面要点</h4><ul class='findings'>{figure_points}</ul>"
        f"<h4>对应结果数据</h4>{data_html}"
        f"<h4>联合结论</h4><ul class='findings'>{conclusion_html}</ul>"
        "</div>"
    )


def render_joint_section(analyses: list[dict], summary: dict) -> str:
    blocks = []
    for analysis in analyses:
        if analysis["kind"] != "figure":
            continue
        companions = companion_data_analyses(Path(analysis["rel"]).name, analyses)
        blocks.append(render_joint_analysis_block(analysis, companions, summary))
    if not blocks:
        return "<p class='muted'>未发现结果图，无法生成联合分析。</p>"
    return "".join(blocks)


def render_markdown_joint(analyses: list[dict], summary: dict) -> str:
    figure_analyses = [a for a in analyses if a["kind"] == "figure"]
    lines = [
        "# 结果图与结果数据联合分析报告",
        "",
        f"- 数据集：{summary.get('dataset', '')}",
        f"- 标题：{summary.get('title', '')}",
        f"- 完成时间：{summary.get('finished_at', '')}",
        "",
    ]
    if not figure_analyses:
        lines.append("未发现结果图，无法生成联合分析。")
        return "\n".join(lines)
    lines.append(
        f"本次共分析 {len(figure_analyses)} 张实际输出的结果图，"
        "并逐图匹配对应结果数据。"
    )
    lines.append("")
    for analysis in figure_analyses:
        rel = analysis["rel"]
        title = analysis.get("title", Path(rel).stem)
        lines.append(f"## {title}")
        lines.append("")
        lines.append(f"![{title}]({rel})")
        lines.append("")
        lines.append(f"- 结果图：`{rel}`，大小 {analysis.get('size', '')}。")
        for item in analysis.get("findings", [])[:4]:
            lines.append(f"- {item}")
        companions = companion_data_analyses(Path(rel).name, analyses)
        if companions:
            lines.append("")
            lines.append("### 对应结果数据")
            lines.append("")
            for companion in companions:
                lines.append(
                    f"- **{Path(companion['rel']).name}**"
                    f"（{companion.get('size', '')}）："
                )
                for item in data_findings_brief(companion, 4):
                    lines.append(f"  - {item}")
        lines.append("")
        lines.append("### 联合结论")
        lines.append("")
        for item in joint_conclusions(Path(rel).name, companions, summary):
            lines.append(f"- {item}")
        lines.append("")
    lines.append("---")
    lines.append(
        "说明：联合结论中的数值来自 results/data 下的实际数据文件；"
        "图面要点仅描述图像内容与检查重点，最终结论请以原始数据为准。"
    )
    return "\n".join(lines)


def build_analysis_json(analyses: list[dict], summary: dict) -> dict:
    figures = []
    for analysis in analyses:
        if analysis["kind"] != "figure":
            continue
        companions = companion_data_analyses(Path(analysis["rel"]).name, analyses)
        figures.append(
            {
                "file": analysis["rel"],
                "title": analysis.get("title", ""),
                "size": analysis.get("size", ""),
                "findings": analysis.get("findings", []),
                "companion_data": [
                    {
                        "file": companion["rel"],
                        "title": companion.get("title", ""),
                        "findings": data_findings_brief(companion, 10),
                    }
                    for companion in companions
                ],
                "conclusions": joint_conclusions(
                    Path(analysis["rel"]).name, companions, summary
                ),
            }
        )
    return {"summary": summary, "figures": figures}


def render_file_analysis_block(analysis: dict) -> str:
    rel = analysis["rel"]
    kind = analysis["kind"]
    badge = "结果图" if kind == "figure" else "结果数据"
    badge_class = "badge-fig" if kind == "figure" else "badge-data"
    findings = "".join(
        f"<li>{esc(item)}</li>" for item in analysis.get("findings", [])
    )
    preview_html = ""
    if analysis.get("preview"):
        preview_html = (
            "<details><summary>数据预览</summary>"
            + render_table(analysis["preview"], MAX_PREVIEW_ROWS)
            + "</details>"
        )
    return (
        "<div class='file-analysis'>"
        f"<div class='file-head'><span class='{badge_class}'>{badge}</span>"
        f"<code>{esc(rel)}</code>"
        f"<span class='file-size'>{esc(analysis.get('size', ''))}</span></div>"
        f"<p class='file-title'>{esc(analysis.get('title', ''))}</p>"
        f"<ul class='findings'>{findings}</ul>"
        f"{preview_html}</div>"
    )


def render_analysis_overview(analyses: list[dict]) -> str:
    rows = []
    for analysis in analyses:
        summary = analysis.get("findings", [""])[0]
        rows.append(
            {
                "文件": analysis["rel"],
                "类型": "结果图" if analysis["kind"] == "figure" else "结果数据",
                "大小": analysis.get("size", ""),
                "分析摘要": summary,
            }
        )
    return render_table(rows, max_rows=300)


def render_overall_conclusion(
    summary: dict,
    analyses: list[dict],
    go_up: list[dict],
    go_down: list[dict],
    kegg_up: list[dict],
    kegg_down: list[dict],
) -> str:
    items = []
    raw = summary.get("n_cells_raw")
    after_qc = summary.get("n_cells_after_qc")
    after_doublet = summary.get("n_cells_after_doublet_removal")
    if isinstance(raw, (int, float)) and isinstance(after_qc, (int, float)):
        items.append(f"QC 保留率：{after_qc} / {raw} = {after_qc / raw * 100:.1f}%。")
    if isinstance(after_qc, (int, float)) and isinstance(after_doublet, (int, float)):
        items.append(
            f"去双细胞后保留：{after_doublet} / {after_qc} = "
            f"{after_doublet / after_qc * 100:.1f}%。"
        )
    if summary.get("n_clusters") is not None:
        items.append(f"共识别 {summary.get('n_clusters')} 个聚类、{summary.get('n_celltypes', 'NA')} 类细胞。")
    if summary.get("deg_total") is not None:
        items.append(
            f"差异表达：共 {summary.get('deg_total')} 个，上调 "
            f"{summary.get('deg_up', 0)} 个，下调 {summary.get('deg_down', 0)} 个。"
        )
    if go_up:
        items.append(f"上调基因 GO 富集 {len(go_up)} 条。")
    if go_down:
        items.append(f"下调基因 GO 富集 {len(go_down)} 条。")
    if kegg_up:
        items.append(f"上调基因 KEGG 富集 {len(kegg_up)} 条。")
    if kegg_down:
        items.append(f"下调基因 KEGG 富集 {len(kegg_down)} 条。")
    figure_count = sum(a["kind"] == "figure" for a in analyses)
    data_count = sum(a["kind"] == "data" for a in analyses)
    items.append(f"逐文件分析覆盖：{figure_count} 张结果图、{data_count} 个结果数据文件。")
    return "".join(f"<li>{esc(item)}</li>" for item in items)


def main() -> int:
    if len(sys.argv) > 1:
        target = Path(sys.argv[1]).expanduser().resolve()
        if target.is_dir():
            if (target / "results").is_dir():
                target = target / "results"
            global RES, FIG, DATA
            RES = target
            FIG = RES / "figures"
            DATA = RES / "data"

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

    go_up = read_named_table("fig_10_enrichment_up_go.csv")
    go_down = read_named_table("fig_11_enrichment_down_go.csv")
    kegg_up = read_named_table("fig_12_enrichment_up_kegg.csv")
    kegg_down = read_named_table("fig_13_enrichment_down_kegg.csv")

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
        data_link("fig_08_deg_all.csv", "全部差异表达基因"),
        data_link("fig_09_deg_significant.csv", "显著差异表达基因"),
        data_link("fig_09_deg_horizontal_violin.csv", "Top DEG 横向小提琴图数据"),
        data_link("fig_01_qc_metrics.csv", "QC 指标"),
        data_link("fig_48_qc_pvalue_comparison.csv", "QC 质控差异度 P 值"),
        data_link("fig_02_doublet_results.csv", "双细胞结果"),
        data_link("fig_05_16_17_cell_annotations.csv", "细胞注释"),
        data_link("fig_07_annotation_confusion.csv", "注释混淆矩阵"),
        data_link("fig_03_04_05_umap_coordinates.csv", "UMAP 坐标"),
        data_link("fig_18_19_30_cluster_composition.csv", "聚类组成"),
        data_link("fig_10_enrichment_up_go.csv", "上调 GO 富集"),
        data_link("fig_11_enrichment_down_go.csv", "下调 GO 富集"),
        data_link("fig_12_enrichment_up_kegg.csv", "上调 KEGG 富集"),
        data_link("fig_13_enrichment_down_kegg.csv", "下调 KEGG 富集"),
        data_link("fig_26_27_cell_cycle_scores.csv", "细胞周期打分"),
        data_link("fig_29_doublet_rate_by_sample.csv", "样本双细胞率"),
        data_link("fig_16_17_31_32_cluster_markers.csv", "聚类 marker 基因"),
        data_link("fig_33_34_35_signature_scores.csv", "功能签名打分"),
        data_link("fig_36_cnv_heatmap.csv", "推断 CNV 热图矩阵"),
        data_link("fig_37_singleR_annotations.csv", "SingleR 注释"),
        data_link("fig_38_singleR_confusion.csv", "SingleR 混淆矩阵"),
        data_link("fig_39_trajectory_pseudotime.csv", "拟时序结果"),
        data_link("fig_40_cellchat_communication.csv", "CellChat 通讯矩阵"),
        data_link("fig_40_cellchat_communication_weight.csv", "CellChat 通讯权重"),
        data_link("fig_42_cellchat_pathways.csv", "CellChat 通路表"),
        data_link("fig_24_ml_feature_importance.csv", "ML 特征重要性"),
        data_link("fig_43_44_45_ml_classification_results.csv", "ML 分类预测"),
        data_link("fig_43_44_45_ml_classification_report.csv", "ML 分类报告"),
        data_link("sample_annotations.csv", "样本分组"),
        data_link("liver_cancer_seurat.rds", "Seurat 对象"),
    ])

    analyses = [
        analyze_figure(path, rel)
        if kind == "figure"
        else analyze_data_file(path, rel)
        for rel, kind, path in iter_result_files()
    ]
    analysis_overview = render_analysis_overview(analyses)
    analysis_blocks = "".join(
        render_file_analysis_block(analysis) for analysis in analyses
    )
    joint_section = render_joint_section(analyses, summary)
    markdown_report = render_markdown_joint(analyses, summary)
    json_report = build_analysis_json(analyses, summary)
    conclusion_html = render_overall_conclusion(
        summary, analyses, go_up, go_down, kegg_up, kegg_down
    )

    html_doc = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>肝癌单细胞分析总报告</title>
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
.file-analysis {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 14px 16px; margin: 12px 0; }}
.file-head {{ display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }}
.file-head code {{ background: #eef2f7; padding: 3px 7px; border-radius: 4px; }}
.file-size {{ color: #7b8794; font-size: 12px; margin-left: auto; }}
.file-title {{ font-weight: 600; margin: 10px 0 6px; }}
.findings {{ margin: 6px 0; padding-left: 20px; color: #334155; }}
.badge-fig, .badge-data {{ display: inline-block; border-radius: 12px; padding: 3px 9px; font-size: 12px; color: #fff; }}
.badge-fig {{ background: #1665c0; }}
.badge-data {{ background: #0e9f6e; }}
.joint-analysis {{ background: #fff; border: 1px solid #d9e2ec; border-left: 4px solid #1665c0; border-radius: 8px; padding: 16px 18px; margin: 16px 0; }}
.joint-analysis h3 {{ margin-top: 0; }}
.joint-analysis h4 {{ margin: 14px 0 6px; color: #334155; }}
.companion-data {{ background: #f7fafc; border: 1px solid #e4e7eb; border-radius: 6px; padding: 10px 12px; margin: 8px 0; }}
.companion-title {{ margin: 0 0 4px; font-weight: 600; }}
.analysis-thumb img {{ max-width: 640px; }}
details {{ margin-top: 8px; }}
summary {{ cursor: pointer; color: #1665c0; font-size: 13px; }}
</style>
</head>
<body>
<div class="wrap">
<h1>单细胞分析总报告</h1>
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
{image_card('fig_48_qc_pvalue_comparison.png', 'QC 质控差异度 P 值')}

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
{image_card('fig_09_deg_horizontal_violin.png', 'Top DEG 横向小提琴图（P 值）')}
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
{image_card('fig_22_go_network.png', 'GO BP 通路网络（Top5 核心 + 5 延伸）')}
{image_card('fig_23_kegg_network.png', 'KEGG 通路网络（Top5 核心 + 5 延伸）')}
{image_card('fig_46_go_top5.png', 'GO BP 筛选后 Top5')}
{image_card('fig_47_kegg_top5.png', 'KEGG 筛选后 Top5')}

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

<h2>11. 结果图与结果数据联合分析</h2>
{joint_section}

<h2>12. 逐文件分析报告</h2>
<h3>文件总览</h3>
{analysis_overview}
<h3>总体结论</h3>
<ul>{conclusion_html}</ul>
{analysis_blocks}

<footer>Pipeline: R / Seurat / scDblFinder / clusterProfiler | 自动生成于 {esc(summary.get('finished_at', ''))}</footer>
</div>
</body>
</html>
"""

    (RES / "result_report.html").write_text(html_doc, encoding="utf-8")
    (RES / "result_analysis_report.md").write_text(markdown_report, encoding="utf-8")
    (RES / "result_analysis.json").write_text(
        json.dumps(json_report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print("HTML report written:", RES / "result_report.html")
    print("Joint analysis report written:", RES / "result_analysis_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
