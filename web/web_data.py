"""Static lookup tables shared by the web console.

Kept in a separate module so ``web_ui.py`` stays focused on request
handling; every name here is a literal table with no runtime state.
"""

from __future__ import annotations


NAV_HTML = (
    '<div class="topnav">'
    '<a href="/full">全自动流水线</a>'
    '<a href="/environment">环境补全</a>'
    '<a href="/">表达分析</a>'
    '<a href="/datasets">数据集搜索</a>'
    '<a href="/dock">虚拟筛选</a>'
    '<a href="/md-simulation">分子动力学</a>'
    '<a href="/molecular-docking">分子对接</a>'
    '<a href="/knockout">虚拟敲除</a>'
    '<a href="/network">网络毒理学</a>'
    '<a href="/faers">FAERS</a>'
    '<a href="/validation">真实数据验证</a>'
    '<a href="/results">结果清单</a>'
    '<a href="/guide">使用教程</a>'
    '<a href="/tasks" class="nav-right">任务进度</a>'
    '</div>'
)
NAV_CSS = (
    ".topnav{position:sticky;top:0;z-index:100;"
    "background:#0f172a;padding:12px 28px;box-shadow:0 2px 8px rgba(15,23,42,.35);"
    "display:flex;gap:18px;align-items:center;flex-wrap:wrap;}"
    ".topnav a{color:#ffffff;text-decoration:none;font-size:15px;"
    "font-weight:600;padding:6px 10px;border-radius:6px;"
    "background:rgba(255,255,255,.08);}"
    ".topnav a:hover,.topnav a.active{background:#1665c0;color:#fff;}"
    ".topnav a .nav-count{display:inline-flex;align-items:center;"
    "justify-content:center;min-width:20px;height:20px;margin-left:6px;"
    "padding:0 6px;border-radius:999px;background:#f59e0b;color:#fff;"
    "font-size:12px;font-weight:700;line-height:1;}"
    ".topnav a .nav-count[hidden]{display:none;}"
    ".topnav .nav-right{margin-left:auto;}"
)

ENV_MODULES = {
    "expression": {
        "title": "表达分析",
        "summary": (
            "单细胞 / bulk RNA-seq / microarray 表达分析依赖，包含 Python 与 R 包"
            "（Seurat、DESeq2、clusterProfiler 等）。"
        ),
        "r_deps": True,
        "dock_tools": False,
        "skills": False,
        "install_bat": "launchers/install_expression_environment.bat",
        "check_bat": "launchers/check_expression_environment.bat",
        "note": "Windows 上未找到 Rscript 时会自动下载安装 R。",
    },
    "datasets": {
        "title": "数据集搜索",
        "summary": "GEO / BioStudies / Expression Atlas 搜索、筛选和下载所需的 Python 依赖。",
        "r_deps": False,
        "dock_tools": False,
        "skills": False,
        "install_bat": "launchers/install_datasets_environment.bat",
        "check_bat": None,
        "note": "自动安装 numpy、pandas、scikit-learn 与 joblib 等搜索/重排依赖。",
    },
    "docking": {
        "title": "虚拟筛选 / 对接",
        "summary": (
            "AutoDock Vina 虚拟筛选环境：RDKit、Meeko、Open Babel、"
            "AutoDockTools、Vina，以及网络毒理学 Cytoscape 推送依赖 "
            "py4cytoscape。"
        ),
        "r_deps": False,
        "dock_tools": True,
        "skills": False,
        "install_bat": "launchers/install_docking_environment.bat",
        "check_bat": "launchers/check_docking_environment.bat",
        "note": (
            "缺失的 AutoDockTools / AutoDock Vina 会下载到 dock/tools/；"
            "Cytoscape 桌面版需单独安装并用 -R 1234 启动 CyREST。"
        ),
    },
    "molecular-docking": {
        "title": "独立分子对接",
        "summary": "独立分子对接板块的 Python 与对接工具依赖，与虚拟筛选工作目录分开。",
        "r_deps": False,
        "dock_tools": True,
        "skills": False,
        "install_bat": "launchers/install_molecular_docking_environment.bat",
        "check_bat": "launchers/check_molecular_docking_environment.bat",
        "note": "复用 dock/tools/ 下的 AutoDockTools 与 Vina。",
    },
    "md": {
        "title": "分子动力学",
        "summary": "GROMACS 分子动力学准备和模拟所需依赖；补齐对接工具，便于处理蛋白与配体。",
        "r_deps": False,
        "dock_tools": True,
        "skills": False,
        "install_bat": "launchers/install_md_environment.bat",
        "check_bat": "launchers/check_md_environment.bat",
        "note": "GROMACS gmx 需要按系统单独安装。",
    },
    "full": {
        "title": "全自动集成流水线",
        "summary": "一次补齐表达分析、虚拟筛选、网页版与项目 Codex Skills，适合完整流程使用。",
        "r_deps": True,
        "dock_tools": True,
        "skills": True,
        "install_bat": "launchers/install_full_environment.bat",
        "check_bat": "launchers/check_full_environment.bat",
        "note": "安装耗时较长；也可以先按单个板块补齐。",
    },
    "web": {
        "title": "网页版",
        "summary": "检查 Python 版本与网页入口；网页本身不自动下载 R 或大型对接工具。",
        "r_deps": False,
        "dock_tools": False,
        "skills": False,
        "install_bat": "launchers/install_web_environment.bat",
        "check_bat": None,
        "note": "主要确认 web/web_ui.py 与本机 Python 可用。",
    },
    "skills": {
        "title": "项目 Codex Skills",
        "summary": "安装 liver-expression-analysis、liver-virtual-screening、liver-full-pipeline 与 liver-dataset-search。",
        "r_deps": False,
        "dock_tools": False,
        "skills": True,
        "install_bat": "launchers/install_codex_skills_environment.bat",
        "check_bat": None,
        "note": "Skill 只指导 Codex 调用项目脚本，不复制核心分析代码。",
    },
}

FIGURES = [
    {"file": "fig_01_qc_raw_violin.png", "label": "QC 小提琴图（原始）"},
    {"file": "fig_01_qc_filtered_violin.png", "label": "QC 小提琴图（过滤后）"},
    {"file": "fig_48_qc_pvalue_comparison.png", "label": "QC 质控差异度 P 值"},
    {"file": "fig_49_qc_umi_feature_correlation.png", "label": "UMI 与基因数关系 QC 图"},
    {"file": "fig_02_doublet_scores.png", "label": "双细胞得分图"},
    {"file": "fig_03_umap_clusters.png", "label": "UMAP 聚类图"},
    {"file": "fig_04_umap_condition.png", "label": "UMAP 分组图"},
    {"file": "fig_05_umap_annotation.png", "label": "UMAP 注释图"},
    {"file": "fig_06_dotplot_markers.png", "label": "Marker 基因 DotPlot"},
    {"file": "fig_07_annotation_confusion_heatmap.png", "label": "注释混淆矩阵热图"},
    {"file": "fig_08_volcano.png", "label": "差异表达图", "styles": ["volcano", "maplot"]},
    {"file": "fig_09_deg_heatmap.png", "label": "Top DEG 热图"},
    {"file": "fig_09_deg_horizontal_violin.png", "label": "Top DEG 横向小提琴图（P 值）"},
    {"file": "fig_10_go_up.png", "label": "GO BP 富集图（上调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_11_go_down.png", "label": "GO BP 富集图（下调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_12_kegg_up.png", "label": "KEGG 富集图（上调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_13_kegg_down.png", "label": "KEGG 富集图（下调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_14_pca.png", "label": "PCA 分组图"},
    {"file": "fig_15_elbow.png", "label": "主成分 Elbow 图"},
    {"file": "fig_16_featureplot_markers.png", "label": "Marker 基因 FeaturePlot"},
    {"file": "fig_17_marker_violin.png", "label": "Marker 基因小提琴图"},
    {"file": "fig_50_marker_ridgeplot.png", "label": "Marker 基因峰峦图"},
    {"file": "fig_51_marker_stacked_violin.png", "label": "Marker 基因堆叠小提琴图"},
    {"file": "fig_18_celltype_proportion.png", "label": "细胞类型比例堆叠图"},
    {"file": "fig_19_condition_proportion.png", "label": "分组构成比例图"},
    {"file": "fig_20_gsea_go.png", "label": "GSEA GO BP 富集图", "styles": ["ridgeplot", "gseaplot2"]},
    {"file": "fig_21_gsea_kegg.png", "label": "GSEA KEGG 富集图", "styles": ["ridgeplot", "gseaplot2"]},
    {"file": "fig_22_go_network.png", "label": "GO BP 通路网络图（Top5 核心 + 5 延伸）", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_23_kegg_network.png", "label": "KEGG 通路网络图（Top5 核心 + 5 延伸）", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_46_go_top5.png", "label": "GO BP 筛选后 Top5", "styles": ["dotplot", "barplot", "cnetplot", "emapplot"]},
    {"file": "fig_47_kegg_top5.png", "label": "KEGG 筛选后 Top5", "styles": ["dotplot", "barplot", "cnetplot", "emapplot"]},
    {"file": "fig_24_ml_feature_importance.png", "label": "ML 特征重要性图"},
    {"file": "fig_25_ml_shap.png", "label": "SHAP 可解释性图"},
    {"file": "fig_26_cellcycle_umap.png", "label": "细胞周期 UMAP"},
    {"file": "fig_27_cellcycle_proportion.png", "label": "细胞周期比例图"},
    {"file": "fig_28_umap_sample.png", "label": "UMAP 按样本"},
    {"file": "fig_29_doublet_rate_sample.png", "label": "样本双细胞率图"},
    {"file": "fig_30_sample_proportion.png", "label": "样本细胞类型比例图"},
    {"file": "fig_31_cluster_marker_heatmap.png", "label": "聚类 Marker 热图"},
    {"file": "fig_32_cluster_marker_dotplot.png", "label": "聚类 Marker DotPlot"},
    {"file": "fig_33_signature_scores_umap.png", "label": "功能签名 UMAP"},
    {"file": "fig_34_signature_scores_boxplot.png", "label": "功能签名箱线图"},
    {"file": "fig_35_celltype_abundance_effect.png", "label": "细胞类型丰度变化图"},
    {"file": "fig_36_cnv_heatmap.png", "label": "推断 CNV 热图"},
    {"file": "fig_37_singler_umap.png", "label": "SingleR 注释 UMAP"},
    {"file": "fig_38_singler_confusion_heatmap.png", "label": "SingleR 混淆矩阵热图"},
    {"file": "fig_39_trajectory_umap.png", "label": "拟时序轨迹图"},
    {"file": "fig_40_cellchat_network.png", "label": "CellChat 通讯网络图"},
    {"file": "fig_41_cellchat_heatmap.png", "label": "CellChat 通讯热图"},
    {"file": "fig_42_cellchat_bubble.png", "label": "CellChat 配体受体气泡图"},
    {"file": "fig_43_ml_confusion_matrix.png", "label": "ML 混淆矩阵"},
    {"file": "fig_44_ml_roc_pr.png", "label": "ML ROC 与 PR 曲线"},
    {"file": "fig_45_ml_cv_scores.png", "label": "ML 交叉验证得分图"},
    {"file": "fig_45_ml_calibration_curve.png", "label": "ML 校准曲线"},
]
FIGURE_NAMES = [item["file"] for item in FIGURES]

STYLE_LABELS = {
    "volcano": "火山图",
    "maplot": "MA 图",
    "dotplot": "气泡图",
    "barplot": "柱状图",
    "cnetplot": "通路网络图",
    "ridgeplot": "峰峦图",
    "gseaplot2": "GSEA 富集曲线",
    "emapplot": "富集关系网络图",
}

SOFTWARE = [
    {"name": "R", "url": "https://www.r-project.org/", "kind": "software"},
    {"name": "Python", "url": "https://www.python.org/", "kind": "software"},
    {"name": "Seurat", "url": "https://satijalab.org/seurat/", "kind": "package", "install": "install.packages('Seurat')"},
    {"name": "scDblFinder", "url": "https://bioconductor.org/packages/scDblFinder/", "kind": "package", "install": "BiocManager::install('scDblFinder')"},
    {"name": "SingleCellExperiment", "url": "https://bioconductor.org/packages/SingleCellExperiment/", "kind": "package", "install": "BiocManager::install('SingleCellExperiment')"},
    {"name": "clusterProfiler", "url": "https://bioconductor.org/packages/clusterProfiler/", "kind": "package", "install": "BiocManager::install('clusterProfiler')"},
    {"name": "enrichplot", "url": "https://bioconductor.org/packages/enrichplot/", "kind": "package", "install": "BiocManager::install('enrichplot')"},
    {"name": "BiocParallel", "url": "https://bioconductor.org/packages/BiocParallel/", "kind": "package", "install": "BiocManager::install('BiocParallel')"},
    {"name": "org.Hs.eg.db", "url": "https://bioconductor.org/packages/org.Hs.eg.db/", "kind": "package", "install": "BiocManager::install('org.Hs.eg.db')"},
    {"name": "org.Mm.eg.db", "url": "https://bioconductor.org/packages/org.Mm.eg.db/", "kind": "package", "install": "BiocManager::install('org.Mm.eg.db')"},
    {"name": "DESeq2", "url": "https://bioconductor.org/packages/DESeq2/", "kind": "package", "install": "BiocManager::install('DESeq2')"},
    {"name": "data.table", "url": "https://rdatatable.gitlab.io/data.table/", "kind": "package", "install": "install.packages('data.table')"},
    {"name": "ggplot2", "url": "https://ggplot2.tidyverse.org/", "kind": "package", "install": "install.packages('ggplot2')"},
]

SINGLE_STAGE_LABELS = {
    "01": "数据加载",
    "02": "QC 过滤",
    "03": "双细胞检测",
    "04": "聚类",
    "05": "细胞注释",
    "06": "差异表达",
    "07": "富集分析",
    "08": "发表级分析",
    "09": "汇总输出",
}

DOCK_STAGE_LABELS = {
    "01": "受体准备",
    "02": "配体准备",
    "03": "分子对接",
    "04": "结果分析",
    "05": "精修重对接",
    "06": "HTML 报告",
}

MOLECULAR_DOCK_STAGE_LABELS = DOCK_STAGE_LABELS

FULL_STAGE_LABELS = {
    "01": "表达分析",
    "02": "关键基因",
    "03": "证据富集",
    "04": "敲除输入",
    "05": "虚拟敲除",
    "06": "分子对接",
    "07": "CADD 下游",
    "08": "网络毒理学",
    "09": "FAERS",
    "10": "细胞反馈",
    "11": "集成报告",
}

RESULT_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".svg"}
RESULT_DATA_SUFFIXES = {".csv", ".xlsx", ".rds", ".xgmml"}
RESULT_FILE_SUFFIXES = RESULT_IMAGE_SUFFIXES | RESULT_DATA_SUFFIXES
