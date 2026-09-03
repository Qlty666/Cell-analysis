# 结果图清单（内容 / 类型 / 用途）

本清单覆盖当前代码中实际存在的结果图文件名，按单细胞分析、虚拟筛选、集成与细胞反馈、网络毒理学分组。仅列结果图，不包含同名 CSV/JSON 数据文件。

可选性说明：`fig_16`、`fig_17` 仅数据集中存在可用 marker 基因时生成；`fig_26`、`fig_27` 受 `LIVER_RUN_CELLCYCLE` 控制（默认 yes，且需成功写入 `Phase`）；`fig_28` 仅多样本（样本数 > 1）时生成；`fig_31`、`fig_32` 受 `LIVER_RUN_CLUSTER_MARKERS` 控制（默认 yes，且需 `FindAllMarkers` 返回非空）；`fig_33`、`fig_34` 受 `LIVER_RUN_SIGNATURES` 控制（默认 yes，且需签名基因充足）；`fig_35` 依赖 `fig_18_19_celltype_proportion_stats.csv` 存在；`fig_36` 受 `LIVER_RUN_CNV` 控制（默认 yes，且需注释库与 CHRLOC 充足）；`fig_37`、`fig_38` 受 `LIVER_RUN_SINGLER` 控制（默认 yes，且需参考数据可用、预测成功）；`fig_39` 受 `LIVER_RUN_TRAJECTORY` 控制（默认 no，且需 `slingshot` 可用）；`fig_40` 至 `fig_42` 受 `LIVER_RUN_CELLCHAT` 控制（默认 no，且需 R/CellChat 可用、存在通讯）；`fig_45_ml_calibration_curve.png` 仅二分类模型成功计算时生成；`fig_24_ml_selected_features.csv` 仅 `LIVER_ML_MODEL=lasso_svm` 时生成；`fig_49_redock_comparison.png` 仅重对接阶段启用且初始/重对接亲和力可合并时生成；`fig_50`/`fig_51` 仅 `ml-train` 阶段生成（`fig_50` 需模型有 `feature_importances_`/`coef_`，`fig_51` 需二分类任务）；`ppi_hub_scores.csv` 仅提供 `--ppi-network-csv` 时生成。

## 一、单细胞分析结果图

对应目录：`<output>/results/figures/`。

### 1.1 质量控制与双细胞检测

01_fig_01_qc_raw_violin.png
内容: 过滤前 nFeature_RNA、nCount_RNA、percent.mt 按条件展示的原始分布。
类型: 小提琴图（多面板）。
用途: 判断原始数据质量，识别低质量细胞尾部、条件间差异和异常样本。

02_fig_01_qc_filtered_violin.png
内容: 过滤后保留细胞的 nFeature_RNA、nCount_RNA、percent.mt 分布。
类型: 小提琴图（多面板）。
用途: 确认过滤阈值是否合理，验证低质量尾部被去除且主要细胞群保留。

03_fig_48_qc_pvalue_comparison.png
内容: 过滤前后各 QC 指标在条件间的 Wilcoxon 检验 -log10 P 值对比，标注 P=0.05 参考线。
类型: 分组条形图（按原始/过滤后面板分面）。
用途: 量化过滤是否改变条件间质量差异，辅助判断 QC 是否引入偏倚。

04_fig_02_doublet_scores.png
内容: scDblFinder 双细胞得分按 singlet/doublet 分类的分布。
类型: 小提琴图叠加箱线图。
用途: 判断双细胞检测的区分度与过滤边界，识别异常高得分的双细胞候选。

### 1.2 聚类、降维与注释

05_fig_14_pca.png
内容: 主成分空间中细胞按条件着色。
类型: PCA 散点图。
用途: 查看主要变异来源、条件/批次分离和离群细胞。

06_fig_15_elbow.png
内容: 各主成分标准差。
类型: 折线图（碎石图）。
用途: 辅助确定下游使用的 PCA 维度数。

07_fig_03_umap_clusters.png
内容: Seurat 聚类结果在 UMAP 上的分布。
类型: UMAP 散点图。
用途: 检查聚类结构、聚类边界和过度/欠聚类。

08_fig_04_umap_condition.png
内容: UMAP 按条件着色。
类型: UMAP 散点图。
用途: 评估条件间分布、批次效应和分组混杂。

09_fig_05_umap_annotation.png
内容: marker 注释与发表注释并排 UMAP。
类型: UMAP 散点图（双面板）。
用途: 比较两种注释的细胞类型分布是否一致。

10_fig_06_dotplot_markers.png
内容: marker 基因在各细胞类型中的表达比例和平均表达量。
类型: 点图（DotPlot）。
用途: 验证细胞类型注释的 marker 特异性。

11_fig_07_annotation_confusion_heatmap.png
内容: marker 注释与发表注释的混淆矩阵。
类型: 热图（带计数）。
用途: 评估注释一致性和主要对应关系。

12_fig_16_featureplot_markers.png
内容: marker 基因在 UMAP 上的表达位置。
类型: UMAP 特征图（FeaturePlot）。
用途: 检查 marker 表达是否集中在预期细胞类型区域。

13_fig_17_marker_violin.png
内容: marker 基因在各细胞类型中的表达分布。
类型: 小提琴图（多面板）。
用途: 验证 marker 在对应细胞类型中表达更高。

14_fig_18_celltype_proportion.png
内容: 各细胞类型在条件中的比例堆叠。
类型: 堆叠条形图。
用途: 比较条件间细胞类型构成差异。

15_fig_19_condition_proportion.png
内容: 每个细胞类型中条件构成比例。
类型: 堆叠条形图（横向）。
用途: 观察每个细胞类型内部的条件组成。

16_fig_28_umap_sample.png
内容: UMAP 按样本着色。
类型: UMAP 散点图。
用途: 检查多样本混合和批次效应。

17_fig_29_doublet_rate_sample.png
内容: 每个样本的双细胞率。
类型: 条形图。
用途: 识别双细胞率异常偏高的样本。

18_fig_30_sample_proportion.png
内容: 每样本各细胞类型比例。
类型: 堆叠条形图。
用途: 比较样本间细胞类型构成和异质性。

19_fig_31_cluster_marker_heatmap.png
内容: 各聚类 marker 基因平均表达热图。
类型: 热图。
用途: 查看聚类 marker 表达模式，辅助定义聚类标签。

20_fig_32_cluster_marker_dotplot.png
内容: 各聚类 marker 基因表达比例与水平。
类型: 点图（DotPlot）。
用途: 检查 marker 是否具有聚类特异性。

21_fig_33_signature_scores_umap.png
内容: 增殖、EMT、缺氧等功能签名得分在 UMAP 上的分布。
类型: UMAP 特征图（多面板）。
用途: 查看功能签名是否集中在特定细胞群。

22_fig_34_signature_scores_boxplot.png
内容: 功能签名得分按条件/细胞类型分布。
类型: 箱线图或小提琴图。
用途: 比较条件间签名活性差异。

23_fig_35_celltype_abundance_effect.png
内容: 细胞类型丰度变化的 log2OR 与 -log10 校正 P 值。
类型: 散点图（类火山图）。
用途: 判断哪些细胞类型在条件间显著增多或减少。

24_fig_36_cnv_heatmap.png
内容: 基于染色体滑动窗口均值的推断 CNV 模式。
类型: 热图。
用途: 查看细胞亚群的染色体拷贝数变化和肿瘤异质性。

25_fig_37_singler_umap.png
内容: SingleR 参考注释在 UMAP 上的分布。
类型: UMAP 散点图。
用途: 查看参考注释与现有注释的空间一致性。

26_fig_38_singler_confusion_heatmap.png
内容: SingleR 注释与现有注释的混淆矩阵。
类型: 热图。
用途: 定量评估注释一致性。

27_fig_39_trajectory_umap.png
内容: slingshot 拟时序轨迹及细胞 UMAP 分布。
类型: UMAP 轨迹图。
用途: 查看分化轨迹、伪时序梯度和谱系结构。

28_fig_40_cellchat_network.png
内容: 细胞类型间的通讯数量网络。
类型: 网络图。
用途: 查看细胞间通讯关系和强度。

29_fig_41_cellchat_heatmap.png
内容: 细胞类型对间的通讯强度。
类型: 热图。
用途: 比较不同细胞类型对之间的通讯强度。

30_fig_42_cellchat_bubble.png
内容: 配体受体对在细胞类型对中的通讯证据。
类型: 气泡图。
用途: 筛选候选配体受体互作。

31_fig_26_cellcycle_umap.png
内容: 细胞周期阶段在 UMAP 上的分布。
类型: UMAP 散点图。
用途: 检查细胞周期是否干扰聚类结构。

32_fig_27_cellcycle_proportion.png
内容: 细胞周期阶段按分组/细胞类型的比例。
类型: 堆叠条形图。
用途: 比较细胞周期构成差异。

### 1.3 差异表达与富集分析

33_fig_08_volcano.png
内容: 差异基因的 log2FC 与 -log10 校正 P 值。
类型: 火山图（可回退为 MA 图）。
用途: 查看差异表达基因数量、方向和显著性。

34_fig_09_deg_horizontal_violin.png
内容: 按校正 P 值排序的 Top DEG 表达分布并标注 P 值。
类型: 横向小提琴图。
用途: 优先查看差异最显著基因的表达差异。

35_fig_09_deg_heatmap.png
内容: Top DEG 在条件/样本中的表达热图。
类型: 热图。
用途: 查看显著基因的整体表达模式。

36_fig_10_go_up.png
内容: 上调基因 GO BP 富集结果。
类型: 富集气泡图/点图。
用途: 查看上调基因富集的生物学过程。

37_fig_11_go_down.png
内容: 下调基因 GO BP 富集结果。
类型: 富集气泡图/点图。
用途: 查看下调基因富集的生物学过程。

38_fig_12_kegg_up.png
内容: 上调基因 KEGG 富集结果。
类型: 富集气泡图/点图。
用途: 查看上调基因相关 KEGG 通路。

39_fig_13_kegg_down.png
内容: 下调基因 KEGG 富集结果。
类型: 富集气泡图/点图。
用途: 查看下调基因相关 KEGG 通路。

40_fig_20_gsea_go.png
内容: GO BP GSEA 富集结果。
类型: GSEA 富集图（折线/条带组合）。
用途: 查看 GO 通路整体富集方向和显著性。

41_fig_21_gsea_kegg.png
内容: KEGG GSEA 富集结果。
类型: GSEA 富集图（折线/条带组合）。
用途: 查看 KEGG 通路整体富集方向和显著性。

42_fig_22_go_network.png
内容: 上调 GO 通路经校正 P 值筛选后的 Top5 核心 + 最佳 5 延伸通路-基因网络。
类型: 网络图（cnetplot）。
用途: 查看通路与基因的关联结构。

43_fig_23_kegg_network.png
内容: 上调 KEGG 通路筛选后的 Top5 核心 + 最佳 5 延伸通路-基因网络。
类型: 网络图（cnetplot）。
用途: 查看 KEGG 通路与基因的关联结构。

44_fig_46_go_top5.png
内容: GO BP 筛选后 Top5 通路富集。
类型: 气泡图。
用途: 快速查看最显著的上调 GO 通路。

45_fig_47_kegg_top5.png
内容: KEGG 筛选后 Top5 通路富集。
类型: 气泡图。
用途: 快速查看最显著的上调 KEGG 通路。

### 1.4 机器学习与模型诊断

46_fig_24_ml_feature_importance.png
内容: 单细胞 ML 分类模型的特征重要性。
类型: 条形图（横向排序）。
用途: 筛选影响样本分类的关键特征。

47_fig_25_ml_shap.png
内容: SHAP 值对预测方向和贡献的解释。
类型: SHAP 蜂群图。
用途: 解释特征如何影响模型预测。

48_fig_43_ml_confusion_matrix.png
内容: 交叉验证分类混淆矩阵。
类型: 热图/混淆矩阵图。
用途: 查看模型错分模式和分类准确度。

49_fig_44_ml_roc_pr.png
内容: 二分类或多分类 ROC 与 PR 曲线。
类型: 折线图（双面板）。
用途: 评估模型区分能力和类别性能。

50_fig_45_ml_cv_scores.png
内容: 交叉验证准确率分布。
类型: 箱线图叠加散点图。
用途: 评估模型稳定性。

51_fig_45_ml_calibration_curve.png
内容: 二分类模型预测概率与观测频率的校准曲线。
类型: 折线图（带完美校准参考线）。
用途: 检查预测概率是否校准，辅助判断概率阈值可靠性。

## 二、虚拟筛选结果图

对应目录：`<workdir>/outputs/run_001/results/`。

52_fig_46_affinity_distribution.png
内容: 全部成功对接配体的亲和力分布及命中阈值线。
类型: 直方图/密度图。
用途: 查看整体对接分数、命中分布和阈值合理性。

53_fig_47_top_hits.png
内容: 排序后的 Top 命中或候选配体亲和力。
类型: 条形图。
用途: 快速选择排名靠前的候选配体。

54_fig_48_diverse_hits.png
内容: Tanimoto 多样性选择后的命中候选。
类型: 条形图。
用途: 减少同一化学骨架重复，提高候选多样性。

55_fig_49_redock_comparison.png
内容: 初始亲和力与重对接亲和力比较及变化量分布。
类型: 散点图 + 直方图（双面板）。
用途: 检查初筛排序和亲和力的稳定性。

56_fig_50_ml_feature_importance.png
内容: 对接 ML 重打分模型的特征重要性。
类型: 条形图。
用途: 解释哪些分子特征影响重打分结果。

57_fig_51_ml_roc.png
内容: 对接 ML 分类模型的 ROC 曲线。
类型: 折线图。
用途: 评估分类模型排序能力。

58_fig_52_knockout_top_candidates.png
内容: Top N 候选基因的 knockout_score。
类型: 条形图。
用途: 快速查看高优先级候选基因。

59_fig_53_knockout_score_distribution.png
内容: 全部基因 knockout_score 分布。
类型: 直方图。
用途: 判断评分是否有区分度和分层能力。

60_md_rmsd_rmsf.png
内容: 生产轨迹的蛋白骨架 RMSD、配体 RMSD 与配体原子 RMSF。
类型: 折线图（三面板）。
用途: 检查蛋白-配体复合物在 MD 中的收敛与配体构象稳定性。

## 三、集成与细胞反馈结果图

对应目录：`<workdir>/outputs/integration/cell_feedback/figures/`。

60_fig_54_feedback_module_umap.png
内容: 筛选靶点模块得分在 UMAP 上的分布。
类型: UMAP 特征图。
用途: 查看候选靶点模块是否集中在特定细胞群。

61_fig_55_feedback_target_expression_umap.png
内容: 候选靶基因在 UMAP 上的表达。
类型: UMAP 特征图（多面板）。
用途: 查看单个候选基因的细胞类型表达位置。

62_fig_56_feedback_celltype_dotplot.png
内容: 候选靶基因在细胞类型中的表达比例和水平。
类型: 点图（DotPlot）。
用途: 比较候选基因的细胞类型特异性。

63_fig_57_feedback_celltype_boxplot.png
内容: 模块得分按细胞类型和条件分布。
类型: 箱线图/小提琴图。
用途: 查看模块活性的细胞类型和条件差异。

64_fig_58_feedback_celltype_heatmap.png
内容: 候选基因平均表达按细胞类型聚类。
类型: 热图。
用途: 查看候选基因表达谱的细胞类型模式。

65_fig_59_feedback_targets_volcano.png
内容: 反馈靶基因在条件间的差异表达。
类型: 火山图。
用途: 查看反馈靶点中显著上调/下调基因及其幅度。

66_fig_60_feedback_condition_violin.png
内容: 反馈靶基因按条件分组的表达分布。
类型: 小提琴图。
用途: 查看反馈靶点在病例/对照间的表达差异。

67_fig_61_feedback_go_network.png
内容: 反馈靶基因 GO BP 富集 Top5 通路-基因网络。
类型: cnetplot 通路-基因网络图。
用途: 直接查看 Top5 GO 通路与反馈靶基因的关联强度。

68_fig_62_feedback_kegg_network.png
内容: 反馈靶基因 KEGG 富集 Top5 通路-基因网络。
类型: cnetplot 通路-基因网络图。
用途: 直接查看 Top5 KEGG 通路与反馈靶基因的关联强度。

## 四、网络毒理学结果图

对应目录：`<workdir>/outputs/run_001/network_toxicology/`。

69_compound_disease_venn.png
内容: 化合物靶点与疾病基因的交集。
类型: Venn 图。
用途: 查看交集规模及来源，确定核心分析基因集。

70_ctpd_network.html
内容: 化合物-靶点-通路-疾病网络可视化。
类型: 交互式网络图（HTML）。
用途: 查看化合物到疾病的多层关联路径。
