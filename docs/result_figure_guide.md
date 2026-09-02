# 现有流程结果图判断与使用指南

本文档用于说明本项目现有流水线实际能够生成哪些结果图、如何判断结果图是否合格、是否可以用于结论，以及每张图的主要用途。

适用范围：

- 单细胞分析：`scripts/run_pipeline.py`，结果图位于 `<output_root>/results/figures/`。
- 虚拟筛选：`scripts/run_docking.py`，结果图位于 `<workdir>/outputs/run_001/results/` 下各阶段目录。
- 全自动集成流水线：`scripts/run_full_pipeline.py`，除复用上述两类结果外，还会在 `<workdir>/outputs/integration/cell_feedback/figures/` 生成细胞反馈结果图。

说明：README 中“48 张”是历史常见口径。按当前代码逐点核对，单细胞流程的图片保存点共有 50 个文件名，其中部分图片由环境变量、数据条件或可选阶段决定是否生成。判断时应以本清单、实际输出目录、阶段状态文件和对应数据文件为准。

每次报告生成时，流水线会读取 `results/figures` 下实际存在的每张结果图，按同编号或同阶段文件自动匹配 `results/data` 中的结果数据，生成 `result_analysis_report.md`（可读版）与 `result_analysis.json`（结构化版）。联合分析中的数值均来自对应数据文件（行数、均值、P 值、差异方向、分类指标、Top 通路等），图面要点只描述图的内容与检查重点；下结论前应回到本指南与原始数据交叉核对。

## 结果图坐标、颜色与图例图解

本节的坐标轴、颜色和图例说明直接对照本仓库绘图代码中实际使用的 `labs()`、`scale_*()`、`color`/`fill` 和 `palette` 定义，属于当前版本结果的官方图解。若后续代码修改轴标签、配色或阈值，应同步更新本节。单细胞图默认均以 `Seurat` 对象中的元数据作图，`UMAP_1`/`UMAP_2` 为该次运行 UMAP 降维后的两个坐标轴，本身无独立生物学单位。

### 质控与双细胞

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_01_qc_raw_violin.png` | 条件（`condition`） | `nFeature_RNA`、`nCount_RNA`、`percent.mt` 的数值 | 小提琴主体默认按条件分组，图例被移除 | 三个并列面板分别对应三个指标，用于查看过滤前分布 |
| `fig_01_qc_filtered_violin.png` | 条件（`condition`） | 过滤后 `nFeature_RNA`、`nCount_RNA`、`percent.mt` 的数值 | 同上 | 三个并列面板，用于查看过滤后分布 |
| `fig_48_qc_pvalue_comparison.png` | QC 指标（`nFeature_RNA`、`nCount_RNA`、`percent.mt`、`percent.ribo`） | `-log10(P value)` | 条形填充色表示不同的条件两两比较（`comparison`） | 灰色虚线为 `P = 0.05`（即 `-log10(0.05)`）；柱顶文字为原始 `P` 值；按 raw/filtered 两个面板分面 |
| `fig_02_doublet_scores.png` | 双细胞判定（`singlet`/`doublet`） | `scDblFinder score`（双细胞得分） | `singlet` 为浅蓝 `#4DBBD5`，`doublet` 为红 `#E64B35` | 小提琴叠加窄箱线图；得分越高越可能为双细胞 |

### 聚类、降维与注释

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_03_umap_clusters.png` | `UMAP_1` | `UMAP_2` | 每个聚类一种颜色，图例为“聚类号 - 主要细胞类型” | 图中直接标注主要细胞类型标签；用于检查聚类边界 |
| `fig_04_umap_condition.png` | `UMAP_1` | `UMAP_2` | 点颜色表示条件（`condition`） | 图例在右侧；用于检查条件分离与批次效应 |
| `fig_05_umap_annotation.png` | `UMAP_1` | `UMAP_2` | 左面板颜色表示 marker 注释，右面板表示发布注释 | 左右双面板并排，均带细胞类型图例 |
| `fig_06_dotplot_markers.png` | marker 基因 | 细胞类型（`celltype_annot`） | 点颜色表示平均表达量（颜色越深表达越高）；点大小表示表达该基因的细胞比例（`pct.exp`） | 为标准 Seurat `DotPlot` 双图例；用于验证注释特异性 |
| `fig_07_annotation_confusion_heatmap.png` | 发布细胞类型 | marker 注释细胞类型 | 热图默认蓝-白渐变，颜色越深计数越高 | 格内数字为细胞计数；用于查看注释对应关系 |
| `fig_14_pca.png` | `PC_1` | `PC_2` | 点颜色表示条件（`condition`） | 主成分散点，用于查看主要变异来源 |
| `fig_15_elbow.png` | 主成分序号（PC 1..N） | 主成分标准差 | 单条折线 | 拐点处建议保留的维度数 |
| `fig_16_featureplot_markers.png` | `UMAP_1` | `UMAP_2` | 点颜色表示 marker 基因表达量（浅灰到深红渐变，默认 `FeaturePlot` 配色） | 每个基因一个面板，最多展示 6 个 marker；输出条件：数据集中存在固定 marker 基因（与表达矩阵交集非空） |
| `fig_17_marker_violin.png` | 细胞类型（`celltype_annot`） | 表达量 | 小提琴按细胞类型分组 | 每个 marker 一个面板；输出条件：同 `fig_16`，需存在可用 marker 基因 |
| `fig_18_celltype_proportion.png` | 细胞类型 | 比例（0-100%，`position="fill"`） | 填充颜色表示条件（`condition`） | 柱高合计为 1，用于比较条件间细胞类型构成 |
| `fig_19_condition_proportion.png` | 条件（`coord_flip` 后显示为纵轴） | 比例（0-100%） | 填充颜色表示细胞类型 | 横向堆叠条形图，用于查看每个细胞类型内部的条件构成 |
| `fig_28_umap_sample.png` | `UMAP_1` | `UMAP_2` | 点颜色表示条件（`condition`），固定为红 `#E64B35` 与蓝 `#4DBBD5` | 标题为 “UMAP by sample”，用于评估样本混合和批次效应 |
| `fig_29_doublet_rate_sample.png` | 样本（短标签） | 双细胞率（百分比） | 条形填充颜色表示条件 | 柱顶标注百分比数字；用于识别异常样本 |
| `fig_30_sample_proportion.png` | 样本 | 比例（0-100%） | 填充颜色表示细胞类型 | 堆叠条形图，用于检查样本间异质性 |
| `fig_31_cluster_marker_heatmap.png` | 基因 | 聚类（`seurat_clusters`） | `viridis` 色阶（深紫到黄，值越大表达越高） | 每聚类 Top marker 表达块，用于辅助定义聚类标签 |
| `fig_32_cluster_marker_dotplot.png` | marker 基因 | 聚类（`seurat_clusters`） | 点颜色表示平均表达量；点大小表示表达比例 | 标准 Seurat `DotPlot` |
| `fig_37_singler_umap.png` | `UMAP_1` | `UMAP_2` | 点颜色表示 SingleR 参考注释 | 图例在右侧，用于查看参考注释与 UMAP 结构的一致性 |
| `fig_38_singler_confusion_heatmap.png` | SingleR 注释 | marker 注释 | 热图默认蓝-白渐变 | 格内数字为细胞计数；用于比较两种注释 |

### 差异表达

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_08_volcano.png` | `avg_log2FC`（平均 log2 差异倍数，标题中标注比较组，如 HCC vs iCCA） | `-log10(p_val_adj)`（校正 P 值的负对数） | 显著上调为红 `#E64B35`，显著下调为蓝 `#4DBBD5`，其余为灰 `grey75` | 灰色竖虚线为 `log2FC = ±0.25`，灰色横虚线为 `P = 0.05`；离群点旁标注基因名 |
| `fig_09_deg_heatmap.png` | 基因（Top DEG） | 样本/细胞（按条件分组） | `viridis` 色阶，深紫到黄表示表达量低到高 | 热图列按条件排列，用于查看差异基因整体表达模式 |
| `fig_09_deg_horizontal_violin.png` | 归一化表达量（`Normalized expression`） | 基因（按校正 P 值排序，标签含 `P = ...`） | 小提琴填充颜色表示条件（`condition`），颜色来自 `hcl.colors("Set 2")` | 基因按校正 P 值升序从上到下排列，用于对比最显著基因 |

### 富集分析

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_10_go_up.png`/`fig_11_go_down.png`/`fig_12_kegg_up.png`/`fig_13_kegg_down.png` | 按条目排序的通路（默认 `barplot`，最多 10 条） | 计数或 P 值指标（由 clusterProfiler 默认绘图函数决定） | 条形/点颜色默认表示 `p.adjust`（颜色越深越显著） | 实际绘图样式可用 `LIVER_FIGURE_STYLES` 切换为 dotplot/cnetplot；无显著条目时显示占位文字 |
| `fig_20_gsea_go.png`/`fig_21_gsea_kegg.png` | 按 log2FC 降序排列的基因排序位置 | 富集分数（GSEA 运行富集得分） | 默认 `ridgeplot` 时按通路着色；`gseaplot2` 时显示富集曲线与基因条带 | NES、p.adjust 等信息见图中标签；无显著条目时显示占位文字 |
| `fig_22_go_network.png`/`fig_23_kegg_network.png` | 无笛卡尔坐标，为网络布局 | 无 | cnetplot：基因节点颜色默认按 `log2FC` 渐变，通路节点按类别着色；emapplot：节点颜色表示 `p.adjust` | 圆形节点为通路，散点节点为基因；默认筛选 `p.adjust <= 0.05` 后保留 Top5 核心通路，并追加最佳 5 个延伸通路 |
| `fig_46_go_top5.png`/`fig_47_kegg_top5.png` | 富集条目（Top5） | 富集条目（按 `GeneRatio` 排序） | 气泡颜色表示 `p.adjust`（越显著越红/越深）；气泡大小表示 `Count`（富集基因数） | 默认 `dotplot`，可切换 barplot/cnetplot/emapplot |

### 机器学习（单细胞）

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_24_ml_feature_importance.png` | 特征重要性（数值） | 特征名 | 统一蓝色 `#4C72B0` | 横向条形图，仅展示重要性 Top15；特征为样本级 QC 均值与细胞类型比例 |
| `fig_25_ml_shap.png` | SHAP 值（对预测的影响） | 特征名 | 点颜色表示该样本/细胞该特征值的大小（红高蓝低，SHAP 默认配色） | SHAP 蜂群图，右侧为特征值颜色条 |
| `fig_43_ml_confusion_matrix.png` | 预测类别（`Predicted`） | 真实类别（`True`） | `Blues` 色阶，颜色越深计数越高；格内数字用白/黑文字对比显示 | 对角线为正确分类；数字为交叉验证预测计数 |
| `fig_44_ml_roc_pr.png` | 左：假阳性率（FPR）；右：召回率（Recall） | 左：真阳性率（TPR）；右：精确率（Precision） | 左/右各一条或多条折线，图例标注类别与 `AUC`/`AP`；灰色虚线为随机基线（ROC） | 二分类时图例显示 `AUC=...`/`AP=...`，多分类时每个类别一条线 |
| `fig_45_ml_cv_scores.png` | 单一类别 `CV accuracy` | 交叉验证准确率（Accuracy） | 箱线图主体默认样式，散点为蓝色 `#4C72B0` | 箱线叠加各折得分点；标题显示折数 `n=...` |
| `fig_45_ml_calibration_curve.png` | 平均预测概率（`Mean predicted probability`） | 观测频率（`Observed frequency`） | 蓝色校准折线与灰色完美校准虚线 | 越接近对角线表示概率校准越好 |

### 高级分析与细胞通讯

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_26_cellcycle_umap.png` | `UMAP_1` | `UMAP_2` | 点颜色表示细胞周期阶段（Seurat `Phase`：G1/S/G2M，默认离散配色） | 用于检查细胞周期是否干扰聚类 |
| `fig_27_cellcycle_proportion.png` | 条件 | 比例（0-100%） | 填充颜色表示细胞周期阶段 | 堆叠条形图，比较组间周期构成 |
| `fig_33_signature_scores_umap.png` | `UMAP_1` | `UMAP_2` | 点颜色表示功能签名得分，浅灰到深红 `grey90`→`#B31B1B` | 每个签名一个面板；签名为增殖、EMT、缺氧、免疫检查点等 |
| `fig_34_signature_scores_boxplot.png` | 条件 | 签名得分（`Signature score`） | 小提琴/箱线按条件填充 | 每个签名一个分面面板 |
| `fig_35_celltype_abundance_effect.png` | `log2 odds ratio`（细胞类型丰度 log2 OR） | `-log10 adjusted p value` | 点颜色表示细胞类型 | 灰色虚线为 `x = 0`；用于查看细胞类型在组间显著增多/减少 |
| `fig_36_cnv_heatmap.png` | 染色体滑动窗口（按染色体顺序） | 细胞（聚类后按相似性排列） | 蓝 `#3B4CC0`-白-红 `#B40426`：蓝为相对低表达（推断拷贝缺失），红为相对高表达（推断拷贝增益） | 顶部注释条：`Condition` 使用红/蓝/绿，`CellType` 使用彩虹色；仅展示每条件最多 750 个细胞 |
| `fig_39_trajectory_umap.png` | `UMAP_1` | `UMAP_2` | 左面板为伪时序（`pseudotime`），颜色使用 `viridis`；右面板点颜色表示聚类，轨迹线为黑色 | 左右双面板：左为拟时序 FeaturePlot，右为聚类 UMAP 叠加 slingshot 谱系曲线 |
| `fig_40_cellchat_network.png` | 无笛卡尔坐标，为环形网络布局 | 无 | 节点颜色/大小表示细胞类型及细胞数；边粗细表示通讯数量 | 弧上标签为通讯数量；用于查看细胞间通讯强度 |
| `fig_41_cellchat_heatmap.png` | 目标细胞类型 | 来源细胞类型 | 红蓝渐变，颜色越红通讯数量/强度越高 | CellChat 默认热图，用于比较细胞类型对之间通讯 |
| `fig_42_cellchat_bubble.png` | 配体受体/信号通路 | 细胞类型对（来源→目标） | 点颜色表示通讯概率（`p` 值相关颜色，默认红蓝渐变）；点大小表示贡献/显著性 | 用于查看关键配体受体对 |

### 虚拟筛选与敲除

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_46_affinity_distribution.png` | 亲和力（`Affinity (kcal/mol)`，越负结合越强） | 配体数量（`Ligand count`） | 直方图统一蓝色 `#4c7bb8`；红色虚线为命中阈值 `cutoff`（默认 `-7.0 kcal/mol`） | 阈值线左侧为命中配体 |
| `fig_47_top_hits.png` | 亲和力（`Affinity (kcal/mol)`） | 配体 ID | 条形统一绿色 `#2e7d32` | 仅展示 Top20；若 `hits=0` 则为按亲和力排序的候选而非“命中” |
| `fig_48_diverse_hits.png` | 亲和力（`Affinity (kcal/mol)`） | 配体 ID | 条形统一紫色 `#8e44ad` | Tanimoto 多样性选择结果，用于减少同一骨架重复 |
| `fig_49_redock_comparison.png` | 左：初始亲和力；右：亲和力变化量 | 左：重对接亲和力；右：配体数量 | 左散点统一绿色 `#2e7d32`，黑色虚线为 `y=x` 对角线；右直方图蓝色 `#4c7bb8`，红色虚线为 0 | 右图围绕 0 分布表示重对接稳定 |
| `fig_50_ml_feature_importance.png` | 特征重要性 | 特征名 | 条形统一蓝色 `#1665c0` | 模型提供 `feature_importances_`/`coef_` 时生成；MLP/torch 无原生重要性时不生成 |
| `fig_51_ml_roc.png` | 假阳性率（`False positive rate`） | 真阳性率（`True positive rate`） | 蓝色 ROC 曲线，灰色虚线为随机基线 | 图例显示 `AUC=...`；用于评估分类模型排序能力 |
| `fig_52_knockout_top_candidates.png` | 敲除优先级得分（`Knockout priority score`，0-1） | 基因名 | 条形统一绿色 `#1d6f42` | 横轴固定为 0-1；分数是统计优先级，不等同于真实敲除表型 |
| `fig_53_knockout_score_distribution.png` | 敲除优先级得分（0-1） | 基因数量（`Gene count`） | 直方图统一蓝色 `#4c7bb8` | 用于判断评分是否有区分度 |
| `compound_disease_venn.png` | 无笛卡尔坐标，为集合圆形布局 | 无 | 各集合圆轮廓分别为蓝 `#4C72B0`、橙 `#DD8452`、绿 `#55A868`，区域数字为基因计数 | 数字表示各交集/独有区域基因数；用于查看化合物靶点与疾病基因交集 |

### 细胞反馈

| 文件 | 横轴 | 纵轴 | 颜色/图例 | 关键符号 |
| --- | --- | --- | --- | --- |
| `fig_54_feedback_module_umap.png` | `UMAP_1` | `UMAP_2` | 点颜色表示筛选靶点模块得分，浅灰 `grey90` 到深红 `#B31B1B` | 用于查看模块是否集中在特定细胞群 |
| `fig_55_feedback_target_expression_umap.png` | `UMAP_1` | `UMAP_2` | 点颜色表示候选靶基因表达量，浅灰 `grey90` 到蓝 `#2E86AB` | 每个基因一个面板 |
| `fig_56_feedback_celltype_dotplot.png` | 候选靶基因 | 细胞类型 | 点颜色表示平均表达量；点大小表示表达比例 | 标准 Seurat `DotPlot` |
| `fig_57_feedback_celltype_boxplot.png` | 细胞类型 | 反馈模块得分（`Feedback module score`） | 小提琴/箱线按条件（`condition`）填充 | 用于比较细胞类型与条件间的模块活性 |
| `fig_58_feedback_celltype_heatmap.png` | 细胞类型 | 候选靶基因 | 白-浅蓝-蓝-深红渐变（`#f7fbff`→`#d6e7ff`→`#2E86AB`→`#7B241C`），颜色越深平均表达越高 | 行和列均聚类，基因按列和均值排序 |
| `fig_59_feedback_targets_volcano.png` | `avg_log2FC` | `-log10(p_val_adj)` | 显著上调为红 `#C0392B`，显著下调为蓝 `#2E86C1`，其余为灰 `grey65` | 灰色虚线为 `P = 0.05` 和 `log2FC = ±0.25`；显著基因标注基因名 |
| `fig_60_feedback_condition_violin.png` | 条件（`condition`） | 归一化表达量 | 小提琴按条件分组（图例被移除） | 每个候选靶基因一个面板 |
| `fig_61_feedback_go_network.png`/`fig_62_feedback_kegg_network.png` | 无笛卡尔坐标，为通路-基因网络布局 | 无 | 基因节点按表达/差异方向着色，通路节点按类别着色（cnetplot 默认配色） | 默认筛选 `p.adjust <= 0.05` 后取 Top5 通路 |

### 逐图详细图解

本节对结果清单中的每一张图给出完整注解：坐标轴定义与单位、颜色/大小/线型编码、阈值与参考线、面板布局、数据来源、读取方法，以及可能出现的回退或占位情况。数据来源列中的“状态文件”指 `results/.stages`、`summary.json`、`pipeline_complete.json`、`ml_model_summary.json`、`cellchat_status.txt` 等过程文件；只有状态、数据、参数和视觉检查都通过，才可作为结论材料。

#### `fig_01_qc_raw_violin.png`

- 坐标轴：横轴为样本条件（Seurat 元数据 `condition`）；纵轴为 QC 指标数值，三个并列面板依次为 `nFeature_RNA`（检测到的基因数）、`nCount_RNA`（总 UMI/总分子数）、`percent.mt`（线粒体基因百分比）。
- 颜色/图例：小提琴按条件分组填充，代码使用 `NoLegend()` 移除图例，因此图中不显示图例；分组由横轴标签区分。
- 读取方法：比较每个条件下小提琴主体位置、宽度和长尾。`nFeature_RNA`/`nCount_RNA` 过低代表低质量或空液滴，`percent.mt` 过高代表细胞状态差或细胞破裂；单条件时只能描述质量分布，不能用于“条件间差异”结论。
- 数据来源：`data/01_qc/fig_01_qc_metrics.csv`（含 raw 前原始指标）与 `results/checkpoints/seurat_raw.rds`。
- 回退/占位：无。若图像空白或全部为 0，应检查输入矩阵与 QC 计算是否成功。

#### `fig_01_qc_filtered_violin.png`

- 坐标轴：横轴为条件；纵轴为过滤后 `nFeature_RNA`、`nCount_RNA`、`percent.mt` 的数值；三个面板同上。
- 颜色/图例：同 `fig_01_qc_raw_violin.png`，图例被移除。
- 读取方法：确认低质量尾部被去除、主要细胞群保留，并检查两个条件是否仍可比较。默认过滤阈值见环境变量 `LIVER_QC_MIN_FEATURES`、`LIVER_QC_MAX_FEATURES`、`LIVER_QC_MIN_COUNTS`、`LIVER_QC_MAX_COUNTS`、`LIVER_QC_MAX_MT`；未设置时由数据分位数自动计算。
- 数据来源：`fig_01_qc_metrics.csv` 的过滤后子集。
- 回退/占位：无。过滤后细胞数过少或分布被过度压缩时，该图不能作为质量合格证据。

#### `fig_48_qc_pvalue_comparison.png`

- 坐标轴：横轴为 QC 指标（`nFeature_RNA`、`nCount_RNA`、`percent.mt`、`percent.ribo`）；纵轴为 `-log10(P value)`，数值越大表示条件间差异越显著。
- 颜色/图例：条形填充颜色表示条件两两比较（`comparison`，格式为“组1 vs 组2”）；图例在底部。
- 阈值与参考线：灰色虚线为 `P = 0.05`（对应 `-log10(0.05) ≈ 1.301`）；柱顶文字直接标注原始 `P` 值。
- 面板布局：按 `raw`/`filtered` 两个阶段分面，纵轴可各自缩放（`scales="free_y"`）。
- 统计方法：Wilcoxon 秩和检验（`wilcox.test`），`padj` 使用 BH FDR 校正，图中显示原始 P 值。
- 读取方法：比较同一指标在过滤前后的 P 值变化，判断 QC 是否改变条件间差异；过滤后 P 值显著变化不代表过滤有误，需结合分布图解释。
- 数据来源：`data/01_qc/fig_48_qc_pvalue_comparison.csv`。
- 回退/占位：条件少于 2 组时显示 “At least two conditions are required”，此时该图不能用于条件差异判断。

#### `fig_02_doublet_scores.png`

- 坐标轴：横轴为双细胞判定（`doublet_call`：`singlet`/`doublet`）；纵轴为 `scDblFinder score`（scDblFinder 双细胞得分）。
- 颜色/图例：`singlet` 为浅蓝 `#4DBBD5`，`doublet` 为红 `#E64B35`；图例为类别本身。
- 图形编码：小提琴（`geom_violin`）叠加窄箱线图（`geom_boxplot`，`outlier.shape=NA`），箱线不显示离群点。
- 读取方法：正常结果中 doublet 得分应整体高于 singlet；若两个类别得分完全重叠或全部为 singlet，需检查 `scDblFinder` 是否成功运行。
- 数据来源：`data/02_doublets/fig_02_doublet_results.csv`（每细胞得分与判定）。
- 回退/占位：`scDblFinder` 失败时代码会把所有细胞标记为 `singlet`、得分设为 0，此时该图不能作为真实双细胞检测结果。

#### `fig_14_pca.png`

- 坐标轴：横轴 `PC_1`，纵轴 `PC_2`，为主成分空间前两个主成分。
- 颜色/图例：点颜色表示条件（`condition`），图例在右侧。
- 读取方法：查看条件/样本在主成分空间中的分离程度、批次混杂和离群点；主成分本身无独立生物学单位。
- 数据来源：Seurat 对象的 `pca` 降维结果；坐标与 `fig_03_04_05_umap_coordinates.csv` 不同，PCA 坐标未单独导出。
- 回退/占位：样本级模式下使用随机占位嵌入（`sample-level`），只能作描述性展示。

#### `fig_15_elbow.png`

- 坐标轴：横轴为主成分序号（PC 1..N，`ndims` 为 `min(30, 细胞数-1)`）；纵轴为主成分标准差。
- 颜色/图例：单条折线，无分类图例。
- 读取方法：找到曲线由陡变缓的“拐点”，结合下游分析所需的可解释方差选择 PCA 维度数；无清晰拐点不代表图片损坏，但应说明维度选择依据。
- 数据来源：Seurat `ElbowPlot` 的标准差向量；样本级模式下为恒定值 1。
- 回退/占位：样本级模式生成水平线，不能用于真实 PCA 维度选择。

#### `fig_03_umap_clusters.png`

- 坐标轴：横轴 `UMAP_1`，纵轴 `UMAP_2`，为该次运行 UMAP 降维后的坐标。
- 颜色/图例：每个 Seurat 聚类一种颜色，图例文本为“聚类号 - 该聚类主要细胞类型”（如 `0 - Hepatocyte`）；图中还直接标注主要细胞类型标签。
- 读取方法：检查聚类边界是否清晰、是否过度聚类（大量碎片）或欠聚类（全部混成一片），并核对图例聚类数与 `summary.json` 的 `n_clusters` 一致。
- 数据来源：`data/03_cluster/fig_03_04_05_umap_coordinates.csv`（UMAP 坐标、聚类、条件、样本）。
- 回退/占位：UMAP 计算失败时使用随机占位嵌入，此时聚类结构无意义。

#### `fig_04_umap_condition.png`

- 坐标轴：`UMAP_1`/`UMAP_2`。
- 颜色/图例：点颜色表示条件，图例在右侧；标题为条件比较（如 `HCC vs iCCA`）。
- 读取方法：评估条件间分布偏移、批次效应和分组混杂。单条件图不能支持条件差异；条件完全分离时需确认是否来自批次而不是真实生物学差异。
- 数据来源：同一 UMAP 坐标表。
- 回退/占位：UMAP 失败时同上。

#### `fig_05_umap_annotation.png`

- 坐标轴：`UMAP_1`/`UMAP_2`；左右双面板并排。
- 颜色/图例：左面板为 marker 注释（`celltype_annot`），右面板为发布注释（`published_type`），各带细胞类型图例。
- 读取方法：比较两种注释在 UMAP 上是否形成连贯且一致的区域；差异大的细胞类型应回到 marker 图与混淆矩阵复核。
- 数据来源：`data/04_annotation/fig_05_16_17_cell_annotations.csv`。
- 回退/占位：无发布注释时右面板可能为空或缺失；占位说明见 `fig_07`。

#### `fig_06_dotplot_markers.png`

- 坐标轴：横轴为 marker 基因，纵轴为细胞类型（`celltype_annot`）。
- 颜色/图例：点颜色表示平均表达量（颜色越深表达越高）；点大小表示表达该基因的细胞比例（`pct.exp`），为 Seurat `DotPlot` 双图例。
- 读取方法：预期 marker 在其对应细胞类型中高表达且高比例，其他类型低表达；若所有点大小/颜色无差异或图像空白，不能验证注释。
- 数据来源：`data/04_annotation/fig_05_16_17_cell_annotations.csv` 与 Seurat 表达矩阵。
- 回退/占位：marker 基因全部未命中时跳过生成。

#### `fig_07_annotation_confusion_heatmap.png`

- 坐标轴：横轴为发布细胞类型，纵轴为 marker 注释细胞类型。
- 颜色/图例：`pheatmap` 默认蓝-白渐变，颜色越深计数越高；格内数字为细胞计数。
- 读取方法：查看对角线是否为高计数主对应关系；出现大片非对角线高计数时说明两种注释存在系统差异。
- 数据来源：`data/04_annotation/fig_07_annotation_confusion.csv`。
- 回退/占位：无发布注释或矩阵过小时显示 “No published annotations”，只能说明缺少发布注释，不能用于一致性判断。

#### `fig_08_volcano.png`

- 坐标轴：横轴 `avg_log2FC`（平均 log2 差异倍数，标题中标注比较组，如 `HCC vs iCCA`）；纵轴 `-log10(p_val_adj)`（BH 校正 P 值的负对数）。
- 颜色/图例：显著上调为红 `#E64B35`，显著下调为蓝 `#4DBBD5`，其余为灰 `grey75`；图例显示 Up/Down/NS。
- 阈值与参考线：灰色竖虚线为 `log2FC = ±0.25`，灰色横虚线为 `P = 0.05`；实际显著性判定使用 `p_val_adj < de_padj` 且 `abs(avg_log2FC) > de_logfc`，默认 `0.05` 与 `0.25`，可用环境变量修改。
- 读取方法：查看显著上/下调基因数量、分布对称性和离群点；Top15 基因（按 P 值排序）标注基因名。
- 回退/占位：可通过 `LIVER_FIGURE_STYLES=fig_08_volcano.png=maplot` 切换为 MA 图（横轴为平均归一化表达，纵轴为 `avg_log2FC`）；无显著 DEG 时仍会生成图，只能说明未检出。

#### `fig_09_deg_heatmap.png`

- 坐标轴：横轴为 Top DEG 基因，纵轴为样本/细胞，按条件分组；列标签旋转 45 度。
- 颜色/图例：`viridis` 色阶（深紫到黄），表示 `ScaleData` 后的表达水平，值越大颜色越黄。
- 读取方法：查看显著基因在条件间的表达块是否清晰；Top30 来自 `deg$gene[seq_len(min(30, nrow(deg)))]`。
- 回退/占位：无显著 DEG 时 Top30 可能来自全部基因排序，需结合 `fig_09_deg_significant.csv` 判断；无基因可绘图时跳过。

#### `fig_09_deg_horizontal_violin.png`

- 坐标轴：横轴为归一化表达量；纵轴为基因，按校正 P 值升序排列，标签为“基因名 + `P = ...`”（P < 0.001 用科学计数法）。
- 颜色/图例：小提琴填充颜色表示条件，颜色来自 `hcl.colors(length(conds), palette="Set 2")`；图例在顶部。
- 图形编码：小提琴叠加窄箱线图和抖动散点（每条件最多 `deg_violin_max_cells` 个细胞，默认 1000，种子固定）。
- 读取方法：对比两组表达分布与 P 值标签，优先看差异最显著的基因；基因数由 `LIVER_DE_VIOLIN_TOP_N` 控制（默认 12）。
- 回退/占位：无显著 DEG 时可能回退到“有限 P 值的全部基因”，不能把这些基因称为显著差异基因。

#### `fig_10_go_up.png`/`fig_11_go_down.png`/`fig_12_kegg_up.png`/`fig_13_kegg_down.png`

- 坐标轴：默认 `barplot` 时横轴为富集基因计数（Count），纵轴为富集条目（Description/ID）；若用 `LIVER_FIGURE_STYLES` 切换为 `dotplot`，横轴为 `GeneRatio`，纵轴为条目。
- 颜色/图例：条形/点颜色默认表示 `p.adjust`（颜色越深越显著），图例为 P 值色条。
- 输入：上调/下调显著 DEG（`p_val_adj < 0.05` 且 `|log2FC| > 0.25`，少于 10 个时放宽为 `p_val_adj < 0.1` 且 `|log2FC| > 0.1`）；GO 使用 `enrichGO(ont="BP", BH, pvalueCutoff=0.1, qvalueCutoff=0.2)`，KEGG 使用 `enrichKEGG`。
- 读取方法：查看显著通路的富集方向、基因数与 P 值；最多展示 10 条。
- 数据来源：`fig_10_enrichment_up_go.csv`、`fig_11_enrichment_down_go.csv`、`fig_12_enrichment_up_kegg.csv`、`fig_13_enrichment_down_kegg.csv`。
- 回退/占位：无显著条目时显示 “No significant enrichment terms”，不能作为阳性富集证据。

#### `fig_20_gsea_go.png`/`fig_21_gsea_kegg.png`

- 坐标轴：默认 `ridgeplot` 时横轴为按 `avg_log2FC` 降序排列的基因排序位置，纵轴为富集分数分布；`gseaplot2` 时横轴为排名位置，纵轴为富集分数（Running Enrichment Score）。
- 颜色/图例：`ridgeplot` 按通路着色；`gseaplot2` 显示富集曲线、基因条带和 NES/p.adjust 标签。
- 输入：全部 DEG 的 `avg_log2FC` 排序向量，经 `bitr` 映射 Entrez ID；`LIVER_SKIP_GSEA=yes` 或映射失败时跳过。
- 读取方法：查看通路整体上调/下调方向（NES 正/负）与显著性；无显著 GSEA 条目时显示 “No significant GSEA terms”。

#### `fig_22_go_network.png`/`fig_23_kegg_network.png`

- 坐标轴：无笛卡尔坐标，为网络布局（cnetplot 默认）。
- 颜色/图例：cnetplot 中基因节点颜色默认按 `avg_log2FC` 渐变，通路节点按类别着色；emapplot 中节点颜色表示 `p.adjust`。
- 输入：上调基因富集结果经 `p.adjust <= 0.05` 筛选后保留 Top5 核心通路，并追加按 p.adjust 排序最佳的 5 个延伸通路；`showCategory=10`。
- 读取方法：圆形节点为通路，散点节点为基因，查看通路-基因关联结构；无显著通路时显示 “No significant pathway network”。

#### `fig_46_go_top5.png`/`fig_47_kegg_top5.png`

- 坐标轴：纵轴为 Top5 富集条目，横轴为 `GeneRatio`（富集基因占输入基因的比例）。
- 颜色/图例：气泡颜色表示 `p.adjust`（越显著越红/越深），气泡大小表示 `Count`（富集基因数）。
- 输入：同 `fig_22/23` 的筛选条件。
- 读取方法：快速查看最显著的上调 GO/KEGG 通路；无筛选后通路时显示 “No significant pathway after filtering”。

#### `fig_24_ml_feature_importance.png`

- 坐标轴：横轴为特征重要性数值，纵轴为特征名。
- 颜色/图例：条形统一蓝色 `#4C72B0`，无分类图例。
- 输入：单细胞 ML 分类模型的 `feature_importances_`/`coef_` 绝对值均值；特征为样本级 QC 均值（`nFeature_RNA`、`nCount_RNA`、`percent.mt`、`percent.ribo`）与细胞类型比例。
- 读取方法：查看影响样本分类的关键特征；仅展示 Top15。
- 数据来源：`data/07_ml/fig_24_ml_feature_importance.csv`、`ml_model_summary.json`（状态为 completed 才可用）。

#### `fig_25_ml_shap.png`

- 坐标轴：横轴为 SHAP 值（特征对预测的贡献，正/负表示方向）；纵轴为特征名。
- 颜色/图例：点颜色表示该样本/细胞该特征值的大小（SHAP 默认红高蓝低），右侧为特征值色条。
- 输入：XGBoost/RF/GBM 使用 `TreeExplainer`，其余使用 `KernelExplainer`；最多展示 15 个特征。
- 回退/占位：SHAP 计算失败时写 `ml_shap_status.txt`（内容为 “SHAP plot skipped”），此时该图不可用。

#### `fig_43_ml_confusion_matrix.png`

- 坐标轴：横轴为预测类别，纵轴为真实类别。
- 颜色/图例：`Blues` 色阶，颜色越深计数越高；格内数字颜色按计数是否超过最大值一半自动选择白/黑。
- 读取方法：查看对角线正确分类与非对角线错分模式；样本量很小或模型状态失败时只能谨慎解释。
- 数据来源：`fig_43_44_45_ml_classification_results.csv` 与 `ml_model_summary.json`。

#### `fig_44_ml_roc_pr.png`

- 坐标轴：左图横轴为假阳性率（FPR），纵轴为真阳性率（TPR）；右图横轴为召回率（Recall），纵轴为精确率（Precision）。
- 颜色/图例：二分类一条曲线，图例标注 `AUC=...`/`AP=...`；多分类每个类别一条曲线，图例标注类别与对应值；ROC 灰色虚线为随机基线。
- 读取方法：AUC/AP 越接近 1 越好；AUC 接近 0.5、样本极少或类别不平衡时不能用于强结论。

#### `fig_45_ml_cv_scores.png`

- 坐标轴：横轴为单一类别 `CV accuracy`，纵轴为交叉验证准确率。
- 颜色/图例：箱线图主体默认样式，散点为蓝色 `#4C72B0`，每个点代表一折。
- 读取方法：查看均值与离散程度；折数 `n=...` 由类别最小样本数决定（最多 5 折）。
- 回退/占位：状态为 `skipped`/`failed` 时不生成或不可用。

#### `fig_45_ml_calibration_curve.png`

- 坐标轴：横轴为平均预测概率（`Mean predicted probability`），纵轴为观测频率（`Observed frequency`）。
- 颜色/图例：蓝色校准折线，灰色虚线为完美校准对角线（标签 Perfect）。
- 读取方法：越接近对角线表示预测概率越可靠；仅二分类且概率可用时生成，样本极少时需谨慎。
- 输出条件：单细胞 ML 阶段中仅当类别数为 2（二分类）且 `calibration_curve` 计算成功时生成；模型状态为 `skipped`/`failed` 时不生成。

#### `fig_26_cellcycle_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`。
- 颜色/图例：点颜色表示细胞周期阶段（Seurat `Phase`：G1/S/G2M，默认离散配色）。
- 读取方法：若 phase 形成强分群，说明需要评估细胞周期回归；图本身不一定损坏。
- 开关：`LIVER_RUN_CELLCYCLE`（默认 `yes`），需至少 5 个 S 期与 5 个 G2M 期基因。
- 输出条件：`LIVER_RUN_CELLCYCLE=yes`（默认）且 `CellCycleScoring` 成功写入 `Phase`（至少 5 个 S 期、5 个 G2M 期基因，或元数据已含 `Phase`）；否则不生成。

#### `fig_27_cellcycle_proportion.png`

- 坐标轴：横轴为条件，纵轴为比例（0-100%，`position="fill"`）。
- 颜色/图例：填充颜色表示细胞周期阶段。
- 读取方法：比较组间细胞周期构成；比例信息过少时需谨慎。
- 输出条件：同 `fig_26`，需 `LIVER_RUN_CELLCYCLE=yes` 且元数据含 `Phase`。

#### `fig_28_umap_sample.png`

- 坐标轴：`UMAP_1`/`UMAP_2`。
- 颜色/图例：点颜色表示条件，固定为红 `#E64B35` 与蓝 `#4DBBD5`，图例标题为 `Condition`。
- 读取方法：评估样本混合和批次效应；仅当样本数 > 1 时生成。
- 回退/占位：单样本时不生成，缺失正常。
- 输出条件：样本数 `length(unique(seurat$sample)) > 1` 时生成，无环境变量开关。

#### `fig_29_doublet_rate_sample.png`

- 坐标轴：横轴为样本短标签（去掉冒号后缀），纵轴为双细胞率（百分比）。
- 颜色/图例：条形填充颜色表示条件；柱顶标注 `%.1f%%` 的百分比。
- 读取方法：识别双细胞率异常偏高的样本；需与 `fig_29_doublet_rate_by_sample.csv` 一致。

#### `fig_30_sample_proportion.png`

- 坐标轴：横轴为样本短标签，纵轴为比例（0-100%）。
- 颜色/图例：填充颜色表示细胞类型。
- 读取方法：检查样本间细胞类型构成与异质性；样本过少或比例全为 0 时不可用。

#### `fig_31_cluster_marker_heatmap.png`

- 坐标轴：横轴为聚类 marker 基因，纵轴为聚类（`seurat_clusters`）。
- 颜色/图例：`viridis` 色阶，深紫到黄表示表达量低到高。
- 输入：`FindAllMarkers`（only.pos=TRUE，`min.pct=0.25`，`logfc.threshold=0.5`），每聚类 Top3 marker。
- 读取方法：查看每聚类 marker 表达块，辅助定义聚类标签；无 marker 时不生成。
- 输出条件：`LIVER_RUN_CLUSTER_MARKERS=yes`（默认）且 `FindAllMarkers` 返回非空、Top 基因存在于表达矩阵；否则不生成。

#### `fig_32_cluster_marker_dotplot.png`

- 坐标轴：横轴为 marker 基因，纵轴为聚类。
- 颜色/图例：点颜色表示平均表达量，点大小表示表达比例。
- 读取方法：marker 应在对应聚类中高表达且高比例；无 marker 或全 0 时不可用。
- 输出条件：同 `fig_31`，需 `LIVER_RUN_CLUSTER_MARKERS=yes` 且存在可绘制的 Top marker 基因。

#### `fig_33_signature_scores_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`，每个签名一个面板。
- 颜色/图例：点颜色表示功能签名得分，浅灰 `grey90` 到深红 `#B31B1B`。
- 输入：`AddModuleScore` 计算的签名（默认展示增殖、EMT、缺氧、免疫检查点；可用 `LIVER_RUN_SIGNATURES=yes` 开启），签名基因不足 3 个时跳过。
- 读取方法：查看签名是否集中在特定细胞群；得分全图均匀时不能解释。
- 输出条件：`LIVER_RUN_SIGNATURES=yes`（默认）且签名基因与表达矩阵交集非空（每个签名至少 3 个基因）。

#### `fig_34_signature_scores_boxplot.png`

- 坐标轴：横轴为条件，纵轴为签名得分。
- 颜色/图例：小提琴/箱线按条件填充（图例隐藏），每个签名一个分面面板，纵轴自由缩放。
- 读取方法：比较条件间签名活性；无差异不代表损坏，但不能当作“有差异”结论。
- 输出条件：同 `fig_33`，需 `LIVER_RUN_SIGNATURES=yes` 且签名基因充足。

#### `fig_35_celltype_abundance_effect.png`

- 坐标轴：横轴为 `log2 odds ratio`（细胞类型丰度 log2 OR），纵轴为 `-log10 adjusted p value`。
- 颜色/图例：点颜色表示细胞类型，点旁标注细胞类型名；灰色虚线为 `x = 0`。
- 输入：`fig_18_19_celltype_proportion_stats.csv` 中 Fisher 精确检验的 `OddsRatio` 与 BH 校正 `Padj`。
- 读取方法：查看哪些细胞类型在条件间显著增多/减少；细胞数过少、OR 不有限或 P 值缺失时需谨慎。

#### `fig_36_cnv_heatmap.png`

- 坐标轴：横轴为染色体滑动窗口（按染色体顺序，标签如 `chr1_w1`），纵轴为细胞（聚类后按相似性排列）。
- 颜色/图例：蓝 `#3B4CC0`-白-红 `#B40426`，蓝为相对低表达（推断拷贝缺失），红为相对高表达（推断拷贝增益），白色为接近中位。
- 注释：顶部注释条 `Condition` 使用红/蓝/绿，`CellType` 使用彩虹色。
- 输入：每条件最多 750 个细胞，基于染色体滑动窗口均值推断；仅展示 `CHRLOC` 注释充足的基因（>500）。
- 读取方法：查看细胞亚群的染色体拷贝数变化和肿瘤异质性；全图白色或窗口过少时不可用。
- 输出条件：`LIVER_RUN_CNV=yes`（默认）且物种注释库（`org.Hs.eg.db`/`org.Mm.eg.db`）可用、`CHRLOC` 非缺失基因 > 500、可映射基因 >= 200；否则不生成。

#### `fig_37_singler_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`。
- 颜色/图例：点颜色表示 SingleR 参考注释（`singleR_label`），图例在右侧。
- 输入：`celldex` 参考（人类 `HumanPrimaryCellAtlasData`、小鼠 `MouseRNAseqData`），最多 20000 个细胞，按聚类填充标签。
- 回退/占位：参考数据不可用或预测失败时不生成，缺失正常。
- 输出条件：`LIVER_RUN_SINGLER=yes`（默认）、参考数据可用且 SingleR 预测成功。

#### `fig_38_singler_confusion_heatmap.png`

- 坐标轴：横轴为 SingleR 注释，纵轴为 marker 注释。
- 颜色/图例：`pheatmap` 默认蓝-白渐变，格内数字为细胞计数。
- 读取方法：查看对角线一致性；无 SingleR 预测时不生成。
- 输出条件：同 `fig_37`，需 `LIVER_RUN_SINGLER=yes` 且 SingleR 预测成功、混淆矩阵行列数均 > 0。

#### `fig_39_trajectory_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`；左右双面板。
- 颜色/图例：左面板点颜色为伪时序（`pseudotime`，`viridis`），右面板点颜色为聚类，黑色线为 slingshot 谱系曲线。
- 输入：`slingshot` 以聚类为起点推断谱系，取第一条谱系伪时序写入 `fig_39_trajectory_pseudotime.csv`。
- 开关：`LIVER_RUN_TRAJECTORY`（默认 `no`），未开启时不生成，缺失正常；轨迹无结构时不能作分化结论。
- 输出条件：`LIVER_RUN_TRAJECTORY=yes`、R 环境安装 `slingshot` 且至少推断出 1 条谱系。

#### `fig_40_cellchat_network.png`

- 坐标轴：无笛卡尔坐标，为环形网络布局。
- 颜色/图例：节点颜色/大小表示细胞类型及该类型细胞数；边粗细表示通讯数量（`cellchat@net$count`，`weight.scale=TRUE`），弧上标签为通讯数。
- 输入：CellChat 标准流程（`computeCommunProb` → `filterCommunication` → `aggregateNet`），`LIVER_RUN_CELLCHAT=yes` 且 R 环境可用时生成。
- 回退/占位：未运行、未安装或无通讯时不生成，缺失正常。
- 输出条件：`LIVER_RUN_CELLCHAT=yes`、R 环境可用、CellChat 已安装且 `filterCommunication` 后存在通讯；否则不生成。

#### `fig_41_cellchat_heatmap.png`

- 坐标轴：横轴为目标细胞类型，纵轴为来源细胞类型。
- 颜色/图例：CellChat 默认热图（红蓝渐变），颜色越红通讯数量/强度越高。
- 读取方法：比较细胞类型对之间的通讯强度；全部为零或空白时不能作为通讯证据。
- 输出条件：同 `fig_40`，需 `LIVER_RUN_CELLCHAT=yes` 且通讯矩阵非空。

#### `fig_42_cellchat_bubble.png`

- 坐标轴：横轴为配体受体/信号通路，纵轴为细胞类型对（来源→目标）。
- 颜色/图例：点颜色表示通讯概率（p 值相关颜色，默认红蓝渐变），点大小表示贡献/显著性；`remove.isolate=TRUE` 移除孤立对。
- 读取方法：筛选候选配体受体互作；无显著对时只能说明未检出。
- 输出条件：同 `fig_40`，需 `LIVER_RUN_CELLCHAT=yes` 且存在通路/配体受体通讯结果。

#### `fig_46_affinity_distribution.png`

- 坐标轴：横轴为亲和力（`Affinity (kcal/mol)`，越负结合越强），纵轴为配体数量（`Ligand count`）。
- 颜色/图例：直方图统一蓝色 `#4c7bb8`，白色描边；红色虚线为命中阈值（`cutoff`，默认 `-7.0 kcal/mol`），阈值线左侧为命中。
- 数据来源：`data/01_analysis/fig_46_47_ranked_results.csv`（仅状态为 ok 且有亲和力的配体）。
- 读取方法：查看整体打分分布与命中数量；`summary.json` 中 `total_docked`、`hits`、`best_affinity` 应与图一致。

#### `fig_47_top_hits.png`

- 坐标轴：横轴为亲和力（`Affinity (kcal/mol)`），纵轴为配体 ID。
- 颜色/图例：条形统一绿色 `#2e7d32`，无分类图例。
- 输入：Top20（`hits` 非空时取命中前 20，否则取排序前 20）。
- 读取方法：查看排名靠前配体；若 `hits=0`，图中只是“Top 排序”而不是“Top 命中”，不能写成命中结果。

#### `fig_48_diverse_hits.png`

- 坐标轴：横轴为亲和力，纵轴为配体 ID。
- 颜色/图例：条形统一紫色 `#8e44ad`。
- 输入：Tanimoto 多样性选择（Morgan 指纹半径 2、2048 位，阈值 0.7 用于减少同一骨架重复）。
- 回退/占位：RDKit 不可用或有效分子少于 2 个时可能不生成；只有一个骨架时不能称为多样。

#### `fig_49_redock_comparison.png`

- 坐标轴：左图横轴为初始亲和力（`Initial affinity`），纵轴为重对接亲和力（`Redock affinity`），单位 kcal/mol；右图横轴为亲和力变化量（`Affinity change`），纵轴为配体数量。
- 颜色/图例：左散点统一绿色 `#2e7d32`，黑色虚线为 `y=x` 对角线；右直方图蓝色 `#4c7bb8`，红色虚线为 0。
- 读取方法：点接近对角线表示排序稳定；右图围绕 0 分布表示无系统性漂移；只有 1 个可比较配体时可能不生成。
- 输出条件：重对接阶段启用（`redock.enabled`，默认 true）、存在可重对接的 Top 命中且至少 1 个配体的初始与重对接亲和力可合并；合格判据建议至少 2 个可比较配体。

#### `fig_50_ml_feature_importance.png`

- 坐标轴：横轴为特征重要性，纵轴为特征名。
- 颜色/图例：条形统一蓝色 `#1665c0`。
- 输入：对接 ML 重打分模型的特征重要性（RF/GBM 的 `feature_importances_`）或系数绝对值（线性模型）；仅 `ml-train` 阶段且模型有重要性时生成。
- 回退/占位：MLP/torch 等无原生重要性时不生成，缺失不等于流程错误。

#### `fig_51_ml_roc.png`

- 坐标轴：横轴为假阳性率（`False positive rate`），纵轴为真阳性率（`True positive rate`）。
- 颜色/图例：蓝色 ROC 曲线，灰色虚线为随机基线，图例显示 `AUC=...`。
- 输入：分类任务测试集预测概率；回归任务不生成。
- 读取方法：AUC 接近 0.5 时模型不能用于优先排序。

#### `fig_52_knockout_top_candidates.png`

- 坐标轴：横轴为敲除优先级得分（`Knockout priority score`，固定 0-1），纵轴为基因名。
- 颜色/图例：条形统一绿色 `#1d6f42`。
- 输入：`fig_52_53_ranked_knockout.csv` 中多维评分（表达差异、增殖共表达、网络 hub、可选 DepMap 等）的加权 `knockout_score`。
- 读取方法：查看高优先级候选；分数是统计优先级，不等同于真实敲除表型。

#### `fig_53_knockout_score_distribution.png`

- 坐标轴：横轴为敲除优先级得分（0-1），纵轴为基因数量（`Gene count`）。
- 颜色/图例：直方图统一蓝色 `#4c7bb8`。
- 读取方法：判断评分是否有区分度和分层能力；评分全部相同或基因数过少时不能用于靶点分层。

#### `compound_disease_venn.png`

- 坐标轴：无笛卡尔坐标，为集合圆形布局。
- 颜色/图例：集合圆轮廓分别为蓝 `#4C72B0`、橙 `#DD8452`、绿 `#55A868`，标签为集合名称，区域数字为基因计数。
- 输入：化合物靶点来源与疾病基因集；交集数与 `compound_disease_overlap.csv` 一致。
- 读取方法：查看交集规模与来源分布；交集为 0 时只能说明无重叠，不能作为阳性证据。

#### `fig_54_feedback_module_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`。
- 颜色/图例：点颜色表示筛选靶点模块得分，浅灰 `grey90` 到深红 `#B31B1B`。
- 输入：`cell_feedback_summary.json` 状态为 completed 且存在模块得分列。
- 读取方法：查看模块是否集中在特定细胞群；全图均匀时不能说明细胞特异性。

#### `fig_55_feedback_target_expression_umap.png`

- 坐标轴：`UMAP_1`/`UMAP_2`，每个候选靶基因一个面板（最多 `max_features` 个，默认 8）。
- 颜色/图例：点颜色表示表达量，浅灰 `grey90` 到蓝 `#2E86AB`。
- 读取方法：查看单个候选基因的细胞类型表达位置；表达全为 0 时不可用。

#### `fig_56_feedback_celltype_dotplot.png`

- 坐标轴：横轴为候选靶基因，纵轴为细胞类型。
- 颜色/图例：点颜色表示平均表达量，点大小表示表达比例。
- 读取方法：比较候选基因的细胞类型特异性；所有基因在所有细胞类型中均为 0 时不可用。

#### `fig_57_feedback_celltype_boxplot.png`

- 坐标轴：横轴为细胞类型，纵轴为反馈模块得分（`Feedback module score`）。
- 颜色/图例：小提琴/箱线按条件（`condition`）填充，图例在右侧。
- 读取方法：比较细胞类型与条件间的模块活性；只有 1 个细胞类型时不生成，分布无差异时只能说明特异性弱。

#### `fig_58_feedback_celltype_heatmap.png`

- 坐标轴：横轴为细胞类型，纵轴为候选靶基因。
- 颜色/图例：白-浅蓝-蓝-深红渐变（`#f7fbff`→`#d6e7ff`→`#2E86AB`→`#7B241C`），颜色越深平均表达越高；行和列均聚类。
- 读取方法：查看候选基因表达谱的细胞类型模式；只有 1 个基因或 `pheatmap` 不可用时不生成，全图同色时不可用。

#### `fig_59_feedback_targets_volcano.png`

- 坐标轴：横轴为 `avg_log2FC`，纵轴为 `-log10(p_val_adj)`。
- 颜色/图例：显著上调为红 `#C0392B`，显著下调为蓝 `#2E86C1`，其余为灰 `grey65`；显著基因标注基因名。
- 阈值与参考线：灰色虚线为 `P = 0.05` 和 `log2FC = ±0.25`。
- 输入：反馈靶基因在条件间的差异表达表 `data/feedback_deg.csv`。

#### `fig_60_feedback_condition_violin.png`

- 坐标轴：横轴为条件，纵轴为归一化表达量。
- 颜色/图例：小提琴按条件分组（图例被移除），每个候选靶基因一个面板。
- 读取方法：查看反馈靶点在病例/对照间的表达差异；仅当条件数 >= 2 时生成。

#### `fig_61_feedback_go_network.png`/`fig_62_feedback_kegg_network.png`

- 坐标轴：无笛卡尔坐标，为通路-基因网络布局。
- 颜色/图例：基因节点按表达/差异方向着色，通路节点按类别着色（cnetplot 默认配色）。
- 输入：反馈靶基因 GO BP/KEGG 富集结果经 `p.adjust <= 0.05` 筛选后取 Top5 通路。
- 回退/占位：无显著通路时显示 “No significant pathway network”。

## 1. 快速判断流程

任何一张结果图都不要只看图片本身。建议按以下顺序判断：

1. 确认运行是否完整：
   - 单细胞：检查 `results/pipeline_complete.json`、`results/summary.json`、`results/result_report.html`。
   - 虚拟筛选：检查 `outputs/run_001/results/01_analysis/summary.json`、`outputs/run_001/results/docking_report.html` 和各阶段 `.done` 标记。
   - 全自动集成：检查 `outputs/integration/run_manifest.json`、`integration_summary.json`，以及 `cell_feedback/cell_feedback_summary.json`。
   - 网络毒理学：检查 `outputs/run_001/network_toxicology/network_toxicology_summary.json` 和 `compound_disease_overlap.csv`。
   - FAERS：检查 `outputs/run_001/faers/faers_summary.json` 和 `faers_signals.csv`。
2. 确认图片文件存在且非空，并能被图片查看器打开。
3. 确认对应状态文件没有 `skipped`、`failed`、`SHAP plot skipped` 等标记。
4. 打开与图片同名的结果数据表，核对样本数、细胞数、P 值、分数、命中数等是否与图片一致。
5. 核对本次运行使用的参数和阈值，例如 `config/docking_config.json` 中的 `analysis.cutoff`、`config/full_pipeline_config.json` 中的 `qc_gate`。
6. 做视觉检查，确认不是空白图、占位图、文字重叠图或明显截断图。
7. 只有状态、数据、参数和视觉检查都通过，才可以把图片作为支持结论的材料。

## 2. 统一判定口径

| 判定 | 含义 | 能否使用 |
| --- | --- | --- |
| 合格 | 文件存在且可打开，状态为完成，对应数据完整，视觉和统计条件符合预期 | 可以作为对应结论的支持证据 |
| 警示 | 文件存在但只是占位图、无显著信号、样本过少、模型性能弱、或存在回退逻辑 | 只能用于说明“未检出/不足”，不能作为阳性证据 |
| 不合格 | 文件缺失、为空、损坏、状态为 failed/skipped、数据对不上、或明显是旧运行残留 | 不应进入结论或交付材料 |

注意：流水线本身主要做文件存在性和非空校验，报告器会补充文件大小、图片尺寸和通用判读提示，但不会自动判断统计显著性、生物学合理性或图片视觉质量。最终是否可用需要人工复核。

## 3. 常见“文件存在但仍不可用”的情况

以下情况由代码行为产生，判断时尤其注意：

| 情况 | 典型表现 | 应如何处理 |
| --- | --- | --- |
| 占位图 | “No significant enrichment terms”“No significant GSEA terms”“No significant pathway network”“No published annotations”“At least two conditions are required” | 只能说明该分析未得到信号，不能当作显著富集、通路或注释一致性证据 |
| 回退图 | 无显著 DEG 时，DEG 图可能使用全部有限 P 值基因而非显著基因 | 不能把这些基因称为“显著差异基因” |
| 状态失败 | `ml_model_summary.json` 为 `failed`，或阶段标记缺失 | 相关 ML 图不可用 |
| 状态跳过 | ML 样本不足、CellChat 未安装、拟时序未开启、单样本不生成样本 UMAP | 缺失是正常结果，不代表图片损坏 |
| 可选图未开启 | `LIVER_RUN_TRAJECTORY=no`、`LIVER_RUN_CELLCHAT` 未设为 `yes` | 对应图不生成，不应要求必须存在 |
| 跳过列表 | 设置 `LIVER_SKIP_FIGURES` 后对应图缺失 | 需确认跳过原因后再判断 |
| 旧运行残留 | 参数已改变但旧图片仍留在目录中 | 应以 `run_manifest.json`、阶段签名、文件修改时间和本轮输出为准 |
| 数据支撑不足 | 只有一个样本、只有一个条件、样本数过少、命中数为 0 | 图可能仍生成，但只能作描述性材料，不能作统计结论 |

## 4. 单细胞分析结果图

### 4.1 输出位置和状态文件

结果图目录：

```text
<output_root>/results/figures/
  01_qc/
  02_doublets/
  03_cluster/
  04_annotation/
  05_deg/
  06_enrichment/
  07_ml/
  08_publication/
  09_cellchat/
```

对应数据目录为 `<output_root>/results/data/`，目录划分与图片一致。关键状态文件包括：

- `results/pipeline_complete.json`：流水线完成标记。
- `results/summary.json`：细胞数、聚类数、DEG 数等汇总。
- `results/data/07_ml/ml_model_summary.json`：ML 是否完成、是否跳过、AUC、CV 得分。
- `results/data/07_ml/ml_shap_status.txt`：SHAP 图是否被跳过。
- `results/data/09_cellchat/cellchat_status.txt`：CellChat 是否运行完成。
- `results/checkpoints/`：Seurat 断点对象，用于复核图背后的数据。

可选开关：

| 环境变量 | 默认 | 影响的图 |
| --- | --- | --- |
| `LIVER_RUN_CELLCYCLE` | `yes` | `fig_26`、`fig_27` |
| `LIVER_RUN_CLUSTER_MARKERS` | `yes` | `fig_31`、`fig_32` |
| `LIVER_RUN_SIGNATURES` | `yes` | `fig_33`、`fig_34`、`fig_35` |
| `LIVER_RUN_CNV` | `yes` | `fig_36` |
| `LIVER_RUN_SINGLER` | `yes` | `fig_37`、`fig_38` |
| `LIVER_RUN_TRAJECTORY` | `no` | `fig_39` |
| `LIVER_RUN_CELLCHAT` | 未默认开启 | `fig_40` 至 `fig_42` |
| `LIVER_SKIP_FIGURES` | 空 | 可跳过指定图片 |

### 4.2 质控阶段

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_01_qc_raw_violin.png` | 过滤前 `nFeature_RNA`、`nCount_RNA`、`percent.mt` 按条件展示的小提琴图，用于判断低质量细胞和条件间质量差异 | 三个指标的小提琴图均可见，细胞数足够，图形不是空白或全部为零 | 只有一个条件、只有一个样本或严重批次混杂时不能用于“条件差异”结论 |
| `fig_01_qc_filtered_violin.png` | 过滤后相同 QC 指标分布，用于确认过滤阈值是否合理 | 低质量尾部被去除，主要细胞群仍保留，两个条件仍可比较 | 过滤后细胞数过少、分布被过度压缩，或图像空白时不可用 |
| `fig_48_qc_pvalue_comparison.png` | 原始/过滤后 QC 指标在条件间的 Wilcoxon P 值对比图，用于量化过滤是否改变条件差异 | 图中显示真实 P 值，`fig_48_qc_pvalue_comparison.csv` 存在且指标、P 值可对应 | 图中出现 “At least two conditions are required” 时，只能说明没有两组条件，不能用于条件差异判断 |
| `fig_02_doublet_scores.png` | `scDblFinder` 双细胞得分按 singlet/doublet 分类展示，用于判断双细胞去除边界 | singlet 与 doublet 得分有区分度，或可明确看到无 doublet；对应 `fig_02_doublet_results.csv` 可核对 | 如果 `scDblFinder` 失败后所有细胞被标记为 singlet，图不能作为真实双细胞检测结果 |

### 4.3 聚类与降维

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_14_pca.png` | PCA 按条件着色，观察主要变异来源和批次/条件分布 | 细胞点可见，图例和坐标轴正确，不是空白图 | 单样本、单条件或点过少时只能作描述性展示 |
| `fig_15_elbow.png` | 主成分标准差 Elbow 图，用于选择 PCA 维度 | 有清晰曲线和拐点，或至少能解释所用维度选择 | 无拐点不代表损坏，但需说明维度选择依据 |
| `fig_03_umap_clusters.png` | Seurat 聚类在 UMAP 上的分布，用于检查聚类结构 | 聚类边界可读，不是大量随机碎片或全部混成一片 | 过度聚类、欠聚类或点太少时不能直接用于细胞类型注释 |
| `fig_04_umap_condition.png` | UMAP 按条件着色，用于评估条件分离与批次效应 | 条件分布与细胞类型结构可解释，不是完全由批次主导 | 单条件图不能支持条件差异；条件完全分离时需确认是否来自批次 |
| `fig_05_umap_annotation.png` | marker 注释与发表注释并排 UMAP，用于验证注释一致性 | 同类型细胞形成较连贯区域，标签清晰可见 | 注释与 UMAP 结构明显冲突时需要回到 marker 图和混淆矩阵复核 |
| `fig_06_dotplot_markers.png` | marker 基因在各细胞类型中的表达比例和表达量，用于验证注释 | 预期 marker 在对应细胞类型中高表达，其他类型低表达 | 所有点大小/颜色无差异、基因缺失或图像空白时不可用 |
| `fig_16_featureplot_markers.png` | marker 基因在 UMAP 上的表达位置 | 表达信号集中在预期细胞类型区域，不是全图均匀灰色 | 数据集中无可用 marker 基因时不生成；基因未匹配或表达全为 0 时不能验证注释 |
| `fig_17_marker_violin.png` | marker 基因在细胞类型中的表达分布 | 预期细胞类型表达明显更高，分布可见 | 同 `fig_16`，无可用 marker 基因时不生成；各细胞类型分布完全一致时不能作为注释支持 |
| `fig_18_celltype_proportion.png` | 细胞类型比例按条件堆叠，用于观察组成变化 | 比例柱可见，条件完整，细胞数足够 | 某些细胞类型细胞数过少时比例不可靠 |
| `fig_19_condition_proportion.png` | 条件构成按细胞类型堆叠，用于观察每个细胞类型中的条件比例 | 条件完整，比例可读，能对应 `fig_18_19_celltype_proportion_stats.csv` | 小样本、单条件或缺失细胞类型时需谨慎 |
| `fig_07_annotation_confusion_heatmap.png` | marker 注释与发表注释混淆矩阵热图 | 对角线计数较高，能看出主要细胞类型对应关系 | 图中出现 “No published annotations” 时只能说明缺少发表注释，不能用于一致性判断 |

### 4.4 差异表达

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_08_volcano.png` | log2FC 与 -log10 校正 P 值火山图，用于查看差异基因数量和方向 | 有 Up/Down/NS 点分布，阈值线、图例和基因标签可读；与 `fig_08_deg_all.csv` 一致 | 无显著基因时仍会生成图，但只能说明未检出；辅助线固定为 0.05 和 ±0.25，若配置阈值不同，应以 CSV 的 `significant` 列为准 |
| `fig_09_deg_horizontal_violin.png` | 按校正 P 值排序的 Top DEG 横向小提琴图，并标注 P 值，用于优先查看差异最显著基因 | 展示基因数量符合 `LIVER_DE_VIOLIN_TOP_N`，P 标签可见，两组表达分布有差异；对应 `fig_09_deg_horizontal_violin.csv` | 无显著 DEG 时可能回退到非显著基因，不能把这些基因称为显著差异基因 |
| `fig_09_deg_heatmap.png` | Top DEG 热图，用于查看差异基因整体表达模式 | 热图能按条件看到表达块，基因和样本标签可读 | 无显著 DEG 时 Top30 可能来自全部基因排序，需结合 `fig_09_deg_significant.csv` 判断 |

### 4.5 富集分析

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_10_go_up.png` | 上调基因 GO BP 富集图 | 有显著条目，P 值、计数、条目名称可读 | 图中出现 “No significant enrichment terms” 或 CSV 为 `note=no significant GO terms` 时不能作为阳性富集证据 |
| `fig_11_go_down.png` | 下调基因 GO BP 富集图 | 同上 | 同上 |
| `fig_12_kegg_up.png` | 上调基因 KEGG 富集图 | 有显著 KEGG 条目，P 值、基因数可读 | 无显著 KEGG 条目时只能说明未检出 |
| `fig_13_kegg_down.png` | 下调基因 KEGG 富集图 | 同上 | 同上 |
| `fig_20_gsea_go.png` | GO BP GSEA 富集图，用于查看通路整体上调/下调方向 | 有 GSEA 条目，NES 或富集方向可读 | “No significant GSEA terms” 时不能作为显著通路证据 |
| `fig_21_gsea_kegg.png` | KEGG GSEA 富集图 | 同上 | 同上 |
| `fig_22_go_network.png` | 上调基因 GO 经 `p.adjust <= 0.05` 筛选后 Top5 核心 + 最佳 5 延伸通路网络图 | 网络节点和基因连接可读，至少一个显著通路 | “No significant pathway network” 时只能说明无显著通路 |
| `fig_23_kegg_network.png` | 上调基因 KEGG 经 `p.adjust <= 0.05` 筛选后 Top5 核心 + 最佳 5 延伸通路网络图 | 同上 | 同上 |
| `fig_46_go_top5.png` | 上调基因 GO 经 `p.adjust <= 0.05` 筛选后的 Top5 富集气泡图 | 最多 5 条通路，P 值、GeneRatio、Count 可读 | “No significant pathway after filtering” 时不能作为阳性富集证据 |
| `fig_47_kegg_top5.png` | 上调基因 KEGG 经 `p.adjust <= 0.05` 筛选后的 Top5 富集气泡图 | 同上 | 同上 |

### 4.6 机器学习分类

使用本组图之前必须先检查 `results/data/07_ml/ml_model_summary.json`。

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_24_ml_feature_importance.png` | ML 模型特征重要性，用于筛选关键样本级特征 | 模型状态为 `completed`，重要性值非负，标签可见 | 状态为 `skipped` 或 `failed` 时不可用；旧运行残留图不可用 |
| `fig_25_ml_shap.png` | SHAP 解释，用于查看特征对预测方向和贡献 | SHAP 图成功生成，`ml_shap_status.txt` 没有 “SHAP plot skipped” | SHAP 被跳过时不可用 |
| `fig_43_ml_confusion_matrix.png` | 交叉验证混淆矩阵，用于查看具体错分 | 类别标签正确，对角线数量清楚，样本数可对应 | 样本量很小或模型状态失败时只能谨慎解释 |
| `fig_44_ml_roc_pr.png` | ROC 与 PR 曲线，用于评估分类能力 | 曲线高于随机线，AUC/AP 可见，不是全部等于随机水平 | AUC 接近 0.5、样本极少或类别不平衡时不能用于强结论 |
| `fig_45_ml_cv_scores.png` | 交叉验证准确率分布，用于评估模型稳定性 | 有 CV 分数点，均值和离散程度可读 | CV 分数波动过大或状态失败时不可用 |
| `fig_45_ml_calibration_curve.png` | 校准曲线，用于检查预测概率是否可靠 | 校准点接近对角线，样本量足够 | 仅二分类模型生成；非二分类、样本极少、曲线无法计算或模型状态失败时不可用 |
| `data/07_ml/fig_24_ml_selected_features.csv` | `lasso_svm` 模型 LASSO 初筛 + SVM-RFE 选定的特征 | 状态为 `completed`，特征数 > 0 | 未选择 `lasso_svm` 时不生成，缺失正常 |

### 4.7 高级分析与发表图

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_26_cellcycle_umap.png` | 细胞周期阶段 UMAP，用于检查细胞周期对聚类的干扰 | 各 phase 分布在 UMAP 上可读 | 未开启 `LIVER_RUN_CELLCYCLE` 或 S/G2M 基因不足时不生成；如果 phase 形成强分群，说明需要评估周期回归，图本身不一定损坏 |
| `fig_27_cellcycle_proportion.png` | 细胞周期阶段比例图 | 分组、细胞类型、phase 比例可读 | 同 `fig_26`，未开启 `LIVER_RUN_CELLCYCLE` 或无 `Phase` 时不生成；比例信息过少时需谨慎 |
| `fig_28_umap_sample.png` | UMAP 按样本着色，用于评估样本混合和批次效应 | 多样本时生成，样本标签可见，同一细胞类型跨样本混合情况可读 | 单样本时不生成，不缺失视为正常 |
| `fig_29_doublet_rate_sample.png` | 每样本双细胞率，用于识别异常样本 | 样本标签和双细胞率可读，与 `fig_29_doublet_rate_by_sample.csv` 一致 | 某样本双细胞率异常高时提示该样本质量风险 |
| `fig_30_sample_proportion.png` | 每样本细胞类型比例，用于检查样本异质性 | 样本和细胞类型比例可读 | 样本过少或比例全为 0 时不可用 |
| `fig_31_cluster_marker_heatmap.png` | 聚类 marker 热图，用于重新定义聚类标签 | 每个聚类有 marker 表达块，热图颜色可读 | 未开启 `LIVER_RUN_CLUSTER_MARKERS` 或无聚类 marker 时不生成，缺失正常 |
| `fig_32_cluster_marker_dotplot.png` | 聚类 marker DotPlot | marker 在对应聚类中表达比例/水平更高 | 同 `fig_31`，未开启 `LIVER_RUN_CLUSTER_MARKERS` 或无 marker 时不生成；全 0 时不可用 |
| `fig_33_signature_scores_umap.png` | 增殖、EMT、缺氧等签名得分 UMAP | 得分梯度可见，不是全图均匀 | 未开启 `LIVER_RUN_SIGNATURES` 或签名基因不足时不生成；签名基因未匹配或得分为常数时不能解释 |
| `fig_34_signature_scores_boxplot.png` | 功能签名得分按条件箱线图/小提琴图 | 条件和签名面板完整，分布可见 | 同 `fig_33`，未开启 `LIVER_RUN_SIGNATURES` 或签名基因不足时不生成；无差异不代表损坏，但不能当作“有差异”结论 |
| `fig_35_celltype_abundance_effect.png` | 细胞类型丰度变化的 log2OR 与 -log10 校正 P 值图 | 点、细胞类型标签和显著性信息可读，数据来自 `fig_18_19_celltype_proportion_stats.csv` | 细胞数过少、OR 不有限或 P 值缺失时需谨慎 |
| `fig_36_cnv_heatmap.png` | 基于染色体滑动窗口均值的推断 CNV 热图 | 细胞按行、染色体窗口按列，红蓝块可见，有条件和细胞类型注释 | 未开启 `LIVER_RUN_CNV`、注释库缺失、CHRLOC 不足或可映射基因 < 200 时不生成；基因注释不足、窗口过少或热图全为白色时不可用 |
| `fig_37_singler_umap.png` | SingleR 参考注释 UMAP | 同类型细胞形成连贯区域，标签可读 | 未开启 `LIVER_RUN_SINGLER`、参考数据不可用或预测失败时不生成，缺失正常 |
| `fig_38_singler_confusion_heatmap.png` | SingleR 与现有注释混淆矩阵热图 | 对角线可读，能看出注释一致性 | 同 `fig_37`，未开启 `LIVER_RUN_SINGLER` 或无 SingleR 预测时不生成；一致性差时应复核注释 |
| `fig_39_trajectory_umap.png` | slingshot 拟时序轨迹图 | 伪时序梯度或轨迹线可见，UMAP 背景可读 | 默认 `LIVER_RUN_TRAJECTORY=no` 时不生成，缺失正常；轨迹无结构时不能作分化结论 |

### 4.8 CellChat 细胞通讯

仅当 `LIVER_RUN_CELLCHAT=yes` 且 R 环境可用时生成。先检查 `results/data/09_cellchat/cellchat_status.txt`。

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_40_cellchat_network.png` | 细胞类型间通讯数量网络图 | 节点、连线、计数可读，有实际通讯 | CellChat 未运行、未安装或无通讯时不生成，缺失正常 |
| `fig_41_cellchat_heatmap.png` | 细胞类型对间通讯强度热图 | 热图有非零颜色块 | 全部为零或空白时不能作为通讯证据 |
| `fig_42_cellchat_bubble.png` | 配体受体对气泡图 | 气泡大小/颜色可读，通路名称可见 | 无显著配体受体对时只能说明未检出 |

## 5. 虚拟筛选结果图

### 5.1 输出位置和状态文件

独立虚拟筛选默认输出：

```text
<workdir>/outputs/run_001/results/
  01_analysis/figures/
  02_redock/figures/
  03_ml/figures/
  04_knockout/figures/
  05_validation/
```

对应状态和数据文件：

- `01_analysis/summary.json`：成功对接数、命中数、Top N、最佳亲和力、中位亲和力。
- `01_analysis/data/fig_46_47_ranked_results.csv`：全部成功对接排序结果。
- `01_analysis/data/fig_47_top_hits.csv`：低于 `analysis.cutoff` 的命中。
- `01_analysis/data/fig_48_diverse_hits.csv`：Tanimoto 多样性选择结果。
- `02_redock/data/fig_49_redock_results.csv` 和 `fig_49_redock_comparison.csv`：重对接结果。
- `03_ml/data/ml_model_info.json`：模型类型、任务类型、模型文件。
- `04_knockout/data/fig_52_53_ranked_knockout.csv`：敲除评分和靶点评分表。
- `docking_report.html`：HTML 报告，会汇总该报告目录下所有 PNG。

默认参数参考 `config/docking_config.json`：命中阈值为 -7.0 kcal/mol，Top N 为 100，多样性 Tanimoto 阈值为 0.7，重对接 Top 20，重对接 exhaustiveness 默认 32。

### 5.2 对接分析与命中

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_46_affinity_distribution.png` | 所有成功对接配体的亲和力分布，叠加命中阈值线，用于查看整体分数和命中分布 | `summary.json` 中 `total_docked > 0`，直方图非空，阈值线可见，亲和力为负值 | 无成功对接、图为空或 `best_affinity` 缺失时不可用；有命中才能说明命中分布 |
| `fig_47_top_hits.png` | 排序后的 Top 命中条形图，用于查看排名靠前配体 | 至少 1 个配体，ID 和亲和力可读，顺序与 CSV 一致 | 若 `hits=0`，图中只是“Top 排序”而不是“Top 命中”，不能写成命中结果 |
| `fig_48_diverse_hits.png` | 多样性选择后的命中条形图，用于减少同一化学骨架重复 | `fig_48_diverse_hits.csv` 非空，ID/SMILES 不重复，与 Top 命中可对应 | RDKit 不可用、有效分子少于 2 个或多样性选择未运行时可能不生成；只有一个骨架时不能称为多样 |

### 5.3 重对接

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_49_redock_comparison.png` | 初始亲和力与重对接亲和力的散点图，以及变化量直方图，用于检查结果稳定性 | 至少 2 个可比较配体，点接近对角线，变化直方图围绕 0 分布，无系统性偏移 | 未启用重对接阶段（`redock.enabled=false`）、重对接全部失败或初始/重对接亲和力无法合并时不生成；只有 1 个配体时可能不生成；出现明显系统性漂移时应视为重对接不稳定，不能直接使用该批排序 |

### 5.4 ML/DL 重打分

仅 `ml-train` 阶段生成。先检查 `03_ml/data/ml_model_info.json`。

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_50_ml_feature_importance.png` | ML 重打分模型的特征重要性，用于解释哪些分子特征影响分数 | 模型有 `feature_importances_`，图非空，特征名可读 | 仅 `ml-train` 阶段且模型提供 `feature_importances_`/`coef_` 时生成；MLP/torch 等没有原生特征重要性时可能不生成；缺失不等于流程错误 |
| `fig_51_ml_roc.png` | ML 重打分分类模型的 ROC 曲线 | 任务为分类，曲线高于随机线，AUC 可见且明显大于 0.5 | 仅 `ml-train` 阶段、分类任务且模型可输出概率时生成；回归任务不生成；AUC 接近 0.5 时模型不能用于优先排序 |

### 5.5 虚拟敲除与靶点优先级

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_52_knockout_top_candidates.png` | Top N 基因的 `knockout_score` 条形图，用于快速查看候选基因 | 至少 5 个基因参与评分，分数在 0 到 1 之间，基因名可读，排序与 CSV 一致 | 基因数不足、分数全部为 0 或图为空时不可用；分数是统计优先级，不是真实敲除表型 |
| `fig_53_knockout_score_distribution.png` | 全部基因 `knockout_score` 分布直方图，用于判断评分是否有区分度 | 直方图非空，评分不是单一值，`genes_scored` 与 CSV 一致 | 评分全部相同或样本数过少时不能用于靶点分层 |

注意：虚拟敲除分数由表达差异、增殖共表达、网络 hub 和可选 DepMap 数据计算，附加维度会进一步生成 `target_score`。该分数是实验前筛选启发式指标，不等同于基因敲除的机制模拟。

## 6. 全自动集成流水线与细胞反馈结果图

### 6.1 输出位置和状态文件

全自动流水线结果：

```text
<workdir>/outputs/integration/
  key_genes.csv
  gene_evidence.csv
  knockout_summary.json
  docking_targets.csv
  cell_feedback/
    cell_feedback_summary.json
    manifest.csv
    figures/
      fig_54_feedback_module_umap.png
      fig_55_feedback_target_expression_umap.png
      fig_56_feedback_celltype_dotplot.png
      fig_57_feedback_celltype_boxplot.png
      fig_58_feedback_celltype_heatmap.png
      fig_59_feedback_targets_volcano.png
      fig_60_feedback_condition_violin.png
      fig_61_feedback_go_network.png
      fig_62_feedback_kegg_network.png
    data/
      feedback_targets.csv
      celltype_summary.csv
      celltype_enrichment.csv
      condition_summary.csv
      feedback_deg.csv
      feedback_enrichment_go.csv
      feedback_enrichment_kegg.csv
  integration_report.html
  integration_summary.json
  run_manifest.json
```

`gene_evidence.csv` 现在除 UniProt/PDB/ChEMBL 外，还包含 STRING 互作伙伴、Reactome 通路、PharmGKB 注释、AlphaFold 结构、Open Targets 靶点关联、KEGG 通路及 `database_sources` 来源覆盖列。

集成流水线中的单细胞图仍位于单细胞输出根目录；独立对接结果图位于每个靶点的工作目录 `<workdir>/work/<gene>/outputs/run_001/results/`。集成报告直接链接 `fig_54`，其余图需到各自结果目录查看。

先检查 `cell_feedback_summary.json`：

- `status=skipped`：可能因 `--skip-cell-feedback`、manifest 为空或找不到 Seurat 对象。
- `status=completed`：继续检查 `genes_matched`、`n_celltypes`、`figures`。

默认参数：`feedback_top_n=12`，`feedback_max_features=8`。

### 6.2 细胞反馈图

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_54_feedback_module_umap.png` | 筛选靶点模块得分在 UMAP 上的分布，用于查看模块是否集中在特定细胞群 | `status=completed`，模块得分列存在，UMAP 有梯度而非全图均匀 | 找不到模块得分或没有 UMAP 时不生成；全图均匀时不能说明细胞特异性 |
| `fig_55_feedback_target_expression_umap.png` | 候选靶基因在 UMAP 上的表达，用于查看单个基因表达位置 | 至少 1 个基因匹配，每基因面板可读，表达不是全为 0 | 基因全部未匹配时不生成；全灰/全 0 时不可用 |
| `fig_56_feedback_celltype_dotplot.png` | 候选靶基因在细胞类型中的表达比例和表达量 | 基因和细胞类型可读，存在非零表达模式 | 所有基因在所有细胞类型中均为 0 时不可用 |
| `fig_57_feedback_celltype_boxplot.png` | 模块得分按细胞类型和条件分布，用于查看细胞类型特异性 | 至少 2 个细胞类型，箱线/小提琴分布可见 | 只有 1 个细胞类型时不生成；分布无差异时只能说明特异性弱 |
| `fig_58_feedback_celltype_heatmap.png` | 候选基因平均表达按细胞类型聚类热图 | 至少 2 个基因，热图有颜色结构，行/列聚类可见 | 只有 1 个基因或 `pheatmap` 不可用时不生成；全图同色时不可用 |
| `fig_59_feedback_targets_volcano.png` | 反馈靶基因在病例/对照条件间的差异表达火山图 | 至少 2 个条件且每组有足够细胞，基因名/方向/P 值阈值可读 | 条件不足时不生成；无显著基因时只能说明反馈靶点在该数据集中没有显著差异 |
| `fig_60_feedback_condition_violin.png` | 反馈靶基因按条件分组的小提琴图 | 基因与条件可读，存在非零表达模式 | 条件不足时不生成；全为 0 时不可用 |
| `fig_61_feedback_go_network.png` | 反馈靶基因 GO BP 富集经 `p.adjust <= 0.05` 筛选后的 Top5 通路-基因网络图 | 节点和基因连接可读，能看出 Top5 通路关联基因 | “No significant pathway network” 时只能说明无显著通路 |
| `fig_62_feedback_kegg_network.png` | 反馈靶基因 KEGG 富集经 `p.adjust <= 0.05` 筛选后的 Top5 通路-基因网络图 | 同上 | 同上 |

细胞反馈阶段的富集 Top5 图使用与 `fig_22_go_network.png` 相同的 `cnetplot` 通路-基因网络类型，不再输出气泡图；对应富集明细表仍保留在 `feedback_enrichment_go.csv` / `feedback_enrichment_kegg.csv` 中。

### 6.3 网络毒理学与 FAERS 信号

这两个分析由虚拟筛选页或独立命令运行，不属于全自动流水线阶段，输出位于：

```text
<workdir>/outputs/run_001/network_toxicology/
  network_toxicology_summary.json
  figures/compound_disease_venn.png
  data/compound_disease_overlap.csv
  data/ppi_hub_scores.csv
  data/ctpd_nodes.csv
  data/ctpd_edges.csv
  data/ctpd_network.html
<workdir>/outputs/run_001/faers/
  faers_summary.json
  data/faers_signals.csv
  data/faers_signals.html
```

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `compound_disease_venn.png` | 化合物靶点与疾病基因的 Venn 图，用于查看交集规模和来源 | 图片可打开，交集基因数与 `compound_disease_overlap.csv` 一致，来源标签清晰 | 只有 1 个基因集或输入文件缺失时不生成；交集为 0 时只能说明无重叠，不能作为阳性证据 |
| `compound_disease_overlap.csv` | 核心交集靶点及来源数据库计数 | 至少 1 个交集基因，`n_sources` 可读 | 交集为 0 或输入表不完整时需复核数据库下载 |
| `ppi_hub_scores.csv` | STRING PPI 的 degree、betweenness、clustering 与 hub 评分 | 提供 PPI 边表时生成，基因名可匹配 | 未提供 `--ppi-network-csv` 时不生成；匹配率过低时 hub 评分代表性不足 |
| `ctpd_network.html` | C-T-P-D 网络可视化页 | HTML 可打开，节点和边数与 CSV 一致 | 输入不完整时不生成 |
| `faers_signals.csv` | ROR/PRR/BCPNN/EBGM 信号表 | 组合数 > 0，计数列被正确汇总，信号判定可解释 | 事件表为空、药物/事件列名错误或计数未识别时不生成 |
| `faers_signals.html` | Top FAERS 信号浏览页 | HTML 可打开，数字与 CSV 一致 | 无信号时只说明当前阈值下未检出 |

注意：FAERS 的 BCPNN IC 与 EBGM 使用常用近似公式，适合筛选，不作为正式药物警戒统计结论。

## 7. 结果图的主要用途

### 7.1 单细胞分析

- QC 图用于决定过滤阈值、发现异常样本/条件质量偏移，并为后续差异表达是否可靠提供依据。
- 聚类和 UMAP 图用于检查分群结构、批次效应和是否过度聚类。
- Marker 图和注释混淆图用于判断细胞类型注释是否可信。
- 火山图、DEG 热图和小提琴图用于产生候选基因、判断差异方向和幅度。
- GO/KEGG/GSEA 图用于把差异基因转化为生物学通路解释。
- ML 图用于判断基于样本组成的分类模型是否可用，以及哪些细胞类型/QC 特征驱动分类。
- 发表级分析图用于检查细胞周期混杂、样本混合、CNV、SingleR、拟时序和细胞通讯。

### 7.2 虚拟筛选

- 亲和力分布图用于判断整体对接分数是否合理、命中阈值是否合适。
- Top 命中图用于快速选择进入精细重对接或 ML 重打分的配体。
- 多样性图用于避免同一化学骨架重复占用候选名额。
- 重对接图用于检查初筛排序是否稳定。
- ML 图用于判断模型是否具备排序能力，而不是直接替代湿实验活性。
- 敲除和靶点优先级图用于从差异表达、网络、DepMap 和多维证据中选出候选靶点。

### 7.3 全自动集成

- 细胞反馈图把虚拟筛选/虚拟敲除结果映射回单细胞图谱，用于判断候选基因表达在哪些细胞类型中，是否具备细胞类型特异性。
- `feedback_targets.csv` 中的 `cell_support_score` 可把筛选优先级和细胞表达特异性合并，用于下一轮靶点收敛。
- 集成报告用于汇总 QC 门控、差异丰度、关键基因、敲除、对接和细胞反馈状态。

### 7.4 网络毒理学与 FAERS

- 网络毒理学 Venn 图用于判断化合物靶点与疾病基因的交集规模，交集表用于确定后续 PPI/富集/对接的核心基因集合。
- PPI hub 评分用于把 STRING 拓扑信息纳入候选靶点排序，配合虚拟敲除可形成“网络毒理学 → 靶点优先级”的闭环。
- FAERS 信号表用于从药物不良反应报告中发现风险信号，不能单独证明因果关系，需要结合机制和临床证据复核。

## 8. 建议验收清单

正式使用或交付结果前，建议逐项确认：

1. 本轮 `run_manifest.json` 存在，且包含输入文件哈希、软件版本和参数。
2. 对应阶段状态文件为完成，没有 skipped/failed。
3. 每个要使用的图片文件存在、非空、可打开，且不是占位图。
4. 图片对应的 CSV/JSON 存在，字段和数量与图片一致。
5. 使用的阈值已记录，例如 DEG 的 `logFC`/`padj`、对接的 `cutoff`、ML 的 AUC。
6. 图片没有明显空白、重叠、截断、全零或全图均匀。
7. 阳性结论必须有对应统计支持，例如显著富集条目、显著 DEG、命中数、AUC、细胞类型表达特异性。
8. 缺失的可选图能由状态文件或参数解释，而不是静默遗漏。
9. 报告中的图片链接指向实际存在的文件。
10. 对不确定或统计支撑不足的图，明确标注为“警示/不可用”，不夸大结论。

## 9. 常见误用

- 把“生成了图片”等同于“结果合格”：图片可能只是占位图或回退图。
- 把无显著富集条目写成显著富集通路。
- 把无显著 DEG 时回退生成的基因图写成显著差异基因。
- 把 AUC 接近 0.5 的 ML 图用于强排序结论。
- 把只有 1 个样本/1 个条件的图写成条件差异证据。
- 把敲除优先级分数解释为真实敲除表型或湿实验验证结果。
- 使用参数或输入已经改变后的旧结果图，而未核对运行标记和 manifest。
