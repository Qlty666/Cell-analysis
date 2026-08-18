# 现有流程结果图判断与使用指南

本文档用于说明本项目现有流水线实际能够生成哪些结果图、如何判断结果图是否合格、是否可以用于结论，以及每张图的主要用途。

适用范围：

- 单细胞分析：`scripts/run_pipeline.py`，结果图位于 `<output_root>/results/figures/`。
- 虚拟筛选：`scripts/run_docking.py`，结果图位于 `<workdir>/outputs/run_001/results/` 下各阶段目录。
- 全自动集成流水线：`scripts/run_full_pipeline.py`，除复用上述两类结果外，还会在 `<workdir>/outputs/integration/cell_feedback/figures/` 生成细胞反馈结果图。

说明：README 中“48 张”是历史常见口径。按当前代码逐点核对，单细胞流程的图片保存点共有 50 个文件名，其中部分图片由环境变量、数据条件或可选阶段决定是否生成。判断时应以本清单、实际输出目录、阶段状态文件和对应数据文件为准。

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
| `fig_16_featureplot_markers.png` | marker 基因在 UMAP 上的表达位置 | 表达信号集中在预期细胞类型区域，不是全图均匀灰色 | 基因未匹配或表达全为 0 时不能验证注释 |
| `fig_17_marker_violin.png` | marker 基因在细胞类型中的表达分布 | 预期细胞类型表达明显更高，分布可见 | 各细胞类型分布完全一致时不能作为注释支持 |
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
| `fig_22_go_network.png` | 上调基因 GO 经 `p.adjust <= 0.05` 筛选后 Top5 通路网络图 | 网络节点和基因连接可读，至少一个显著通路 | “No significant pathway network” 时只能说明无显著通路 |
| `fig_23_kegg_network.png` | 上调基因 KEGG 经 `p.adjust <= 0.05` 筛选后 Top5 通路网络图 | 同上 | 同上 |
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

### 4.7 高级分析与发表图

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_26_cellcycle_umap.png` | 细胞周期阶段 UMAP，用于检查细胞周期对聚类的干扰 | 各 phase 分布在 UMAP 上可读 | 如果 phase 形成强分群，说明需要评估周期回归，图本身不一定损坏 |
| `fig_27_cellcycle_proportion.png` | 细胞周期阶段比例图 | 分组、细胞类型、phase 比例可读 | 比例信息过少时需谨慎 |
| `fig_28_umap_sample.png` | UMAP 按样本着色，用于评估样本混合和批次效应 | 多样本时生成，样本标签可见，同一细胞类型跨样本混合情况可读 | 单样本时不生成，不缺失视为正常 |
| `fig_29_doublet_rate_sample.png` | 每样本双细胞率，用于识别异常样本 | 样本标签和双细胞率可读，与 `fig_29_doublet_rate_by_sample.csv` 一致 | 某样本双细胞率异常高时提示该样本质量风险 |
| `fig_30_sample_proportion.png` | 每样本细胞类型比例，用于检查样本异质性 | 样本和细胞类型比例可读 | 样本过少或比例全为 0 时不可用 |
| `fig_31_cluster_marker_heatmap.png` | 聚类 marker 热图，用于重新定义聚类标签 | 每个聚类有 marker 表达块，热图颜色可读 | 无聚类 marker 时不生成，缺失正常 |
| `fig_32_cluster_marker_dotplot.png` | 聚类 marker DotPlot | marker 在对应聚类中表达比例/水平更高 | 无 marker 或全 0 时不可用 |
| `fig_33_signature_scores_umap.png` | 增殖、EMT、缺氧等签名得分 UMAP | 得分梯度可见，不是全图均匀 | 签名基因未匹配或得分为常数时不能解释 |
| `fig_34_signature_scores_boxplot.png` | 功能签名得分按条件箱线图/小提琴图 | 条件和签名面板完整，分布可见 | 无差异不代表损坏，但不能当作“有差异”结论 |
| `fig_35_celltype_abundance_effect.png` | 细胞类型丰度变化的 log2OR 与 -log10 校正 P 值图 | 点、细胞类型标签和显著性信息可读，数据来自 `fig_18_19_celltype_proportion_stats.csv` | 细胞数过少、OR 不有限或 P 值缺失时需谨慎 |
| `fig_36_cnv_heatmap.png` | 基于染色体滑动窗口均值的推断 CNV 热图 | 细胞按行、染色体窗口按列，红蓝块可见，有条件和细胞类型注释 | 基因注释不足、窗口过少或热图全为白色时不可用 |
| `fig_37_singler_umap.png` | SingleR 参考注释 UMAP | 同类型细胞形成连贯区域，标签可读 | 参考数据不可用或预测失败时不生成，缺失正常 |
| `fig_38_singler_confusion_heatmap.png` | SingleR 与现有注释混淆矩阵热图 | 对角线可读，能看出注释一致性 | 无 SingleR 预测时不生成；一致性差时应复核注释 |
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
| `fig_49_redock_comparison.png` | 初始亲和力与重对接亲和力的散点图，以及变化量直方图，用于检查结果稳定性 | 至少 2 个可比较配体，点接近对角线，变化直方图围绕 0 分布，无系统性偏移 | 只有 1 个配体、重对接全部失败或 CSV 无法合并时不生成；出现明显系统性漂移时应视为重对接不稳定，不能直接使用该批排序 |

### 5.4 ML/DL 重打分

仅 `ml-train` 阶段生成。先检查 `03_ml/data/ml_model_info.json`。

| 文件 | 内容与用途 | 合格判据 | 不可用或警示 |
| --- | --- | --- | --- |
| `fig_50_ml_feature_importance.png` | ML 重打分模型的特征重要性，用于解释哪些分子特征影响分数 | 模型有 `feature_importances_`，图非空，特征名可读 | MLP/torch 等没有原生特征重要性时可能不生成；缺失不等于流程错误 |
| `fig_51_ml_roc.png` | ML 重打分分类模型的 ROC 曲线 | 任务为分类，曲线高于随机线，AUC 可见且明显大于 0.5 | 回归任务不生成；AUC 接近 0.5 时模型不能用于优先排序 |

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
    data/
      feedback_targets.csv
      celltype_summary.csv
      celltype_enrichment.csv
      condition_summary.csv
  integration_report.html
  integration_summary.json
  run_manifest.json
```

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
