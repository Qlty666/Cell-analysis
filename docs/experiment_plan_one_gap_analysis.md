# 实验方案一脚本覆盖度评估

## 评估范围

依据 `实验方案一 .docx` 的 Figure 1-5，共 42 个有效 Panel：

- Figure 1：8 个 Panel
- Figure 2：8 个 Panel
- Figure 3：8 个 Panel
- Figure 4：10 个 Panel
- Figure 5：8 个 Panel

本评估使用 `src/experiment_plan_one/coverage.py` 的加权覆盖矩阵，权重反映每个
Panel 对方案闭环的重要性，不把外部数据库授权、湿实验和计算资源当作代码缺陷。

## 原脚本水平

已有 `experiment_plan_one` 专用流水线可以完成数据下载、化合物表征、开放疾病
靶点、证据中心、PPI、bulk 差异分析、11 模型比较、SHAP、单细胞处理、对接和
GROMACS 输入准备。对 `y1` 现有结果的审计显示：

- 代码实现覆盖：96.26%
- 当前已执行结果覆盖：75.61%（旧 `y1` 未运行 MD，且未重建靶点 Panel）
- 配置必要外部数据和工具后的预计覆盖：98.33%
- 旧版性能目标计数：39/42；新评估把 Figure 3c/3d 列为“终点不同，不按
  健康/NAFLD 阈值评价”，Figure 3e 仍作为可评价目标单独报告。

当前未达到 90% 结果覆盖的主要原因是：

1. GeneCards/OMIM/TTD 授权导出缺失，疾病 Panel 1e 无法形成三数据库 Venn。
2. 旧 `y1` 中 STITCH 对 6PPD-Q 无记录，ChEMBL 相似性没有映射到人类基因，
   化合物 Panel 1d 只有 SwissTargetPrediction 一个可用来源；新代码已加入
   SEA 连接器和多本地预测表入口。
3. GSE49541 是纤维化分期，GSE164441 是 HCC 肿瘤/癌旁，不能替代健康/NAFLD
   终点；GSE135251 的 AUC 未达到 0.80。
4. Hosmer-Lemeshow 校准在当前 63 例训练集上未通过。
5. 原 Figure 5d-h 在旧的 `y1` 结果中没有真实 100 ns 轨迹；新代码已经能从
   GROMACS XVG/MM-PBSA 输出生成最终 Figure 5d-h，但需要设置 `md.run=true`。
6. Figure 5a/5b 只有 Cα 散点和近邻残基放射图，缺少相互作用类型和发布级
   3D 结构呈现。

## 已完成的结构性调整

- 新增 `plan_coverage` 加权覆盖模块，分别报告代码实现覆盖、当前结果覆盖、
  外部前提后的预计覆盖和性能目标未达标 Panel。
- 新增本地化合物靶点预测表入口，可把 SEA、SuperPred、TargetNet 或内部预测
  结果作为独立来源加入 Panel 1d Venn。
- 疾病阶段可补充 GWAS Catalog 和 ClinVar 开放证据，并在授权导出缺失时形成
  可审计的多来源交集，而不是伪造 GeneCards/OMIM/TTD。
- CellChat-like 分析改为细胞类型内条件标签置换、p 值和 BH-FDR，不再只报告
  简单的 ligand/receptor 最小值。
- 对接结果新增氢键、疏水、盐桥、芳香接触和 vdW 分类，并改进 3D 结合口袋
  与 ligand 呈现。
- 新增 `md_figures.py`，从真实 GROMACS XVG 和 MM-PBSA 输出生成 Figure
  5d-h；没有轨迹时明确生成 not-run 面板，不制造数值。
- Figure 4g 的 Top10 变化图改为从虚拟敲除源数据按 600 dpi 重绘。
- 新增 SEA 化合物靶点预测连接器（开放 API、缓存、超时和失败状态），并补充
  多本地预测表入口，使 ChEMBL/STITCH/SwissTargetPrediction 之外可利用
  SEA 或 PharmMapper 导出，形成可审计的真实多来源 Venn。
- 新增参考文献常用的 WGCNA 风格共表达模块筛查：方差/MAD 预筛、软阈值、
  模块-疾病关联、kME 和 hub 输出。该方法明确标记为本地可复现的近似实现，
  不冒充 R WGCNA 的 TOM 和动态树切割。
- 机器学习从单次 5 折 CV 扩展为可配置重复分层 CV，并新增决策曲线、重复
  区间和每个外部队列的校准图。GSE49541 与 GSE164441 因终点不同，不再
  套用健康/NAFLD 的 AUC 阈值，而是单独报告实际终点；GSE135251 作为
  补充 NAFLD 终点保留阈值评价。
- 实验参数按参考文献收紧：AutoDock Vina 默认 exhaustiveness=100、
  num_modes=20、3 个重复种子；MD 明确使用 100 ns、310 K、1 bar、
  TIP3P、0.15 M ions，并将实验方案一默认蛋白力场改为 `amber14sb`。

## 仍然无法由代码单独保证的部分

- GeneCards、OMIM、TTD、DrugBank 等授权数据必须由用户提供本地导出。
- 外部 AUC 和校准是否达到方案阈值属于模型性能，不能由代码保证。GSE49541
  和 GSE164441 只评价各自的纤维化/HCC 终点，不能据其宣称健康/NAFLD
  分类性能；真正可评价的补充 NAFLD 终点目前只有 GSE135251。如果最终
  必须达到 AUC >=0.80 且 H-L p>0.05，需要扩大训练样本、增加独立 NAFLD
  队列、改进终点定义，或在方案中预先规定替代性能标准。
- 100 ns GROMACS 和 gmx_MMPBSA 需要本机算力与外部工具；默认只准备输入。
- CellChat、Discovery Studio 和湿实验不能由 Python 脚本完全替代，只能提供
  经过审计的近似分析或可导入的外部工具交接文件。

## 当前机器环境差异

本机自动检查到的版本与方案文字存在差异，必须在论文方法中记录或更换版本：

- AutoDock Vina：本机为 `1.2.7`，方案写的是 `1.2.3`。
- GROMACS：本机为 `2020.6-MODIFIED`，方案写的是 `GROMACS 2022`。
- `gmx_MMPBSA`：当前未安装，Figure 5h 无法完成真实自由能分解。
- R `CellChat`：当前未安装，Figure 4f 使用本仓库的置换近似；如要求原版
  CellChat，需安装 R 包并保留完整参数。
- PoseBusters：当前 Python 环境未安装，结构合理性独立复核需要使用其他环境。

## 结论

调整后的脚本在代码实现层面达到 96.26%，并在授权数据、NAFLD 外部队列和
GROMACS/gmx_MMPBSA 前提满足时达到 98.33% 的预计 Panel 覆盖。旧 `y1` 的
当前结果覆盖为 75.61%，不能代表新代码重新执行后的状态。若要求“当前已有
输出立即达到 90%”，必须重新运行 targets/disease、bulk/coexpression、ML、
single-cell、docking、MD、classify 和 figure audit；其中 100 ns MD 和
gmx_MMPBSA 是不可省略的前提。性能目标现在单独报告“可评价/未评价”，避免把
纤维化分期或 HCC 终点错误地当成 NAFLD 分类达标或未达标。
