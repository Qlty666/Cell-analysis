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

- 代码实现覆盖：95.77%
- 当前已执行结果覆盖：75.61%（旧 `y1` 未运行 MD，且未重建靶点 Panel）
- 配置必要外部数据和工具后的预计覆盖：97.85%
- 性能目标达标 Panel：39/42

当前未达到 90% 结果覆盖的主要原因是：

1. GeneCards/OMIM/TTD 授权导出缺失，疾病 Panel 1e 无法形成三数据库 Venn。
2. STITCH 对 6PPD-Q 无记录，ChEMBL 相似性没有映射到人类基因，化合物 Panel
   1d 只有 SwissTargetPrediction 一个可用来源。
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

## 仍然无法由代码单独保证的部分

- GeneCards、OMIM、TTD、DrugBank 等授权数据必须由用户提供本地导出。
- 外部 AUC 和校准是否达到方案阈值属于模型性能，不能由代码保证。如果最终
  必须达到 AUC >=0.80 且 H-L p>0.05，需要扩大训练样本、增加独立 NAFLD
  队列、改进终点定义，或在方案中预先规定替代性能标准。
- 100 ns GROMACS 和 gmx_MMPBSA 需要本机算力与外部工具；默认只准备输入。
- CellChat、Discovery Studio 和湿实验不能由 Python 脚本完全替代，只能提供
  经过审计的近似分析或可导入的外部工具交接文件。

## 结论

调整后的脚本在代码实现层面达到 95.77%，并在授权数据、NAFLD 外部队列和
GROMACS/gmx_MMPBSA 前提满足时达到 97.85% 的预计 Panel 覆盖。旧 `y1` 的
当前结果覆盖为 75.61%，不能代表新代码重新执行后的状态。若要求“当前已有
输出立即达到 90%”，必须重新运行 targets/disease、ML、single-cell、docking、
MD、classify 和 figure audit；其中 100 ns MD 和 gmx_MMPBSA 是不可省略的前提。
