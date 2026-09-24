# 实验方案一脚本覆盖度评估

> 说明：本文件保留为历史覆盖度评估。权重百分比不能表示科学结论完整或投稿级，当前脚本冻结、交付分级和状态口径以
> `docs/experiment_plan_one_script_freeze_20260923.md` 及运行产物
> `10_reports/governance/delivery_readiness.json` 为准。

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
- Figure 4f 默认恢复 R CellChat，将细胞级概率与独立生物单位统计分开；旧显式
  ligand/receptor 评分仅保留为显式后备方法，不再冒充 CellChat。
- 对接结果优先调用 PLIP 3.0.1 生成验证后的相互作用表；缺少 PLIP 结果时才
  回到距离/几何候选层级，并保留证据边界。
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
- 100 ns GROMACS 和 gmx_MMPBSA 需要可用算力、兼容版本和真实轨迹；当前工具
  已发现，但只应生成准备状态和 `not_run` 面板，不能由代码补造。
- PLIP 可以验证计算构象中的相互作用，但不能替代湿实验直接结合证据；
  Discovery Studio 只属于可选复核工具。

## 当前机器环境差异

本轮重新核查后的状态如下；工具可用不等于真实结果已经生成：

- AutoDock Vina：本机为 `1.2.7`，方案写的是 `1.2.3`。版本差异须在 Methods
  中如实记录，不能把 1.2.7 写成 1.2.3。
- GROMACS：环境审计发现 `E:\BaiduNetdiskDownload\gmx2020.6_GPU\bin\gmx.EXE`
  为 `2020.6-MODIFIED`，方案文字写的是 2022；版本差异须记录。由于尚未执行
  Figure 5d–5g 仍为 `not_run`，不能由占位图或旧结果替代。
- `gmx_MMPBSA`：可执行文件已安装在
  `D:\AAA Liver cancer\envs\experiment_plan_one\Scripts`，版本为 1.7.0；缺少
  真实 GROMACS 轨迹和拓扑，Figure 5h 仍不能标记完成。
- R `MCL`：已安装 `MCL 1.0`，Figure 2b 默认恢复原方案 MCL，并记录
  `inflation`。
- R `CellChat`：已安装 `CellChat 1.6.1`，Figure 4f 默认调用原方案 R 实现。
  当前数据只有 1 个 NCD 文库和 2 个 HFD 文库，细胞级概率不能当成独立
  生物学重复的显著组间结果。
- R `scTenifoldKnk`：已安装 `scTenifoldKnk 1.1` 和 `scTenifoldNet 1.4`，
  Figure 4g/4h 默认调用 R 原版实现；工具缺失或运行失败时会标记
  `blocked/failed`，不再静默回退局部 GRN。
- PLIP：已安装 `PLIP 3.0.1`，Figure 5b 在存在真实 docking pose 时运行 PLIP
  并保存机器可读结果；当前尚未基于冻结版本完成全流程对接重跑。
- PoseBusters：当前 Python 环境未安装，结构合理性独立复核需要使用其他环境。

## 本轮真实小规模验证

- Cohort 冻结：基于旧 `y1` 已处理矩阵构建并冻结了 434 条样本/文库记录。GSE164441 的 20 个文库正确归并为 10 名患者，GSE202379 的 59 个文库归并为 47 名 donor；GSE270583 的 4 个文库保持 `unresolved_library_level_only`，只允许描述性分析。Cohort 和 donor mapping 已达到 `frozen_complete`，但 GeneCards、OMIM、TTD 授权凭证仍缺失，所以总状态为 `partially_frozen_authorization_blocked`。
- 独立同终点队列：GSE126848 已从 GEO 下载原始 counts、SOFT 和 series matrix 到 `D:\AAA Liver cancer\external_validation\GSE126848`，SHA-256 校验通过。条件是 14 个正常体重对照、12 个肥胖对照、15 个 NAFL、16 个 NASH，共 57 名独立患者。当前冻结状态为 `reserved_ready_not_evaluated`，不得在模型冻结前用于训练、筛选或调参。
- 动物重复：公开 GEO 检索未找到满足至少两个独立 NCD 和两个独立 HFD 动物且可用于单细胞复现的队列。`animal_replication_registry.frozen.json` 保持 blocked，GSE270583 仍只能做描述性分析。
- 湿实验规划：已冻结 3 项 `protocol_only` 实验记录，覆盖暴露表型、靶点干预/救援和直接结合/靶点占位；伦理批准、剂量、样本量、真实数据和结果文件均为空，T3 继续 blocked。
- R MCL：在旧 `y1` 的 60 节点、389 边 STRING 网络上实际运行，输出 1 个 MCL 模块；这证明调用链可用，不代表最终 42 Panel 已验收。
- R CellChat：在旧 `y1` 的 31,922 个细胞上实际运行。首次运行暴露小鼠基因名大小写不一致问题；修复后先映射到 1,057 个 CellChatDB.mouse 相互作用基因，再按条件/细胞类型抽样 4,811 个细胞，仅分析 NCD/HFD，JQF 明确排除。导出矩阵由约 1.26 GB 降到约 34 MB。当前只输出细胞级概率，不作独立生物单位组间显著性。
- R scTenifoldKnk：在旧数据抽取的 300 个基因、500 个细胞上实际运行，返回 300 行 `diffRegulation`，并记录随机种子、网络数和细胞数。
- PLIP：在旧 `y1` 的 APP、EGFR、GPAT3、PPP2R2A、STAT3 五个 docking pose 上实际运行，分别得到 21、43、27、37、22 个相互作用。当前环境缺少 Open Babel InChIKey 插件，代码使用进程内兼容写回，并在状态中记录 `inchikey_writer_patched=true`。

## 结论

当前结论不再使用 96.26% 或 98.33% 作为投稿成熟度：这些数值只反映旧版本
代码入口和文件覆盖，不能证明原方案方法已经运行、统计单位正确或结论可复现。
本轮已完成的关键修复是恢复原方案方法、删除静默回退、强化状态和溯源，并补充
真实小规模工具验证；尚未完成一次冻结版本的真实全流程重跑。

因此当前应明确分为：

1. **代码和方法对齐已改善**：MCL、R CellChat、R scTenifoldKnk 和 PLIP 已接入，
   环境审计能够识别缺失工具。
2. **真实结果仍未全部生成**：旧 `y1` 结果属于历史运行，不能与新代码混合；
   必须使用新 `run_id` 重跑 00–10 阶段。
3. **科学前提仍然缺失**：GeneCards/OMIM/TTD 授权导出、未参与调参的同终点
   外部队列、GSE270583 独立动物重复、100 ns GROMACS 轨迹和湿实验因果验证
   不能由脚本补造。

在上述缺口关闭前，实验方案一只能声明为可复现计算流程持续修复中；不能声明
42 个 Panel 已完成或达到投稿级。
