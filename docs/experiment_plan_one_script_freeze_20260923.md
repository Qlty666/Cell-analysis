# 实验方案一：脚本冻结与交付分级说明

日期：2026-09-23

本文件只约束脚本实现、运行状态和证据交付。原《实验方案一》的研究问题、终点、队列用途、Figure 1–5 和 42 个 Panel 不因本文件而改变。

## 冻结范围

- 冻结标识：`experiment-plan-one-script-freeze-20260923`
- 冻结基线：`e6c0370b3b7210ac557ae9c5432e72af74f14b35`
- 实验方案一包版本：`0.2.0`
- 冻结状态：`frozen`
- 当前实现不再新增分析模块、模型族、终点、数据集角色或必须交付 Panel。

## 为什么多轮优化后仍反复出现问题

前几轮调整主要在处理表面完成率和结果展示，但没有同时完成下面五件事，因此问题会不断“换名字”后重新出现：

1. 覆盖度百分比把代码入口、文件存在和方法有效混在一起，导致实现率看起来很高，但原方案方法未必真正运行。
2. 工具缺失或运行失败时曾存在静默回退，图注仍可能沿用原方法名称，结果与方法边界不一致。
3. 同一 Panel 的“代码已接入”“真实结果已生成”“科学结论可成立”没有分别记录，状态会被上层流程过早标记为完成。
4. R、Python、GROMACS、PLIP 等外部工具版本和调用路径没有统一进入环境审计，换机器后结果不可复现。
5. 修复后没有先执行一次冻结版本全流程重跑，旧结果、旧缓存和新代码混在一起，导致同一问题在不同阶段再次被审计发现。

因此本轮只做冻结范围内允许的修复：恢复原方案方法、删除静默回退、修正状态口径、补齐版本和输入输出溯源；没有增加分析模块，也没有改变实验方案一的终点和 42 个 Panel。

## 原始方法恢复记录

本轮已安装并接入以下原方案或原方案指定工具，替代方法不再作为默认路径：

| Panel | 原方案方法 | 当前实际调用 | 失败处置 |
|---|---|---|---|
| 2b | MCL | R `MCL` 1.0，记录 `inflation` | 工具不可用时 `blocked`，不改成 Louvain 后沿用 MCL 名称 |
| 4f | CellChat | R `CellChat` 1.6.1，使用物种数据库和原始计数 | 工具不可用时 Figure 4f `blocked`，不把显式 LR 分数冒充 CellChat |
| 4g/4h | scTenifoldKnk | R `scTenifoldKnk` 1.1，保存网络数、细胞数、种子和 `diffRegulation` | 工具不可用或失败时 `blocked/failed`，不静默回退局部 GRN |
| 5b | PLIP 或经验证几何分析 | `PLIP 3.0.1`，输出机器可读相互作用表 | 没有真实 pose 时 `not_run`；几何候选保持候选层级 |

## 允许的后续修改

只有以下四类修改可以在冻结版本上继续提交：

1. 代码缺陷；
2. 数据与元数据对齐；
3. 运行状态、Panel 状态和错误状态判定；
4. 输入输出可追溯性、校验和、manifest 和缓存失效。

以下事项必须新建方案修订或独立任务，不能作为脚本优化继续推进：

1. 新增分析模块或模型族；
2. 修改实验终点或数据集用途；
3. 用替代算法冒充原方案算法；
4. 为了让 AUC、P 值或覆盖率变好而调参；
5. 在缺少真实轨迹或实验数据时补造结果。

## 交付分级

| 级别 | 名称 | 允许的结论 | 必备证据 |
|---|---|---|---|
| T1 | 可复现计算流程 | 脚本可运行、状态真实、输入输出可追溯 | 冻结策略、run_id、输入输出校验和、测试与小规模真实数据核查 |
| T2 | 完整 42 Panel 证据 | 真实方法、真实输出、完整限制和 Panel 证据 | 冻结版本全流程重跑、实际方法或已登记替代、Panel 证据账本、源数据 |
| T3 | 湿实验因果验证 | 可支持暴露、靶点干预、救援或直接结合主张 | 真实实验数据、伦理信息、实验与计算证据对齐 |

T1 通过不代表 T2 或 T3 通过。T2 未完成时，综合风险、Panel 数量、文件数量或百分比不能替代科学验收。T3 未完成时，计算结论必须限定为候选机制、预测或探索性关联。

T3 记录模板见 `config/experiment_plan_one_experimental_validation.example.json`。只有真实完成并核验的实验记录才可把状态改为 `completed_verified`。

## Panel 状态规则

缺少工具、授权数据、独立重复、轨迹或实验条件时，统一使用：

- `blocked`：外部必要条件缺失；
- `not_run`：方法或运行没有开始；
- `failed`：已尝试但失败；
- `exploratory`：结果可展示但不能作主结论；
- `valid_negative`：方法有效且结果为阴性；
- `valid_positive`：方法有效且结果支持预先定义的结论。

文件存在不等于方法完成。一个 Panel 完成不自动代表同一 Figure 或同一分析模块完成。

## Cohort、donor 与授权数据冻结

- `00_plan/cohort_manifest.tsv`：运行期间更新；计划阶段的占位行只能形成 draft。
- `00_plan/donor_mapping.tsv`：按 library 记录 donor/animal 映射、生物单位、映射证据、状态和冲突。
- `00_plan/cohort_freeze.json`：全流程结束时生成；只有没有 planned 行、患者/动物映射没有 unresolved、授权数据没有错误时才是 `frozen_complete`。
- `00_plan/cohort_manifest.frozen.tsv`：仅在 `frozen_complete` 时生成，并保存 SHA-256。
- `00_plan/source_authorizations.frozen.json`：冻结 GeneCards、OMIM、TTD 等授权来源的导出文件、授权引用和文件哈希。缺少合法导出或授权引用时保持 `blocked_incomplete`。
- `00_plan/external_validation_registry.frozen.json`：冻结独立同终点候选队列、文件哈希、条件计数和 `not_evaluated` 状态。GSE126848 已作为保留验证队列登记，模型冻结前不得用于训练、筛选或调参。
- `00_plan/animal_replication_registry.frozen.json`：记录独立动物重复要求和公开数据检索结果。没有满足 NCD/HFD 独立动物门槛的数据集时保持 blocked，GSE270583 只能描述性展示。
- `00_plan/experimental_validation_manifest.json`：记录暴露、靶点干预/救援和直接结合实验的计划、伦理与结果路径。没有真实 `completed_verified` 记录时 T3 必须保持 blocked。
- `config/experiment_plan_one_source_authorizations.json`：只填写文件路径和授权引用，不把密钥或个人隐私写入仓库。

Cohort/donor 已满足冻结条件但授权凭证缺失时，状态为 `partially_frozen_authorization_blocked`；此时 cohort 可以冻结，但 T2 仍不能通过。

如果 `cohort_freeze.json` 不是 `frozen_complete`，T2 必须保持 `not_run`，不得用旧 cohort manifest、GSM 数量或文件存在状态替代样本冻结。

## 方法替代登记

实际方法与原方案不一致时，必须同时记录：

- 原方法；
- 实际执行方法；
- 替代原因；
- 影响的 Figure/Panel；
- 是否算法等价；
- 是否需要方案修订；
- 方法名称、版本和证据文件。

登记替代方法不等于该方法与原方案等价。当前默认路径已恢复 MCL、CellChat 和 R scTenifoldKnk；若未来重新启用 Louvain、Python scTenifoldpy 或局部 GRN，必须创建新的替代登记并把 Figure 图注和状态同步改名。几何候选接触也不能写成已证实氢键。

## 投稿级门槛

只有满足以下全部条件，才允许重新讨论投稿级结论：

1. 冻结脚本完成一次真实全流程重跑；
2. 42 个 Panel 均有方法、数据、输出、证据和限制记录；
3. 所有 `blocked` 和 `not_run` 项都有明确的补足、降级或删除处置；
4. 所有方法替代已登记，且没有冒充原方法；
5. 论文主张与计算证据、外部验证和湿实验证据逐项匹配。

当前状态由运行产物中的 `10_reports/governance/delivery_readiness.json` 自动计算；该文件不等于任何期刊或审稿人的最终认定。
