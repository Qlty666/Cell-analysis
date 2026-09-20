# 全仓库评估与优化建议

评估日期：2026-09-16。范围：当前工作区的 Liver Cancer Bioinformatics Workflow，README 标示版本 1.7.1。用户未指定单个文件，因此按整个分析系统评估。此次仅评估与生成报告，没有修改业务代码。

## 总体判断

项目已经具备较完整的科研工作流骨架：统一 CLI、Web 入口、表达分析、证据中心、虚拟筛选、独立对接、MD、报告与测试。现阶段最需要投入的是结果正确性、断点恢复一致性与验证质量，其后才是提速和继续扩展功能。

可以作为探索性研究和分析自动化工具使用，但目前不能仅凭“流水线完成”或内部 publication_grade 标签认定结果具备论文级可信度。已发现并复现旧对接结果复用、MD 完成状态误判和下载完整性检查缺口；高级分析代码还存在交叉验证前按全体样本标签筛选基因的问题。

没有进行真实数据全集重跑、真实 Vina/GROMACS 长任务、联网数据库全量可用性验证或生产规模性能测试，因此不提供虚构的准确率、覆盖率或提速百分比。

## 实际检查与测试结果

| 项目 | 结果 | 解释 |
| --- | --- | --- |
| Git 跟踪的 Python 文件语法解析 | 188 个通过 | AST 解析，不代表运行正确 |
| Git 跟踪的 JSON 文件解析 | 11 个通过 | 不代表业务配置校验通过 |
| src/scripts 下 R 脚本语法解析 | 19 个通过 | 不等于 Seurat 等依赖可用 |
| Git 跟踪测试目录静态统计 | 55 个 Python 文件，485 个 test_ 函数 | 不是实际执行数量或覆盖率 |
| 默认 Python 的 pytest 启动 | 失败，缺 pygments | 默认解释器本身还缺 numpy |
| Anaconda 环境完整测试收集 | 2 个模块收集失败 | test_experiment_plan_one.py 与 test_figure_optimization.py 缺 anndata |
| 排除上述两个模块后的测试 | 455 通过、3 失败、4 跳过，27 条警告 | pytest 报告用时 41.91 秒 |
| 3 个失败的直接原因 | 1 个缺 RDKit；2 个 R 驱动加载缺 Seurat | 不能将这些环境失败直接归为算法缺陷，也不能报告全部测试通过 |
| 最小复现 | 7 项观察输出 | 使用模拟工具/响应与临时目录，没有真实对接、MD 或网络请求 |

可执行测试命令：

```powershell
& 'D:\anaconda\python.exe' -m pytest -q -p no:cacheprovider --disable-warnings --durations=15 --ignore=tests/test_experiment_plan_one.py --ignore=tests/test_figure_optimization.py --basetemp validation_output/review_20260916_available
```

最小复现脚本位于 [reproduce.py](<D:/AAA Liver cancer/Script/validation_output/review_20260916/reproduce.py>)，输出位于 [reproduction_results.json](<D:/AAA Liver cancer/Script/validation_output/review_20260916/reproduction_results.json>)。两者属于被 Git 忽略的本地评估产物。

## 优先修复的问题

以下 P1 表示可能改变结果或错误报告成功，宜优先修复；P2 表示工程保证或分析解释存在缺口。

### 1. P1：对接断点恢复不验证输入，--force 也未彻底强制重算

位置：[docking/pipeline.py:53](<D:/AAA Liver cancer/Script/src/docking/pipeline.py:53>)、[docking/docking.py:91](<D:/AAA Liver cancer/Script/src/docking/docking.py:91>)、[integration.py:1731](<D:/AAA Liver cancer/Script/src/pipeline/integration.py:1731>)。

阶段完成标记只存时间戳；配体恢复键只有 id 和 replicate，未纳入受体、配体内容、seed、盒参数、exhaustiveness 或工具版本，也不核验姿势文件是否存在。--force 删除阶段标记，但未关闭内部 docking.resume，因此重新进入 dock 阶段仍可能跳过已有成功行。独立分子对接复用同一 runner；集成流程虽有上层签名，进入靶点子流程时仍受这个问题影响。

复现结果：首次 seed=42 运行后改为 999，模拟对接函数累计只调用一次；随后 pipeline(force=True) 新增调用数仍为 0。模拟成功记录指向不存在的 pose 文件，也照样被恢复逻辑接受。

建议：阶段与每个配体采用统一内容指纹，记录输入 SHA256、规范化配置、工具版本和输出校验；输入变更失效当前阶段及下游。明确将 force 传至配体恢复层和 redock 层。结果用每个运行独立目录或可更新的唯一任务键存储，历史结果保留为单独运行，避免混进本次 CSV。

验收：更换 seed、受体、同名配体内容或盒参数必须触发重算；删除 pose 必须重算；参数未变时只补失败项；force 必须实际调用对接程序。

### 2. P1：MD 只有输入文件，也可能被计为完成

位置：[md_simulation.py:372](<D:/AAA Liver cancer/Script/src/docking/md_simulation.py:372>)、[md_simulation.py:833](<D:/AAA Liver cancer/Script/src/docking/md_simulation.py:833>)、[md_simulation.py:132](<D:/AAA Liver cancer/Script/src/docking/md_simulation.py:132>)。

_run_gromacs_hit 发现 md.tpr 就走“复用已有模拟”。缺少 md_nojump.xtc 时，分析函数返回空指标，不抛失败；上层随后无条件设置 completed。

复现结果：只创建占位 md.tpr，不提供轨迹，得到 completed=1、failed=0。准备步骤在复现中被 mock；完成判断路径使用真实实现，未运行 GROMACS。

建议：明确 prepared、running、interrupted、completed、analysis_failed 等状态；用真实轨迹最后时间、目标生产时长、正常结束信息和有效输出共同判定完成。存在 checkpoint 时支持继续 mdrun，不能把 tpr 当完成证据。输入与 MDP 指纹变化时禁止直接复用旧轨迹。

验收：只有 tpr、零长度轨迹、半程轨迹、损坏轨迹、已完成但分析失败分别得到正确状态；完成率只计真实满足条件的运行。

### 3. P1：高级 ML 分析存在交叉验证前的标签泄漏

位置：[advanced_analysis.py:1330](<D:/AAA Liver cancer/Script/src/analysis/advanced_analysis.py:1330>)、[advanced_analysis.py:685](<D:/AAA Liver cancer/Script/src/analysis/advanced_analysis.py:685>)、[advanced_analysis.py:813](<D:/AAA Liver cancer/Script/src/analysis/advanced_analysis.py:813>)。

程序先用主队列全体样本做差异分析，按 padj/top_genes 构造 candidate_genes，再对筛后的矩阵交叉验证。后面的 Pipeline 虽将 SelectKBest 放进训练折，仍无法消除前面全队列 DEG 筛选已带入的验证标签信息。此外，选取 CV 指标最好的模型后再展示其同队列 OOF 指标，不能替代独立的外层模型选择评估。

这是代码路径确认的问题；本次没有测量其对具体队列 AUC 的抬高幅度，也不否定独立外部验证本身的价值。scikit-learn 官方说明要求包括特征选择在内的数据拟合只发生在训练部分，并说明嵌套交叉验证的作用：[数据泄漏](https://scikit-learn.org/stable/common_pitfalls.html)、[嵌套验证](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html)。

建议：将 DEG 筛选纳入训练折，或使用分析前固定、来源独立的候选基因集合；外层评估泛化性能，内层选择基因数和模型。对同患者、多部位或多队列样本增加 patient/cohort 分组切分。批次校正及归一化也需检查是否借用了验证集信息。可参考现有 experiment_plan_one/ml.py 的嵌套评估实现，但不能只复制 CV 外壳。

验收：每折候选基因可追踪；验证标签改变不影响该折训练特征；随机标签多次验证不持续高于机会水平；独立队列保持完全隔离。

### 4. P1：maxwarn=0 实际被转换为 20

位置：[md_simulation.py:670](<D:/AAA Liver cancer/Script/src/docking/md_simulation.py:670>)、[md_simulation.py:756](<D:/AAA Liver cancer/Script/src/docking/md_simulation.py:756>)。

两处命令构建采用 int(value or 20)，数值 0 会被当成缺省值。当前配置默认也是 20，因此用户即使显式设置 0，也无法按预期禁止 grompp 警告放行。

复现配置表达式得到 configured=0、command_value=20。GROMACS 官方提醒应理解警告后才考虑覆盖，不能普遍视为无害：[gmx grompp 文档](https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html)。

建议：仅在 None/缺键时使用缺省值；默认采用 0；允许的特殊警告需记录原因与完整日志。

验收：数值 0 在两处命令中均保留 0；非法负数或非整数在运行前报错。

### 5. P1：下载器未确认完整长度和正确续传范围

位置：[common/http.py:137](<D:/AAA Liver cancer/Script/src/common/http.py:137>)、[common/http.py:153](<D:/AAA Liver cancer/Script/src/common/http.py:153>)。

代码计算 expected 后只用于最大大小限制，没有在流结束后验证实际长度。续传只检查 HTTP 206，没有校验 Content-Range 起点及资源身份。并非所有短响应都会被底层库自动拒绝，不能依赖这个行为作为完整性保证。

模拟响应复现：声明 100 字节、实际只有 3 字节仍返回成功；已有 abc、请求从字节 3 开始续传，但服务器报告起点 0、返回 def，函数仍接受并拼成 abcdef。

建议：使用 .part 文件，校验响应长度和 Content-Range，保存 ETag/Last-Modified 并配合 If-Range；有官方校验和时验证校验和，完成后原子重命名。对 416 区分“已完整”与“本地状态无效”，限制可重试错误类型。

验收：截断、错误范围、远端文件替换、忽略 Range、正常续传均有模拟响应测试；不完整下载不得进入完成缓存。

### 6. P2：机制性扰动质量门控只检查文件存在

位置：[integrated_report.py:202](<D:/AAA Liver cancer/Script/src/pipeline/integrated_report.py:202>)。

mechanistic_perturbation 直接取 insilico_summary.json.exists()，未解析状态、数据内容或本次运行归属。复现中零字节文件使该项门控为 True。其他门控仍然需要通过，所以此结果不意味着一个空文件单独就能把整个报告升级为 publication_grade。

建议：解析结果 schema，检查成功状态、有效结果数量、输入与运行指纹；将计算扰动证据与湿实验机制证据分别展示。质量标签应表达“哪些检查通过”，避免被误读为发表认证。

验收：空文件、损坏 JSON、失败状态及其他运行遗留文件都不能使该项通过。

### 7. P2：单细胞统计路径需要更明确的分析层级

位置：[analysis_pipeline.R:1863](<D:/AAA Liver cancer/Script/src/analysis/analysis_pipeline.R:1863>)、[analysis_pipeline.R:1984](<D:/AAA Liver cancer/Script/src/analysis/analysis_pipeline.R:1984>)。

当前主 pseudobulk 按 sample+condition 聚合全部细胞，得到的是混合细胞群的样本级差异；它不是细胞类型内的差异状态检验。细胞组成变化可能影响该结果。生物学重复不足或部分计算失败时会退回细胞级 Wilcoxon，虽写 warning，但应确保下游不会把这种探索性结果与样本级检验等同处理。

建议：保留总体差异，同时增加 cell_type×sample pseudobulk；根据设计支持 batch、配对 patient 等协变量。生物学重复不足时明确标记“仅探索性”，将这一状态传至目标筛选和最终报告。muscat 的官方实现展示了按 cluster_id 与 sample_id 聚合并检验的设计：[muscat 官方代码与说明](https://code.bioconductor.org/browse/muscat/)。

这属于方法适用性改进，不等于现有总体差异计算在所有场景都错误。当前 differential.py 已按样本比例做 Welch 检验，README 仍描述细胞计数卡方，应同步文档，不能据旧文档认定代码仍在使用卡方。

### 8. P2：依赖安装失败可以被 CI 的非锁定回退掩盖

位置：[tests.yml:31](<D:/AAA Liver cancer/Script/.github/workflows/tests.yml:31>)、[requirements.txt](<D:/AAA Liver cancer/Script/requirements.txt>)、[launchers/_common.bat](<D:/AAA Liver cancer/Script/launchers/_common.bat>)。

CI 在 requirements.txt 安装失败后改装无版本约束的最小依赖，对接依赖允许失败，R 语法检查也是可选。这样单元测试成功并不证明声明的完整环境能安装或主要科研工具可运行。本机还存在默认 Python 与具备科学计算包的 Anaconda 环境不一致。

建议：核心测试与完整环境安装验证分开；后者必须严格按锁定依赖安装。建立 core、expression、docking、md 的清晰环境入口，启动器支持显式解释器并显示版本。R 设置固定包库与锁定文件；不在正常分析运行中隐式升级依赖。增加必需的 R 小数据运行检查和真实工具小规模 smoke 作业。

验收：干净机器按文档完成一个小规模全流程；每个作业记录解释器与工具版本；安装失败不能被另一个未锁定环境的成功掩盖。

## 性能、结构与运行管理建议

这些建议依据代码结构提出，尚无真实规模基准，不承诺具体提速。

| 方向 | 当前证据 | 可行改进 | 验收指标 |
| --- | --- | --- | --- |
| ML 重复训练 | ml_analysis.py:270 起分别调用一次 cross_val_score 和两次 cross_val_predict | 单次折循环保存标签、概率与分数，再统一计算指标；保留一次全量最终拟合 | 每折拟合次数、总耗时、输出一致性 |
| ML 尺度处理 | 基础 ml_analysis.py 将细胞比例与原始 QC 均值一起输入，SVM/MLP 路径未放缩放器 | 折内 StandardScaler；训练数据质量与每类最少样本检查 | 同尺度模型比较、极端量纲数据验证 |
| 对接资源预算 | 默认 4 workers×4 Vina CPU；Web 各领域各自排队 | 全局 CPU/内存/GPU 预算，按任务实际资源申请；大库限制在途 futures 数量 | 峰值线程与内存、吞吐、取消响应时间 |
| H5 转换 | h5_converter.py:140 附近对稠密 Dataset 整体载入；MTX 非零元逐项 Python 写出 | 按块读写，减少 CSR/COO 和浮点验证副本；评估批量/编译实现 | 百万细胞或大矩阵峰值 RSS、文件一致性 |
| 长任务恢复 | web_state.py 的进程注册表和队列在内存，JSON 主要保存历史 | SQLite 任务状态表、运行目录锁、进程身份/心跳；重启时核对孤儿任务 | Web 重启后能识别运行/中断任务，避免重复启动 |
| 下载与证据请求 | 已有共享重试和 SQLite 证据存储 | 在修好完整性后，增加按来源限流、Retry-After、TTL/版本缓存；批量 API 优先 | 请求数、缓存命中率、429 率、失败可解释性 |
| 模块体积 | integration.py 4227 行，web_ui.py 4187 行，analysis_pipeline.R 3556 行 | 先定义阶段输入/输出与状态协议，再逐个提取领域服务、配置校验、路由；保留 CLI 兼容层 | 单次变更影响范围、契约测试、原有命令兼容 |
| 数据与报告契约 | 一些路径仅验证文件存在或非空 | 统一 schema、状态、单位、输入指纹、实际方法、降级原因；报告消费同一结构 | 缺失不等于 0，失败不等于完成，前后端状态一致 |

表达流程自身也有仅凭 pipeline_complete.json 跳过重算的路径（orchestrator.py:511）。重构恢复协议时应一起覆盖 accession、species、环境变量 QC 参数、原始数据哈希与代码版本，不能只修上层 integration.py。

## 应保留的现有设计

- CLI、Web、分析实现和报告已分层，独立对接能够复用通用 runner，有渐进重构的基础。
- 集成流程已有阶段签名、输出检查及可复现 manifest；证据中心已有来源记录、覆盖率、来源消融和外部验证设计。应修正薄弱环节并统一协议。
- H5 转换对 counts 层优先级与整数计数做了检查，能够拒绝部分归一化矩阵误当原始计数的情形。
- Web 默认绑定回环地址，非回环启用令牌；具备来源检查、请求体限制、路径约束。压缩包提取也已有路径越界与链接检查。本次未做完整渗透测试，不将这些检查等同于全面安全认证。
- 已有大量功能测试和模拟工具测试。后续重点应补充故障注入、状态与统计设计的测试，不能单纯追求测试数量。

## 建议实施顺序

| 阶段 | 交付内容 | 估算投入 | 完成标准 |
| --- | --- | --- | --- |
| 第一批 | 修复对接 force/恢复、MD 完成判定、maxwarn、下载完整性、机制门控 | 约 3–5 个开发日 | 本报告最小复现全部变成能阻止旧行为的回归测试 |
| 第二批 | 高级 ML 折内筛选与嵌套验证、按患者/队列分组、单细胞统计层级和报告降级 | 约 5–10 个开发日 | 小型已知真值数据、置换标签、独立队列三类验证 |
| 第三批 | 环境锁定、严格 CI、小型真实工具 smoke、完整运行说明 | 约 2–4 个开发日 | 干净环境可重现，失败可定位 |
| 第四批 | 全局资源预算、H5/ML 性能基准、任务持久化与渐进拆分 | 约 5–10 个开发日 | 在固定数据与硬件上比较耗时、峰值内存和恢复行为 |

投入为单人开发粗估，不含真实大数据计算、数据库配额等待和 100 ns MD 等计算成本。先修正确性可以避免把旧结果、半成品或乐观验证指标通过提速进一步放大。

此次工作位于 codex/comprehensive-review-20260916 分支；没有合并、推送、修改版本号或发布 Release。
