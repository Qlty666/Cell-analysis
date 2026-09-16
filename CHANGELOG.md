# 更新日志


### v1.7.1

- 补齐实验方案一使用的 `beautifulsoup4`、`anndata`、`scanpy` 运行依赖，并同步 GitHub Actions 最小依赖回退列表，修复 Linux/Windows CI 测试收集阶段的模块缺失问题。
- 全量测试通过：427 个测试用例 + 80 个 subtests；GitHub Actions 的 Ubuntu 与 Windows 任务均通过。

### v1.7.0

- 新增 `src/evidence/` 多数据库证据中心：统一证据模型、SQLite 证据库、来源运行记录、实体和关系证据长表、覆盖感知评分与 leave-one-source-out 消融分析。
- 新增 Open Targets、ChEMBL、BindingDB、PubChem BioAssay、GWAS Catalog、GTEx、Human Protein Atlas 和 DepMap/本地快照连接器；CTD、Tox21/ToxCast、LINCS、DisGeNET 及授权数据库支持本地表接入。
- 实验方案一新增 `evidence` 阶段，位于 `disease` 与 `ppi` 之间，输出 `02b_evidence/target_priority.csv`、证据矩阵、来源消融和完整 SQLite 证据库。
- 机器学习改为外层折内完成特征集和模型选择，并使用训练折内 sigmoid 校准后的概率进行外部验证；新增嵌套模型选择记录。
- 修复样本级表达数据使用随机 PCA/UMAP 占位图的问题：改为真实 PCA，样本数不足时明确跳过 UMAP，单细胞 PCA/UMAP 失败时不再伪造降维结果。
- 对接重复种子不足时自动生成互不相同的确定性随机种子，避免把同一 seed 的重复运行误报为独立重复。
- 默认 GROMACS 生产步数调整为 100 ns，并增加证据阶段签名、阶段输出校验和虚拟敲除输入签名，配置或输入变化时不再静默复用旧结果。
- 全量测试通过：426 个测试用例 + 80 个 subtests。
- 将 v1.6.0 新增的高级分析、MR/共定位、分析工作区导出、多队列验证、网络多来源靶点、重复对接、阳性对照、PCA/FEL、MM/PBSA 和 KEGG GSEA 状态同步到网页端。
- 新增网页版“高级分析”页，支持任务日志、状态轮询、历史记录和结果文件下载；全自动流水线页可配置 Advanced/MR 并衔接高级分析优先级表。
- 真实数据验证页新增随机全流程、多队列靶点、真实 PDB 证据和随机证据/对接盒四类验证入口。
- 拆分网页任务层：高级分析与 MR/导出任务、真实数据验证、结果文件发现和通用参数解析分别进入 `web_analysis.py`、`web_validation.py`、`web_files.py` 和 `web_utils.py`。
- 拆分集成流水线公共能力：共享异常、QC 门控、样本级差异丰度与关键基因排序分别进入 `errors.py`、`qc.py`、`differential.py` 和 `key_targets.py`；`integration.py` 与 `web_ui.py` 保留兼容导出和原有运行行为。

### v1.6.0

- 新增 Advanced 多队列分析与本机 MR/共定位流水线，覆盖批量表达验证、WGCNA、免疫浸润、生存分析、因果分析和结果导出。
- 新增本地 Codex 分析工作区导出功能，支持数据集编号、项目根目录登记、增量同步和结果清单刷新。
- 修复 KEGG GSEA 与当前 `clusterProfiler` 的参数兼容问题，并新增 KEGG 状态文件、失败门控和回归测试。
- 加固真实数据验证脚本：PDB 下载自动重试，验证输出状态可检测，单数据集 smoke test 与完整 20 队列验证使用独立门槛。
- 修复 CADD、QC、CLI、环境检查、网页结果同步和文件导出中的多项运行缺陷，并补充 R 模块快照和语法测试。
- 强化仓库隐私与发布安全：清理本机绝对路径示例，扩充凭据和临时文件忽略规则，并同步重写公开 Git 历史中的作者信息。
- 全量测试通过：387 个测试用例 + 50 个 subtests。

### v1.5.1

- 依据只读代码审查修复 P0/P1 问题：MD 默认模拟时长与完成信息、可复现随机种子、HETATM 丢弃告警、信号检测零细胞修正、虚拟敲除输入校验、配置跨字段约束、`docking.cpu` / tanimoto 去重阈值接线、HTTP 超时重试、HTML 转义与失败原因记录。
- 收敛 docking 与 molecular_docking 共享逻辑：参数化保存配置、流水线阶段/日志名与 CLI 公共参数，删除重复的盒子检测、路径解析、读取/转义/JSON 写入等私有副本。
- 拆分巨型文件：`analysis_pipeline.R` 的 R 分析函数抽到 `src/analysis/R/*.R`；`web_ui.py` 抽成 `web/web_handler.py`、`web/web_data.py`、`web/web_state.py`、`web/web_results.py`；`integration.py` / `generate_report.py` 的阶段路径、HTML 报告与指南表抽到独立模块。
- 补齐 R 模块快照与冒烟测试：orchestrator 运行时会同步 `src/analysis/R` 并设置 `LIVER_R_MODULES_DIR`；新增 `tests/test_r_pipeline_syntax.py`，逐文件解析 R 并验证模块可独立加载。
- 清理静默异常与仓库卫生：网页/脚本路径读取失败改为带原因日志；离线环境检查、公共 `.bat` Python 探测、`.gitignore`、`pytest.ini` / `tests/conftest.py` 与 GitHub Actions Windows/Linux 测试工作流同步。
- 全量测试通过：339 个测试用例 + 45 个 subtests。

### v1.5.0

- 把现有 CADD 功能全部接入全自动集成流水线：新增 `07_cadd_downstream`（GROMACS MD 准备/运行、对接 ML 重打分、Amber/GROMACS 与 UniDock-Pro/HDOCK/HADDOCK 工具导出）、`08_network`（网络毒理学）和 `09_faers`（FAERS 不相称性信号）阶段；细胞反馈顺延为 `10_cell_feedback`，集成报告顺延为 `11_report`。
- 每个成功对接靶点在独立工作目录自动完成 MD/导出；MD 无 GROMACS 时默认只准备输入，不会阻断流水线。提供 `--md-mode auto`、`--md-top-n` 后可自动运行 GROMACS 模拟。
- 对接 ML 重打分支持全流程训练/预测：`--docking-ml-training-csv`、`--docking-ml-label-column`、`--docking-ml-model`；已存在 `ml_model_info.json` 时自动直接重打分。
- 网络毒理学与 FAERS 作为输入可选的流水线阶段接入：提供 `network_toxicology.compound_targets_csv` / `target_sources` 和 `faers.input_csv` 后自动运行，缺输入时写 `skipped` 汇总而不是中断或伪造结果。
- 网页全自动流水线页同步增加 MD/ML/导出、网络毒理学与 FAERS 的输入控件、运行开关、阶段显示和结果表；结果清单文件扫描加入 MD 导出目录。
- 集成报告与 `integration_summary.json` 新增 CADD 下游、网络毒理学和 FAERS 汇总；阶段标记升级后旧运行目录会自然重建新增阶段。
- 审计修复：网络毒理学/FAERS 用户显式提供输入但运行失败时不再吞错，流水线返回失败；MD 逐靶点记录 `md_requested/md_failed`，全部失败时阶段状态为 failed 并中断。
- 审计修复：运行前自动清理旧版 `07_cell_feedback` / `08_report` 等不再属于当前阶段表的遗留标记，避免网页进度虚高或阶段显示错误；网络毒理学 Venn 默认改为 true，与独立命令和文档一致。
- 审计修复：未配置 ML 训练 CSV 时不再把空路径当作训练文件尝试读取；修复后会如实跳过 ML 重打分。
- 全量测试通过：271 个测试用例 + 18 个 subtests（含新增 4 个审计修复回归测试）。

### v1.4.0

- `virtual-knockout` 默认检测原始计数矩阵并调用官方 scTenifoldpy/scTenifoldKnk 引擎，同时保留 CellOracle 风格的 GRN 传播；新增 scTenifold QC/网络/流形配置与 DrugReflector checkpoint 化合物排序。
- 对接结果增加 strong / moderate / weak 亲和力分级并写入 CSV、图片与报告；GROMACS MD 扩展 Rg、SASA、蛋白-配体氢键、结合口袋残基 RMSF 与 last-half 稳定性标签，网页端与结果清单同步展示。
- 单细胞 QC 新增 UMI-基因数 log-log 关系图和统计表，marker 图新增 RidgePlot 与堆叠小提琴视图；报告、结果清单和结果图指南同步更新。
- GO 富集气泡图按 BP / CC / MF 分面展示，并优化 colorbar 与图例布局。
- 报告主函数改为显式接收输出目录，避免 pytest 参数被误当成报告目录，并补充回归测试。
- 网页版把原“虚拟筛选”页拆成独立功能页面：虚拟筛选、分子动力学、虚拟敲除、网络毒理学、FAERS 和真实数据验证；导航全站统一，CADD 页面脚本抽到 `web/static/dock_app.js`，不再在单个页面堆叠全部模块。
- 各功能板块使用优化：参数分组折叠、设置保存/恢复/重置与启动自动保存；虚拟敲除支持建模基因数/细胞数/传播轮数/DrugReflector 等高级参数；网络毒理学可指定疾病基因列并控制 Venn；真实数据验证可设置数据集数与种子并实时查看运行状态。
- 补充 `docs/project_structure.md` 代码结构说明和全站网页模板、测试与文档同步。
- 全量测试通过：252 个测试用例 + 10 个 subtests。

### v1.3.0

- 新增血红蛋白 QC 跟踪与可配置污染过滤：表达分析 QC 统计 `percent.hb`，对明显血红蛋白高占比细胞默认做 99% 分位数上限过滤，并可用 `LIVER_QC_MAX_HB` / `LIVER_QC_MAX_RIBO` 覆盖阈值；血红蛋白模式同步纳入小鼠 `Hbq1b`。
- 新增 `liverbio` 统一软件入口：`scripts/liverbio.py` 把表达分析、虚拟筛选、全自动流水线、数据集搜索、网页端和环境检查集中到同一命令，支持 `help`、`version`、`expression`、`docking`、`full`、`datasets`、`web` 和 `doctor` 子命令。
- 新增数据集搜索、表达分析、全自动流水线、虚拟筛选四个功能 Codex skill 源文件，配套 `scripts/install_codex_skills.py` 与 `docs/software_guide.md`；skill 只指导 Codex 调用项目现有脚本，不复制核心分析代码。
- 新增 GROMACS 分子动力学模拟模块：支持准备 GROMACS 输入、运行 MD、容错缺失分析文件并自动填充生产时间、按蛋白骨架对齐计算配体 RMSD；模块已接入虚拟筛选 CLI 与网页端。
- 新增独立分子对接板块：`src/molecular_docking/`、`scripts/run_molecular_docking.py`、`config/molecular_docking_config.json` 与网页“分子对接”页面，提供独立工作目录、任务日志、暂停/继续、自动检测对接盒、结果表、图库、HTML 报告和 PDBQT 构象下载。
- 网页版整体体验优化：统一页头与快捷入口、表单分组折叠、设置保存/恢复/重置、结果统计卡片、任务数量统计、结果清单实时筛选，并新增导航高亮与进行中任务数量徽标。
- 新增 `tests/test_md_simulation.py`、`tests/test_molecular_docking.py`、`tests/test_liverbio_cli.py`，并同步补充网页端模板与运行接口测试。
- 新增 `docs/project_structure.md`，按入口层、实现层、网页层、配置技能层、测试文档层整理全部代码文件；同步更新 README 的脚本文件一览与目录结构。
- 全量测试通过：236 个测试用例 + 10 个 subtests。

### v1.2.0

- `virtual-knockout` 新增单细胞调控网络虚拟敲除扩展：传入 `--insilico-gene` 后，基于 CellOracle 思路的 KNN 平滑、稀疏 GRN 与迭代信号传播模拟基因敲除，并把结果整理成 UMAP 命运偏转矢量场、TF-靶基因调控网络、WT/KO 表达变化、Top 15 定量表、GO/KEGG 富集气泡图和中文 HTML 报告。
- 输出目录统一为 `04_knockout/in_silico/`，数据与图编号为 `fig_63` 至 `fig_68`，支持 `--insilico-embedding-csv`、`--insilico-regulators-csv`、`--insilico-species` 和 `--insilico-photo-dir`。
- 新增 `insilico_enrichment.R`，复用项目现有 `org.Hs.eg.db` / `org.Mm.eg.db` / `clusterProfiler` 环境完成 GO（BP/CC/MF）与 KEGG 富集。
- 新增 `tests/test_insilico_knockout.py` 单元测试与 `requirements.txt` 中的 `umap-learn` 依赖。
- 真实数据验证后修复密集 UMAP 箭头重叠问题：`fig_66_ko_shift_umap.png` 在细胞数较多时自动改为网格聚合箭头，避免箭头成片遮挡散点；新增 `src/docking/export_single_cell_insilico.R` 从 Seurat 对象导出虚拟敲除输入。
- 图表细节优化：`fig_63_ko_target_expression_table.png` 增加中文标题与中文表头；`fig_65_ko_regulatory_network.png` 补全中心基因标签、连接线变细并放大节点；GO/KEGG 富集气泡图增加独立颜色条区域，避免右侧标签压入绘图区或被裁剪。

### v1.1.0

- 富集分析网络图拓展：`fig_22_go_network.png` / `fig_23_kegg_network.png` 在保留 Top5 核心通路的基础上，按 `p.adjust` 追加最佳 5 个延伸通路（`showCategory=10`），保持 `cnetplot` 通路-基因网络样式；报告、网页图开关、结果图指南与结果清单同步更新，并使用真实数据完成验证。
- 网页版加固：任务运行/排队/暂停期间不再因页面空闲自动退出；POST 接口增加同源与跨站校验、请求体大小限制；任务历史与队列写入加锁，避免并发重复记录；非回环地址启动时输出安全警告。
- 网页版结果清单搜索优化：支持按结果图名或完整本地路径直接定位对应清单，搜索后只显示命中的结果清单，不再展开指南章节或关联文件。
- 修复 GEO manifest 中 `single_cell_hint` 被 bulk 文件名规则覆盖的问题。
- 全量测试通过：208 个测试用例 + 4 个 subtests。

### v1.0.1

- 真实数据集全流程兼容性修复：GEO 下载器对同时包含 bulk 与单细胞文件的系列改为优先选择 bulk 计数表，避免把样本级矩阵与单细胞矩阵混载导致下游失败；新增 `.h5ad.gz` 自动解压转换。
- 修复 10x 多文库数据合并时重复 barcode 报 `duplicate row.names` 的问题：每个文库按文件名生成唯一样本后缀，跨文库细胞名不再冲突。
- 修复 GEO series matrix 中带引号逗号标题解析错误，并兼容不同 platform 的 series matrix 列不一致（`rbind.fill` 合并）。
- 样本级（bulk/microarray）数据集差异表达改为直接对原始样本计数运行 DESeq2，不再先聚合成伪 bulk，解决小样本量下 DEG 全为空的问题；保留单细胞伪 bulk 逻辑。
- 条件自动推断增强：支持按 `geo_accession` 匹配样本元数据，并把多类别条件（如 ICC/HCC、Normal/Tumor）自动归并为两组，避免无法推断分组导致流水线中断。
- 大细胞数单细胞数据集 `FindAllMarkers` 增加 `max.cells.per.ident` 限制，避免聚类 marker 计算长时间停滞。
- 全流程证据阶段增加空 `key_genes.csv` 保护，DEG 为空时跳过证据收集并生成空表，不再因缺列崩溃。
- 新增 `scripts/run_web_full_new_datasets.py`：通过网页端 `/full/start` 自动提交并监控多个真实数据集的全流程任务。
- 项目规则新增第 13 条：禁止“准备做”式中间回复，任务执行完成或遇到必须由用户处理的阻塞前不发送文字回复。

### v1.0.0

- 正式发布 1.0.0：表达分析、虚拟筛选、全自动集成流水线、网页端、数据集搜索与验证脚本形成完整可交付链路。
- 修复验证脚本“全部失败仍报成功”的问题：`validate_real_evidence.py` / `validate_real_random.py` 新增最少成功靶点、最少配体记录和对接盒成功数阈值，未达标时返回非零退出码；技能脚本缺失时返回失败而不是崩溃。
- 修复长时间验证脚本未捕获超时的问题：`validate_pipeline.py` / `validate_real_pipeline.py` 超时后明确报错并继续/退出，不再让单个数据集拖垮整个验证。
- 真实数据验证兼容 ML 阶段跳过：`validate_real_pipeline.py` 与合成验证一致，ML 被跳过时不再强制要求 ML 图。
- 数据集检索验证不再把同义词启发式标签当作真实相关性：新增 `--manual-labels` 人工标签、`label_source` 字段和 `manual_rounds` / `manual_found_rate` 汇总；提供人工标签后 `found` / `found_rate` 以人工标签为准。
- ML/DL 相关性排序修复：训练要求同时包含相关与不相关样本，重排校验二分类模型；`--eval` 默认要求人工标签，`--allow-heuristic-eval` 仅用于冒烟。
- 修复 `validate_random_real_full_pipeline.py --only` 未校验 GEO 编号导致的路径逃逸风险，新增 `GSE\d+` 白名单与单元测试。
- 修复 docking 环境检查/安装假成功：AutoDockTools 源码与 zip 均缺失时安装返回失败；`check-dock-env` 与 `check_dock_environment.py` 按真实检查结果返回退出码。
- 收敛 Rscript 探测到 `src/common/env.py`，统一 Windows PATH、Program Files 与用户目录探测，消除多处重复实现。
- 修复真实 TCGA 生存分析 p 值固定为 1.0 的问题：改用 Cox 似然比检验输出真实 p 值；GSE165816 文件名解析增加回退规则。
- 更新数据集搜索 `run_supported` 表述：搜索结果只标记“可自动运行候选”，实际文件可用性在下载时校验。
- 全量测试通过：201 个测试用例 + 4 个 subtests。

### v0.9.0

- 数据集搜索从单一 NCBI GEO 扩展为多数据库检索：新增 EBI ArrayExpress/BioStudies 与 Expression Atlas，支持按数据库选择范围，并给每条结果增加来源数据库、质量分和 `run_supported` 标记。
- 新增 `src/data/biostudies_downloader.py`：自动选择 BioStudies/Expression Atlas 的 processed 表达文件，生成与 GEO 相同的 manifest，可直接进入表达分析与全自动流水线；`E-GEOD-xxxxx` 自动映射到 GEO `GSExxxxx`。
- 流水线入口、网页端和 CLI 支持 `GSE125449`、`E-MTAB-1234`、`S-BSST123` 等数据集编号；搜索页和全自动流水线页文案同步更新。
- 新增/更新数据集搜索、BioStudies 下载、网页参数和 accession 校验测试。
- 修复 ArrayExpress 等样本级数据集缺少 condition 元数据时无法自动推断分组的问题：R 流水线会按样本名 token 自动推断条件组，非单细胞数据可继续跑差异表达。
- 修复富集为空时 stage 09 汇总读取 `go_up` 占位表导致 `undefined columns selected` 的问题，样本级数据集可正常生成 `summary.json` 并完成流水线。
- 细胞反馈阶段新增反馈靶基因差异表达：把虚拟敲除/虚拟筛选得到的反馈基因放回 Seurat 对象，输出 `feedback_deg.csv`、`fig_59_feedback_targets_volcano.png` 和 `fig_60_feedback_condition_violin.png`。
- 细胞反馈阶段新增 GO/KEGG 富集分析：对反馈靶基因运行 `enrichGO` / `enrichKEGG`，输出 `feedback_enrichment_go.csv` / `feedback_enrichment_kegg.csv`，富集 Top5 使用与 `fig_22_go_network.png` 相同的 `cnetplot` 通路-基因网络图（`fig_61_feedback_go_network.png` / `fig_62_feedback_kegg_network.png`），不再使用气泡图。
- KEGG 富集在线注释超时放宽到 180 秒，并用 `setReadable` 把 KEGG 网络和结果表中的基因 ID 转换为基因符号，便于直接查看 Top5 通路关联基因。
- 集成报告、网页全自动流水线结果区和结果图指南同步展示反馈差异表达与 GO/KEGG 富集结果；独立 `cell-feedback` 命令新增 `--feedback-species hs/mm`。

### v0.8.2

- 修复 stage 08 发表分析阶段给 Seurat 对象写入 `sample_label` 时携带样本名而非细胞条形码名，导致 `AddMetaData` 报 `No cell overlap between new meta data and Seurat object` 的问题：写入前移除向量 names，多样本数据集可正常完成发表分析。

### v0.8.1

- 修复 R 表达流水线在读取 `dataset_mode.txt` 时调用尚未定义的 `ckpt_path()` 导致启动即报错的问题：将 `ckpt_path()` 定义提前到读取数据集模式之前，单细胞与样本级数据均可正常启动流水线。

### v0.8.0

- 取消“全自动流水线仅支持单细胞”的限制，全面支持 bulk RNA-seq、microarray 和样本级表达矩阵：GEO 下载器不再拒绝非单细胞数据集，这些数据可直接进入表达分析与全自动流水线。
- 修复逐样本单细胞 CSV（如 GSE165816 的 `GSMxxxx_counts.csv.gz`）被误判为 bulk 的问题：下载器会结合 Series Matrix 标题和矩阵表头中的细胞 barcode 自动识别为单细胞。
- R 表达流水线新增 `dataset_mode` 模式：单细胞数据保留 QC、双细胞、聚类、注释和细胞级反馈；bulk RNA-seq、microarray 等非单细胞数据按样本级表达矩阵运行，自动跳过双细胞检测、细胞级发表分析和细胞反馈，保留差异表达、富集、虚拟敲除和虚拟筛选链路。
- 全自动流水线对非单细胞样本级数据集自动跳过细胞类型差异丰度检验和细胞反馈阶段，并在 `differential_abundance_summary.json` / `cell_feedback_summary.json` 中注明原因。
- 网页数据集搜索与全自动流水线页移除 bulk 选择限制，`single-cell`、`bulk`、`other` 均可带入并运行；网页文案统一改为“表达分析”。
- 更新 GEO 数据集搜索、manifest 分类、流水线编排和全流程相关测试。

### v0.7.2

- 彻底修复 h5ad/loom 转换：正确解码 h5ad 中 bytes 类型的 obs/var 索引，避免 barcodes 为空或基因名变成 `b'...'`。
- 支持通过 obs/var 的 `_index` 属性定位真实细胞/基因列；索引确实缺失时自动生成 `Cell1...CellN` / `Gene1...GeneN`，不会生成空 barcodes 文件。
- h5ad 优先读取 `layers/counts` 原始计数层，X 为归一化数据时不再丢失真实计数；矩阵值非数值时自动置零。
- 多文件 h5ad/loom 转换改为并行执行，大样本 GEO 数据集的转换耗时显著下降。
- R 流水线增加防御性回退：barcodes/genes 文件为空或损坏时自动补位，不再因单个空文件导致整个 stage 01 反复失败。

### v0.7.1

- 修复 GEO 下载器未识别 `.h5ad`/`.loom` 单细胞补充文件的问题，h5ad 数据集（如 GSE315928）可自动转换为 10x MTX 后进入流水线。
- 下载器增加已有文件跳过下载、已解压归档跳过重复解压逻辑，断点续跑不会重复下载大文件。
- 修复 h5ad/loom 转换结果丢失 `_extracted/` 子目录的问题，manifest 现在引用实际生成的矩阵文件。

### v0.7.0

- 扩大证据数据库覆盖范围：虚拟筛选证据收集新增 STRING 互作网络、Open Targets 靶点关联、Reactome 通路、PharmGKB 药物基因组、AlphaFold 结构和 KEGG 通路来源。
- 证据报告升级为数据库覆盖表，逐项展示每个数据库的状态与记录数；`check-cadd` 自动检查新增 skill 脚本。
- 全自动流水线 stage 03 同步接入扩展数据库：`gene_evidence.csv` 新增 STRING/Reactome/PharmGKB/AlphaFold/Open Targets/KEGG 的计数、ID 和 `database_sources`，集成报告与结果清单同步展示。
- 更新 `VIRTUAL_SCREENING_REQUIREMENTS.md` 环境清单和网页端证据收集文案。

### v0.6.1

- 网页端整合近期全部功能：单细胞分析页新增 ML 模型选择，支持 `xgb` / `rf` / `gbm` / `mlp` / `lasso_svm`。
- 全自动流水线页新增单细胞 ML 模型、DepMap 依赖表和 PPI 网络边表输入，后端同步追加 `--ml-model`、`--depmap-csv`、`--ppi-network-csv`。
- 虚拟筛选页新增重打分 ML 模型/训练 CSV/标签列，以及虚拟敲除 PPI 网络边表输入。
- 结果清单和结果图指南补充 ML 校准曲线与 `lasso_svm` 选定特征表。
- `run_full_pipeline.py` 新增 `--ml-model`，使全自动流水线可传递单细胞 ML 模型选择。
- 补充网页控件、命令行参数和结果清单测试。

### v0.6.0

- 数据集搜索新增更多筛选条件：数据类型、最小/最大样本数、起始/结束日期、平台 GPL/名称、数据集类型。
- 命令行 `scripts/search_datasets.py` 同步支持 `--data-type`、`--min-samples`、`--max-samples`、`--start-date`、`--end-date`、`--platform`、`--dataset-type`。
- 网页数据集搜索页的“过滤与下载选项”和搜索结果过滤栏同步新增对应条件，可在已加载结果中即时筛选。
- 搜索结果接口返回 `filters` 字段，便于网页端展示当前筛选状态。
- 补充筛选函数、日期解析、网页参数传递和模板字段测试。

### v0.5.1

- 网页版虚拟筛选页新增“网络毒理学分析”和“FAERS 不相称性信号检测”卡片，支持填写工作目录/输入文件后直接运行，并展示交集表、PPI hub 评分、Venn 图、C-T-P-D 网络和 FAERS 信号表。
- 新增 `/dock/network`、`/dock/faers` 同步运行接口及对应文件下载接口。
- 结果清单页新增“网络毒理学与 FAERS 信号”信息板块；结果图指南同步补充输出路径、用途和判读标准。
- 全自动流水线结果文件扫描与下载权限加入 `network_toxicology`、`faers` 输出目录。
- 补充网页运行接口与文件下载安全测试。

### v0.5.0

- 新增网络毒理学命令 `python scripts\run_docking.py network`：支持化合物-疾病靶点交集、多数据库来源计数、STRING PPI hub 评分、Venn 图和 C-T-P-D 节点/边导出。
- 新增 FAERS 风格信号检测命令 `python scripts\run_docking.py faers`：支持 ROR、PRR、BCPNN IC、EBGM 四种不相称性指标及组合信号判定。
- 虚拟敲除新增 PPI hub 维度：`virtual-knockout --ppi-network-csv` 可把 STRING 边表的 degree/betweenness/clustering 合入 `target_score`；全自动流水线同步支持 `--ppi-network-csv` 和 `full_pipeline_config.json` 中的 `ppi_network_csv`。
- 单细胞 ML 扩展：新增 `--ml-model` / `LIVER_ML_MODEL`，支持 `xgb`、`rf`、`gbm`、`mlp`、`lasso_svm`；`lasso_svm` 使用 LASSO 初筛 + SVM-RFE 特征选择，并输出选定特征表。
- 单细胞 ML 新增校准曲线图 `fig_45_ml_calibration_curve.png`，与 ROC/PR、SHAP 一起构成更完整的模型诊断。
- 同步更新 `config/docking_config.json`、`config/full_pipeline_config.json`，并新增网络毒理学与 FAERS 单元测试。

### v0.4.2

- 网页版同步 bulk 数据集检测：数据集搜索页新增数据类型过滤与说明，bulk 数据集只可下载、不能直接运行全自动流水线。
- 全自动流水线页内嵌搜索新增数据类型列，bulk 数据集不可选择；带入链接携带 `data_type` 参数，启动前校验可拦截 bulk 数据集。
- 单细胞分析页标注仅支持单细胞数据集。
- 补充网页端 bulk 链接与搜索数据类型的单元测试。
- 单细胞富集分析优化 `fig_22_go_network.png` / `fig_23_kegg_network.png`：网络图先按 `p.adjust <= 0.05` 筛选，再取前 5 个通路，保留原图并降低网络杂乱度。
- 新增 `fig_46_go_top5.png` / `fig_47_kegg_top5.png`：展示上调基因 GO BP / KEGG 筛选后 Top5 富集结果。
- 新增 `fig_48_qc_pvalue_comparison.png` / `fig_48_qc_pvalue_comparison.csv`：比较过滤前后 QC 指标在条件间的差异程度（Wilcoxon 秩和检验 + BH FDR）。
- 同步报告、网页图开关、结果清单和流水线输出校验。

### v0.4.1

- 修复 bulk RNA-seq 数据集（如 GSE299321）下载时报 `No count matrix files found` 的误导性错误：下载器现可识别 RAW tar 中的逐样本计数表并标记为 bulk 数据集。
- 单细胞与全自动流水线新增 bulk 数据集检测：遇到 bulk 数据集时明确提示当前流水线仅支持单细胞数据，原始计数文件和 manifest 仍保留在缓存中。
- GEO 数据集搜索结果新增 `data_type` 字段（`single-cell` / `bulk` / `other`），网页端显示数据类型徽标，bulk 数据集不再提供“全自动流水线”入口。
- 补充 bulk 检测、manifest 拒绝和搜索分类单元测试。

### v0.4.0

- 网页版整体布局优化：新增共享样式 `web/static/app.css`，统一卡片、按钮、表单、表格、状态卡片和移动端响应式细节；各页面统一页头、副标题与快捷入口。
- 全自动流水线页表单按“基础设置 / 分析参数 / 反馈与起始阶段 / 运行开关”分组折叠，支持保存设置、恢复设置、恢复默认；结果区新增统计卡片，并保留 QC 门控与细胞类型差异丰度表。
- 单细胞分析页改为分组折叠布局，支持表单保存/恢复/重置，结果图开关可折叠。
- 数据集搜索页将搜索条件与过滤/下载选项分组，页面信息层级更清晰。
- 任务进度页新增运行中、排队中、已暂停、进行中合计统计卡片。
- 结果清单页新增按文件名、内容或用途实时筛选的工具条。
- 修复 `/static/app.css` 的 MIME 类型，浏览器可正确应用共享样式。
- 补充网页模板布局与静态样式测试；全量单元/集成测试 99 个全部通过。
- 全自动流水线标记升级为“配置/输入指纹”标记：每个阶段记录签名，参数、物种、标签、配置或输入表变化时自动使当前及下游阶段失效，修复“改参数后仍跳过旧结果”的静默错误。
- 新增阶段输出校验：每个阶段按必需输出清单验证，缺失或空输出不会写入完成标记，重跑时自动重建。
- 新增 QC 门控：单细胞阶段后汇总 `qc_metrics.json`，按 `qc_gate` 阈值输出 `pass/warn/fail`，可配置最小细胞数、最小基因数、最大双细胞率和伪 bulk 强制要求。
- 新增细胞类型组成差异检验：基于单细胞注释表计算条件间细胞比例变化（2×2 卡方 + BH FDR），输出 `differential_abundance.csv`，与差异表达配对防止漏掉组成偏移。
- 虚拟筛选加固：对接盒中心/尺寸校验、非法盒自动跳过并记录原因；PDB 下载失败自动重试 3 次。
- 全自动流水线新增 `--dry-run`、`--skip-qc-gate`、`--skip-differential-abundance`。
- 集成报告新增 QC 门控和差异丰度结果表，汇总信息同步写入 `integration_summary.json`。
- 补充阶段签名失效、QC 门控、差异丰度、对接盒校验、PDB 重试和 dry-run 单元测试。
- 新增细胞反馈闭环：全自动流水线新增 `07_cell_feedback` 阶段，把虚拟敲除评分和虚拟筛选命中重新写回 Seurat 单细胞对象，计算每细胞靶基因表达、筛选靶点模块评分、细胞类型表达汇总与富集检验，并输出 UMAP/DotPlot/箱线图/热图和 `feedback_targets.csv`。
- 新增独立命令 `python scripts\run_docking.py cell-feedback --single-cell-root <单细胞结果目录>`，可对已有虚拟敲除/虚拟筛选结果单独执行细胞反馈分析。
- 全自动流水线支持 `--feedback-top-n`、`--feedback-max-features`、`--feedback-timeout` 和 `--skip-cell-feedback`。
- 网页版全自动流水线新增细胞反馈参数、阶段显示和结果表；结果清单页新增细胞反馈输出说明。
- 加固单细胞 R 流水线进程管理：运行前快照 `analysis_pipeline.R`，暂停或停滞时终止完整 R 进程树，避免残留子进程。
- 优化停滞判断：日志长时间无更新但 R 进程仍在计算时延长等待，不再误杀活跃任务。
- 补充流水线编排器单元测试，覆盖 R 脚本快照、CPU 采样和进程树终止。
- `AGENTS.md` 明确任务分支隔离：不同任务必须创建不同功能分支，提交时只包含当前任务相关文件。

### v0.3.0

- 网页版统一界面新增 GEO 数据集搜索页：按疾病、研究方向或原始查询搜索，支持 ML/DL 相关性重排序、CSV/JSON 结果下载与批量下载。
- 单细胞报告升级为总报告：对 `results/figures` 和 `results/data` 下每个结果文件生成独立分析，包括文件说明、表格规模、关键字段统计、P 值/差异方向/样本分组等结论，并汇总到 `result_report.html`；DOCX/PDF 导出同步加入结果文件清单与数量统计。
- 单细胞分析新增 `fig_09_deg_horizontal_violin.png`：按校正 P 值排序的差异最显著基因横向小提琴图，并在图中标注 P 值；新增 `LIVER_DE_VIOLIN_TOP_N` 和 `LIVER_DE_VIOLIN_MAX_CELLS` 环境变量。
- 网页版单细胞分析完成后的结果报告入口：可直接打开包含逐文件分析的 `result_report.html`。
- UMAP 聚类图直接标注细胞类型名称。
- 结果清单页优化为仅展示结果图和结果数据，并按阶段目录整理全流程输出。
- 网页任务进度页支持历史记录、一键清空与完成/中断弹窗提醒。
- 统一结果文件名与阶段输出目录，增强路径处理和网页端安全校验。
- 自动 GEO 数据集搜索脚本 `scripts/search_datasets.py` 与 ML/DL 重排序脚本 `scripts/dataset_search_ml.py` 上线。
