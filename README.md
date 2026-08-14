# Liver Cancer Bioinformatics Workflow

> 当前版本：0.4.2

面向肝癌研究的本地生信自动化工作流，整合三条可实际运行的流水线：

- 单细胞转录组分析：GEO 数据下载、QC、双细胞检测、聚类注释、差异表达、富集分析和 ML 可解释性分析。
- 虚拟筛选（CADD）：靶点证据收集、受体/配体准备、AutoDock Vina 并行对接、命中排序、精细重对接、ML/DL 重打分和 MD/外部工具交接。
- 全自动集成流水线：从单细胞分析直接筛选关键基因/蛋白，再自动完成证据富集、虚拟敲除和虚拟筛选，最终输出集成报告和湿实验验证方案。

## 1. 项目解决什么问题

单细胞分析和虚拟筛选通常依赖多个分散工具和手工流程，结果难以复现，也难以从“细胞图谱”推进到“可干预靶点”。本项目把以下环节串成一条可检查、可断点续跑、可交付的本地链路：

- GEO 单细胞数据从下载、QC、注释到差异表达与富集分析的完整分析链。
- 从显著差异基因中自动筛选关键基因/蛋白，并补充 UniProt、PDB、ChEMBL、BindingDB、PubChem、ChEBI 证据。
- 对有 PDB 结构的靶点自动准备受体、收集已知配体并运行 AutoDock Vina 对接、命中排序和精细重对接。
- 通过虚拟敲除对候选靶点做多维评分：表达差异、增殖共表达、共表达网络 hub、DepMap 依赖、疾病逆转、通路控制、细胞类型特异性、预后、成药性和脱靶风险。
- 将排序后的候选靶点导出为分阶段湿实验验证方案（细胞系、类器官、药物剂量反应、动物模型、PDX）。
- 每次运行写入 `run_manifest.json`，记录配置哈希、输入文件 SHA256、软件版本和本次参数，保证结果可复现、可追溯。

## 2. 主要功能

### 2.1 单细胞分析流水线

- 输入 GSE 编号，自动下载 GEO 数据并识别常见格式（10x MTX、CSV、Series Matrix）。
- 内置 GSE125449 适配；支持 `--species hs/mm/auto`。
- QC、双细胞检测（`scDblFinder`）、PCA/UMAP 聚类和细胞注释。
- 差异表达（DESeq2 pseudobulk 或 Seurat Wilcoxon 回退）与 GO/KEGG/GSEA 富集分析。
- 按校正后 P 值排序输出差异最显著基因的横向小提琴图，并在图中标注 P 值。
- 报告生成器会对每个结果图和结果数据文件生成逐文件分析，并汇总为一份总报告。
- 可选 CellChat 细胞通讯分析（设置 `LIVER_RUN_CELLCHAT=yes`）。
- ML 疾病分类（XGBoost 或 RandomForest）、特征重要性和 SHAP 解释。
- 发表级分析：细胞周期打分与回归校正、聚类 marker 发现、功能签名打分、推断 CNV、SingleR 自动注释，以及可选 slingshot 拟时序。
- 输出 48 张分析图、HTML 报告，支持断点续跑、暂停/继续、停滞自动重启。
- GO/KEGG 通路网络图按校正 P 值筛选后展示 Top5，并新增对应的 Top5 富集气泡图。

新增分析可通过环境变量控制：

| 环境变量 | 默认 | 说明 |
| --- | --- | --- |
| `LIVER_RUN_CELLCYCLE` | `yes` | 细胞周期打分与 UMAP/比例图 |
| `LIVER_REGRESS_CELLCYCLE` | `no` | 在 ScaleData 时回归 S/G2M 得分 |
| `LIVER_RUN_CLUSTER_MARKERS` | `yes` | `FindAllMarkers` 与聚类 marker 热图/DotPlot |
| `LIVER_RUN_SIGNATURES` | `yes` | 增殖/EMT/缺氧/免疫检查点等签名打分 |
| `LIVER_RUN_CNV` | `yes` | 基于染色体窗口均值的推断 CNV 热图 |
| `LIVER_RUN_SINGLER` | `yes` | SingleR 参考注释与混淆矩阵 |
| `LIVER_RUN_TRAJECTORY` | `no` | slingshot 拟时序轨迹（需安装 `slingshot`） |
| `LIVER_DE_VIOLIN_TOP_N` | `12` | Top DEG 横向小提琴图展示的基因数 |
| `LIVER_DE_VIOLIN_MAX_CELLS` | `1000` | 每个条件下用于该图的抽样细胞数 |

### 2.2 虚拟筛选（CADD）流水线

| 命令 | 说明 |
| --- | --- |
| `init` | 创建 docking 工作目录骨架和配置文件 |
| `evidence` | 调用数据库 skill 收集靶点和已知配体证据 |
| `prepare-receptor` | 通过 Meeko/Open Babel/MGLTools 准备受体 PDBQT |
| `prepare-ligands` | RDKit 标准化、3D 构象生成并输出 PDBQT |
| `dock` | AutoDock Vina 并行对接，按配体写入结果，支持断点续跑 |
| `analyze` | 按亲和力排序、阈值筛选、Tanimoto 多样性选择 |
| `redock` | 对 Top 命中用更高 exhaustiveness 精细重对接 |
| `ml-train` / `ml-predict` | 随机森林、GBDT、MLP 或 PyTorch MLP 重打分 |
| `export-md` / `export-external` | 导出 Amber/GROMACS 和 UniDock-Pro/HDOCK/HADDOCK 模板 |
| `report` | 生成 HTML 汇总报告 |
| `virtual-knockout` | 基因敲除优先级和多维靶点评分 |
| `export-validation` | 把排序后的靶点导出为湿实验验证方案 |
| `cell-feedback` | 把虚拟敲除/对接结果返回 Seurat 单细胞对象做细胞级分析 |
| `detect-box` | 从共晶配体自动检测对接盒并写回配置 |
| `check-env` / `check-cadd` | 检查软件、库和数据库 skill 环境 |
| `pipeline` | 按阶段执行准备、对接、分析、重对接、报告，支持断点续跑 |

### 2.3 全自动集成流水线

`scripts/run_full_pipeline.py` 把单细胞分析、关键基因筛选、证据富集、虚拟敲除和虚拟筛选串成一条流水线：

```text
01 single_cell           GEO 单细胞分析（下载、QC、注释、差异表达、富集）
02 key_targets           从显著 DEG 中筛选并排序关键基因/蛋白
03 evidence              UniProt / PDB / ChEMBL 证据富集（带本地缓存）
04 knockout_inputs       导出样本级伪 bulk 表达矩阵并生成敲除输入
05 knockout              虚拟敲除 + 多维靶点评分 + 湿实验验证方案
06 docking               对有 PDB 结构的靶点自动收集已知配体并跑 Vina 对接
07 cell_feedback         把虚拟敲除/对接结果返回 Seurat 做细胞级反馈分析
08 report                生成集成 HTML 报告和 run_manifest.json
```

每一阶段写标记文件，重跑时自动断点续跑；`--start-stage` 可从任意阶段开始。标记文件不再只是时间戳：每个阶段会记录配置和输入指纹（`signature`），当 `top_genes`、物种、标签、证据/对接配置、关键基因表或证据表发生变化时，会自动使当前阶段及下游阶段失效，避免“参数改了但结果仍是旧值”的静默错误。每个阶段完成后还会按 `STAGE_OUTPUTS` 校验必需输出，缺失或空文件不会写入完成标记。

新增 `--dry-run`，不执行任何阶段，只打印每个阶段会 `RUN` 还是 `DONE` 及原因；新增 `--skip-qc-gate` 和 `--skip-differential-abundance` 可分别关闭 QC 门控和细胞组成差异检验。

单细胞阶段完成后会汇总 `qc_metrics.json`，包括细胞数、基因数、双细胞率、伪 bulk 使用情况和下游阶段统计，并按 `config/full_pipeline_config.json` 中的 `qc_gate` 阈值给出 `pass/warn/fail` 门控结果；默认阈值不强制拦截，需要拦截时在配置中填写阈值即可。差异表达阶段同时补做细胞类型组成差异检验（2×2 卡方 + Benjamini-Hochberg FDR），输出 `differential_abundance.csv` 并写入集成报告，避免“表达没变但比例变了”的组成偏移被漏掉。

虚拟筛选阶段增加对接盒有效性校验：中心/尺寸非有限值或尺寸非正数时跳过该靶点并写明原因；PDB 下载失败会自动重试 3 次，避免单次网络抖动直接丢弃有结构靶点。

细胞反馈阶段会把虚拟敲除评分和虚拟筛选命中合并成反馈清单，重新读取单细胞 Seurat 对象，为每个候选基因写入细胞级表达、计算筛选靶点模块评分，并输出细胞类型表达汇总、模块富集检验、条件×细胞类型汇总和 UMAP/DotPlot/热图等结果；同时生成 `feedback_targets.csv`，把筛选优先级与细胞表达特异性合并为 `cell_support_score`，用于下一轮靶点收敛。

### 2.4 虚拟敲除与多维靶点评分

核心评分 `knockout_score` 由表达差异、增殖共表达、共表达网络 hub 程度和 DepMap CRISPR 依赖加权得到。提供附加数据后还会计算：

- `reversal_score`：表达变化方向与疾病模块的一致性。
- `pathway_score`：与细胞周期、凋亡、EMT、p53、PI3K-AKT 等通路基因集的平均相关性。
- `specificity_score`：基于注释细胞类型的表达特异性，降低“误伤友军”风险。
- `prognosis_score`：方向感知的风险比评分。
- `druggability_score`：已知配体、蛋白结构和生物活性计数。
- `safety_concern` / `off_target_paralogs`：脱靶和旁系同源标记。

以上维度按 `target_weights` 加权合并为 `target_score`，并把候选基因分为 `core_driver`、`microenvironment_regulator`、`biomarker`、`high_priority`、`low_priority`。

### 2.5 网页版统一界面

`web/web_ui.py` 提供本地网页端，顶部导航默认顺序为全自动流水线、单细胞分析、数据集搜索、虚拟筛选、结果清单，“任务进度”固定在右上角：

- 全自动流水线：`/full`
- 单细胞分析：`/`
- 数据集搜索：`/datasets`
- 虚拟筛选：`/dock`
- 结果清单：`/results`
- 任务进度：`/tasks`

网页端支持任务启动、实时日志、暂停/继续、结果表和文件下载、环境检查与自动补全。任务完成或中断时会弹窗提醒，中断提示会显示运行到的阶段和原因；“任务进度”页集中显示流水线页面启动的排队、运行和暂停任务，保存已完成任务的历史记录并支持一键清空，同时可直接跳转到对应页面继续查看。切换顶部导航不会清除正在运行的任务记录，返回原页面后会自动恢复日志轮询；刷新或关闭网页时才会清除会话中的任务记录。单细胞页的 GEO 数据集编号需手动输入，不再预填示例；全自动流水线页支持直接填写工作目录加载已有结果，工作目录需填到 `outputs` 的上一层目录，未找到结果时会明确提示。数据集搜索页支持疾病、研究方向或原始查询搜索 GEO 数据集，可选 ML/DL 模型重排序、CSV/JSON 结果下载与批量下载；搜索结果可直接带入全自动流水线页，全自动流水线页也可内嵌搜索并选择数据集自动填入 GSE 编号，启动后自动执行 GEO 数据下载；单细胞分析完成后可在结果报告区直接打开包含逐文件分析的 `result_report.html`。
全自动流水线页新增细胞反馈阶段的基因数、展示基因数和跳过选项，并在流程结果中显示 `feedback_targets.csv` 的细胞支持度排序表。
全自动流水线页同步新增 QC 门控与差异丰度检验开关、`dry-run` 仅预演选项；流程结果区新增 QC 门控表和细胞类型差异丰度表，结果清单页补充 `qc_metrics.json` 与 `differential_abundance.csv` 的说明。
网页版整体布局优化：各页面统一页头与快捷入口、表单按“基础/分析/运行”分组折叠、全自动流水线与单细胞表单支持设置保存/恢复/重置、结果区增加统计卡片、任务页增加数量统计、结果清单页支持按文件名/用途筛选。

### 2.6 真实数据验证与可复现性

- `scripts/validate_new_features.py` 使用 20 个 TCGA PanCancer Atlas 队列和 GSE165816 真实单细胞数据运行虚拟敲除与验证方案导出。
- 其余 `scripts/validate_*.py` 分别验证合成数据流水线、对接流水线、真实 GEO 数据、证据收集和随机真实数据。
- 每次评分/导出写入 `run_manifest.json`，记录配置、输入哈希、软件版本和参数。

## 3. 安装方法

### 环境要求

- Python 3.10+（推荐 3.11）。
- R 4.5+（仅单细胞分析和伪 bulk 导出需要）。
- AutoDock Vina（虚拟筛选需要，可放在 `dock/tools/vina.exe` 或加入 PATH）。
- 可选 Codex skills（仅证据收集需要）：`uniprot-skill`、`rcsb-pdb-skill`、`chembl-skill`、`bindingdb-skill`、`pubchem-pug-skill`、`chebi-skill`。

### 安装步骤

复制整个项目文件夹后，在项目根目录执行：

```text
launchers\check_pipeline_environment.bat
launchers\install_pipeline_dependencies.bat
```

检查并安装虚拟筛选 Python 依赖（RDKit、Meeko、Open Babel、AutoDockTools 等）：

```text
launchers\check_dock_environment.bat
launchers\install_dock_dependencies.bat
```

也可以用 conda 直接创建虚拟筛选环境：

```text
conda env create -f environment_dock.yml
```

网页端“软件环境”区域也可以填写安装地址后点击“自动补全环境”。完整软件/插件清单见 `VIRTUAL_SCREENING_REQUIREMENTS.md`。

## 4. 使用方法

### 4.1 单细胞分析命令行

```bash
python scripts\run_pipeline.py GSE125449 --output ../liver_cancer --species auto
```

常用参数：

```bash
python scripts\run_pipeline.py GSE125449 --output ../liver_cancer --species hs --force
python scripts\run_pipeline.py GSE125449 --output ../liver_cancer --skip-download
python scripts\run_pipeline.py GSE125449 --output ../liver_cancer --skip-deps
```

Windows 下也可以直接使用：

```text
launchers\run_GSE125449.bat
launchers\run_pipeline_prompt.bat
```

### 4.2 虚拟筛选命令行

先初始化工作目录：

```bash
python scripts\run_docking.py init
```

运行完整对接流程（准备受体、准备配体、对接、分析、精细重对接、HTML 报告，支持断点续跑）：

```bash
python scripts\run_docking.py pipeline --config config/docking_config.json
```

分阶段运行：

```bash
python scripts\run_docking.py evidence --uniprot P00533 --pdb 1M17 --target-name EGFR
python scripts\run_docking.py prepare-receptor
python scripts\run_docking.py prepare-ligands
python scripts\run_docking.py dock
python scripts\run_docking.py analyze
python scripts\run_docking.py redock
python scripts\run_docking.py report
```

ML/DL 重打分：

```bash
python scripts\run_docking.py ml-train --training-csv data/ml/training.csv --model rf
python scripts\run_docking.py ml-predict
```

导出 MD/外部交接模板：

```bash
python scripts\run_docking.py export-md
python scripts\run_docking.py export-external
```

虚拟敲除（基础评分）：

```bash
python scripts\run_docking.py virtual-knockout \
  --expression-csv data/knockout/expression.csv \
  --metadata-csv data/knockout/metadata.csv \
  --depmap-csv data/knockout/depmap_gene_effect.csv \
  --case-label Tumor --normal-label Normal
```

多维评分与验证方案导出：

```bash
python scripts\run_docking.py virtual-knockout \
  --expression-csv data/knockout/expression.csv \
  --metadata-csv data/knockout/metadata.csv \
  --prognosis-csv data/knockout/prognosis.csv \
  --druggability-csv data/knockout/druggability.csv \
  --off-target-csv data/knockout/off_target.csv \
  --cell-type-column cell_type \
  --case-label Tumor --normal-label Normal

python scripts\run_docking.py export-validation --validation-top-n 10
```

把已有虚拟敲除/虚拟筛选结果返回单细胞分析：

```bash
python scripts\run_docking.py cell-feedback \
  --workdir y3 \
  --single-cell-root ../liver_cancer \
  --feedback-top-n 12 \
  --feedback-max-features 8
```

环境检查：

```bash
python scripts\run_docking.py check-env
python scripts\run_docking.py check-cadd
```

### 4.3 全自动集成流水线

一键启动：

```text
launchers\run_full_pipeline.bat
```

或直接运行：

```bash
python scripts\run_full_pipeline.py \
  --accession GSE125449 \
  --output ../liver_cancer \
  --workdir y3 \
  --top-genes 50 \
  --docking-targets 3
```

常用参数：

- `--skip-scrna`：复用已完成的单细胞结果，直接从关键基因筛选开始。
- `--skip-docking`：只跑虚拟敲除和验证方案，跳过对接。
- `--skip-evidence-fetch`：不联网，使用已有证据缓存或置零。
- `--skip-download` / `--skip-deps` / `--skip-pseudobulk` / `--skip-knockout` / `--skip-cell-feedback`。
- `--top-genes`：关键基因数量，默认 50。
- `--docking-targets`：参与对接的靶点数量，默认 3。
- `--feedback-top-n`：进入细胞反馈的基因数，默认 12。
- `--feedback-max-features`：细胞反馈图中展示的基因数，默认 8。
- `--feedback-timeout`：细胞反馈 R 分析超时秒数，默认 3600。
- `--ligand-library`：自定义配体库（`.smi` / `.sdf` / `.csv`），也可放到 `dock/data/ligands/`。
- `--case-label` / `--normal-label`：虚拟敲除的病例/正常分组标签。
- `--start-stage 08`：从指定阶段继续，之前阶段自动标记为跳过。
- `--dry-run`：不执行任何阶段，只打印每个阶段会运行还是跳过及原因。
- `--skip-qc-gate` / `--skip-differential-abundance`：分别关闭 QC 门控和细胞组成差异检验。

查看阶段清单：

```bash
python scripts\run_full_pipeline.py --list-stages
```

### 4.4 网页版

```text
launchers\run_web_ui.bat
```

浏览器默认打开 `http://127.0.0.1:8000/full`（全自动流水线页）。直接打开指定页面：

```text
launchers\run_web_ui.bat --page dock
launchers\run_web_ui.bat --page full
launchers\run_web_ui.bat --page results
launchers\run_web_ui.bat --page tasks
```

全自动流水线页的“单细胞结果目录”和“工作目录”均为必填项，须手动填写；工作目录需填到包含 `outputs` 的上一层目录，目录中没有结果时页面会显示错误提示。单细胞页需要手动填写 GEO 数据集编号和结果保存地址。“结果清单”页展示 `scripts/run_full_pipeline.py` 成功且完整运行后应输出的图片、数据、报告、断点和溯源文件清单；“任务进度”页支持查看任务进度并跳转到对应任务页面。

关闭所有网页标签后，本地网页服务会在数秒内自动退出并释放端口；正常退出时启动窗口也会自动关闭。再次启动时，如果检测到旧网页服务仍占用端口，会自动关闭旧实例后再启动；若端口被其他非网页程序占用，窗口会保留错误信息等待确认后关闭。

### 4.5 验证脚本与测试

```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

各验证脚本用途：

| 脚本 | 用途 |
| --- | --- |
| `scripts/validate_pipeline.py` | 用 10 个合成单细胞数据集跑通完整单细胞流水线 |
| `scripts/validate_dock_pipeline.py` | 用假 Vina 可执行文件验证对接流水线 |
| `scripts/validate_real_pipeline.py` | 用 10 个真实肝病 GEO 数据集跑通完整流水线 |
| `scripts/validate_real_evidence.py` | 用 10 个真实 PDB 结构验证证据收集 |
| `scripts/validate_real_random.py` | 随机真实数据验证证据收集和对接盒检测 |
| `scripts/validate_new_features.py` | 用 20 个 TCGA PanCancer 队列 + GSE165816 验证靶点评分 |
| `scripts/validate_random_real_full_pipeline.py` | 随机真实 GSE 数据集跑通全自动流水线 |
| `scripts/validate_dataset_search.py` | 用 50 轮随机疾病+研究方向组合验证 GEO 数据集搜索召回 |

### 4.6 自动搜索数据集

`scripts/search_datasets.py` 通过 NCBI E-utilities 自动搜索 GEO 数据集，并可按需下载命中的 GSE 系列：

```bash
python scripts\search_datasets.py \
  --query "hepatocellular carcinoma single cell" \
  --max-results 20 \
  --organism "Homo sapiens"
```

也可以直接指定疾病和研究方向，脚本会自动组合成搜索词：

```bash
python scripts\search_datasets.py \
  --disease "liver cancer" \
  --research-direction "single cell RNA-seq" \
  --max-results 20 \
  --organism "Homo sapiens"
```

搜索结束后会在 `data_cache/dataset_search/` 写出 `dataset_search_results.csv` 和 `dataset_search_results.json`，每行包含 `data_type` 字段（`single-cell` / `bulk` / `other`），网页端也会显示数据类型徽标。需要下载时追加下载参数：

```bash
python scripts\search_datasets.py \
  --query "liver cancer scRNA-seq" \
  --keyword "HCC" \
  --download GSE125449 \
  --download-root ../liver_cancer
```

也可以用 `--download-top N` 下载搜索结果的前 N 个数据集；下载状态写入 `data_cache/dataset_search/download_results.json`。

网页端搜索结果每行提供“全自动流水线”入口，点击后自动带入 GSE 编号；全自动流水线页也可直接搜索并选择数据集，填入后启动即可复用流水线的自动下载流程。注意：全自动流水线只支持单细胞数据，`bulk` 类型数据集会被识别并明确提示，不会再用“No count matrix files found”这类误导性错误中断。

GEO 下载器会把 RAW tar 中的逐样本 bulk 计数表（如 `GSMxxxx_GCxxxx.txt.gz`）识别为 bulk 数据集：数据文件与 manifest 仍会缓存到 `data_cache/<GSE>` 供手工批量分析，但单细胞与全自动流水线会拒绝继续运行，提示改用单细胞 GSE 编号。

网页版已同步该限制：数据集搜索页提供数据类型过滤（`single-cell` / `bulk` / `other`）并显示数据类型徽标，bulk 数据集只保留下载按钮；全自动流水线页的内嵌搜索、直接带入和启动前校验都会拦截 bulk 数据集并给出明确提示；单细胞分析页也标注仅支持单细胞数据集。

批量随机验证搜索是否命中“疾病+研究方向”数据集：

```bash
python scripts\validate_dataset_search.py --rounds 50 --seed 20260812
```

验证结果写入 `data_cache/dataset_search/validation_50_rounds.csv` 和 `validation_50_rounds.json`；未命中时会自动用疾病名或研究方向名扩大搜索范围。

当前 50 轮随机验证（`--seed 20260812`）命中率 100%（50/50），其中 4 轮通过扩大搜索范围命中。

### 4.7 数据集搜索 ML/DL 相关性排序

`scripts/dataset_search_ml.py` 用 TF-IDF 特征训练机器学习/深度学习模型，对搜索结果按“疾病 + 研究方向”相关性重新排序：

```bash
python scripts\validate_dataset_search.py --rounds 50 --seed 20260812
python scripts\dataset_search_ml.py \
  --train \
  --samples data_cache/dataset_search/training_samples.csv \
  --output data_cache/dataset_search/relevance_model.joblib \
  --model-type mlp
```

训练完成后，普通搜索和批量验证都可以传入模型：

```bash
python scripts\search_datasets.py \
  --disease "liver cancer" \
  --research-direction "single cell RNA-seq" \
  --model data_cache/dataset_search/relevance_model.joblib

python scripts\validate_dataset_search.py \
  --rounds 50 \
  --seed 20260812 \
  --model data_cache/dataset_search/relevance_model.joblib \
  --rerank-top 5
```

支持 `lr`、`rf`、`gbm`、`mlp` 四种模型；其中 `mlp` 为多层感知机。可用 `--eval` 对标注样本做交叉验证，例如当前 285 条样本上 MLP 的 ROC AUC 为 0.81。

## 5. 输入输出示例

### 5.1 单细胞分析

输入：

- GSE 编号，例如 `GSE125449`。
- 输出目录，例如 `../liver_cancer`。
- 物种：`hs` / `mm` / `auto`。

输出（以 `../liver_cancer` 为例）：

- `results/figures/`：48 张结果图，按 `01_qc`、`02_doublets`、`03_cluster`、`04_annotation`、`05_deg`、`06_enrichment`、`07_ml`、`08_publication`、`09_cellchat` 阶段分子目录。
- `results/figures/06_enrichment/fig_46_go_top5.png` 与 `fig_47_kegg_top5.png`：上调基因 GO/KEGG 经 `p.adjust <= 0.05` 筛选后的前 5 条通路富集气泡图。
- `results/data/`：QC、双细胞、注释、差异表达、富集、ML 分类和可选 CellChat 表格，按与 `figures/` 相同的阶段分子目录存放；数据文件名与对应结果图编号一致，例如 `data/01_qc/fig_01_qc_metrics.csv`、`data/02_doublets/fig_02_doublet_results.csv`、`data/05_deg/fig_09_deg_significant.csv`、`data/05_deg/fig_09_deg_horizontal_violin.csv`、`data/07_ml/fig_24_ml_feature_importance.csv`、`data/08_publication/fig_36_cnv_heatmap.csv`、`data/08_publication/fig_37_singleR_annotations.csv`、`data/08_publication/fig_39_trajectory_pseudotime.csv`、`data/09_cellchat/fig_40_cellchat_communication.csv`、`data/07_ml/fig_43_44_45_ml_classification_report.csv`。
- `results/data/05_deg/fig_09_deg_significant.csv`：显著差异基因表。
- `results/figures/05_deg/fig_09_deg_horizontal_violin.png`：按校正 P 值排序的差异最显著基因横向小提琴图。
- `results/data/05_deg/fig_09_deg_horizontal_violin.csv`：该图对应的基因与 P 值数据。
- `results/checkpoints/`：Seurat 断点对象。
- `results/result_report.html`：最终 HTML 报告。

### 5.2 虚拟筛选

输入：

- 受体文件：`dock/data/receptors/receptor.pdb`。
- 配体库：`dock/data/ligands/library.sdf`（也支持 `.smi` 和 CSV）。
- 对接盒中心/尺寸：`config/docking_config.json`，或使用 `detect-box` 自动检测。
- 可选训练标签：`dock/data/ml/training.csv`（包含 `smiles` 和 `active` 或 `affinity`）。
- 可选敲除数据：`expression.csv`、`metadata.csv`、`depmap_gene_effect.csv`、`prognosis.csv`、`druggability.csv`、`off_target.csv`。

输出（工作目录默认 `dock`）：

- `data/receptors/receptor.pdbqt`：受体。
- `data/ligands/prepared/`：配体准备结果和 `manifest.csv`。
- `outputs/run_001/docked/results.csv`：对接结果。
- `outputs/run_001/results/`：按阶段汇总的结果目录，`01_analysis`、`02_redock`、`03_ml`、`04_knockout`、`05_validation` 下分别用 `figures/` 和 `data/` 区分结果图与结果数据。
- `outputs/run_001/results/01_analysis/data/fig_46_47_ranked_results.csv`、`fig_47_top_hits.csv`、`fig_48_diverse_hits.csv`、`docking_results.xlsx` 和编号 `fig_46` 至 `fig_49` 的结果图（含 `fig_49_redock_comparison.png` 与 `fig_49_redock_comparison.csv`）。
- `outputs/run_001/results/docking_report.html`：HTML 报告。
- `outputs/run_001/results/02_redock/data/fig_49_redock_results.csv`：精细重对接结果。
- `outputs/run_001/results/04_knockout/data/fig_52_53_ranked_knockout.csv`、`fig_52_target_candidates.csv`、`target_report.md`，以及编号 `fig_52`、`fig_53` 的敲除结果图。
- `outputs/run_001/results/05_validation/data/validation_candidates.csv`、`validation_plan.md`。
- `evidence/evidence_report.md`、`known_ligands.csv`。

附加数据文件格式：

| 文件 | 必需列 | 说明 |
| --- | --- | --- |
| expression.csv | `gene` + 样本列，或 `gene/sample/value` | 宽矩阵或长表 |
| metadata.csv | `sample` + `condition`，可加 `cell_type` | 分组和细胞类型注释 |
| depmap.csv | 宽表 `ModelID + 基因列`，或长表 `gene/effect` | CRISPR 基因效应 |
| prognosis.csv | `gene` + `hr` | 也支持 `hazard_ratio`、`cox_hr` |
| druggability.csv | `gene` + `known_ligands`、`pdb_structures`、`chembl_bioactivities` | 成药性计数 |
| off_target.csv | `gene` + `off_target_paralogs`、`safety_concern` | 脱靶风险 |

### 5.3 全自动集成流水线

输入：

- GSE 编号（默认 `GSE125449`）。
- 单细胞输出目录（必填，例如 `y2`）。
- 工作目录（必填，例如 `y3`）。
- 可选配体库、病例/正常标签、DepMap CSV。

输出（`<workdir>/outputs/integration/`）：

- `key_genes.csv`：关键基因排序表。
- `gene_evidence.csv`：每个基因的 UniProt、PDB、ChEMBL 证据。
- `knockout_summary.json`：虚拟敲除与验证方案汇总。
- `docking_targets.csv`：每个靶点的对接状态、命中数和最佳亲和力。
- `cell_feedback/`：细胞反馈阶段输出，包括 `data/cell_scores.csv`、`data/feedback_targets.csv`、`data/celltype_summary.csv`、`data/celltype_enrichment.csv`、`data/condition_summary.csv`，以及 `fig_54` 至 `fig_58` 的结果图。
- `integration_report.html`：全流程集成报告。
- `integration_summary.json` / `run_manifest.json`：本次运行的汇总和溯源信息。

每个靶点的对接在独立目录 `<workdir>/work/<gene>/` 下运行，支持单独断点续跑；配体优先使用 ChEMBL/BindingDB 已知活性分子，无数据库配体时自动提取共晶配体作为对照，最后回退到用户提供的配体库。

## 6. 脚本文件一览

| 脚本 | 说明 |
| --- | --- |
| `scripts/run_pipeline.py` | 单细胞分析 CLI 入口 |
| `scripts/run_docking.py` | 虚拟筛选 CLI 入口 |
| `scripts/run_full_pipeline.py` | 全自动集成流水线 CLI 入口 |
| `scripts/search_datasets.py` | GEO 数据集搜索与下载 |
| `scripts/dataset_search_ml.py` | 数据集搜索 ML/DL 相关性重排序 |
| `scripts/validate_pipeline.py` | 合成数据单细胞流水线验证 |
| `scripts/validate_dock_pipeline.py` | 假 Vina 对接流水线验证 |
| `scripts/validate_real_pipeline.py` | 真实 GEO 数据流水线验证 |
| `scripts/validate_real_evidence.py` | 真实 PDB 证据收集验证 |
| `scripts/validate_real_random.py` | 随机真实数据验证 |
| `scripts/validate_new_features.py` | 真实数据靶点评分/验证方案验证 |
| `scripts/validate_random_real_full_pipeline.py` | 随机真实 GSE 全流程验证 |
| `scripts/validate_dataset_search.py` | GEO 数据集搜索随机验证 |
| `launchers/check_*.bat/.py` | 环境检查 |
| `launchers/install_*.bat/.py` | 环境自动补全 |
| `launchers/run_web_ui.bat` | 启动网页端 |
| `launchers/run_docking.bat` | 虚拟筛选快捷入口 |
| `launchers/run_full_pipeline.bat` | 全自动流水线快捷入口 |
| `launchers/run_GSE125449.bat` | GSE125449 单细胞快捷入口 |
| `launchers/run_pipeline_prompt.bat` | 交互式单细胞入口 |
| `src/analysis/*` | R/Python 分析实现（QC、聚类、DEG、富集、CellChat、ML） |
| `src/data/*` | GEO 下载、格式转换、合成数据生成 |
| `src/docking/*` | 对接、证据、敲除、验证、报告等实现 |
| `src/pipeline/orchestrator.py` | 单细胞流水线编排 |
| `src/pipeline/integration.py` | 全自动集成流水线编排 |
| `src/pipeline/cell_feedback.py` / `cell_feedback.R` | 虚拟敲除/对接结果返回单细胞的闭环分析 |
| `src/pipeline/export_pseudobulk.R` | 伪 bulk 表达矩阵导出 |
| `src/report/*` | HTML/Word 报告生成 |
| `web/web_ui.py` | 本地网页服务 |
| `web/templates/*` | 三个页面的 HTML 模板 |
| `config/*.json` | 单细胞、对接和全流程配置 |
| `tests/test_*.py` | 单元/集成测试 |

## 7. 目录结构

```text
Script/
├── scripts/
│   ├── run_pipeline.py
│   ├── run_docking.py
│   ├── run_full_pipeline.py
│   ├── search_datasets.py
│   ├── dataset_search_ml.py
│   └── validate_*.py
├── README.md
├── AGENTS.md            # Codex 项目执行规则
├── VIRTUAL_SCREENING_REQUIREMENTS.md
├── requirements.txt
├── requirements_dock.txt
├── environment_dock.yml
├── config/
│   ├── project_config.json
│   ├── docking_config.json
│   └── full_pipeline_config.json
├── src/
│   ├── analysis/
│   ├── data/
│   ├── docking/
│   ├── pipeline/
│   └── report/
├── web/
│   ├── web_ui.py
│   ├── static/
│   └── templates/
├── launchers/
├── tests/
├── data_cache/          # 运行时下载缓存（gitignore）
├── dock/                # 虚拟筛选工作目录（产物 gitignore）
│   ├── config/
│   ├── data/
│   ├── outputs/
│   ├── evidence/
│   ├── validation_real/
│   └── tools/
└── results/             # 单细胞结果（gitignore）
```

`dock/tools/`、`dock/outputs/`、`dock/logs/`、`dock/evidence/`、`dock/validation_real/`、`dock/work/`、`data_cache/` 等运行产物和二进制文件默认被 `.gitignore` 排除，不会上传 GitHub。

## 8. 数据来源与许可

GSE125449: Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. PMID: 31588021

GSE165816 和 TCGA PanCancer Atlas 仅用于真实数据验证。

MIT License. See `LICENSE` for details.

## 9. 更新日志

### v0.4.2

- 网页版同步 bulk 数据集检测：数据集搜索页新增数据类型过滤与说明，bulk 数据集只可下载、不能直接运行全自动流水线。
- 全自动流水线页内嵌搜索新增数据类型列，bulk 数据集不可选择；带入链接携带 `data_type` 参数，启动前校验可拦截 bulk 数据集。
- 单细胞分析页标注仅支持单细胞数据集。
- 补充网页端 bulk 链接与搜索数据类型的单元测试。
- 单细胞富集分析优化 `fig_22_go_network.png` / `fig_23_kegg_network.png`：网络图先按 `p.adjust <= 0.05` 筛选，再取前 5 个通路，保留原图并降低网络杂乱度。
- 新增 `fig_46_go_top5.png` / `fig_47_kegg_top5.png`：展示上调基因 GO BP / KEGG 筛选后 Top5 富集结果。
- 同步报告、网页图开关和流水线输出校验。

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
