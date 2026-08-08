# Liver Cancer Bioinformatics Workflow

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
- 可选 CellChat 细胞通讯分析（设置 `LIVER_RUN_CELLCHAT=yes`）。
- ML 疾病分类（XGBoost 或 RandomForest）、特征重要性和 SHAP 解释。
- 输出 25 张分析图、HTML 报告，支持断点续跑、暂停/继续、停滞自动重启。

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
| `detect-box` | 从共晶配体自动检测对接盒并写回配置 |
| `check-env` / `check-cadd` | 检查软件、库和数据库 skill 环境 |
| `pipeline` | 按阶段执行准备、对接、分析、重对接、报告，支持断点续跑 |

### 2.3 全自动集成流水线

`run_full_pipeline.py` 把单细胞分析、关键基因筛选、证据富集、虚拟敲除和虚拟筛选串成一条流水线：

```text
01 single_cell           GEO 单细胞分析（下载、QC、注释、差异表达、富集）
02 key_targets           从显著 DEG 中筛选并排序关键基因/蛋白
03 evidence              UniProt / PDB / ChEMBL 证据富集（带本地缓存）
04 knockout_inputs       导出样本级伪 bulk 表达矩阵并生成敲除输入
05 knockout              虚拟敲除 + 多维靶点评分 + 湿实验验证方案
06 docking               对有 PDB 结构的靶点自动收集已知配体并跑 Vina 对接
07 report                生成集成 HTML 报告和 run_manifest.json
```

每一阶段写标记文件，重跑时自动断点续跑；`--start-stage` 可从任意阶段开始。

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

`web/web_ui.py` 提供本地网页端，顶部导航包含三个页面：

- 单细胞分析：`/`
- 虚拟筛选：`/dock`
- 全自动流水线：`/full`

网页端支持任务启动、实时日志、暂停/继续、结果表和文件下载、环境检查与自动补全；全自动流水线页还支持直接填写工作目录加载已有结果。

### 2.6 真实数据验证与可复现性

- `validate_new_features.py` 使用 20 个 TCGA PanCancer Atlas 队列和 GSE165816 真实单细胞数据运行虚拟敲除与验证方案导出。
- 其余 `validate_*.py` 分别验证合成数据流水线、对接流水线、真实 GEO 数据、证据收集和随机真实数据。
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
python run_pipeline.py GSE125449 --output ../liver_cancer --species auto
```

常用参数：

```bash
python run_pipeline.py GSE125449 --output ../liver_cancer --species hs --force
python run_pipeline.py GSE125449 --output ../liver_cancer --skip-download
python run_pipeline.py GSE125449 --output ../liver_cancer --skip-deps
```

Windows 下也可以直接使用：

```text
launchers\run_GSE125449.bat
launchers\run_pipeline_prompt.bat
```

### 4.2 虚拟筛选命令行

先初始化工作目录：

```bash
python run_docking.py init
```

运行完整对接流程（准备受体、准备配体、对接、分析、精细重对接、HTML 报告，支持断点续跑）：

```bash
python run_docking.py pipeline --config config/docking_config.json
```

分阶段运行：

```bash
python run_docking.py evidence --uniprot P00533 --pdb 1M17 --target-name EGFR
python run_docking.py prepare-receptor
python run_docking.py prepare-ligands
python run_docking.py dock
python run_docking.py analyze
python run_docking.py redock
python run_docking.py report
```

ML/DL 重打分：

```bash
python run_docking.py ml-train --training-csv data/ml/training.csv --model rf
python run_docking.py ml-predict
```

导出 MD/外部交接模板：

```bash
python run_docking.py export-md
python run_docking.py export-external
```

虚拟敲除（基础评分）：

```bash
python run_docking.py virtual-knockout \
  --expression-csv data/knockout/expression.csv \
  --metadata-csv data/knockout/metadata.csv \
  --depmap-csv data/knockout/depmap_gene_effect.csv \
  --case-label Tumor --normal-label Normal
```

多维评分与验证方案导出：

```bash
python run_docking.py virtual-knockout \
  --expression-csv data/knockout/expression.csv \
  --metadata-csv data/knockout/metadata.csv \
  --prognosis-csv data/knockout/prognosis.csv \
  --druggability-csv data/knockout/druggability.csv \
  --off-target-csv data/knockout/off_target.csv \
  --cell-type-column cell_type \
  --case-label Tumor --normal-label Normal

python run_docking.py export-validation --validation-top-n 10
```

环境检查：

```bash
python run_docking.py check-env
python run_docking.py check-cadd
```

### 4.3 全自动集成流水线

一键启动：

```text
launchers\run_full_pipeline.bat
```

或直接运行：

```bash
python run_full_pipeline.py \
  --accession GSE125449 \
  --output ../liver_cancer \
  --workdir dock \
  --top-genes 50 \
  --docking-targets 3
```

常用参数：

- `--skip-scrna`：复用已完成的单细胞结果，直接从关键基因筛选开始。
- `--skip-docking`：只跑虚拟敲除和验证方案，跳过对接。
- `--skip-evidence-fetch`：不联网，使用已有证据缓存或置零。
- `--skip-download` / `--skip-deps` / `--skip-pseudobulk` / `--skip-knockout`。
- `--top-genes`：关键基因数量，默认 50。
- `--docking-targets`：参与对接的靶点数量，默认 3。
- `--ligand-library`：自定义配体库（`.smi` / `.sdf` / `.csv`），也可放到 `dock/data/ligands/`。
- `--case-label` / `--normal-label`：虚拟敲除的病例/正常分组标签。
- `--start-stage 07`：从指定阶段继续，之前阶段自动标记为跳过。

查看阶段清单：

```bash
python run_full_pipeline.py --list-stages
```

### 4.4 网页版

```text
launchers\run_web_ui.bat
```

浏览器打开 `http://127.0.0.1:8000`。直接打开指定页面：

```text
launchers\run_web_ui.bat --page dock
launchers\run_web_ui.bat --page full
```

### 4.5 验证脚本与测试

```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

各验证脚本用途：

| 脚本 | 用途 |
| --- | --- |
| `validate_pipeline.py` | 用 10 个合成单细胞数据集跑通完整单细胞流水线 |
| `validate_dock_pipeline.py` | 用假 Vina 可执行文件验证对接流水线 |
| `validate_real_pipeline.py` | 用 10 个真实肝病 GEO 数据集跑通完整流水线 |
| `validate_real_evidence.py` | 用 10 个真实 PDB 结构验证证据收集 |
| `validate_real_random.py` | 随机真实数据验证证据收集和对接盒检测 |
| `validate_new_features.py` | 用 20 个 TCGA PanCancer 队列 + GSE165816 验证靶点评分 |

## 5. 输入输出示例

### 5.1 单细胞分析

输入：

- GSE 编号，例如 `GSE125449`。
- 输出目录，例如 `../liver_cancer`。
- 物种：`hs` / `mm` / `auto`。

输出（以 `../liver_cancer` 为例）：

- `results/figures/`：25 张结果图（R 图 + ML 解释图）。
- `results/data/`：QC、双细胞、注释、差异表达、富集、ML 分类和可选 CellChat 表格。
- `results/data/deg_significant.csv`：显著差异基因表。
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
- `outputs/run_001/reports/ranked_results.csv`、`hits.csv`、`diverse_hits.csv`、`docking_results.xlsx` 和图。
- `outputs/run_001/reports/docking_report.html`：HTML 报告。
- `outputs/run_001/redock/results.csv`：精细重对接结果。
- `outputs/run_001/knockout/ranked_knockout.csv`、`target_candidates.csv`、`target_report.md`。
- `outputs/run_001/validation/validation_candidates.csv`、`validation_plan.md`。
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
- 单细胞输出目录（默认 `../liver_cancer`）。
- 工作目录（默认 `dock`）。
- 可选配体库、病例/正常标签、DepMap CSV。

输出（`dock/outputs/integration/`）：

- `key_genes.csv`：关键基因排序表。
- `gene_evidence.csv`：每个基因的 UniProt、PDB、ChEMBL 证据。
- `knockout_summary.json`：虚拟敲除与验证方案汇总。
- `docking_targets.csv`：每个靶点的对接状态、命中数和最佳亲和力。
- `integration_report.html`：全流程集成报告。
- `integration_summary.json` / `run_manifest.json`：本次运行的汇总和溯源信息。

每个靶点的对接在独立目录 `dock/work/<gene>/` 下运行，支持单独断点续跑；配体优先使用 ChEMBL/BindingDB 已知活性分子，无数据库配体时自动提取共晶配体作为对照，最后回退到用户提供的配体库。

## 6. 脚本文件一览

| 脚本 | 说明 |
| --- | --- |
| `run_pipeline.py` | 单细胞分析 CLI 入口 |
| `run_docking.py` | 虚拟筛选 CLI 入口 |
| `run_full_pipeline.py` | 全自动集成流水线 CLI 入口 |
| `validate_pipeline.py` | 合成数据单细胞流水线验证 |
| `validate_dock_pipeline.py` | 假 Vina 对接流水线验证 |
| `validate_real_pipeline.py` | 真实 GEO 数据流水线验证 |
| `validate_real_evidence.py` | 真实 PDB 证据收集验证 |
| `validate_real_random.py` | 随机真实数据验证 |
| `validate_new_features.py` | 真实数据靶点评分/验证方案验证 |
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
| `src/pipeline/export_pseudobulk.R` | 伪 bulk 表达矩阵导出 |
| `src/report/*` | HTML/Word 报告生成 |
| `web/web_ui.py` | 本地网页服务 |
| `web/templates/*` | 三个页面的 HTML 模板 |
| `config/*.json` | 单细胞、对接和全流程配置 |
| `tests/test_*.py` | 单元/集成测试 |

## 7. 目录结构

```text
Script/
├── run_pipeline.py
├── run_docking.py
├── run_full_pipeline.py
├── validate_*.py
├── README.md
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
