# Liver Cancer Bioinformatics Workflow

面向肝癌研究的本地生信工作流，整合两条流水线：

- 单细胞转录组分析：输入 GEO 单细胞数据集编号，自动完成下载、质控、双细胞检测、聚类、注释、差异表达、富集分析和可选的细胞通讯分析。
- CADD 虚拟筛选：从靶点证据收集开始，完成受体/配体准备、AutoDock Vina 并行对接、结果分析、ML/DL 重打分、精细重对接，并可导出分子动力学和外部对接交接模板。

两条流水线共用同一个本地网页端，可通过顶部导航切换。

## 1. 项目解决什么问题

单细胞分析和虚拟筛选通常依赖多个分散工具和手工路径，结果难以复现，也难以从“细胞图谱”推进到“可干预靶点”。本项目把以下环节串成一条可检查、可断点续跑、可交付的本地链路：

- GEO 单细胞数据从下载、质控、注释到差异表达与富集分析的完整分析链。
- 靶点证据（UniProt、RCSB PDB、ChEMBL、BindingDB、PubChem、ChEBI）收集。
- 受体/配体准备、Vina 并行对接、命中排序、精细重对接与 ML/DL 重打分。
- 虚拟敲除的多维靶点评分：表达差异、增殖共表达、共表达网络 hub、DepMap 依赖，以及可选的疾病逆转、通路控制、细胞类型特异性、临床预后、成药性和脱靶风险。
- 把排序后的候选靶点导出为分阶段湿实验验证方案（细胞系、类器官、药物剂量反应、动物模型、PDX）。
- 每次运行写入 `run_manifest.json`，记录配置哈希、输入文件 SHA256、软件版本和本次参数，保证结果可复现、可追溯。
- 内置真实数据验证脚本，可用 20 个 TCGA PanCancer 队列和本地缓存单细胞数据自动跑通多维靶点评分、运行清单和验证方案导出，并生成汇总报告。
- 网页端实时日志、分析历史、软件环境检测与自动补全、任务暂停/继续。

## 2. 主要功能

### 2.1 单细胞分析

- GEO 数据自动下载与格式识别，内置 GSE125449 适配。
- QC、双细胞检测（`scDblFinder`）、PCA/UMAP 聚类和细胞注释。
- 差异表达（DESeq2 pseudobulk 或 Wilcoxon 回退）与 GO/KEGG/GSEA 富集分析。
- 23 张 R 分析结果图加 2 张 ML 解释图（共 25 张），逐张可开关，火山图/MA 图、气泡图/柱状图/cnet 图等样式可切换。
- 可选 CellChat 细胞通讯分析（设置环境变量 `LIVER_RUN_CELLCHAT=yes` 后运行）。
- ML 疾病分类（XGBoost 或 RandomForest 回退）、特征重要性和 SHAP 解释；样本或分组不足时自动输出 skipped 状态。
- HTML 报告、断点续跑、暂停/继续、停滞检测自动重启。

### 2.2 虚拟筛选（CADD）

| 命令 | 说明 |
| --- | --- |
| `init` | 创建 docking 工作目录骨架和配置文件 |
| `evidence` | 调用数据库 skill 收集靶点和已知配体证据 |
| `prepare-receptor` | 通过 Meeko/Open Babel/MGLTools 准备受体 PDBQT |
| `prepare-ligands` | RDKit 标准化、去盐、3D 构象生成并输出 PDBQT |
| `dock` | AutoDock Vina 并行对接，按配体写入结果，支持断点续跑 |
| `analyze` | 按亲和力排序、阈值筛选、Tanimoto 多样性选择，输出 CSV/Excel/图片 |
| `redock` | 对 Top 命中用更高 exhaustiveness 精细重对接 |
| `detect-box` | 从共晶配体自动检测对接盒并写回配置 |
| `ml-train` / `ml-predict` | 随机森林、GBDT、MLP 或 PyTorch MLP 重打分 |
| `export-md` / `export-external` | 导出 Amber/GROMACS 和 UniDock-Pro/HDOCK/HADDOCK 模板 |
| `report` | 生成 HTML 汇总报告 |
| `virtual-knockout` | 基因敲除优先级和多维靶点评分 |
| `export-validation` | 把排序后的靶点导出为分阶段湿实验验证方案 |
| `check-env` / `check-cadd` | 检查软件、库和数据库 skill 环境 |
| `pipeline` | 按阶段执行准备、对接、分析、重对接、报告，支持断点续跑 |

`virtual-knockout` 的核心评分 `knockout_score` 由表达差异、增殖共表达、共表达网络 hub 程度和 DepMap CRISPR 依赖加权得到。提供附加数据后，还会计算：

- `reversal_score`：表达变化方向与疾病模块的一致性，衡量“逆转疾病特征”的潜力。
- `pathway_score`：与细胞周期、凋亡、EMT、p53、PI3K-AKT 等疾病通路基因集的平均相关性。
- `specificity_score`：基于注释细胞类型的熵值表达特异性，降低“误伤友军”风险。
- `prognosis_score`：方向感知的风险比评分。
- `druggability_score`：已知配体、蛋白结构和生物活性计数。
- `safety_concern` 和 `off_target_paralogs`：脱靶/旁系同源标记，可对最终得分施加罚分。

以上维度按 `target_weights` 加权合并为 `target_score`，并将候选基因分为 `core_driver`、`microenvironment_regulator`、`biomarker`、`high_priority`、`low_priority`。

### 2.3 网页端

- 单细胞分析和虚拟筛选共用 `web_ui.py` 服务，通过顶部导航切换。
- 虚拟筛选页支持完整流程、分阶段选择、受体/配体路径、对接盒参数、强制重跑。
- 支持证据收集、ML/DL 重打分、精细重对接、MD/外部交接、HTML 报告等阶段。
- 内置虚拟敲除表单，可填写表达矩阵、样本元数据、DepMap 和分组标签。
- 虚拟敲除结果表直接展示 `target_score`、`target_class`、逆转分、通路分、特异性分、预后分、成药性分和安全标记。
- 提供“导出湿实验验证方案”按钮，一键生成验证方案文件。
- 提供“真实数据验证报告”卡片，可查看 20 个 TCGA 队列 + GSE165816 的验证汇总，也可从网页后台重新运行验证。
- 显示分析历史和软件环境；提供“检查环境”和“自动补全环境”。
- 任务支持暂停/继续，暂停后可从断点恢复。

### 2.4 真实数据验证

`validate_new_features.py` 会使用真实公共数据对新增功能做端到端验证：

- 通过 cBioPortal 公共 API 拉取 20 个 TCGA PanCancer 队列的真实 mRNA 表达和 OS 生存数据。
- 通过 GEO 缓存使用 GSE165816 真实单细胞数据（50 个样本）。
- 通过 ChEMBL、RCSB PDB 公共 API 抓取每个基因的真实活性计数和结构计数作为成药性证据。
- 对每个数据集运行 `virtual-knockout` 和 `export-validation`，并校验 `run_manifest.json` 与验证方案文件是否生成。
- 输出 `validation_report.md` 和 `validation_summary.csv`，判定真实数据集成功数是否达到 20 个。

说明：TCGA PanCancer 队列均为肿瘤样本，因此 `reversal_score` 与 `specificity_score` 主要依靠含分组的 GSE165816 验证；`off_target_paralogs` 在外部同源数据库不可用时按 0 记录。

## 3. 安装方法

### 环境要求

- Python 3.10+（推荐 3.11）。
- R 4.5+（仅单细胞分析需要）。
- AutoDock Vina（虚拟筛选需要，可放在 `dock/tools/vina.exe` 或加入 PATH）。
- Codex skills（仅证据收集需要）：`uniprot-skill`、`rcsb-pdb-skill`、`chembl-skill`、`bindingdb-skill`、`pubchem-pug-skill`、`chebi-skill`。
- 可选：Amber/AmberTools、GROMACS、UniDock-Pro、HDOCK、HADDOCK、xgboost、shap、torch。

### 安装步骤

复制整个项目文件夹后，在项目根目录执行环境检查：

```text
launchers\check_pipeline_environment.bat
```

环境不满足时安装单细胞依赖（R 包和基础 Python 包）：

```text
launchers\install_pipeline_dependencies.bat
```

检查并安装虚拟筛选 Python 依赖（RDKit、Meeko、Open Babel、AutoDockTools 等）：

```text
launchers\check_dock_environment.bat
launchers\install_dock_dependencies.bat
```

网页端也可以在“软件环境”区域填写安装地址后点击“自动补全环境”。

完整软件/插件清单见 `VIRTUAL_SCREENING_REQUIREMENTS.md`。

## 4. 使用方法

### 网页启动

```text
launchers\run_web_ui.bat
```

浏览器打开 `http://127.0.0.1:8000`。直接打开虚拟筛选页：

```text
launchers\run_web_ui.bat --page dock
```

### 单细胞命令行

```bash
python run_pipeline.py GSE125449 --output results/GSE125449 --species hs
```

`--species` 可选 `hs`（人）、`mm`（小鼠）或 `auto`（默认）。其他参数：

```bash
python run_pipeline.py GSE125449 --output results/GSE125449 --species hs --force
python run_pipeline.py GSE125449 --output results/GSE125449 --skip-download
python run_pipeline.py GSE125449 --output results/GSE125449 --skip-deps
```

### 虚拟筛选命令行

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

导出交接模板：

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
  --depmap-csv data/knockout/depmap_gene_effect.csv \
  --prognosis-csv data/knockout/prognosis.csv \
  --druggability-csv data/knockout/druggability.csv \
  --off-target-csv data/knockout/off_target.csv \
  --cell-type-column cell_type \
  --case-label Tumor --normal-label Normal

python run_docking.py export-validation --validation-top-n 10
```

真实数据验证（自动拉取 20 个 TCGA PanCancer 队列并运行全部 21 个数据集）：

```bash
python validate_new_features.py --max-studies 20
```

已在本地建好缓存后，可以只复用已有数据集并跳过网络构建：

```bash
python validate_new_features.py --max-studies 20 --skip-build
```

环境检查：

```bash
python run_docking.py check-env
python run_docking.py check-cadd
```

Windows 下也可以使用包装脚本：

```text
launchers\run_docking.bat pipeline
launchers\run_GSE125449.bat
launchers\run_pipeline_prompt.bat
```

## 5. 输入输出示例

### 单细胞分析

输入：

- GSE 编号，例如 `GSE125449`。
- 输出目录，例如 `results/GSE125449`。
- 物种，`hs` / `mm` / `auto`。
- 可选参数：`--force`、`--skip-download`、`--skip-deps`。

输出（以 `results/GSE125449` 为例）：

- `results/GSE125449/results/figures/`：结果图（23 张 R 图 + 2 张 ML 图）。
- `results/GSE125449/results/data/`：QC、双细胞、注释、差异表达、富集、ML 分类和 CellChat 表格。
- `results/GSE125449/results/checkpoints/`：Seurat 断点对象。
- `results/GSE125449/results/result_report.html`：最终 HTML 报告。

### 虚拟筛选

输入：

- 受体文件：`dock/data/receptors/receptor.pdb`。
- 配体库：`dock/data/ligands/library.sdf`（也支持 `.smi` 和 CSV）。
- 对接盒中心和尺寸：`config/docking_config.json`，或用 `detect-box` 自动检测。
- 训练标签（可选）：`dock/data/ml/training.csv`，包含 `smiles` 和 `active` 或 `affinity` 列。
- 敲除表达矩阵（可选）：`dock/data/knockout/expression.csv`。
- 样本元数据（可选）：`dock/data/knockout/metadata.csv`。
- DepMap 依赖表（可选）：`dock/data/knockout/depmap_gene_effect.csv`。
- 预后表（可选）：`dock/data/knockout/prognosis.csv`。
- 成药性表（可选）：`dock/data/knockout/druggability.csv`。
- 脱靶表（可选）：`dock/data/knockout/off_target.csv`。

附加数据文件格式：

| 文件 | 必需列 | 说明 |
| --- | --- | --- |
| expression.csv | `gene` + 样本列，或 `gene/sample/value` | 宽矩阵或长表 |
| metadata.csv | `sample` + `condition`，可加 `cell_type` | 分组和细胞类型注释 |
| depmap.csv | 宽表 `ModelID + 基因列`，或长表 `gene/effect` | CRISPR 基因效应 |
| prognosis.csv | `gene` + `hr` | 也支持 `hazard_ratio`、`cox_hr` 等列名 |
| druggability.csv | `gene` + `known_ligands`、`pdb_structures`、`chembl_bioactivities` | 成药性计数 |
| off_target.csv | `gene` + `off_target_paralogs`、`safety_concern` | 脱靶风险 |

输出（工作目录默认为 `dock`）：

- 受体 PDBQT：`dock/data/receptors/receptor.pdbqt`。
- 配体准备：`dock/data/ligands/prepared/` 和 `manifest.csv`。
- 对接结果：`dock/outputs/run_001/docked/results.csv`。
- 分析报告：`dock/outputs/run_001/reports/ranked_results.csv`、`hits.csv`、`diverse_hits.csv`、`docking_results.xlsx` 和图片。
- HTML 汇总报告：`dock/outputs/run_001/reports/docking_report.html`。
- ML 重打分：`dock/outputs/run_001/reports/ml_ranked_results.csv`、`ml_feature_importance.csv`。
- 精细重对接：`dock/outputs/run_001/redock/results.csv`。
- 证据报告：`dock/evidence/evidence_report.md` 和 `known_ligands.csv`。
- MD 交接：`dock/outputs/md/`。
- 虚拟敲除：`dock/outputs/run_001/knockout/ranked_knockout.csv`、`target_candidates.csv`、`target_report.md`、`run_manifest.json` 和结果图。
- 验证方案：`dock/outputs/run_001/validation/validation_candidates.csv`、`validation_plan.md`、`run_manifest.json`。
- 真实数据验证：`dock/validation_real/pan_cancer_20/validation_report.md`、`validation_summary.csv`，每个数据集下有独立的表达/元数据/预后输入和 `work/outputs/run_001/` 输出。

## 6. 目录结构

```text
Script/
├── README.md
├── requirements.txt
├── requirements_dock.txt
├── environment_dock.yml
├── VIRTUAL_SCREENING_REQUIREMENTS.md
├── run_pipeline.py
├── run_docking.py
├── validate_pipeline.py
├── validate_dock_pipeline.py
├── validate_real_pipeline.py
├── validate_real_evidence.py
├── validate_real_random.py
├── validate_new_features.py
├── config/
│   ├── project_config.json
│   └── docking_config.json
├── src/
│   ├── pipeline/
│   ├── data/
│   ├── analysis/
│   ├── report/
│   └── docking/
│       ├── config.py
│       ├── cli.py
│       ├── pipeline.py
│       ├── evidence.py
│       ├── knockout.py
│       ├── provenance.py
│       ├── validation.py
│       └── ...
├── web/
│   ├── web_ui.py
│   ├── static/
│   └── templates/
├── launchers/
├── tests/
│   ├── test_docking.py
│   ├── test_knockout.py
│   └── test_target_scoring.py
└── dock/
    ├── config/
    ├── data/
    │   ├── receptors/
    │   ├── ligands/
    │   └── knockout/
    ├── outputs/
    ├── evidence/
    ├── validation_real/
    └── tools/
```

`dock/tools/`、`dock/outputs/`、`dock/logs/`、`dock/evidence/`、`data_cache/` 等运行产物和二进制文件默认被 `.gitignore` 排除，不上传 GitHub。

## 7. 验证与测试

```bash
python -m unittest discover -s tests -p "test_*.py" -v
python validate_pipeline.py
python validate_dock_pipeline.py
python validate_real_pipeline.py
python validate_real_evidence.py
python validate_real_random.py
python validate_new_features.py --max-studies 20
```

## 8. 数据来源

GSE125449: Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer.

PMID: 31588021

## License

MIT License. See `LICENSE` for details.
