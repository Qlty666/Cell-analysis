# Liver Cancer Bioinformatics Workflow

面向肝癌研究的本地生信工作流，整合两套流水线：

- 单细胞转录组分析：输入任意 GEO 单细胞数据集编号，自动完成下载、质控、双细胞检测、聚类、注释、差异表达和富集分析。
- CADD 虚拟筛选：从靶点/结构证据收集开始，完成受体准备、配体库准备、AutoDock Vina 并行对接、结果分析、ML/DL 重打分，并可导出 Amber/GROMACS 和外部对接交接模板。

两套流水线共用同一个本地网页端，可通过顶部导航切换。

## 1. 项目解决什么问题

单细胞分析和虚拟筛选通常依赖多个分散工具和手工路径，结果难以复现。本项目把以下环节串成可检查、可断点续跑、可交付的本地流程：

- GEO 单细胞数据从下载到富集分析的完整分析链。
- 靶点证据（UniProt、RCSB PDB、ChEMBL、BindingDB、PubChem、ChEBI）收集。
- 受体/配体准备、Vina 并行对接、命中排序。
- 随机森林、GBDT、MLP、PyTorch MLP 的 ML/DL 重打分。
- Amber/GROMACS 分子动力学交接和 UniDock-Pro/HDOCK/HADDOCK 外部对接模板。
- 网页端实时日志、分析历史、软件环境检测与自动补全、任务暂停/继续。

## 2. 主要功能

### 单细胞分析

- GEO 数据自动下载与格式识别，内置 GSE125449 适配。
- QC、双细胞检测（`scDblFinder`）、PCA/UMAP 聚类和细胞注释。
- 差异表达（DESeq2 pseudobulk 或 Wilcoxon 回退）与 GO/KEGG/GSEA 富集。
- 24 张论文常用结果图，逐张可开关。
- HTML 报告、断点续跑、暂停/继续、停滞检测自动重启。

### 虚拟筛选

- `evidence`：调用数据库 skill 收集靶点和已知配体证据。
- `prepare-receptor`：通过 Meeko/Open Babel/MGLTools 准备受体 PDBQT。
- `prepare-ligands`：RDKit 标准化、去盐、3D 构象生成，Meeko/Open Babel 输出 PDBQT。
- `dock`：AutoDock Vina 并行对接，按配体写入结果，支持断点续跑。
- `analyze`：按亲和力排序、阈值筛选、Tanimoto 多样性选择，输出 CSV/Excel/图片。
- `ml-train` / `ml-predict`：随机森林、GBDT、MLP 或 PyTorch MLP 重打分。
- `export-md` / `export-external`：导出 Amber/GROMACS 和 UniDock-Pro/HDOCK/HADDOCK 模板。
- `check-env` / `check-cadd`：检查软件、库和数据库 skill 环境。

### 网页端

- 单细胞分析和虚拟筛选共用 `web_ui.py` 服务。
- 虚拟筛选页支持证据收集、ML 重打分、MD 交接等阶段选择。
- 显示分析历史和软件环境；提供“检查环境”和“自动补全环境”。
- 任务支持暂停/继续，暂停后可从断点恢复。

## 3. 安装方法

### 环境要求

- Python 3.10+（推荐 3.11）
- R 4.5+（仅单细胞分析需要）
- AutoDock Vina（虚拟筛选需要，可放在 `dock/tools/vina.exe` 或加入 PATH）
- Codex skills（仅证据收集需要）：`uniprot-skill`、`rcsb-pdb-skill`、`chembl-skill`、`bindingdb-skill`、`pubchem-pug-skill`、`chebi-skill`

### 安装步骤

复制整个项目文件夹后，在项目根目录执行：

```text
launchers\check_pipeline_environment.bat
```

环境不满足时安装单细胞依赖：

```text
launchers\install_pipeline_dependencies.bat
```

安装虚拟筛选 Python 依赖（RDKit、Meeko、Open Babel、AutoDockTools 等）：

```text
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

强制从头重跑：

```bash
python run_pipeline.py GSE125449 --output results/GSE125449 --species hs --force
```

### 虚拟筛选命令行

先初始化工作目录：

```bash
python run_docking.py init
```

运行完整流程：

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

环境检查：

```bash
python run_docking.py check-env
python run_docking.py check-cadd
```

## 5. 输入输出示例

### 单细胞分析

输入：

- GSE 编号，例如 `GSE125449`
- 输出目录，例如 `results/GSE125449`
- 物种，例如 `hs`

输出：

- `results/GSE125449/results/figures/`：结果图
- `results/GSE125449/results/data/`：QC、双细胞、注释、差异表达和富集表格
- `results/GSE125449/results/result_report.html`：最终报告

### 虚拟筛选

输入：

- 受体文件：`dock/data/receptors/receptor.pdb`
- 配体库：`dock/data/ligands/library.sdf`（也支持 `.smi` 和 CSV）
- 对接盒中心和尺寸：`config/docking_config.json`
- 训练标签（可选）：`dock/data/ml/training.csv`，包含 `smiles` 和 `active` 或 `affinity` 列

输出：

- 受体 PDBQT：`dock/data/receptors/receptor.pdbqt`
- 配体准备：`dock/data/ligands/prepared/` 和 `manifest.csv`
- 对接结果：`dock/outputs/run_001/docked/results.csv`
- 分析报告：`dock/outputs/run_001/reports/`
- ML 重打分：`dock/outputs/run_001/reports/ml_ranked_results.csv`
- 证据报告：`dock/evidence/evidence_report.md`
- MD 交接：`dock/outputs/md/`

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
├── validate_real_evidence.py
├── config/
│   ├── project_config.json
│   └── docking_config.json
├── src/
│   ├── pipeline/
│   ├── data/
│   ├── analysis/
│   ├── report/
│   └── docking/
├── web/
│   ├── web_ui.py
│   └── templates/
├── launchers/
├── tests/
└── dock/
    ├── config/
    ├── data/
    │   ├── receptors/
    │   └── ligands/
    ├── outputs/
    └── tools/
```

`dock/tools/`、`dock/outputs/`、`dock/logs/`、`dock/evidence/` 等运行产物和二进制文件默认被 `.gitignore` 排除，不上传 GitHub。

## 7. 验证与测试

```bash
python -m unittest discover -s tests -p "test_*.py" -v
python validate_dock_pipeline.py
python validate_real_evidence.py
```

## 8. 数据来源

GSE125449: Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer.

PMID: 31588021

## License

MIT License. See `LICENSE` for details.
