# 环境要求（按功能板块）

本页用于在新电脑上按需准备运行环境。准备前先确认要运行哪个功能板块，再安装该板块列出的软件和依赖；不需要一次性装完全部可选工具。

整仓复制到新电脑后，`dock/tools/`、`data_cache/`、`results/` 等运行产物和二进制文件通常不会被一起复制，需要重新准备第三方软件与 Python/R 依赖。

## 速查表

| 功能板块 | 运行入口 | 基础软件 | Python 依赖 | R 环境 | 检查 / 安装 |
| --- | --- | --- | --- | --- | --- |
| 表达分析 | `scripts/run_pipeline.py`、`liverbio expression` | Python 3.10+、R 4.5+ | `requirements.txt` | 需要（Seurat 等） | `check_pipeline_environment.bat` / `install_pipeline_dependencies.bat` |
| 数据集搜索 | `scripts/search_datasets.py`、`liverbio datasets` | Python 3.10+ | 标准库为主；ML 重排需 `requirements.txt` | 不需要 | 无独立检查脚本 |
| 虚拟筛选 | `scripts/run_docking.py`、`liverbio docking` | Python 3.10+、AutoDock Vina | `requirements.txt` + `requirements_dock.txt` + AutoDockTools | 仅部分子命令需要 | `check_dock_environment.bat` / `install_dock_dependencies.bat` |
| 独立分子对接 | `scripts/run_molecular_docking.py` | Python 3.10+、AutoDock Vina | `requirements.txt` + `requirements_dock.txt` + AutoDockTools | 不需要 | 同上（docking 检查脚本） |
| 分子动力学 | `scripts/run_docking.py md-simulation` | Python、Vina、GROMACS | 上述对接依赖；可选 `torch` | 不需要 | docking 检查 + MD 运行时检查 |
| 全自动集成流水线 | `scripts/run_full_pipeline.py`、`liverbio full` | Python、R、Vina | 两套 requirements + AutoDockTools | 需要 | 两个检查 / 安装脚本都要执行 |
| 网页版 | `web/web_ui.py`、`liverbio web` | Python | 按实际使用的页面功能安装 | 按实际使用的页面功能安装 | 页面内“环境状态与自动补全” |
| Codex Skills 证据收集 | `scripts/run_docking.py evidence`、全自动 stage 03 | Python + 对应 Codex skills | `requirements.txt` 基础包 | 不需要 | `python scripts/install_codex_skills.py` + `check-cadd` |

## 1. 表达分析（单细胞 / bulk RNA-seq / microarray）

**运行入口：**

- `scripts/run_pipeline.py`
- `launchers/run_GSE125449.bat`
- `launchers/run_pipeline_prompt.bat`
- `liverbio expression ...`
- 网页端“表达分析”页

**软件：**

- Python 3.10+，推荐 3.11。
- R 4.5+，且 `Rscript.exe` 已加入 `PATH`。若 R 安装在自定义目录，需要把对应的 `bin` 目录加入 `PATH`。
- 能访问 NCBI GEO / EBI 的网络连接（下载数据用）。

**Python 包：**

```bat
python -m pip install -r requirements.txt
```

基础检查覆盖：`numpy`、`pandas`、`matplotlib`、`pyyaml`、`openpyxl`、`h5py`、`scipy`、`scikit-learn`、`fpdf2`。`requirements.txt` 还包含 `umap-learn`、`scTenifoldpy>=0.3.0`，供虚拟敲除 / 全自动流水线等扩展模块使用。

**R 包：**

```bat
launchers\install_pipeline_dependencies.bat
```

实际安装清单以 `src/analysis/install_deps.R` 为准，主要包括：

- CRAN：`Seurat`、`dplyr`、`ggplot2`、`patchwork`、`Matrix`、`data.table`、`jsonlite`、`ggrepel`、`pheatmap`、`RColorBrewer`、`harmony`、`R.utils`、`WGCNA`、`survival`、`survminer`、`timeROC`、`glmnet`
- Bioconductor：`BiocManager`、`scDblFinder`、`SingleCellExperiment`、`clusterProfiler`、`org.Hs.eg.db`、`org.Mm.eg.db`、`enrichplot`、`BiocParallel`、`SingleR`、`celldex`、`DESeq2`、`hdf5r`、`limma`、`edgeR`、`sva`、`GSVA`、`celda`

可选功能：CellChat（`LIVER_RUN_CELLCHAT=yes`）、slingshot 拟时序（`LIVER_RUN_TRAJECTORY=yes`）、decontX（`LIVER_DECONTX=yes`）和指定细胞类型再聚类（`LIVER_SUBCLUSTER_CELLTYPES=...`）需要额外单独安装对应 R 包。

多队列、WGCNA、基因级 ML、免疫浸润和生存分析入口：

```bat
python scripts\run_advanced_analysis.py --config config\advanced_analysis.json --output <OUTPUT_DIR>\advanced
liverbio advanced --config config\advanced_analysis.json --output <OUTPUT_DIR>\advanced
```

本地 MR/共定位入口：

```bat
python scripts\run_mr_coloc.py --config config\mr_coloc.json --output <OUTPUT_DIR>\mr
liverbio mr --config config\mr_coloc.json --output <OUTPUT_DIR>\mr
```

MR 的 R 后端可选依赖 `MendelianRandomization` 和 `coloc`；未安装时 Python 会输出 IVW、weighted median、MR-Egger、异质性和 leave-one-out 结果。配置 `clump.bfile` 与 `clump.plink_executable` 时需要本地 PLINK 1.9 和 LD 参考数据；否则使用透明标记的距离剪枝，并自动剔除链方向不确定的回文 SNP。

**检查：**

```bat
launchers\check_pipeline_environment.bat
```

## 2. 数据集搜索

**运行入口：**

- `scripts/search_datasets.py`
- `scripts/dataset_search_ml.py`
- `liverbio datasets ...`
- 网页端“数据集搜索”页

**依赖：**

- Python 3.10+。
- 能访问 NCBI GEO、EBI BioStudies/ArrayExpress、Expression Atlas 的网络连接。
- 基础搜索不依赖第三方包；ML/DL 相关性重排需要 `requirements.txt` 中的 `scikit-learn` 等依赖和训练标签数据。

无独立 R 环境要求，也没有独立的环境检查脚本；可用 `liverbio doctor pipeline` 检查公共 Python 基础包。

## 3. 虚拟筛选（CADD）

**运行入口：**

- `scripts/run_docking.py`，子命令见 `python scripts\run_docking.py --help`
- `launchers/run_docking.bat`
- `liverbio docking ...`
- 网页端“虚拟筛选”页

**软件：**

- Python 3.10+，推荐 3.11。
- AutoDock Vina 1.2.x：放在 `dock/tools/vina.exe`，或把 `vina.exe` 所在目录加入 `PATH`。
- Meeko、Open Babel、AutoDockTools（由 pip / 安装脚本提供）。

**Python 包：**

```bat
python -m pip install -r requirements.txt
python -m pip install -r requirements_dock.txt
launchers\install_dock_dependencies.bat
```

`install_dock_dependencies.bat` 会安装 `requirements_dock.txt` 并安装 `dock/tools/AutoDockTools_py3`，但不会自动安装 Vina，也不会自动安装 `torch`。

`check_dock_environment.bat` 会检查 `sklearn`、`joblib`、`torch`：

- 只做对接主流程可不使用 `torch`，但该检查项会显示未安装。
- 需要 ML/DL 重打分 `ml-train` / `ml-predict` 时执行 `python -m pip install scikit-learn joblib torch`。

**R 环境：**

- 核心对接流程不需要 R。
- 虚拟敲除后的 GO/KEGG 富集、`cell-feedback` 返回单细胞对象时，需要第 1 节中的 R 4.5+ 和表达分析 R 包。

**可选 Codex skills（`evidence` / 全自动证据阶段）：**

```bat
python scripts\install_codex_skills.py
python scripts\run_docking.py check-cadd
```

证据收集使用 `%USERPROFILE%\.codex\skills\` 下的外部 skill 脚本：`uniprot-skill`、`rcsb-pdb-skill`、`chembl-skill`、`bindingdb-skill`、`pubchem-pug-skill`、`chebi-skill`、`string-skill`、`reactome-skill`、`pharmgkb-skill`、`alphafold-skill`、`opentargets-skill`。KEGG 使用公开 REST API，不需要额外 skill。

**虚拟敲除子板块：**

- `scTenifoldKnk` 引擎：需要 `requirements.txt` 中的 `scTenifoldpy>=0.3.0`。
- UMAP 命运偏转图：需要 `umap-learn`。
- GO/KEGG 富集与细胞反馈：需要第 1 节 R 环境。

## 4. 独立分子对接

**运行入口：**

- `scripts/run_molecular_docking.py`
- `launchers/run_molecular_docking.bat`
- 网页端“分子对接”页

**依赖：**

- Python 3.10+。
- 与虚拟筛选核心相同的对接依赖：`requirements_dock.txt`、AutoDockTools、AutoDock Vina。
- 不需要证据收集 Codex skills，也不需要 R。

```bat
python -m pip install -r requirements.txt
python -m pip install -r requirements_dock.txt
launchers\install_dock_dependencies.bat
```

检查：

```bat
python scripts\run_molecular_docking.py check-env
```

## 5. 分子动力学（GROMACS MD）

**运行入口：**

- `scripts/run_docking.py md-simulation ...`
- 网页端“分子动力学”页

**依赖：**

- 第 3 节虚拟筛选基础环境。
- GROMACS `gmx`：加入 `PATH`，或通过配置 `md_simulation.executable` / `--md-executable` 指定。
- 若 GROMACS 数据目录不在可执行文件旁，配置 `md_simulation.gmx_data_dir` 为包含 `share/gromacs/top` 的安装根目录。
- 配体参数化二选一：
  - ACPYPE + AmberTools：`conda install -c conda-forge gromacs ambertools`，再 `python -m pip install acpype`
  - 预生成 `<id>.itp` 与 `<id>.gro` 放入 `md_simulation.topology_dir`

可选结合自由能与构象分析：

- `md_simulation.mmpbsa_command`：配置本地 gmx_MMPBSA 或 Amber MM/PBSA 命令，支持 `{run_dir}`、`{tpr}`、`{trajectory}`、`{topology}`、`{index}` 占位符。
- 蛋白 PCA/FEL 会在 GROMACS 命令可用时自动尝试，并输出 `pca_fel.csv` 与对应指标。

## 6. 全自动集成流水线

**运行入口：**

- `scripts/run_full_pipeline.py`
- `launchers/run_full_pipeline.bat`
- `liverbio full ...`
- 网页端“全自动流水线”页

**依赖：**

全自动流水线会按阶段调用表达分析、证据收集、虚拟敲除、虚拟筛选和细胞反馈，实际启用哪些阶段取决于参数。运行完整流程前需要同时具备：

- Python 3.10+ 与 `requirements.txt`
- R 4.5+ 与表达分析 R 包
- `requirements_dock.txt`、AutoDockTools、AutoDock Vina
- 证据收集所需 Codex skills（除非使用 `--skip-evidence-fetch`）
- 细胞反馈 / GO/KEGG 富集所需 R 环境

只运行部分阶段时可减少依赖：

- `--skip-docking`：不需要 Vina 与 docking Python 依赖。
- `--skip-evidence-fetch`：不需要证据 skills，但会使用已有缓存或输出空证据。
- `--skip-cell-feedback`：不要求反馈阶段 R 脚本。

## 7. 网页版

**运行入口：**

- `web/web_ui.py`
- `launchers/run_web_ui.bat`
- `liverbio web ...`

网页服务本身使用 Python 标准库 HTTP 服务。页面能否实际启动分析，取决于该页面对应功能板块的环境；各页面底部提供“环境状态与自动补全”区域，可分别检查/补全单细胞环境和虚拟筛选环境。

## 8. 新电脑最小安装流程

在项目根目录执行：

```bat
cd /d "项目根目录"

python --version
Rscript --version

python -m pip install -r requirements.txt
launchers\install_pipeline_dependencies.bat

python -m pip install -r requirements_dock.txt
launchers\install_dock_dependencies.bat

REM AutoDock Vina 放到 dock\tools\vina.exe，或把 vina.exe 目录加入 PATH

python scripts\install_codex_skills.py

liverbio doctor
```

只使用表达分析时，不需要执行 `requirements_dock.txt`、`install_dock_dependencies.bat` 和 Vina 安装；只使用独立分子对接时，不需要执行 R 包安装和 skills 安装。最终以各板块的 `check_*.bat` 输出为准。
