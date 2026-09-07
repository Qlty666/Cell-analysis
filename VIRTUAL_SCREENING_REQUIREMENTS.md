# 虚拟筛选环境清单

整个 `Script` 目录满足以下环境后，可以整体复制到新电脑直接使用。
新电脑一键安装与干净源码包说明见 `NEW_COMPUTER_SETUP.md`。

## Python 核心

- Python 3.10+（推荐 3.11）
- numpy
- pandas
- pyyaml
- matplotlib
- openpyxl

安装：

```bash
python -m pip install numpy pandas pyyaml matplotlib openpyxl
```

## 配体/结构准备

- RDKit：分子标准化、3D 构象、描述符和指纹
- Meeko：`mk_prepare_ligand.py` / `mk_prepare_receptor.py`
- gemmi：Meeko 依赖
- openbabel-wheel：Open Babel `obabel`
- versioneer：AutoDockTools_py3 构建依赖
- AutoDockTools_py3：`prepare_receptor4.py`

安装：

```bash
python -m pip install rdkit meeko gemmi openbabel-wheel versioneer
python -m pip install ./dock/tools/AutoDockTools_py3
```

## 对接引擎

- AutoDock Vina 1.2.x：`dock/tools/vina.exe` 或加入 PATH
- 可选：UniDock-Pro、HDOCK、HADDOCK（外部对接模板）
- 可选：Amber/AmberTools、GROMACS（对接后 MD 模拟）

## 分子动力学模拟

`md-simulation` 命令会把 Top 对接姿态转成 GROMACS 蛋白-配体复合物，并在 `auto`
模式下执行能量最小化、NVT/NPT 平衡与生产模拟。运行前需要：

- GROMACS `gmx`（加入 PATH，或通过 `md_simulation.executable` 指定）
- 如果 GROMACS 数据目录不在可执行文件旁，配置 `md_simulation.gmx_data_dir`
  为包含 `share/gromacs/top` 的安装根目录
- 配体参数化二选一：
  - ACPYPE + AmberTools：`acpype -a gaff2 -c bcc` 自动生成 GAFF2 拓扑
  - 把预先生成好的 `<id>.itp` 与 `<id>.gro` 放入
    `md_simulation.topology_dir`

Linux/WSL/conda 环境示例：

```bash
conda install -c conda-forge gromacs ambertools
python -m pip install acpype
```

## ML/DL 重打分

- scikit-learn：随机森林、GBDT、MLP
- joblib：模型保存
- torch：PyTorch MLP 深度学习重打分

安装：

```bash
python -m pip install scikit-learn joblib torch
```

## Cytoscape 网络可视化（网络毒理学可选）

- Cytoscape 3.9+ 桌面版（Cytoscape 3.10 已预装 CyREST）
- `py4cytoscape>=1.13.0`（已写入 `requirements_dock.txt`）

启动 Cytoscape 并启用 REST 服务后，网络毒理学流程会在 `cytoscape: auto`
模式下自动推送 C-T-P-D 网络、应用样式并导出 PNG；未启动服务时仍会输出
可直接导入的 `ctpd_network.xgmml`。

```bash
# Cytoscape 安装目录示例（Windows）
"<CYTOSCAPE_INSTALL_DIR>\cytoscape.bat" -R 1234
python -m pip install py4cytoscape
```

## 数据库证据收集（Codex skills）

需要在 `%USERPROFILE%\.codex\skills\` 下存在：

- uniprot-skill
- rcsb-pdb-skill
- chembl-skill
- bindingdb-skill
- pubchem-pug-skill
- chebi-skill
- string-skill
- reactome-skill
- pharmgkb-skill
- alphafold-skill
- opentargets-skill

`kegg` 使用公开 REST API（`rest.kegg.jp`），不需要额外 skill。

## 单细胞分析（原有系统）

- R 4.5+
- Seurat、scDblFinder、SingleCellExperiment、clusterProfiler、DESeq2 等 R 包

详情见 `README.md` 和 `launchers/install_pipeline_dependencies.py`。

- 细胞反馈阶段还需要 Seurat 对象（`results/checkpoints/seurat_annotated.rds` 或 `results/data/*.rds`）。

## 一键自动补全

虚拟筛选：

```text
launchers\install_docking_environment.bat
launchers\install_dock_dependencies.bat
```

`install_docking_environment.bat` 会先安装基础 Python 依赖，再自动补全
AutoDockTools；若 `dock/tools/vina.exe` 缺失，会自动下载 Windows 版
AutoDock Vina 到该目录。需要 ML/DL 重打分时使用
`launchers\install_docking_environment.bat --with-ml`。

网页端：虚拟筛选页 → 软件环境 → 自动补全环境。
