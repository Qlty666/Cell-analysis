# 虚拟筛选环境清单

整个 `Script` 目录满足以下环境后，可以整体复制到新电脑直接使用。

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
- 可选：Amber/AmberTools、GROMACS（对接后 MD 交接）

## ML/DL 重打分

- scikit-learn：随机森林、GBDT、MLP
- joblib：模型保存
- torch：PyTorch MLP 深度学习重打分

安装：

```bash
python -m pip install scikit-learn joblib torch
```

## 数据库证据收集（Codex skills）

需要在 `%USERPROFILE%\.codex\skills\` 下存在：

- uniprot-skill
- rcsb-pdb-skill
- chembl-skill
- bindingdb-skill
- pubchem-pug-skill
- chebi-skill

## 单细胞分析（原有系统）

- R 4.5+
- Seurat、scDblFinder、SingleCellExperiment、clusterProfiler、DESeq2 等 R 包

详情见 `README.md` 和 `launchers/install_pipeline_dependencies.py`。

- 细胞反馈阶段还需要 Seurat 对象（`results/checkpoints/seurat_annotated.rds` 或 `results/data/*.rds`）。

## 一键自动补全

虚拟筛选：

```text
launchers\install_dock_dependencies.bat
```

网页端：虚拟筛选页 → 软件环境 → 自动补全环境。
