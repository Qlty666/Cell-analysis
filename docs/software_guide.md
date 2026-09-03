# 肝癌生信软件使用指南

本目录不仅提供命令行脚本，也可以通过一个统一的入口使用。`liverbio` 只负责转发参数，不会改变任何原有脚本的实际行为；原有 `launchers`、网页端和 `scripts/run_*.py` 仍可直接使用。

## 统一入口

Windows 下在项目根目录打开命令提示符或 PowerShell：

```text
liverbio help
liverbio version
```

`liverbio.bat` 位于项目根目录。若希望在任何目录直接使用，可把项目根目录加入 `PATH`。

主要命令：

```text
liverbio expression GSE125449 --output ../liver_cancer --species auto
liverbio full --accession GSE125449 --output ../liver_cancer --workdir ../liver_cancer_full
liverbio docking pipeline --config config/docking_config.json
liverbio datasets --disease "liver cancer" --max-results 20
liverbio web --page full
liverbio doctor
```

每条命令后面的参数与原脚本完全一致。对某个功能需要查看详细参数时，先查看原脚本帮助：

```text
python scripts\run_pipeline.py --help
python scripts\run_docking.py <子命令> --help
python scripts\run_full_pipeline.py --help
python scripts\search_datasets.py --help
```

## 常用工作流

### 打开网页软件

```text
liverbio web
```

默认打开 `http://127.0.0.1:8000/full`。可使用 `liverbio web --page dock`、`--page datasets`、`--page results`、`--page tasks` 等打开对应页面。

### 检查运行环境

```text
liverbio doctor
liverbio doctor pipeline
liverbio doctor docking
```

`liverbio doctor` 按顺序执行表达分析和虚拟筛选环境检查；缺少组件时按输出提示运行 `launchers` 下的安装脚本。

## Codex Skills

`skills/` 目录保存按功能拆分的 Codex skill 定义：

- `liver-expression-analysis`：表达分析。
- `liver-virtual-screening`：虚拟筛选。
- `liver-full-pipeline`：全自动集成流水线。
- `liver-dataset-search`：数据集搜索与下载。

安装到当前用户的 Codex skill 目录：

```text
python scripts\install_codex_skills.py
python scripts\install_codex_skills.py --force
```

Skill 只是指导 Codex 调用本项目现有脚本，不复制或替代核心分析代码。

## 目录结构

```text
scripts\         命令行入口
src\             可复用 Python 模块
  common\        环境与通用工具
  data\          数据下载、转换与校验
  pipeline\      表达分析、全流程编排
  docking\       虚拟筛选、靶点评分
  report\        结果报告
  liverbio_suite\ 统一命令入口
web\             本地网页软件
skills\          Codex skill 源文件
launchers\       Windows 快捷启动与安装脚本
config\          默认配置
tests\           单元与集成测试
```

优化原则：优先复用现有已测试模块；不通过复制代码方式“打包”，而是用统一入口、环境检查和文档把现有功能组织成软件。
