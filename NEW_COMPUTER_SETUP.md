# 新电脑快速部署

本目录是源码项目，不携带本机结果、缓存或下载工具。克隆仓库或解压
`portable/` 下生成的干净源码包后，在新电脑上按下面步骤部署。

## Windows 最快路径

1. 安装或确认已有 Python 3.10+（推荐 3.11）。
2. 双击或命令行运行根目录的 `setup_new_computer.bat`：

```bat
setup_new_computer.bat
```

该命令等价于：

```bat
python launchers\install_environment.py install full --with-ml
```

脚本会按需安装：

- `requirements.txt` 与 `requirements_dock.txt` 中的 Python 依赖；
- 表达分析 R 包；Windows 上找不到 R 时会从 CRAN 自动安装；
- AutoDockTools 与 AutoDock Vina 到 `dock/tools/`；
- 项目 Codex skills；
- `joblib` 与 `torch`（ML/DL 可选依赖）。

3. 检查环境：

```bat
check_new_computer.bat
```

4. 使用统一入口：

```bat
liverbio help
liverbio version
liverbio web
```

## Linux / macOS

R 4.5+、GROMACS、AutoDock Vina 等外部程序需要先按系统方式安装；Python 部分可运行：

```bash
python3 launchers/install_environment.py install full --with-ml
python3 launchers/install_environment.py check full
./setup_new_computer.sh
```

## 只安装某一板块

不想安装完整环境时，可只补当前使用的板块：

```bash
python launchers/install_environment.py list
python launchers/install_environment.py install expression
python launchers/install_environment.py install docking --with-ml
python launchers/install_environment.py check full
```

网页版“环境补全”页也提供相同的按板块安装/检查入口。

## 生成干净源码包

在保留 git 仓库的电脑上生成不含本机运行数据、缓存、日志和下载工具的 zip：

```bat
liverbio package
```

等价命令：

```bash
python launchers/package_portable.py
```

输出位于 `portable/Cell-analysis-portable_v<版本号>.zip`。把该 zip 拷贝到新电脑，
解压后从“Windows 最快路径”继续即可。

## 输出目录

分析结果默认写到源码目录之外，避免污染可重新部署的源码：

- 表达分析/全流程常用 `--output ../liver_cancer`，可改成任意路径；
- 全自动流水线用 `--workdir <目录>`；
- 独立对接默认在 `molecular_docking/`，可用 `--workdir` 改变；
- 网页端验证任务结果默认写入 `data_cache/validation_runs/`，可通过
  `LIVER_VALIDATION_ROOT` 环境变量修改。

## 仍需手动安装的可选软件

- Cytoscape 桌面版：网络毒理学自动推送需要；启动后执行
  `cytoscape.bat -R 1234`；
- GROMACS：`md-simulation auto` 模式需要；
- ACPYPE/AmberTools：MD 配体拓扑自动生成需要，也可预先放入
  `md_simulation.topology_dir`；
- 外部 Codex skills：证据收集数据库 skill 需要安装在当前用户的 Codex skill
  目录中，项目功能 skills 可由安装脚本自动复制。
