# 项目代码结构说明

本文把当前代码按“入口层、实现层、资源与配置层、测试与文档层”整理，不改变任何现有脚本的调用方式。

## 1. 总体结构

```text
Script/
├── liverbio.bat                     # 统一 CLI 根入口
├── setup_new_computer.bat           # 新电脑一键安装（Windows）
├── check_new_computer.bat           # 新电脑环境检查（Windows）
├── setup_new_computer.sh            # 新电脑一键安装（Linux/macOS）
├── check_new_computer.sh            # 新电脑环境检查（Linux/macOS）
├── package_for_new_computer.bat     # 生成干净源码包（Windows）
├── package_for_new_computer.sh      # 生成干净源码包（Linux/macOS）
├── NEW_COMPUTER_SETUP.md            # 新电脑部署说明
├── README.md                        # 总说明、使用方法和更新日志
├── AGENTS.md                        # Codex 项目执行规则
├── VIRTUAL_SCREENING_REQUIREMENTS.md
├── requirements.txt / requirements_dock.txt
├── environment_dock.yml
├── pytest.ini                       # pytest 配置（testpaths = tests）
├── .github/workflows/tests.yml      # GitHub Actions 测试工作流
├── scripts/                         # Python 命令行入口与验证脚本
├── launchers/                       # Windows 快捷启动、环境检查与安装（_common.bat 为公共 Python 探测）
├── src/                             # 可复用实现模块
│   ├── analysis/                    # 表达分析 R/Python 实现（analysis_pipeline.R 驱动 + R/ 模块）
│   ├── common/                      # 通用环境探测、HTTP 与 HTML 工具
│   ├── data/                        # 数据集下载、转换与合成数据
│   ├── docking/                     # 虚拟筛选、虚拟敲除、MD 等 CADD 实现
│   ├── molecular_docking/           # 独立分子对接板块
│   ├── liverbio_suite/              # liverbio 统一入口实现
│   ├── pipeline/                    # 流水线编排、阶段路径、集成报告与细胞反馈
│   └── report/                      # HTML/Word 报告生成与结果图指南
├── web/                             # 本地网页服务、请求处理、模板与静态资源
├── skills/                          # Codex skill 源文件
├── config/                          # 默认 JSON 配置
├── docs/                            # 使用与结构文档
├── tests/                           # 单元与集成测试
├── data_cache/                      # 运行时下载缓存（gitignore）
├── dock/                            # 虚拟筛选运行产物（gitignore）
└── results/                         # 表达分析结果（gitignore）
```

## 2. 入口层

### 2.1 根入口

| 文件 | 作用 |
| --- | --- |
| `liverbio.bat` | 根目录统一入口，把参数转发给 `scripts/liverbio.py` |
| `setup_new_computer.bat` / `check_new_computer.bat` | 新电脑完整安装与环境检查 |
| `package_for_new_computer.bat` | 生成可拷贝的干净源码 zip |

### 2.2 scripts 运行入口

| 文件 | 作用 |
| --- | --- |
| `scripts/liverbio.py` | 统一 CLI，转发到各运行/验证入口 |
| `scripts/run_pipeline.py` | 表达分析 CLI |
| `scripts/run_docking.py` | 虚拟筛选 CLI |
| `scripts/run_full_pipeline.py` | 全自动集成流水线 CLI |
| `scripts/run_molecular_docking.py` | 独立分子对接 CLI |
| `scripts/search_datasets.py` | 多数据库数据集搜索与下载 |
| `scripts/dataset_search_ml.py` | 搜索结果 ML/DL 相关性排序 |
| `scripts/run_web_full_new_datasets.py` | 批量提交真实数据集全流程 |
| `scripts/install_codex_skills.py` | 安装 `skills/` 下的 Codex skill |

### 2.3 scripts 验证入口

| 文件 | 验证范围 |
| --- | --- |
| `scripts/validate_pipeline.py` | 合成数据表达流水线 |
| `scripts/validate_dock_pipeline.py` | 假 Vina 对接流水线 |
| `scripts/validate_real_pipeline.py` | 真实 GEO 表达流水线 |
| `scripts/validate_real_evidence.py` | 真实 PDB 证据收集 |
| `scripts/validate_real_random.py` | 随机真实数据证据收集与对接盒 |
| `scripts/validate_new_features.py` | 真实数据靶点评分与验证方案 |
| `scripts/validate_random_real_full_pipeline.py` | 随机真实 GSE 全流程 |
| `scripts/validate_dataset_search.py` | 多数据库数据集搜索召回 |

### 2.4 launchers

| 分组 | 文件 | 作用 |
| --- | --- | --- |
| 运行 | `run_web_ui.bat` | 启动网页端 |
| 运行 | `run_docking.bat` | 虚拟筛选快捷入口 |
| 运行 | `run_molecular_docking.bat` | 独立分子对接快捷入口 |
| 运行 | `run_full_pipeline.bat` | 全自动流水线快捷入口 |
| 运行 | `run_GSE125449.bat` | GSE125449 快捷入口 |
| 运行 | `run_pipeline_prompt.bat` | 交互式表达分析入口 |
| 环境 | `check_dock_environment.bat/.py` | 对接环境检查 |
| 环境 | `check_pipeline_environment.bat/.py` | 流水线环境检查（`.py` 支持 `--offline` 跳过 NCBI GEO 网络探测） |
| 环境 | `install_environment.py` | 按功能板块安装/检查环境的统一入口 |
| 环境 | `check_*_environment.bat` | 各功能板块环境检查快捷入口 |
| 安装 | `install_dock_dependencies.bat/.py` | 对接依赖安装 |
| 安装 | `install_pipeline_dependencies.bat/.py` | 流水线依赖安装 |
| 安装 | `install_*_environment.bat` / `setup_environment.bat` | 各功能板块一键补全与选择式安装 |
| 公共 | `_common.bat` | 各 `.bat` 共用的 Python 3 探测（设置 `PYTHON_EXE`） |
| 打包 | `package_portable.py` | 用 git 跟踪文件生成干净源码包 |

## 3. 实现层

### 3.1 `src/analysis`

表达分析实现，R 为主、Python 为辅：

| 文件 | 作用 |
| --- | --- |
| `analysis_pipeline.R` | 表达分析主流程 |
| `cellchat_analysis.R` | CellChat 细胞通讯 |
| `install_deps.R` | R 依赖安装 |
| `ml_analysis.py` | 表达谱 ML 分类与可解释性 |

### 3.2 `src/common`

| 文件 | 作用 |
| --- | --- |
| `env.py` | Rscript/工具路径与环境探测 |
| `http.py` | 通用 HTTP 请求封装（重试、退避、大小限制与明确错误） |
| `html_utils.py` | 报告生成共用的 HTML 转义工具 |

### 3.3 `src/data`

| 文件 | 作用 |
| --- | --- |
| `download_data.py` | 通用数据下载 |
| `geo_downloader.py` | GEO 系列下载与 manifest 分类 |
| `biostudies_downloader.py` | ArrayExpress/BioStudies/Expression Atlas 下载 |
| `h5_converter.py` | h5ad/loom 转 10x MTX |
| `validation_generator.py` | 合成验证数据生成 |

### 3.4 `src/docking`

虚拟筛选与靶点分析核心：

| 文件 | 作用 |
| --- | --- |
| `cli.py` | 虚拟筛选 CLI 路由 |
| `config.py` | 配置加载与默认值 |
| `docking.py` | AutoDock Vina 对接主逻辑 |
| `analysis.py` / `box.py` / `receptor.py` / `ligands.py` / `redock.py` | 对接分析、盒检测、受体/配体处理 |
| `evidence.py` | 多数据库靶点证据收集 |
| `knockout.py` / `insilico.py` / `export_single_cell_insilico.R` / `insilico_enrichment.R` | 虚拟敲除与单细胞 GRN 模拟 |
| `md_simulation.py` | GROMACS 分子动力学模拟 |
| `network_toxicology.py` / `cytoscape_network.py` / `signal_detection.py` | 网络毒理学、Cytoscape XGMML/推送与 FAERS 信号 |
| `pipeline.py` / `provenance.py` / `validation.py` | 对接流水线、溯源与校验 |
| `report.py` | 对接/验证报告 |
| `ml.py` / `handoff.py` / `environment.py` / `utils.py` | 重打分、工具交接、环境与通用工具 |
| `http.py` / `html_utils.py` | 证据收集 HTTP 下载（超时与重试）与报告 HTML 转义工具 |

### 3.5 `src/molecular_docking`

独立于虚拟筛选的分子对接板块：

| 文件 | 作用 |
| --- | --- |
| `cli.py` | 独立对接 CLI |
| `config.py` | 独立配置 |
| `pipeline.py` | 独立流水线 |
| `report.py` | HTML 报告 |

### 3.6 `src/liverbio_suite`

| 文件 | 作用 |
| --- | --- |
| `cli.py` | `liverbio` 子命令路由 |

### 3.7 `src/pipeline`

| 文件 | 作用 |
| --- | --- |
| `orchestrator.py` | 表达流水线进程编排 |
| `integration.py` | 全自动集成流水线编排 |
| `stage_paths.py` | 全自动流水线阶段目录、完成标记与输出清单 |
| `integrated_report.py` | 集成 HTML 报告生成 |
| `cell_feedback.py` / `cell_feedback.R` | 结果写回单细胞对象的反馈分析 |
| `export_pseudobulk.R` | 伪 bulk 导出 |

### 3.8 `src/report`

| 文件 | 作用 |
| --- | --- |
| `generate_report.py` | 图/数据联合分析与总报告 |
| `export_report.py` | DOCX/PDF 等导出 |
| `guides.py` | 结果图/数据指南表（从 `generate_report.py` 抽出的纯数据） |

## 4. 网页层

| 路径 | 作用 |
| --- | --- |
| `web/web_ui.py` | 本地网页服务、路由注册与任务调度 |
| `web/web_handler.py` | HTTP 请求处理器（页面、运行接口、文件下载与安全校验） |
| `web/web_data.py` | 静态数据表（导航、环境模块、图目录、软件清单与阶段标签） |
| `web/web_state.py` | 进程内运行时状态（任务注册表、队列、锁与历史读写） |
| `web/web_results.py` | 任务状态、历史记录与结果清单读取辅助 |
| `web/templates/*.html` | 全流程、表达分析、数据集、虚拟筛选、分子对接、结果清单、任务进度等页面模板 |
| `web/static/app.css` | 共享样式 |
| `web/static/nav.js` | 导航与任务徽标脚本 |
| `web/static/result_details.json` | 结果清单/图指南数据 |

## 5. 配置、技能与测试

| 路径 | 作用 |
| --- | --- |
| `config/*.json` | 表达分析、虚拟筛选、独立分子对接与项目默认配置 |
| `skills/liver-*/SKILL.md` | Codex skill 定义 |
| `skills/liver-*/agents/openai.yaml` | skill agent 配置 |
| `pytest.ini` | pytest 配置（`testpaths = tests`） |
| `tests/conftest.py` | pytest 共享 `sys.path` 引导（仓库根、`src`、`web`、`launchers`、`scripts`） |
| `tests/test_*.py` | 与 `src`/`scripts`/`web` 对应的单元与集成测试 |
| `tests/test_r_pipeline_syntax.py` | R 脚本解析、驱动加载与 `R/` 模块独立加载冒烟测试 |
| `tests/test_web_security.py` | 网页端同源/跨站校验、请求体限制与路径安全测试 |
| `tests/test_web_server_smoke.py` | 在回环端口启动真实 `ThreadingHTTPServer` 的路由与安全冒烟测试 |
| `tests/test_web_script_mode.py` | `python web/web_ui.py` 脚本模式下状态不被重复导入的回归测试 |
| `.github/workflows/tests.yml` | GitHub Actions：push/PR 在 Windows 与 Linux + Python 3.11 上安装依赖并运行 pytest |

## 6. 代码组织原则

- `scripts/` 只放命令行入口与验证驱动，具体分析逻辑放到 `src/`。
- `launchers/` 中的 `.bat` 负责双击启动，`launchers/*.py` 负责同名的环境检查/安装实现。
- `src/` 按领域分包：表达分析、数据下载、虚拟筛选/CADD、独立分子对接、全流程编排、报告生成。
- `web/` 通过调用 `scripts/` 与 `src/` 复用分析能力，不复制核心算法。
- `skills/` 只描述调用方式，不复制分析代码。
- `tests/` 与实现模块同名对应，运行 `python -m pytest` 可整仓验证。

更细的功能说明见 `README.md` 与 `docs/software_guide.md`。
