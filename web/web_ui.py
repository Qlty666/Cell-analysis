#!/usr/bin/env python3
"""Local HTML web UI for the single-cell pipeline."""

import argparse
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlparse

WEB_DIR = Path(__file__).resolve().parent
APP_ROOT = WEB_DIR.parent
SCRIPTS_DIR = APP_ROOT / "scripts"
JOBS = {}
QUEUE = []
QUEUE_LOCK = threading.Lock()
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"
INDEX_PATH = TEMPLATE_DIR / "index.html"
PAGE_TEMPLATE_PATH = TEMPLATE_DIR / "web_page_template.html"
HISTORY_PATH = WEB_DIR / "history.json"
INSTALL_LOG = WEB_DIR / "install_log.txt"
INSTALL_JOB = {}
FINISHED_NOTIFICATIONS: list[dict] = []
NOTIFY_LOCK = threading.Lock()
TASK_HISTORY_PATH = WEB_DIR / "task_history.json"
TASK_HISTORY_LOCK = threading.Lock()

HEARTBEAT_CLIENTS: dict[str, float] = {}
HEARTBEAT_LAST_SEEN_AT: float | None = None
HEARTBEAT_LOCK = threading.Lock()
HEARTBEAT_INTERVAL_SECONDS = 5
HEARTBEAT_IDLE_TIMEOUT_SECONDS = 120
HEARTBEAT_START_GRACE_SECONDS = 20
HEARTBEAT_SHUTDOWN_GRACE_SECONDS = 5

SRC_DIR = APP_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

DOCK_TEMPLATE_PATH = TEMPLATE_DIR / "dock_page_template.html"
DOCK_JOBS = {}
DOCK_QUEUE = []
DOCK_QUEUE_LOCK = threading.Lock()
DOCK_HISTORY_PATH = WEB_DIR / "dock_history.json"
FULL_TEMPLATE_PATH = TEMPLATE_DIR / "full_page_template.html"
RESULTS_TEMPLATE_PATH = TEMPLATE_DIR / "results_manifest_optimized.html"
RESULT_GUIDE_PATH = APP_ROOT / "docs" / "result_figure_guide.md"
RESULT_DETAILS_PATH = STATIC_DIR / "result_details.json"
TASKS_TEMPLATE_PATH = TEMPLATE_DIR / "tasks_template.html"
DATASET_TEMPLATE_PATH = TEMPLATE_DIR / "datasets_template.html"
DATASET_SEARCH_DIR = APP_ROOT / "data_cache" / "dataset_search"
DATASET_DATABASES = ("geo", "biostudies", "atlas")
DATASET_DOWNLOAD_JOBS = {}
DATASET_DOWNLOAD_LOCK = threading.Lock()
FULL_JOBS = {}
FULL_QUEUE = []
FULL_QUEUE_LOCK = threading.Lock()
VALIDATION_REPORT_DIR = APP_ROOT.parent / "y3" / "validation"
VALIDATION_REPORT_PATH = VALIDATION_REPORT_DIR / "validation_summary.json"
VALIDATION_LOG = WEB_DIR / "validation_run.log"
VALIDATION_JOB = {"proc": None, "log": None, "handle": None, "started": None}
NAV_HTML = (
    '<div class="topnav">'
    '<a href="/full">全自动流水线</a>'
    '<a href="/">表达分析</a>'
    '<a href="/datasets">数据集搜索</a>'
    '<a href="/dock">虚拟筛选</a>'
    '<a href="/results">结果清单</a>'
    '<a href="/tasks" class="nav-right">任务进度</a>'
    '</div>'
)
NAV_CSS = (
    ".topnav{position:sticky;top:0;z-index:100;"
    "background:#0f172a;padding:12px 28px;box-shadow:0 2px 8px rgba(15,23,42,.35);"
    "display:flex;gap:18px;align-items:center;flex-wrap:wrap;}"
    ".topnav a{color:#ffffff;text-decoration:none;font-size:15px;"
    "font-weight:600;padding:6px 10px;border-radius:6px;"
    "background:rgba(255,255,255,.08);}"
    ".topnav a:hover,.topnav a.active{background:#1665c0;color:#fff;}"
    ".topnav .nav-right{margin-left:auto;}"
)

HEARTBEAT_SCRIPT = """
<script>
(function () {
  var clientId = "";
  try {
    clientId = crypto.randomUUID ? crypto.randomUUID() : "c" + Date.now().toString(36) + Math.random().toString(36).slice(2);
  } catch (e) {
    clientId = "c" + Date.now().toString(36) + Math.random().toString(36).slice(2);
  }
  var heartbeatUrl = "/heartbeat?client=" + encodeURIComponent(clientId);
  function sendHeartbeat(left) {
    var url = left ? heartbeatUrl + "&left=1" : heartbeatUrl;
    try {
      if (navigator.sendBeacon) {
        navigator.sendBeacon(url);
      } else {
        fetch(url, {method: "GET", keepalive: true, credentials: "same-origin"}).catch(function () {});
      }
    } catch (e) {}
  }
  sendHeartbeat(false);
  window.setInterval(function () { sendHeartbeat(false); }, __HEARTBEAT_MS__);
  window.addEventListener("pagehide", function () { sendHeartbeat(true); });
  document.addEventListener("visibilitychange", function () {
    if (!document.hidden) sendHeartbeat(false);
  });
})();
</script>
""".replace("__HEARTBEAT_MS__", str(HEARTBEAT_INTERVAL_SECONDS * 1000)).encode("utf-8")


def register_heartbeat(client_id: str, now: float | None = None) -> None:
    global HEARTBEAT_LAST_SEEN_AT
    client_id = (client_id or "").strip()
    if not client_id or len(client_id) > 128:
        return
    with HEARTBEAT_LOCK:
        seen_at = time.monotonic() if now is None else now
        HEARTBEAT_CLIENTS[client_id] = seen_at
        HEARTBEAT_LAST_SEEN_AT = seen_at


def unregister_heartbeat(client_id: str) -> None:
    client_id = (client_id or "").strip()
    if not client_id:
        return
    with HEARTBEAT_LOCK:
        HEARTBEAT_CLIENTS.pop(client_id, None)


def _heartbeat_client_ids() -> set[str]:
    with HEARTBEAT_LOCK:
        return set(HEARTBEAT_CLIENTS)


def _heartbeat_last_seen_at() -> float | None:
    with HEARTBEAT_LOCK:
        return HEARTBEAT_LAST_SEEN_AT


def _purge_stale_heartbeats(
    now: float,
    timeout: float = HEARTBEAT_IDLE_TIMEOUT_SECONDS,
) -> None:
    cutoff = now - timeout
    with HEARTBEAT_LOCK:
        for client_id, seen in list(HEARTBEAT_CLIENTS.items()):
            if seen < cutoff:
                HEARTBEAT_CLIENTS.pop(client_id, None)


def _inject_heartbeat_script(body: bytes) -> bytes:
    if b"</body>" not in body:
        return body
    return body.replace(b"</body>", HEARTBEAT_SCRIPT + b"</body>", 1)


def _run_idle_shutdown_monitor(
    server: ThreadingHTTPServer,
    started_at: float,
) -> None:
    idle_since: float | None = None
    while True:
        now = time.monotonic()
        _purge_stale_heartbeats(now)
        if _heartbeat_client_ids():
            idle_since = None
        else:
            if idle_since is None:
                idle_since = now
            if (
                _heartbeat_last_seen_at() is not None
                or now - started_at >= HEARTBEAT_START_GRACE_SECONDS
            ) and now - idle_since >= HEARTBEAT_SHUTDOWN_GRACE_SECONDS:
                print("No web page connected; stopping web UI")
                server.shutdown()
                return
        time.sleep(1.0)


def _port_is_listening(host: str, port: int) -> bool:
    target = "127.0.0.1" if host in ("0.0.0.0", "::") else host
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(0.3)
            return sock.connect_ex((target, port)) == 0
    except OSError:
        return False


def _stop_stale_web_ui(host: str, port: int) -> bool:
    if os.name != "nt":
        return True
    current_pid = os.getpid()
    script = f"""
$ErrorActionPreference = 'SilentlyContinue'
$conns = Get-NetTCPConnection -LocalPort {port} -State Listen
$ids = @($conns | Select-Object -ExpandProperty OwningProcess -Unique)
if ($ids.Count -eq 0) {{ exit 0 }}
$procs = Get-CimInstance Win32_Process | Where-Object {{
  $_.ProcessId -in $ids -and
  $_.CommandLine -like '*web_ui.py*' -and
  $_.ProcessId -ne {current_pid}
}}
foreach ($proc in $procs) {{
  Stop-Process -Id $proc.ProcessId -Force
}}
"""
    try:
        subprocess.run(
            ["powershell", "-NoProfile", "-NonInteractive", "-Command", script],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=20,
        )
    except Exception:
        return False
    deadline = time.monotonic() + 6.0
    while time.monotonic() < deadline:
        if not _port_is_listening(host, port):
            return True
        time.sleep(0.2)
    return False


def _cleanup_stale_web_ui(host: str, port: int) -> bool:
    if not _port_is_listening(host, port):
        return True
    print(f"Port {port} is in use; checking for stale web UI process...")
    return _stop_stale_web_ui(host, port)


FIGURES = [
    {"file": "fig_01_qc_raw_violin.png", "label": "QC 小提琴图（原始）"},
    {"file": "fig_01_qc_filtered_violin.png", "label": "QC 小提琴图（过滤后）"},
    {"file": "fig_48_qc_pvalue_comparison.png", "label": "QC 质控差异度 P 值"},
    {"file": "fig_02_doublet_scores.png", "label": "双细胞得分图"},
    {"file": "fig_03_umap_clusters.png", "label": "UMAP 聚类图"},
    {"file": "fig_04_umap_condition.png", "label": "UMAP 分组图"},
    {"file": "fig_05_umap_annotation.png", "label": "UMAP 注释图"},
    {"file": "fig_06_dotplot_markers.png", "label": "Marker 基因 DotPlot"},
    {"file": "fig_07_annotation_confusion_heatmap.png", "label": "注释混淆矩阵热图"},
    {"file": "fig_08_volcano.png", "label": "差异表达图", "styles": ["volcano", "maplot"]},
    {"file": "fig_09_deg_heatmap.png", "label": "Top DEG 热图"},
    {"file": "fig_09_deg_horizontal_violin.png", "label": "Top DEG 横向小提琴图（P 值）"},
    {"file": "fig_10_go_up.png", "label": "GO BP 富集图（上调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_11_go_down.png", "label": "GO BP 富集图（下调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_12_kegg_up.png", "label": "KEGG 富集图（上调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_13_kegg_down.png", "label": "KEGG 富集图（下调）", "styles": ["dotplot", "barplot", "cnetplot"]},
    {"file": "fig_14_pca.png", "label": "PCA 分组图"},
    {"file": "fig_15_elbow.png", "label": "主成分 Elbow 图"},
    {"file": "fig_16_featureplot_markers.png", "label": "Marker 基因 FeaturePlot"},
    {"file": "fig_17_marker_violin.png", "label": "Marker 基因小提琴图"},
    {"file": "fig_18_celltype_proportion.png", "label": "细胞类型比例堆叠图"},
    {"file": "fig_19_condition_proportion.png", "label": "分组构成比例图"},
    {"file": "fig_20_gsea_go.png", "label": "GSEA GO BP 富集图", "styles": ["ridgeplot", "gseaplot2"]},
    {"file": "fig_21_gsea_kegg.png", "label": "GSEA KEGG 富集图", "styles": ["ridgeplot", "gseaplot2"]},
    {"file": "fig_22_go_network.png", "label": "GO BP 通路网络图（筛选后 Top5）", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_23_kegg_network.png", "label": "KEGG 通路网络图（筛选后 Top5）", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_46_go_top5.png", "label": "GO BP 筛选后 Top5", "styles": ["dotplot", "barplot", "cnetplot", "emapplot"]},
    {"file": "fig_47_kegg_top5.png", "label": "KEGG 筛选后 Top5", "styles": ["dotplot", "barplot", "cnetplot", "emapplot"]},
    {"file": "fig_24_ml_feature_importance.png", "label": "ML 特征重要性图"},
    {"file": "fig_25_ml_shap.png", "label": "SHAP 可解释性图"},
    {"file": "fig_26_cellcycle_umap.png", "label": "细胞周期 UMAP"},
    {"file": "fig_27_cellcycle_proportion.png", "label": "细胞周期比例图"},
    {"file": "fig_28_umap_sample.png", "label": "UMAP 按样本"},
    {"file": "fig_29_doublet_rate_sample.png", "label": "样本双细胞率图"},
    {"file": "fig_30_sample_proportion.png", "label": "样本细胞类型比例图"},
    {"file": "fig_31_cluster_marker_heatmap.png", "label": "聚类 Marker 热图"},
    {"file": "fig_32_cluster_marker_dotplot.png", "label": "聚类 Marker DotPlot"},
    {"file": "fig_33_signature_scores_umap.png", "label": "功能签名 UMAP"},
    {"file": "fig_34_signature_scores_boxplot.png", "label": "功能签名箱线图"},
    {"file": "fig_35_celltype_abundance_effect.png", "label": "细胞类型丰度变化图"},
    {"file": "fig_36_cnv_heatmap.png", "label": "推断 CNV 热图"},
    {"file": "fig_37_singler_umap.png", "label": "SingleR 注释 UMAP"},
    {"file": "fig_38_singler_confusion_heatmap.png", "label": "SingleR 混淆矩阵热图"},
    {"file": "fig_39_trajectory_umap.png", "label": "拟时序轨迹图"},
    {"file": "fig_40_cellchat_network.png", "label": "CellChat 通讯网络图"},
    {"file": "fig_41_cellchat_heatmap.png", "label": "CellChat 通讯热图"},
    {"file": "fig_42_cellchat_bubble.png", "label": "CellChat 配体受体气泡图"},
    {"file": "fig_43_ml_confusion_matrix.png", "label": "ML 混淆矩阵"},
    {"file": "fig_44_ml_roc_pr.png", "label": "ML ROC 与 PR 曲线"},
    {"file": "fig_45_ml_cv_scores.png", "label": "ML 交叉验证得分图"},
    {"file": "fig_45_ml_calibration_curve.png", "label": "ML 校准曲线"},
]
FIGURE_NAMES = [item["file"] for item in FIGURES]

STYLE_LABELS = {
    "volcano": "火山图",
    "maplot": "MA 图",
    "dotplot": "气泡图",
    "barplot": "柱状图",
    "cnetplot": "通路网络图",
    "ridgeplot": "峰峦图",
    "gseaplot2": "GSEA 富集曲线",
    "emapplot": "富集关系网络图",
}

SOFTWARE = [
    {"name": "R", "url": "https://www.r-project.org/", "kind": "software"},
    {"name": "Python", "url": "https://www.python.org/", "kind": "software"},
    {"name": "Seurat", "url": "https://satijalab.org/seurat/", "kind": "package", "install": "install.packages('Seurat')"},
    {"name": "scDblFinder", "url": "https://bioconductor.org/packages/scDblFinder/", "kind": "package", "install": "BiocManager::install('scDblFinder')"},
    {"name": "SingleCellExperiment", "url": "https://bioconductor.org/packages/SingleCellExperiment/", "kind": "package", "install": "BiocManager::install('SingleCellExperiment')"},
    {"name": "clusterProfiler", "url": "https://bioconductor.org/packages/clusterProfiler/", "kind": "package", "install": "BiocManager::install('clusterProfiler')"},
    {"name": "enrichplot", "url": "https://bioconductor.org/packages/enrichplot/", "kind": "package", "install": "BiocManager::install('enrichplot')"},
    {"name": "BiocParallel", "url": "https://bioconductor.org/packages/BiocParallel/", "kind": "package", "install": "BiocManager::install('BiocParallel')"},
    {"name": "org.Hs.eg.db", "url": "https://bioconductor.org/packages/org.Hs.eg.db/", "kind": "package", "install": "BiocManager::install('org.Hs.eg.db')"},
    {"name": "org.Mm.eg.db", "url": "https://bioconductor.org/packages/org.Mm.eg.db/", "kind": "package", "install": "BiocManager::install('org.Mm.eg.db')"},
    {"name": "DESeq2", "url": "https://bioconductor.org/packages/DESeq2/", "kind": "package", "install": "BiocManager::install('DESeq2')"},
    {"name": "data.table", "url": "https://rdatatable.gitlab.io/data.table/", "kind": "package", "install": "install.packages('data.table')"},
    {"name": "ggplot2", "url": "https://ggplot2.tidyverse.org/", "kind": "package", "install": "install.packages('ggplot2')"},
]


def validate_accession(accession: str) -> str:
    acc = accession.strip().upper()
    match = re.fullmatch(r"E-GEOD-(\d+)", acc)
    if match:
        acc = "GSE" + match.group(1)
    if not re.fullmatch(
        r"(?:GSE\d+|E-[A-Z0-9]+-\d+|S-BSST\d+)",
        acc,
    ):
        raise ValueError(
            "数据集编号格式不正确，支持 GSE125449、E-MTAB-1234、S-BSST123"
        )
    return acc


def find_rscript() -> str:
    found = shutil.which("Rscript")
    if found:
        return found
    for base in [
        Path(r"C:\Program Files\R"),
        Path(r"C:\Program Files\Microsoft\R Open"),
    ]:
        if base.exists():
            candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
            if candidates:
                return str(candidates[0])
    raise RuntimeError("Rscript not found")


def get_versions() -> list[dict]:
    r_packages = [
        item["name"] for item in SOFTWARE
        if item["kind"] == "package"
    ]
    script = (
        "pkgs <- c(" + ",".join(f'"{p}"' for p in r_packages) + "); "
        "ip <- installed.packages(); "
        "cat(sapply(pkgs, function(p) if(p %in% rownames(ip)) ip[p,'Version'] else 'not installed'), sep='|')"
    )
    r_version = "unknown"
    versions = {}
    try:
        result = subprocess.run(
            [find_rscript(), "-e", script],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        values = result.stdout.strip().split("|")
        versions = dict(zip(r_packages, values))
        ver_res = subprocess.run(
            [find_rscript(), "-e", "cat(as.character(getRversion()))"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
        )
        r_version = ver_res.stdout.strip() or "unknown"
    except Exception:
        versions = {}

    out = []
    for item in SOFTWARE:
        if item["name"] == "R":
            version = r_version
        elif item["name"] == "Python":
            version = sys.version.split()[0]
        elif item["kind"] == "package":
            version = versions.get(item["name"], "not installed")
        else:
            version = "unknown"
        out.append({**item, "version": version})
    return out


def start_job(
    accession: str,
    output: str,
    species: str,
    skip_figs: list[str],
    figure_styles: dict[str, str],
    params: dict[str, str] | None = None,
) -> dict:
    acc = validate_accession(accession)
    out = Path(output).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    job_id = uuid.uuid4().hex[:8]
    log_dir = out / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"web_{job_id}.log"

    cmd = [
        sys.executable,
        str(SCRIPTS_DIR / "run_pipeline.py"),
        acc,
        "--output",
        str(out),
        "--species",
        species,
    ]
    env = os.environ.copy()
    env["LIVER_SKIP_FIGURES"] = ",".join(skip_figs)
    if figure_styles:
        env["LIVER_FIGURE_STYLES"] = ",".join(
            f"{name}={style}" for name, style in figure_styles.items()
        )
    for key, value in (params or {}).items():
        env[key] = str(value)

    JOBS[job_id] = {
        "job_id": job_id,
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "out": out,
        "accession": acc,
        "species": species,
        "skip_figs": skip_figs,
        "figure_styles": figure_styles,
        "cmd": cmd,
        "env": env,
        "queued": True,
        "recorded": False,
        "paused": False,
        "notified": False,
    }
    QUEUE.append(JOBS[job_id])
    _drain_queue()
    return {
        "job": job_id,
        "log_url": f"/log?job={job_id}",
        "status_url": f"/status?job={job_id}",
    }


def _start_process(info: dict) -> None:
    log_handle = info["log"].open("w", encoding="utf-8", errors="replace")
    proc = subprocess.Popen(
        info["cmd"],
        cwd=APP_ROOT,
        env=info["env"],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    info["proc"] = proc
    info["queued"] = False
    info["started"] = time.time()


def _drain_queue() -> None:
    with QUEUE_LOCK:
        for info in QUEUE:
            if info.get("proc") is None:
                _start_process(info)
                break


def load_history() -> list[dict]:
    if not HISTORY_PATH.exists():
        return []
    try:
        return json.loads(HISTORY_PATH.read_text(encoding="utf-8"))
    except Exception:
        return []


def save_history(records: list[dict]) -> None:
    HISTORY_PATH.write_text(
        json.dumps(records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def record_job(info: dict, ok: bool) -> None:
    fig_dir = info["out"] / "results" / "figures"
    figure_count = 0
    if fig_dir.exists():
        figure_count = len(list(fig_dir.rglob("*.png")))
    record = {
        "job": info["log"].stem.replace("web_", ""),
        "accession": info["accession"],
        "output": str(info["out"]),
        "species": info["species"],
        "status": "success" if ok else "failed",
        "started": info["started"],
        "finished": time.time(),
        "figures": figure_count,
        "report": str(info["out"] / "results" / "result_report.html"),
    }
    records = load_history()
    records.insert(0, record)
    save_history(records)


def resume_job(job_id: str) -> dict:
    info = JOBS.get(job_id)
    if not info:
        raise ValueError("job not found")
    if info.get("proc") is None:
        _drain_queue()
        return {"job": job_id, "queued": True}
    if info["proc"].poll() is None:
        return {"job": job_id, "already_running": True}

    pause_path = info["out"] / "pause_request.flag"
    if pause_path.exists():
        pause_path.unlink()

    log_path = info["log"]
    log_handle = log_path.open("a", encoding="utf-8", errors="replace")
    log_handle.write("\n[pipeline] resume requested\n")
    log_handle.flush()
    proc = subprocess.Popen(
        info["cmd"],
        cwd=APP_ROOT,
        env=info.get("env", os.environ.copy()),
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    info["proc"] = proc
    info["notified"] = False
    info["started"] = time.time()
    return {"job": job_id, "resumed": True}


PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>表达分析流水线</title>
<style>
body { font-family: "Segoe UI", Arial, sans-serif; background: #f5f7fa; margin: 0; color: #1f2933; }
.wrap { max-width: 1100px; margin: 0 auto; padding: 28px; }
h1 { font-size: 26px; }
h2 { font-size: 18px; }
.card { background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 20px; margin-bottom: 18px; }
label { display: block; margin: 12px 0 4px; font-weight: 600; }
label.toggle { display: inline-block; width: 48%; margin: 5px 0; font-weight: normal; }
label.toggle select { width: 58px; margin-left: 6px; padding: 4px; display: inline-block; }
input, select { width: 100%; padding: 10px; border: 1px solid #cbd2d9; border-radius: 6px; box-sizing: border-box; }
button { margin-top: 18px; padding: 10px 18px; border: 0; border-radius: 6px; background: #1665c0; color: #fff; font-size: 15px; cursor: pointer; }
button:disabled { background: #9aa5b1; cursor: not-allowed; }
pre { background: #0f172a; color: #dbeafe; padding: 14px; border-radius: 8px; height: 420px; overflow: auto; font-size: 12px; white-space: pre-wrap; }
.error { color: #b91c1c; margin-top: 10px; }
.ok { color: #047857; margin-top: 10px; }
.gallery { display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 16px; }
.gallery figure { margin: 0; }
.gallery img { width: 100%; border: 1px solid #e4e7eb; border-radius: 6px; background: #fff; }
.gallery figcaption { color: #52606d; font-size: 12px; margin-top: 5px; }
.job-toast {
  position: fixed;
  top: 76px;
  right: 24px;
  z-index: 999;
  max-width: 440px;
  background: #ffffff;
  border: 1px solid #d7dde4;
  border-left: 6px solid #1665c0;
  border-radius: 8px;
  box-shadow: 0 12px 32px rgba(15, 23, 42, 0.28);
  padding: 14px 16px;
}
.job-toast.ok { border-left-color: #047857; }
.job-toast.error { border-left-color: #b91c1c; }
.job-toast.paused { border-left-color: #b45309; }
.job-toast-title { font-weight: 700; font-size: 15px; margin-bottom: 6px; }
.job-toast-body {
  font-size: 13px;
  color: #334155;
  white-space: pre-wrap;
  word-break: break-all;
  max-height: 190px;
  overflow: auto;
}
.job-toast-close {
  margin-top: 10px;
  padding: 6px 14px;
  border: 0;
  border-radius: 6px;
  background: #1665c0;
  color: #fff;
  cursor: pointer;
}
</style>
</head>
<body>
<div class="wrap">
<h1>表达分析流水线</h1>
<div class="card">
  <form id="form">
    <label for="acc">GSE 数据集编号</label>
    <input id="acc" name="accession" placeholder="GSE125449" required>

    <label for="out">结果保存路径</label>
    <input id="out" name="output" placeholder="请输入结果保存地址" required>

    <label for="sp">物种</label>
    <select id="sp" name="species">
      <option value="hs">hs - Human</option>
      <option value="mm">mm - Mouse</option>
    </select>

    <div id="figureToggles" style="margin-top: 16px;">
      <h2>结果图开关</h2>
      <p class="muted">每个结果图默认 yes；选择 no 则跳过该图。</p>
<!--FIGURE_TOGGLES-->
    </div>

    <button id="startBtn" type="button" onclick="startRun()">开始分析</button>
  </form>
  <div id="message"></div>
</div>
<div class="card">
  <h2>实时进度</h2>
  <pre id="log">等待开始...</pre>
</div>
<div class="card">
  <h2>结果图</h2>
  <div id="gallery" class="gallery"><p class="muted">分析完成后在此显示。</p></div>
</div>
</div>

<script>
function esc(value) {
  return String(value == null ? '' : value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}
let currentJob = null;
let pollTimer = null;
let jobAlerted = false;
const SINGLE_JOB_KEY = 'liver_ui_single_job';

function saveJobRecord() {
  try { sessionStorage.setItem(SINGLE_JOB_KEY, currentJob); } catch (e) {}
}

function clearJobRecord() {
  try { sessionStorage.removeItem(SINGLE_JOB_KEY); } catch (e) {}
}

function restoreJobRecord() {
  let saved = null;
  let target = null;
  try { saved = sessionStorage.getItem(SINGLE_JOB_KEY); } catch (e) {}
  try { target = new URLSearchParams(window.location.search).get('job'); } catch (e) {}
  const job = target || saved;
  if (!job) return;
  currentJob = job;
  saveJobRecord();
  if (target) {
    try {
      const url = new URL(window.location.href);
      url.searchParams.delete('job');
      window.history.replaceState(null, '', url.toString());
    } catch (e) {}
  }
  const msg = document.getElementById('message');
  msg.className = 'ok';
  msg.textContent = '已恢复正在运行的分析任务，日志继续刷新。';
  if (pollTimer) clearInterval(pollTimer);
  pollTimer = setInterval(pollLog, 1000);
}

function showJobToast(kind, title, body) {
  const toast = document.getElementById('jobToast');
  if (!toast) return;
  toast.className = 'job-toast ' + (kind || 'info');
  document.getElementById('jobToastTitle').textContent = title;
  document.getElementById('jobToastBody').textContent = body || '';
  toast.hidden = false;
  jobAlerted = true;
}

function closeJobToast() {
  const toast = document.getElementById('jobToast');
  if (toast) toast.hidden = true;
}

async function startRun() {
  const form = document.getElementById('form');
  const data = new URLSearchParams(new FormData(form));
  const btn = document.getElementById('startBtn');
  const msg = document.getElementById('message');
  btn.disabled = true;
  jobAlerted = false;
  closeJobToast();
  msg.className = '';
  msg.textContent = '正在启动...';

  try {
    const resp = await fetch('/start', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '启动失败');
    currentJob = result.job;
    saveJobRecord();
    document.getElementById('log').textContent = '';
    document.getElementById('gallery').innerHTML = '<p class="muted">分析完成后在此显示。</p>';
    msg.className = 'ok';
    msg.textContent = '任务已启动，日志会实时刷新。';
    if (pollTimer) clearInterval(pollTimer);
    pollTimer = setInterval(pollLog, 1000);
  } catch (e) {
    msg.className = 'error';
    msg.textContent = String(e.message || e);
  } finally {
    btn.disabled = false;
  }
}

async function pollLog() {
  if (!currentJob) return;
  try {
    const resp = await fetch('/log?job=' + currentJob);
    const text = await resp.text();
    const box = document.getElementById('log');
    box.textContent = text;
    box.scrollTop = box.scrollHeight;
    const statusResp = await fetch('/status?job=' + currentJob);
    const status = await statusResp.json();
    if (status.queued) {
      const msg = document.getElementById('message');
      msg.className = 'ok';
      msg.textContent = '任务已进入队列，等待前面的分析完成...';
      return;
    }
    if (!status.running) {
      clearInterval(pollTimer);
      const msg = document.getElementById('message');
      if (status.ok) {
        msg.className = 'ok';
        msg.textContent = '流水线已完成';
        if (!jobAlerted) showJobToast('ok', '任务已完成', '表达分析任务已全部完成。');
      } else {
        msg.className = 'error';
        msg.textContent = '流水线运行失败，请查看日志';
        if (!jobAlerted) {
          showJobToast('error', '任务中断', '运行到阶段：' + (status.stage || '未知') + '\n原因：' + (status.error || '请查看日志'));
        }
      }
      if (status.ok) await loadResults(currentJob);
      currentJob = null;
      clearJobRecord();
    }
  } catch (e) {
    // keep polling on transient errors
  }
}

async function loadResults(job) {
  const gallery = document.getElementById('gallery');
  try {
    const resp = await fetch('/files?job=' + job);
    const data = await resp.json();
    gallery.innerHTML = '';
    if (!data.figures || data.figures.length === 0) {
      gallery.innerHTML = '<p class="muted">未生成结果图。</p>';
      return;
    }
    data.figures.forEach(name => {
      const figure = document.createElement('figure');
      const img = document.createElement('img');
      img.src = '/figure?job=' + job + '&name=' + encodeURIComponent(name);
      const caption = document.createElement('figcaption');
      caption.textContent = name;
      figure.appendChild(img);
      figure.appendChild(caption);
      gallery.appendChild(figure);
    });
  } catch (e) {
    gallery.innerHTML = '<p class="error">结果图加载失败：' + esc(e.message || e) + '</p>';
  }
}

restoreJobRecord();
</script>
<div id="jobToast" class="job-toast" hidden>
  <div class="job-toast-title" id="jobToastTitle"></div>
  <div class="job-toast-body" id="jobToastBody"></div>
  <button type="button" class="job-toast-close" onclick="closeJobToast()">知道了</button>
</div>
</body>
</html>
"""


def render_page() -> str:
    template = (
        PAGE_TEMPLATE_PATH.read_text(encoding="utf-8")
        if PAGE_TEMPLATE_PATH.exists()
        else PAGE_TEMPLATE
    )
    rows = []
    for item in FIGURES:
        options = (
            '<select name="fig_' + item["file"] + '">'
            '<option value="yes" selected>yes</option>'
            '<option value="no">no</option>'
            '</select>'
        )
        if item.get("styles"):
            style_options = "".join(
                (
                    f'<option value="{style}">'
                    f'{STYLE_LABELS.get(style, style)} ({style})'
                    f'</option>'
                )
                for style in item["styles"]
            )
            options += (
                '<select class="style-select" name="style_' +
                item["file"] + '">' + style_options + '</select>'
            )
        rows.append(
            f'<div class="figrow" id="figrow_{item["file"]}">'
            f'<div class="fighead">'
            f'<span>{item["label"]}</span>'
            f'<div class="fig-opts">{options}</div>'
            f'</div>'
            f'<div class="figslot" data-name="{item["file"]}"></div>'
            f'</div>'
        )
    html = template.replace("<!--FIGURE_TOGGLES-->", "\n".join(rows))
    html = html.replace("</style>", NAV_CSS + "</style>", 1)
    html = html.replace('<div class="wrap">', NAV_HTML + '\n<div class="wrap">', 1)
    return html


def get_page() -> str:
    return render_page()


def render_dock_page() -> str:
    if DOCK_TEMPLATE_PATH.exists():
        return DOCK_TEMPLATE_PATH.read_text(encoding="utf-8")
    return (
        "<html><body><h1>dock template missing</h1>"
        "<p>web/templates/dock_page_template.html not found</p></body></html>"
    )


def render_full_page() -> str:
    if FULL_TEMPLATE_PATH.exists():
        return FULL_TEMPLATE_PATH.read_text(encoding="utf-8")
    return (
        "<html><body><h1>full pipeline template missing</h1>"
        "<p>web/templates/full_page_template.html not found</p></body></html>"
    )


def render_results_page() -> str:
    if RESULTS_TEMPLATE_PATH.exists():
        return RESULTS_TEMPLATE_PATH.read_text(encoding="utf-8")
    return (
        "<html><body><h1>results template missing</h1>"
        "<p>web/templates/results_manifest_optimized.html not found</p></body></html>"
    )


def result_guide_data() -> dict:
    """Return the result-figure guide as searchable sections."""
    if not RESULT_GUIDE_PATH.exists():
        return {
            "available": False,
            "source": str(RESULT_GUIDE_PATH),
            "sections": [],
            "files": [],
        }
    text = RESULT_GUIDE_PATH.read_text(encoding="utf-8", errors="replace")
    sections: list[dict] = []
    title = "总览"
    lines: list[str] = []

    def flush_section() -> None:
        body = "\n".join(lines).strip()
        if body:
            sections.append(
                {
                    "id": f"guide-{len(sections) + 1}",
                    "title": title,
                    "text": body,
                    "files": sorted(
                        set(
                            re.findall(
                                r"\b(fig_[A-Za-z0-9_.-]+\.(?:png|csv|json|rds|txt))\b",
                                body,
                                flags=re.IGNORECASE,
                            )
                        )
                    ),
                }
            )
        lines.clear()

    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if line.startswith("## "):
            flush_section()
            title = line[3:].strip()
        else:
            lines.append(line)
    flush_section()
    return {
        "available": True,
        "source": str(RESULT_GUIDE_PATH),
        "sections": sections,
        "files": sorted(
            set(
                re.findall(
                    r"\b(fig_[A-Za-z0-9_.-]+\.(?:png|csv|json|rds|txt))\b",
                    text,
                    flags=re.IGNORECASE,
                )
            )
        ),
    }


def result_details_data() -> dict:
    """Return detailed content, type and use descriptions for result files."""
    if not RESULT_DETAILS_PATH.exists():
        return {"available": False, "entries": []}
    try:
        data = json.loads(RESULT_DETAILS_PATH.read_text(encoding="utf-8"))
        data["available"] = True
        return data
    except (OSError, ValueError):
        return {"available": False, "entries": []}


def render_tasks_page() -> str:
    if TASKS_TEMPLATE_PATH.exists():
        return TASKS_TEMPLATE_PATH.read_text(encoding="utf-8")
    return (
        "<html><body><h1>tasks template missing</h1>"
        "<p>web/templates/tasks_template.html not found</p></body></html>"
    )


def render_datasets_page() -> str:
    if DATASET_TEMPLATE_PATH.exists():
        return DATASET_TEMPLATE_PATH.read_text(encoding="utf-8")
    return (
        "<html><body><h1>datasets template missing</h1>"
        "<p>web/templates/datasets_template.html not found</p></body></html>"
    )


def _first(data: dict, key: str, default: str = "") -> str:
    values = data.get(key)
    if not values:
        return default
    return str(values[0])


def _float3(data: dict, prefix: str):
    vals = [
        _first(data, f"{prefix}_x"),
        _first(data, f"{prefix}_y"),
        _first(data, f"{prefix}_z"),
    ]
    if any(str(v).strip() == "" for v in vals):
        return None
    return [float(v) for v in vals]


def _int_field(data: dict, key: str):
    value = _first(data, key, "")
    if str(value).strip() == "":
        return None
    return int(float(value))


def _float_field(data: dict, key: str):
    value = _first(data, key, "")
    if str(value).strip() == "":
        return None
    return float(value)


def dataset_search_request(data: dict) -> dict:
    disease = _first(data, "disease", "").strip()
    research_direction = _first(data, "research_direction", "").strip()
    query = _first(data, "query", "").strip()
    if not disease and not query:
        raise ValueError("疾病名称或原始查询至少填写一项")
    max_results = 20
    try:
        max_results = max(
            1,
            min(100, int(float(_first(data, "max_results", "20") or 20))),
        )
    except ValueError:
        pass
    raw_organism = _first(data, "organism", "").strip()
    organism_aliases = {
        "hs": "Homo sapiens",
        "human": "Homo sapiens",
        "homo sapiens": "Homo sapiens",
        "mm": "Mus musculus",
        "mouse": "Mus musculus",
        "mus musculus": "Mus musculus",
        "auto": "",
        "all": "",
    }
    organism = organism_aliases.get(raw_organism.lower(), raw_organism) or None
    keyword = _first(data, "keyword", "").strip() or None
    data_type = _first(data, "data_type", "").strip() or None
    min_samples = _int_field(data, "min_samples")
    max_samples = _int_field(data, "max_samples")
    start_date = _first(data, "start_date", "").strip() or None
    end_date = _first(data, "end_date", "").strip() or None
    platform = _first(data, "platform", "").strip() or None
    dataset_type = _first(data, "dataset_type", "").strip() or None
    raw_databases = data.get("databases") or []
    if isinstance(raw_databases, str):
        raw_databases = [raw_databases]
    databases = [
        str(item).strip().lower()
        for item in raw_databases
        if str(item).strip()
    ]
    if not databases:
        databases = None
    model_value = _first(data, "model", "").strip()

    import search_datasets as sd

    query_text = sd.build_query(disease, research_direction, query or None)
    rows = sd.search_datasets(
        query_text,
        max_results=max_results,
        organism=organism,
        keyword=keyword,
        data_type=data_type,
        min_samples=min_samples,
        max_samples=max_samples,
        start_date=start_date,
        end_date=end_date,
        platform=platform,
        dataset_type=dataset_type,
        disease=disease,
        research_direction=research_direction,
        databases=databases,
    )
    model_applied = False
    model_path = ""
    if model_value:
        model_file = Path(model_value).expanduser()
        if not model_file.is_file():
            raise ValueError(f"模型文件不存在：{model_file}")
        from dataset_search_ml import load_model, rerank

        rows = rerank(
            rows,
            disease,
            research_direction,
            model=load_model(model_file),
        )
        model_applied = True
        model_path = str(model_file)

    DATASET_SEARCH_DIR.mkdir(parents=True, exist_ok=True)
    csv_path, json_path = sd.write_outputs(rows, DATASET_SEARCH_DIR)
    result_rows = [
        dict(row, full_pipeline_url=dataset_full_pipeline_url(row))
        for row in rows
    ]
    return {
        "query": query_text,
        "disease": disease,
        "research_direction": research_direction,
        "model_applied": model_applied,
        "model_path": model_path,
        "filters": {
            "organism": organism or "",
            "keyword": keyword or "",
            "data_type": data_type or "",
            "min_samples": min_samples,
            "max_samples": max_samples,
            "start_date": start_date or "",
            "end_date": end_date or "",
            "platform": platform or "",
            "dataset_type": dataset_type or "",
                "databases": (
                    databases
                    if databases
                    else list(DATASET_DATABASES)
                ),
        },
        "count": len(rows),
        "results": result_rows,
        "output_dir": str(DATASET_SEARCH_DIR),
        "csv_url": f"/datasets/file?name={csv_path.name}",
        "json_url": f"/datasets/file?name={json_path.name}",
    }


def start_dataset_download(data: dict) -> dict:
    accessions: list[str] = []
    for value in data.get("accessions", []):
        accessions.extend(
            acc.strip().upper()
            for acc in str(value).split(",")
            if acc.strip()
        )
    if not accessions:
        for value in data.get("accession", []):
            accessions.extend(
                acc.strip().upper()
                for acc in str(value).split(",")
                if acc.strip()
            )
    if not accessions:
        raise ValueError("至少选择一个数据集")
    accessions = list(dict.fromkeys(validate_accession(acc) for acc in accessions))
    download_root = Path(
        _first(data, "download_root", str(APP_ROOT.parent / "liver_cancer"))
    ).expanduser().resolve()
    DATASET_SEARCH_DIR.mkdir(parents=True, exist_ok=True)
    job_id = uuid.uuid4().hex[:8]
    log_path = DATASET_SEARCH_DIR / f"download_{job_id}.log"
    info = {
        "job_id": job_id,
        "log": log_path,
        "results": {},
        "error": "",
        "running": True,
        "started": time.time(),
    }
    with DATASET_DOWNLOAD_LOCK:
        DATASET_DOWNLOAD_JOBS[job_id] = info
    threading.Thread(
        target=_run_dataset_download,
        args=(info, accessions, download_root),
        daemon=True,
    ).start()
    return {
        "job": job_id,
        "status_url": f"/datasets/download/status?job={job_id}",
    }


def _run_dataset_download(
    info: dict,
    accessions: list[str],
    download_root: Path,
) -> None:
    import search_datasets as sd

    log_path = info["log"]
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_path.open("w", encoding="utf-8", errors="replace") as log_handle:
            log_handle.write(f"starting download: {', '.join(accessions)}\n")
            log_handle.flush()
            def log(*args):
                log_handle.write(
                    " ".join(str(arg) for arg in args) + "\n"
                )
                log_handle.flush()
            try:
                results = sd.download_accessions(
                    accessions,
                    download_root,
                    log=log,
                )
                info["results"] = results
                (DATASET_SEARCH_DIR / "download_results.json").write_text(
                    json.dumps(results, ensure_ascii=False, indent=2),
                    encoding="utf-8",
                )
                log_handle.write("download results:\n")
                log_handle.write(
                    json.dumps(results, ensure_ascii=False, indent=2)
                )
                log_handle.write("\n")
            except Exception as exc:  # noqa: BLE001
                info["error"] = str(exc)
                log_handle.write(f"download error: {exc}\n")
    finally:
        info["running"] = False


def dataset_download_status(job_id: str) -> dict:
    with DATASET_DOWNLOAD_LOCK:
        info = DATASET_DOWNLOAD_JOBS.get(job_id)
    if not info:
        raise ValueError("下载任务不存在")
    running = bool(info.get("running"))
    log_text = (
        info["log"].read_text(encoding="utf-8", errors="replace")
        if info["log"].exists()
        else ""
    )
    return {
        "job": job_id,
        "running": running,
        "ok": not running and not info.get("error"),
        "error": info.get("error", ""),
        "results": info.get("results", {}),
        "log": log_text[-8000:],
    }


def dataset_file_path(name: str) -> Path | None:
    target = (DATASET_SEARCH_DIR / Path(name).name).resolve()
    if target.parent != DATASET_SEARCH_DIR.resolve() or not target.is_file():
        return None
    return target


def dataset_full_pipeline_url(row: dict) -> str:
    """Build the full-pipeline page URL prefilled with a searched dataset."""
    if row.get("run_supported") is False:
        return ""
    accession = validate_accession(str(row.get("accession") or ""))
    params = {"accession": accession}
    title = str(row.get("title") or "").strip()
    if title:
        params["dataset_title"] = title
    data_type = str(row.get("data_type") or "").strip().lower()
    if data_type in ("single-cell", "bulk", "other"):
        params["data_type"] = data_type
    organism = str(row.get("organism") or "").lower()
    species = "auto"
    if "homo sapiens" in organism or "human" in organism:
        species = "hs"
    elif "mus musculus" in organism or "mouse" in organism:
        species = "mm"
    if species != "auto":
        params["species"] = species
    return "/full?" + urlencode(params)


def _single_report_path(info: dict | None = None, output: str = "") -> Path | None:
    if info is not None:
        out = Path(info["out"]).resolve()
    elif output:
        out = Path(output).expanduser().resolve()
    else:
        return None
    target = (out / "results" / "result_report.html").resolve()
    if target.is_file() and target.is_relative_to(out):
        return target
    return None


def start_dock_job(data: dict) -> dict:
    from docking.config import load_config, save_config

    workdir = Path(
        _first(data, "workdir", str(APP_ROOT / "dock"))
    ).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    stage = _first(data, "stage", "pipeline")
    job_id = uuid.uuid4().hex[:8]
    cfg_path = workdir / "config" / f"docking_web_{job_id}.json"

    overrides = {
        "workdir": str(workdir),
        "receptor": _first(data, "receptor", "") or None,
        "ligand": _first(data, "ligand", "") or None,
        "center": _float3(data, "center"),
        "size": _float3(data, "size"),
        "exhaustiveness": _int_field(data, "exhaustiveness"),
        "num_modes": _int_field(data, "num_modes"),
        "energy_range": _float_field(data, "energy_range"),
        "max_workers": _int_field(data, "max_workers"),
        "cutoff": _float_field(data, "cutoff"),
        "top_n": _int_field(data, "top_n"),
        "model": _first(data, "model", "") or None,
        "training_csv": _first(data, "training_csv", "") or None,
        "label_column": _first(data, "label_column", "") or None,
    }
    cfg = load_config(APP_ROOT / "config" / "docking_config.json", overrides)
    save_config(cfg, cfg_path)

    log_path = workdir / "logs" / f"web_dock_{job_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPTS_DIR / "run_docking.py"),
        stage,
        "--config",
        str(cfg_path),
    ]
    if _first(data, "force", "") in ("1", "true", "on", "yes"):
        cmd.append("--force")

    env = os.environ.copy()
    env["DOCK_WORKDIR"] = str(workdir)
    DOCK_JOBS[job_id] = {
        "job_id": job_id,
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "workdir": workdir,
        "output_dir": cfg.output_dir,
        "stage": stage,
        "cmd": cmd,
        "env": env,
        "queued": True,
        "recorded": False,
        "notified": False,
    }
    DOCK_QUEUE.append(DOCK_JOBS[job_id])
    _drain_dock_queue()
    return {
        "job": job_id,
        "log_url": f"/dock/log?job={job_id}",
        "status_url": f"/dock/status?job={job_id}",
    }


def _start_dock_process(info: dict) -> None:
    log_handle = info["log"].open("w", encoding="utf-8", errors="replace")
    proc = subprocess.Popen(
        info["cmd"],
        cwd=APP_ROOT,
        env=info["env"],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    info["proc"] = proc
    info["queued"] = False
    info["started"] = time.time()


def _drain_dock_queue() -> None:
    with DOCK_QUEUE_LOCK:
        for info in DOCK_QUEUE:
            if info.get("proc") is None:
                _start_dock_process(info)
                break


def _dock_status(info: dict) -> dict:
    queued = info.get("proc") is None
    if queued:
        return {"running": False, "ok": False, "queued": True, "paused": False, "stage": "", "error": ""}
    running = info["proc"].poll() is None
    paused = bool(info.get("paused"))
    ok = False
    if not running:
        marker_dir = info["output_dir"] / ".stages"
        log_paths = [info["log"]]
        if paused:
            stage, error = _finished_info(marker_dir, DOCK_STAGE_LABELS, log_paths, True)
            _notify_finished(
                info,
                "dock",
                "虚拟筛选",
                f"{info.get('stage', 'pipeline')} 虚拟筛选",
                "paused",
                stage,
                error,
                exit_code=info["proc"].returncode,
            )
            return {
                "running": False,
                "ok": False,
                "queued": False,
                "paused": True,
                "stage": stage,
                "error": error,
            }
        ok = info["proc"].returncode == 0
        if not info.get("recorded"):
            record_dock_job(info, ok)
            info["recorded"] = True
        _drain_dock_queue()
        stage, error = _finished_info(marker_dir, DOCK_STAGE_LABELS, log_paths, False)
        status = "completed" if ok else "interrupted"
        if status == "completed":
            error = ""
        _notify_finished(
            info,
            "dock",
            "虚拟筛选",
            f"{info.get('stage', 'pipeline')} 虚拟筛选",
            status,
            stage,
            error,
            exit_code=info["proc"].returncode,
        )
        return {
            "running": False,
            "ok": ok,
            "queued": False,
            "paused": False,
            "stage": stage,
            "error": error,
        }
    return {"running": True, "ok": False, "queued": False, "paused": False, "stage": "", "error": ""}


def start_full_job(data: dict) -> dict:
    output_value = _first(data, "output", "").strip()
    if not output_value:
        raise ValueError("表达分析结果目录不能为空，请输入结果保存地址")
    workdir_value = _first(data, "workdir", "").strip()
    if not workdir_value:
        raise ValueError("\u5de5\u4f5c\u76ee\u5f55\u4e0d\u80fd\u4e3a\u7a7a\uff0c\u8bf7\u624b\u52a8\u8f93\u5165")
    accession = _first(data, "accession", "").strip()
    if accession:
        accession = validate_accession(accession)
    workdir = Path(workdir_value).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    job_id = uuid.uuid4().hex[:8]
    log_path = workdir / "logs" / f"web_full_{job_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(SCRIPTS_DIR / "run_full_pipeline.py"),
        "--config",
        str(APP_ROOT / "config" / "full_pipeline_config.json"),
        "--docking-config",
        str(APP_ROOT / "config" / "docking_config.json"),
    ]
    cmd += ["--workdir", str(workdir)]
    if accession:
        cmd += ["--accession", accession]
    output = _first(data, "output", "").strip()
    if output:
        cmd += ["--output", output]
    species = _first(data, "species", "").strip()
    if species:
        cmd += ["--species", species]
    ml_model = _first(data, "ml_model", "").strip()
    if ml_model:
        cmd += ["--ml-model", ml_model]
    top_genes = _int_field(data, "top_genes")
    if top_genes:
        cmd += ["--top-genes", str(top_genes)]
    docking_targets = _int_field(data, "docking_targets")
    if docking_targets is not None:
        cmd += ["--docking-targets", str(docking_targets)]
    ko_top_n = _int_field(data, "ko_top_n")
    if ko_top_n:
        cmd += ["--ko-top-n", str(ko_top_n)]
    feedback_top_n = _int_field(data, "feedback_top_n")
    if feedback_top_n:
        cmd += ["--feedback-top-n", str(feedback_top_n)]
    feedback_max_features = _int_field(data, "feedback_max_features")
    if feedback_max_features:
        cmd += ["--feedback-max-features", str(feedback_max_features)]
    ligand_library = _first(data, "ligand_library", "").strip()
    if ligand_library:
        cmd += ["--ligand-library", ligand_library]
    depmap_csv = _first(data, "depmap_csv", "").strip()
    if depmap_csv:
        cmd += ["--depmap-csv", depmap_csv]
    ppi_network_csv = _first(data, "ppi_network_csv", "").strip()
    if ppi_network_csv:
        cmd += ["--ppi-network-csv", ppi_network_csv]
    case_label = _first(data, "case_label", "").strip()
    if case_label:
        cmd += ["--case-label", case_label]
    normal_label = _first(data, "normal_label", "").strip()
    if normal_label:
        cmd += ["--normal-label", normal_label]
    start_stage = _first(data, "start_stage", "").strip()
    if start_stage:
        cmd += ["--start-stage", start_stage]
    for flag in [
        "skip_scrna",
        "skip_download",
        "skip_deps",
        "skip_evidence_fetch",
        "skip_pseudobulk",
        "skip_knockout",
        "skip_docking",
        "skip_cell_feedback",
        "skip_qc_gate",
        "skip_differential_abundance",
        "keep_all_genes",
        "force",
        "dry_run",
    ]:
        if _first(data, flag, "") in ("1", "true", "on", "yes"):
            cmd.append("--" + flag.replace("_", "-"))

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    FULL_JOBS[job_id] = {
        "job_id": job_id,
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "workdir": workdir,
        "accession": accession,
        "output": output,
        "cmd": cmd,
        "env": env,
        "queued": True,
        "notified": False,
    }
    FULL_QUEUE.append(FULL_JOBS[job_id])
    _drain_full_queue()
    return {
        "job": job_id,
        "log_url": f"/full/log?job={job_id}",
        "status_url": f"/full/status?job={job_id}",
    }


def _start_full_process(info: dict) -> None:
    log_handle = info["log"].open("w", encoding="utf-8", errors="replace")
    proc = subprocess.Popen(
        info["cmd"],
        cwd=APP_ROOT,
        env=info["env"],
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    info["proc"] = proc
    info["queued"] = False
    info["started"] = time.time()


def _drain_full_queue() -> None:
    with FULL_QUEUE_LOCK:
        for info in FULL_QUEUE:
            if info.get("proc") is None:
                _start_full_process(info)
                break


def _full_status(info: dict) -> dict:
    queued = info.get("proc") is None
    if queued:
        return {"running": False, "ok": False, "queued": True, "paused": False, "stage": "", "error": ""}
    if info.get("paused"):
        marker_dir = info["workdir"] / "outputs" / "integration" / ".stages"
        stage, error = _finished_info(
            marker_dir,
            FULL_STAGE_LABELS,
            _full_log_paths(info),
            True,
            SINGLE_STAGE_LABELS,
        )
        _notify_finished(
            info,
            "full",
            "全自动流水线",
            f"{info.get('accession') or '全自动流水线'} 全自动流水线",
            "paused",
            stage,
            error,
            exit_code=info["proc"].returncode if info.get("proc") else None,
        )
        return {
            "running": False,
            "ok": False,
            "queued": False,
            "paused": True,
            "stage": stage,
            "error": error,
        }
    running = info["proc"].poll() is None
    if not running:
        paused = bool(info.get("paused")) or info["proc"].returncode == 98
        if paused:
            marker_dir = info["workdir"] / "outputs" / "integration" / ".stages"
            stage, error = _finished_info(
                marker_dir,
                FULL_STAGE_LABELS,
                _full_log_paths(info),
                True,
                SINGLE_STAGE_LABELS,
            )
            _notify_finished(
                info,
                "full",
                "全自动流水线",
                f"{info.get('accession') or '全自动流水线'} 全自动流水线",
                "paused",
                stage,
                error,
                exit_code=info["proc"].returncode,
            )
            return {
                "running": False,
                "ok": False,
                "queued": False,
                "paused": True,
                "stage": stage,
                "error": error,
            }
        ok = info["proc"].returncode == 0
        _drain_full_queue()
        status = "completed" if ok else "interrupted"
        if ok:
            stage = "\u5168\u6d41\u7a0b\u5b8c\u6210"
            error = ""
        else:
            marker_dir = info["workdir"] / "outputs" / "integration" / ".stages"
            stage, error = _finished_info(
                marker_dir,
                FULL_STAGE_LABELS,
                _full_log_paths(info),
                False,
                SINGLE_STAGE_LABELS,
            )
        _notify_finished(
            info,
            "full",
            "全自动流水线",
            f"{info.get('accession') or '全自动流水线'} 全自动流水线",
            status,
            stage,
            error,
            exit_code=info["proc"].returncode,
        )
        return {
            "running": False,
            "ok": ok,
            "queued": False,
            "paused": False,
            "stage": stage,
            "error": error,
        }
    return {"running": True, "ok": False, "queued": False, "paused": False, "stage": "", "error": ""}


def _single_status(info: dict) -> dict:
    queued = info.get("proc") is None
    if queued:
        return {"running": False, "ok": False, "queued": True, "paused": False, "stage": "", "error": ""}
    running = info["proc"].poll() is None
    ok = False
    paused = False
    if not running:
        ok = info["proc"].returncode == 0
        paused = info["proc"].returncode == 98
        if not info.get("recorded") and not paused:
            record_job(info, ok)
            info["recorded"] = True
        _drain_queue()
        marker_dir = info["out"] / "results" / ".stages"
        log_paths = [info["log"], info["out"] / "logs" / "pipeline_r.log"]
        stage, error = _finished_info(marker_dir, SINGLE_STAGE_LABELS, log_paths, paused)
        status = "paused" if paused else "completed" if ok else "interrupted"
        if status == "completed":
            error = ""
        _notify_finished(
            info,
            "single",
            "表达分析",
            f"{info.get('accession', 'GSE')} 表达分析",
            status,
            stage,
            error,
            exit_code=info["proc"].returncode,
        )
        return {
            "running": False,
            "ok": ok,
            "queued": False,
            "paused": paused,
            "stage": stage,
            "error": error,
        }
    return {"running": True, "ok": False, "queued": False, "paused": False, "stage": "", "error": ""}


SINGLE_STAGE_LABELS = {
    "01": "数据加载",
    "02": "QC 过滤",
    "03": "双细胞检测",
    "04": "聚类",
    "05": "细胞注释",
    "06": "差异表达",
    "07": "富集分析",
    "08": "发表级分析",
    "09": "汇总输出",
}

DOCK_STAGE_LABELS = {
    "01": "受体准备",
    "02": "配体准备",
    "03": "分子对接",
    "04": "结果分析",
    "05": "精修重对接",
    "06": "HTML 报告",
}

FULL_STAGE_LABELS = {
    "01": "表达分析",
    "02": "关键基因",
    "03": "证据富集",
    "04": "敲除输入",
    "05": "虚拟敲除",
    "06": "分子对接",
    "07": "细胞反馈",
    "08": "集成报告",
}


def _marker_progress(marker_dir: Path, total: int) -> int:
    if not marker_dir.exists() or total <= 0:
        return 0
    try:
        done = len(list(marker_dir.glob("*.done")))
    except OSError:
        return 0
    return max(0, min(100, int(round(done * 100 / total))))


def _current_stage(marker_dir: Path, labels: dict) -> str:
    if not marker_dir.exists():
        return ""
    try:
        done = [path.stem for path in marker_dir.glob("*.done")]
    except OSError:
        return ""
    if not done:
        return ""
    codes = sorted(code[:2] for code in done)
    return labels.get(codes[-1], "")


def _log_tail(path: Path, limit: int = 1200) -> str:
    if not path.exists():
        return ""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        return text[-limit:]
    except OSError:
        return ""


def _stage_from_log(log_text: str, labels: dict) -> str:
    if not log_text:
        return ""
    matches = []
    for pattern in [
        r"start stage:\s*(\d+)[_\s-][A-Za-z0-9_-]+",
        r"=== stage (\d+) [A-Za-z0-9_-]+ ===",
    ]:
        matches.extend(re.finditer(pattern, log_text))
    if not matches:
        return ""
    return labels.get(matches[-1].group(1), "")


def _full_stage_from_log(
    log_text: str,
    full_labels: dict,
    single_labels: dict,
) -> str:
    full_matches = re.findall(
        r"=== stage (\d+) [A-Za-z0-9_-]+ ===",
        log_text,
    )
    if full_matches:
        return full_labels.get(full_matches[-1], "")
    single_matches = re.findall(
        r"start stage:\s*(\d+)[_\s-][A-Za-z0-9_-]+",
        log_text,
    )
    if single_matches:
        single_label = single_labels.get(single_matches[-1], "")
        return "表达分析" if single_label else "表达分析"
    return ""


def _full_log_paths(info: dict) -> list[Path]:
    paths = [Path(info["log"])]
    roots: list[Path] = []
    output = info.get("output")
    if output:
        roots.append(Path(output).expanduser().resolve())
    context_root = _single_cell_root_from_workdir(Path(info["workdir"]))
    if context_root is not None:
        roots.append(context_root)
    for root in roots:
        paths.append(root / "logs" / "pipeline_r.log")
    return paths


def _full_stage_label(info: dict, marker_dir: Path) -> str:
    log_text = "\n".join(
        _log_tail(path, 20000) for path in _full_log_paths(info)
    ).strip()
    stage = _full_stage_from_log(
        log_text,
        FULL_STAGE_LABELS,
        SINGLE_STAGE_LABELS,
    )
    if stage:
        return stage
    return _current_stage(marker_dir, FULL_STAGE_LABELS) or ""


def _extract_error(log_text: str, limit: int = 700) -> str:
    lines = [line.strip() for line in log_text.splitlines() if line.strip()]
    if not lines:
        return ""
    keywords = ("error", "traceback", "exception", "failed", "fatal", "cannot", "missing")
    error_lines = [
        line for line in lines
        if any(keyword in line.lower() for keyword in keywords)
    ]
    tail = "\n".join(error_lines[-6:]) if error_lines else "\n".join(lines[-8:])
    return tail[-limit:]


def _finished_info(
    marker_dir: Path,
    labels: dict,
    log_paths: list[Path],
    paused: bool,
    single_cell_labels: dict | None = None,
) -> tuple[str, str]:
    log_text = "\n".join(_log_tail(path, 20000) for path in log_paths).strip()
    stage = ""
    if single_cell_labels is not None:
        stage = _full_stage_from_log(log_text, labels, single_cell_labels)
        if not stage:
            stage = _stage_from_log(log_text, single_cell_labels)
            if stage:
                stage = "表达分析"
    else:
        stage = _stage_from_log(log_text, labels)
    stage = stage or _current_stage(marker_dir, labels) or "未知"
    if paused:
        error = "任务已暂停"
    else:
        error = _extract_error(log_text) or "进程已退出，请查看日志"
    return stage, error


def _load_task_history() -> list[dict]:
    if not TASK_HISTORY_PATH.exists():
        return []
    try:
        data = json.loads(TASK_HISTORY_PATH.read_text(encoding="utf-8"))
        return data if isinstance(data, list) else []
    except Exception:
        return []


def _save_task_history(records: list[dict]) -> None:
    TASK_HISTORY_PATH.write_text(
        json.dumps(records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _append_task_history(item: dict) -> None:
    with TASK_HISTORY_LOCK:
        records = _load_task_history()
        records.insert(0, item)
        if len(records) > 100:
            del records[100:]
        _save_task_history(records)


def task_history_data() -> dict:
    with TASK_HISTORY_LOCK:
        records = _load_task_history()
    return {"history": records, "count": len(records)}


def clear_task_history() -> dict:
    with TASK_HISTORY_LOCK:
        try:
            _save_task_history([])
            return {"cleared": True, "count": 0}
        except Exception as exc:
            return {"cleared": False, "error": str(exc)}


def _notify_finished(
    info: dict,
    page: str,
    page_label: str,
    title: str,
    status: str,
    stage: str,
    error: str,
    exit_code: int | None = None,
) -> None:
    if info.get("notified"):
        return
    info["notified"] = True
    item = {
        "page": page,
        "page_label": page_label,
        "job": info.get("job_id", ""),
        "title": title,
        "status": status,
        "stage": stage,
        "error": error,
        "exit_code": exit_code,
        "finished_at": time.time(),
    }
    started_at = float(info.get("started") or time.time())
    item["started_at"] = started_at
    item["elapsed"] = max(0, int(time.time() - started_at))
    with NOTIFY_LOCK:
        FINISHED_NOTIFICATIONS.append(item)
        if len(FINISHED_NOTIFICATIONS) > 50:
            del FINISHED_NOTIFICATIONS[:-50]
    try:
        _append_task_history(item)
    except Exception:
        pass


def running_tasks_data() -> dict:
    tasks = []
    now = time.time()
    for job_id, info in list(JOBS.items()):
        status = _single_status(info)
        if not (status.get("running") or status.get("queued") or status.get("paused")):
            continue
        state = (
            "queued" if status.get("queued")
            else "paused" if status.get("paused")
            else "running"
        )
        marker_dir = info["out"] / "results" / ".stages"
        progress = 0 if state == "queued" else _marker_progress(marker_dir, 9)
        stage_label = (
            "排队中" if state == "queued"
            else "已暂停" if state == "paused"
            else _current_stage(marker_dir, SINGLE_STAGE_LABELS) or "准备中"
        )
        log_text = "\n".join(
            [
                _log_tail(info["log"]),
                _log_tail(info["out"] / "logs" / "pipeline_r.log"),
            ]
        ).strip()
        started = float(info.get("started") or now)
        tasks.append(
            {
                "page": "single",
                "page_label": "表达分析",
                "job": job_id,
                "url": f"/?job={job_id}",
                "title": (
                    f"{info.get('accession', 'GSE')} · "
                    f"{info.get('species', '')} 表达分析"
                ),
                "detail": str(info["out"]),
                "status": state,
                "started": started,
                "elapsed": int(now - started),
                "progress": progress,
                "stage_label": stage_label,
                "log_tail": log_text,
            }
        )

    for job_id, info in list(DOCK_JOBS.items()):
        status = _dock_status(info)
        if not (status.get("running") or status.get("queued") or status.get("paused")):
            continue
        state = (
            "queued" if status.get("queued")
            else "paused" if status.get("paused")
            else "running"
        )
        marker_dir = info["output_dir"] / ".stages"
        progress = 0 if state == "queued" else _marker_progress(marker_dir, 6)
        stage_label = (
            "排队中" if state == "queued"
            else "已暂停" if state == "paused"
            else _current_stage(marker_dir, DOCK_STAGE_LABELS) or "准备中"
        )
        started = float(info.get("started") or now)
        tasks.append(
            {
                "page": "dock",
                "page_label": "虚拟筛选",
                "job": job_id,
                "url": f"/dock?job={job_id}",
                "title": f"{info.get('stage', 'pipeline')} · 虚拟筛选",
                "detail": str(info["workdir"]),
                "status": state,
                "started": started,
                "elapsed": int(now - started),
                "progress": progress,
                "stage_label": stage_label,
                "log_tail": _log_tail(info["log"]),
            }
        )

    for job_id, info in list(FULL_JOBS.items()):
        status = _full_status(info)
        if not (status.get("running") or status.get("queued") or status.get("paused")):
            continue
        state = (
            "queued" if status.get("queued")
            else "paused" if status.get("paused")
            else "running"
        )
        marker_dir = info["workdir"] / "outputs" / "integration" / ".stages"
        progress = 0 if state == "queued" else _marker_progress(marker_dir, 7)
        stage_label = (
            "排队中" if state == "queued"
            else "已暂停" if state == "paused"
            else _full_stage_label(info, marker_dir) or "准备中"
        )
        started = float(info.get("started") or now)
        output = info.get("output") or "-"
        tasks.append(
            {
                "page": "full",
                "page_label": "全自动流水线",
                "job": job_id,
                "url": f"/full?job={job_id}",
                "title": (
                    f"{info.get('accession') or '全自动流水线'} · "
                    "全自动流水线"
                ),
                "detail": f"{info['workdir']} | output: {output}",
                "status": state,
                "started": started,
                "elapsed": int(now - started),
                "progress": progress,
                "stage_label": stage_label,
                "log_tail": "\n".join(
                    _log_tail(path, 4000) for path in _full_log_paths(info)
                ).strip(),
            }
        )

    tasks.sort(key=lambda item: float(item.get("started") or 0))
    return {"tasks": tasks}


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


RESULT_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".svg"}
RESULT_DATA_SUFFIXES = {".csv", ".xlsx", ".rds"}
RESULT_FILE_SUFFIXES = RESULT_IMAGE_SUFFIXES | RESULT_DATA_SUFFIXES


def _is_result_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in RESULT_FILE_SUFFIXES


def _list_result_files(root: Path) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    return sorted(
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if _is_result_file(p)
    )


def _list_result_images(root: Path) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    return sorted(
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and p.suffix.lower() in RESULT_IMAGE_SUFFIXES
    )


def _single_cell_root_from_workdir(workdir: Path) -> Path | None:
    context = _read_json(
        workdir / "outputs" / "integration" / ".stages" / "run_context.json"
    )
    value = context.get("single_cell_root")
    if not value:
        return None
    return Path(str(value)).expanduser().resolve()


def _full_result_files(workdir: Path) -> list[str]:
    files: list[str] = []
    roots: list[tuple[Path, str]] = [
        (workdir / "outputs" / "integration", ""),
        (workdir / "outputs" / "run_001" / "results", "outputs/run_001/results"),
        (workdir / "outputs" / "run_001" / "docked", "outputs/run_001/docked"),
        (
            workdir / "outputs" / "run_001" / "network_toxicology",
            "outputs/run_001/network_toxicology",
        ),
        (
            workdir / "outputs" / "run_001" / "faers",
            "outputs/run_001/faers",
        ),
    ]
    single_cell_root = _single_cell_root_from_workdir(workdir)
    if single_cell_root:
        roots.append(
            (
                single_cell_root / "results" / "figures",
                "single_cell/results/figures",
            )
        )
        roots.append(
            (
                single_cell_root / "results" / "data",
                "single_cell/results/data",
            )
        )
    targets = workdir / "work"
    if targets.exists():
        for gene_dir in sorted(
            p for p in targets.iterdir() if p.is_dir()
        ):
            base = f"work/{gene_dir.name}/outputs/run_001"
            roots.append(
                (gene_dir / "outputs" / "run_001" / "docked", f"{base}/docked")
            )
            roots.append(
                (gene_dir / "outputs" / "run_001" / "results", f"{base}/results")
            )
    for root, prefix in roots:
        for rel in _list_result_files(root):
            files.append(f"{prefix}/{rel}" if prefix else rel)
    return sorted(set(files))


def full_results(workdir: Path) -> dict:
    workdir = Path(workdir).expanduser().resolve()
    out = workdir / "outputs" / "integration"
    result = {
        "workdir": str(workdir),
        "exists": out.exists(),
        "files": [],
        "summary": {},
        "qc_metrics": {},
        "differential_abundance": [],
        "key_genes": [],
        "knockout": [],
        "docking": [],
        "cell_feedback": [],
        "evidence": [],
    }
    if not out.exists():
        return result
    result["summary"] = _read_json(out / "integration_summary.json")
    result["qc_metrics"] = _read_json(out / "qc_metrics.json")
    result["files"] = _full_result_files(workdir)
    try:
        result["key_genes"] = json.loads(
            pd_read_csv(out / "key_genes.csv").to_json(orient="records")
        )
    except Exception:
        result["key_genes"] = []
    try:
        result["knockout"] = json.loads(
            pd_read_csv(
                workdir
                / "outputs"
                / "run_001"
                / "results"
                / "04_knockout"
                / "data"
                / "fig_52_53_ranked_knockout.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["knockout"] = []
    try:
        result["docking"] = json.loads(
            pd_read_csv(out / "docking_targets.csv").to_json(orient="records")
        )
    except Exception:
        result["docking"] = []
    try:
        result["cell_feedback"] = json.loads(
            pd_read_csv(
                out / "cell_feedback" / "data" / "feedback_targets.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["cell_feedback"] = []
    try:
        result["evidence"] = json.loads(
            pd_read_csv(out / "gene_evidence.csv").to_json(orient="records")
        )
    except Exception:
        result["evidence"] = []
    try:
        result["differential_abundance"] = json.loads(
            pd_read_csv(out / "differential_abundance.csv").to_json(
                orient="records"
            )
        )
    except Exception:
        result["differential_abundance"] = []
    return result


def pd_read_csv(path: Path):
    import pandas as pd

    return pd.read_csv(path)


def _full_file_path(workdir: Path, name: str) -> Path | None:
    workdir = Path(workdir).expanduser().resolve()
    single_cell_root = _single_cell_root_from_workdir(workdir)
    allowed_roots = [
        (workdir / "outputs" / "integration").resolve(),
        (workdir / "outputs" / "run_001" / "results").resolve(),
        (workdir / "outputs" / "run_001" / "docked").resolve(),
        (workdir / "outputs" / "run_001" / "network_toxicology").resolve(),
        (workdir / "outputs" / "run_001" / "faers").resolve(),
        (workdir / "work").resolve(),
        (workdir / "data" / "knockout").resolve(),
    ]
    if single_cell_root:
        allowed_roots.append(single_cell_root.resolve())
    name_path = Path(name)
    if name_path.is_absolute():
        return None
    if name.startswith("single_cell/") and single_cell_root:
        target = single_cell_root.joinpath(*name_path.parts[1:]).resolve()
    elif name.startswith(("outputs/", "work/", "data/")):
        target = (workdir / name_path).resolve()
    else:
        integration = (workdir / "outputs" / "integration").resolve()
        target = integration.joinpath(*name_path.parts).resolve()
        if not target.is_file() and (workdir / name_path).is_file():
            target = (workdir / name_path).resolve()
    if not _is_result_file(target):
        return None
    if not any(target.is_relative_to(root) for root in allowed_roots):
        return None
    return target


def validation_report_text() -> str:
    if VALIDATION_REPORT_PATH.exists():
        try:
            data = json.loads(
                VALIDATION_REPORT_PATH.read_text(encoding="utf-8", errors="replace")
            )
        except Exception:
            data = {}
        lines = [
            "# 随机真实 GSE 全流程验证汇总",
            "",
            f"- 请求数量：{data.get('requested', 0)}",
            f"- 通过数量：{data.get('passed', 0)}",
            f"- 随机种子：{data.get('seed', '')}",
            f"- 完成时间：{data.get('finished_at', '')}",
            "",
        ]
        for record in data.get("results", []):
            lines.append(
                f"### {record.get('accession', '')}：{record.get('status', '')}"
            )
            lines.append(f"- 耗时：{record.get('elapsed_seconds', '')} 秒")
            lines.append(f"- 工作目录：{record.get('workdir', '')}")
            lines.append("")
        return "\n".join(lines)
    return (
        "报告尚未生成。\n"
        "请在网页版点击“重新运行随机真实 GSE 验证”，"
        "或在命令行运行：python scripts\\validate_random_real_full_pipeline.py --count 10"
    )


def start_validation_job() -> dict:
    proc = VALIDATION_JOB.get("proc")
    if proc is not None and proc.poll() is None:
        return {"running": True, "message": "验证任务已在运行"}
    VALIDATION_LOG.parent.mkdir(parents=True, exist_ok=True)
    handle = VALIDATION_LOG.open("w", encoding="utf-8", errors="replace")
    proc = subprocess.Popen(
        [
            sys.executable,
            str(SCRIPTS_DIR / "validate_random_real_full_pipeline.py"),
            "--result-root",
            str(APP_ROOT.parent / "y3"),
            "--count",
            "10",
            "--seed",
            "20260813",
        ],
        cwd=APP_ROOT,
        stdout=handle,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    VALIDATION_JOB.update(
        {
            "proc": proc,
            "log": VALIDATION_LOG,
            "handle": handle,
            "started": time.time(),
        }
    )
    return {"running": True, "job": "validation"}


def validation_job_status() -> dict:
    proc = VALIDATION_JOB.get("proc")
    if proc is None:
        return {"running": False, "ok": False, "started": False, "log": ""}
    running = proc.poll() is None
    log_text = ""
    if VALIDATION_LOG.exists():
        log_text = VALIDATION_LOG.read_text(
            encoding="utf-8", errors="replace"
        )[-6000:]
    return {
        "running": running,
        "ok": not running and proc.returncode == 0,
        "started": True,
        "log": log_text,
    }


def load_dock_history() -> list[dict]:
    if not DOCK_HISTORY_PATH.exists():
        return []
    try:
        return json.loads(DOCK_HISTORY_PATH.read_text(encoding="utf-8"))
    except Exception:
        return []


def save_dock_history(records: list[dict]) -> None:
    DOCK_HISTORY_PATH.write_text(
        json.dumps(records, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def record_dock_job(info: dict, ok: bool) -> None:
    records = load_dock_history()
    records.insert(
        0,
        {
            "job": info["log"].stem.replace("web_dock_", ""),
            "stage": info.get("stage", ""),
            "workdir": str(info["workdir"]),
            "output": str(info["output_dir"]),
            "status": "success" if ok else "failed",
            "started": info["started"],
            "finished": time.time(),
        },
    )
    save_dock_history(records)


def dock_results(info: dict) -> dict:
    import csv as csv_module

    reports = info["output_dir"] / "results"
    summary = {}
    summary_path = reports / "01_analysis" / "summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            summary = {}
    rows: list[dict] = []
    ranked = reports / "01_analysis" / "data" / "fig_46_47_ranked_results.csv"
    if ranked.exists():
        with ranked.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= 200:
                    break
                rows.append(row)
    figures = _list_result_images(reports)
    files = _list_result_files(reports)
    docked_files = _list_result_files(info["output_dir"] / "docked")
    files.extend(f"docked/{rel}" for rel in docked_files)
    files = sorted(set(files))
    return {
        "summary": summary,
        "rows": rows,
        "figures": figures,
        "files": files,
        "output_dir": str(reports),
        "stage": info.get("stage", ""),
    }


def _dock_file_path(info: dict, name: str):
    out = info["output_dir"].resolve()
    results_root = (out / "results").resolve()
    docked_root = (out / "docked").resolve()
    target = (results_root / Path(name)).resolve()
    if _is_result_file(target) and results_root in target.parents:
        return target
    base = Path(name).name
    target = (docked_root / base).resolve()
    if _is_result_file(target) and target.parent == docked_root:
        return target
    return None


def run_knockout_request(data: dict) -> dict:
    from docking.config import load_config
    from docking.knockout import run_knockout

    workdir = Path(
        _first(data, "ko_workdir", str(APP_ROOT / "dock"))
    ).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    overrides = {
        "workdir": str(workdir),
        "expression_csv": _first(data, "ko_expression", "") or None,
        "metadata_csv": _first(data, "ko_metadata", "") or None,
        "depmap_csv": _first(data, "ko_depmap", "") or None,
        "ppi_network_csv": _first(data, "ko_ppi", "") or None,
        "case_label": _first(data, "ko_case", "") or None,
        "normal_label": _first(data, "ko_normal", "") or None,
        "ko_top_n": _int_field(data, "ko_top_n"),
    }
    cfg = load_config(APP_ROOT / "config" / "docking_config.json", overrides)
    import logging

    summary = run_knockout(cfg, logging.getLogger("docking.web_knockout"))
    ko_dir = cfg.knockout_dir()
    ranked = ko_dir / "data" / "fig_52_53_ranked_knockout.csv"
    rows: list[dict] = []
    if ranked.exists():
        import csv as csv_module

        limit = int(cfg.get("knockout", "top_n", 50))
        with ranked.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= limit:
                    break
                rows.append(row)
    figures = _list_result_images(ko_dir)
    files = _list_result_files(ko_dir)
    return {
        "summary": summary,
        "rows": rows,
        "figures": figures,
        "files": files,
        "output_dir": str(ko_dir),
        "workdir": str(workdir),
    }


def run_validation_request(data: dict) -> dict:
    from docking.config import load_config
    from docking.validation import export_validation

    workdir = Path(
        _first(data, "ko_workdir", str(APP_ROOT / "dock"))
    ).expanduser().resolve()
    overrides = {
        "workdir": str(workdir),
        "validation_top_n": _int_field(data, "validation_top_n"),
    }
    cfg = load_config(APP_ROOT / "config" / "docking_config.json", overrides)
    import logging

    summary = export_validation(cfg, logging.getLogger("docking.web_validation"))
    val_dir = cfg.validation_dir()
    files = _list_result_files(val_dir)
    return {
        "summary": summary,
        "files": files,
        "output_dir": str(val_dir),
        "workdir": str(workdir),
    }


def run_network_request(data: dict) -> dict:
    from docking.config import load_config
    from docking.network_toxicology import run_network_toxicology

    workdir = Path(
        _first(data, "net_workdir", str(APP_ROOT / "dock"))
    ).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    overrides = {
        "workdir": str(workdir),
        "compound_name": _first(data, "net_compound_name", "") or None,
        "disease_name": _first(data, "net_disease_name", "") or None,
        "compound_targets_csv": (
            _first(data, "net_compound_targets", "") or None
        ),
        "disease_genes_csv": _first(data, "net_disease_genes", "") or None,
        "ppi_network_csv": _first(data, "net_ppi", "") or None,
    }
    cfg = load_config(
        APP_ROOT / "config" / "docking_config.json",
        overrides,
    )
    import logging

    summary = run_network_toxicology(
        cfg,
        logging.getLogger("docking.web_network"),
    )
    out_dir = cfg._resolve(
        cfg.data.get("network_toxicology", {}).get("output_dir")
        or "outputs/run_001/network_toxicology",
        cfg.workdir,
    )
    rows: list[dict] = []
    overlap = out_dir / "data" / "compound_disease_overlap.csv"
    if overlap.exists():
        import csv as csv_module

        with overlap.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= 200:
                    break
                rows.append(row)
    return {
        "summary": summary,
        "rows": rows,
        "figures": _analysis_files(out_dir, images_only=True),
        "files": _analysis_files(out_dir),
        "output_dir": str(out_dir),
        "workdir": str(workdir),
    }


def run_faers_request(data: dict) -> dict:
    from docking.config import load_config
    from docking.signal_detection import run_faers

    workdir = Path(
        _first(data, "faers_workdir", str(APP_ROOT / "dock"))
    ).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    overrides = {
        "workdir": str(workdir),
        "faers_input": _first(data, "faers_input", "") or None,
        "faers_drug_column": _first(data, "faers_drug_column", "drug"),
        "faers_event_column": _first(data, "faers_event_column", "event"),
        "faers_count_column": _first(data, "faers_count_column", "") or None,
        "faers_min_count": _int_field(data, "faers_min_count") or 3,
    }
    cfg = load_config(
        APP_ROOT / "config" / "docking_config.json",
        overrides,
    )
    import logging

    summary = run_faers(
        cfg,
        logging.getLogger("docking.web_faers"),
    )
    out_dir = cfg._resolve(
        cfg.data.get("faers", {}).get("output_dir")
        or "outputs/run_001/faers",
        cfg.workdir,
    )
    rows: list[dict] = []
    signals = out_dir / "data" / "faers_signals.csv"
    if signals.exists():
        import csv as csv_module

        with signals.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= 200:
                    break
                rows.append(row)
    return {
        "summary": summary,
        "rows": rows,
        "figures": _analysis_files(out_dir, images_only=True),
        "files": _analysis_files(out_dir),
        "output_dir": str(out_dir),
        "workdir": str(workdir),
    }


def _analysis_files(root: Path, images_only: bool = False) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    suffixes = (
        RESULT_IMAGE_SUFFIXES
        if images_only
        else RESULT_IMAGE_SUFFIXES
        | {".csv", ".html", ".json", ".md", ".xlsx"}
    )
    return sorted(
        p.relative_to(root).as_posix()
        for p in root.rglob("*")
        if p.is_file() and p.suffix.lower() in suffixes
    )


def _analysis_file_path(workdir: str, name: str, kind: str):
    from docking.config import load_config

    work = Path(workdir).expanduser().resolve()
    cfg = load_config(
        APP_ROOT / "config" / "docking_config.json",
        {"workdir": str(work)},
    )
    if kind == "network":
        section = cfg.data.get("network_toxicology", {})
        default = "outputs/run_001/network_toxicology"
    else:
        section = cfg.data.get("faers", {})
        default = "outputs/run_001/faers"
    root = cfg._resolve(
        section.get("output_dir") or default,
        cfg.workdir,
    ).resolve()
    target = (root / Path(name)).resolve()
    allowed_suffixes = RESULT_IMAGE_SUFFIXES | {
        ".csv",
        ".html",
        ".json",
        ".md",
        ".xlsx",
    }
    if (
        target.is_file()
        and target.suffix.lower() in allowed_suffixes
        and root in target.parents
    ):
        return target
    return None


def _ko_file_path(workdir: str, name: str):
    from docking.config import load_config

    work = Path(workdir).expanduser().resolve()
    cfg = load_config(
        APP_ROOT / "config" / "docking_config.json",
        {"workdir": str(work)},
    )
    for folder in (cfg.knockout_dir(), cfg.validation_dir()):
        folder = folder.resolve()
        target = (folder / Path(name)).resolve()
        if _is_result_file(target) and folder in target.parents:
            return target
    return None


def _detect_and_save_box(workdir: str, receptor: str) -> dict:
    from docking.box import detect_box_data
    from docking.config import load_config, save_config

    workdir = Path(workdir).expanduser().resolve()
    receptor_path = Path(receptor).expanduser()
    if not receptor_path.is_absolute():
        receptor_path = workdir / receptor_path
    receptor_path = receptor_path.resolve()
    if not receptor_path.is_file():
        raise ValueError(f"receptor file not found: {receptor_path}")
    center, size, mode = detect_box_data(receptor_path)
    config_path = workdir / "config" / "docking_config.json"
    if config_path.exists():
        cfg = load_config(config_path)
    else:
        cfg = load_config(
            APP_ROOT / "config" / "docking_config.json",
            {"workdir": str(workdir)},
        )
    rel = os.path.relpath(receptor_path, workdir)
    rel_path = rel if not rel.startswith("..") else str(receptor_path)
    cfg.data["receptor"]["detect_input"] = rel_path
    cfg.data["receptor"]["input"] = rel_path
    cfg.data["receptor"]["center"] = center
    cfg.data["receptor"]["size"] = size
    save_config(cfg, config_path)
    return {
        "center": center,
        "size": size,
        "mode": mode,
        "config": str(config_path),
    }


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass

    def _send(self, code: int, body: bytes, content_type: str) -> None:
        if content_type.startswith("text/html"):
            body = _inject_heartbeat_script(body)
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        self.end_headers()
        self.wfile.write(body)

    def _handle_heartbeat(self, params: dict) -> None:
        client_id = _first(params, "client")
        if _first(params, "left", "0") == "1":
            register_heartbeat(client_id)
            unregister_heartbeat(client_id)
        else:
            register_heartbeat(client_id)
        self._send(200, b'{"ok": true}', "application/json")

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/heartbeat":
            self._handle_heartbeat(parse_qs(parsed.query))
            return
        if parsed.path == "/":
            self._send(200, get_page().encode("utf-8"), "text/html; charset=utf-8")
            return
        if parsed.path == "/datasets":
            self._send(
                200,
                render_datasets_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/dock":
            self._send(200, render_dock_page().encode("utf-8"), "text/html; charset=utf-8")
            return
        if parsed.path == "/log":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            parts = []
            if info["log"].exists():
                parts.append(info["log"].read_text(encoding="utf-8", errors="replace"))
            r_log = info["out"] / "logs" / "pipeline_r.log"
            if r_log.exists():
                parts.append(r_log.read_text(encoding="utf-8", errors="replace"))
            text = "\n".join(parts)
            self._send(200, text.encode("utf-8"), "text/plain; charset=utf-8")
            return
        if parsed.path == "/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(_single_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/files":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            fig_dir = info["out"] / "results" / "figures"
            figures = sorted(
                p.name for p in fig_dir.rglob("*.png")
            ) if fig_dir.exists() else []
            body = json.dumps({"figures": figures}).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/figure":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            name = Path(query.get("name", [""])[0]).name
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            fig_dir = (info["out"] / "results" / "figures").resolve()
            target = next(
                (p for p in fig_dir.rglob(name) if p.is_file()),
                None,
            )
            if target is None:
                self._send(404, b"figure not found", "text/plain; charset=utf-8")
                return
            suffix = target.suffix.lower()
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".pdf": "application/pdf",
            }.get(suffix, "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/report":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            output = query.get("output", [""])[0]
            info = JOBS.get(job) if job else None
            target = (
                _single_report_path(info)
                if info is not None
                else _single_report_path(None, output)
            )
            if not target:
                self._send(404, b"report not found", "text/plain; charset=utf-8")
                return
            self._send(200, target.read_bytes(), "text/html; charset=utf-8")
            return
        if parsed.path == "/datasets/file":
            query = parse_qs(parsed.query)
            name = query.get("name", [""])[0]
            target = dataset_file_path(name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".csv": "text/csv; charset=utf-8",
                ".json": "application/json",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/datasets/download/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            try:
                body = json.dumps(
                    dataset_download_status(job),
                    ensure_ascii=False,
                ).encode("utf-8")
                self._send(200, body, "application/json")
            except ValueError as exc:
                self._send(
                    404,
                    json.dumps({"error": str(exc)}).encode("utf-8"),
                    "application/json",
                )
            return
        if parsed.path == "/dock/log":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            text = (
                info["log"].read_text(encoding="utf-8", errors="replace")
                if info["log"].exists()
                else ""
            )
            self._send(200, text.encode("utf-8"), "text/plain; charset=utf-8")
            return
        if parsed.path == "/dock/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(_dock_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full":
            self._send(
                200,
                render_full_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/results":
            self._send(
                200,
                render_results_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/result-guide":
            self._send(
                200,
                json.dumps(result_guide_data(), ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
            )
            return
        if parsed.path == "/result-details":
            self._send(
                200,
                json.dumps(result_details_data(), ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
            )
            return
        if parsed.path == "/tasks":
            self._send(
                200,
                render_tasks_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/tasks/data":
            body = json.dumps(
                running_tasks_data(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/tasks/notifications":
            with NOTIFY_LOCK:
                items = FINISHED_NOTIFICATIONS[:]
                FINISHED_NOTIFICATIONS.clear()
            body = json.dumps(
                {"notifications": items},
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/tasks/history":
            body = json.dumps(
                task_history_data(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full/log":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = FULL_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            text = (
                info["log"].read_text(encoding="utf-8", errors="replace")
                if info["log"].exists()
                else ""
            )
            self._send(200, text.encode("utf-8"), "text/plain; charset=utf-8")
            return
        if parsed.path == "/full/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = FULL_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(_full_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full/results":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            workdir = query.get("workdir", [""])[0]
            if job and job in FULL_JOBS:
                workdir = str(FULL_JOBS[job]["workdir"])
            if not workdir:
                self._send(400, b"job or workdir required", "application/json")
                return
            body = json.dumps(
                full_results(workdir),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full/file":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            name = query.get("name", [""])[0]
            workdir = query.get("workdir", [""])[0]
            if job and job in FULL_JOBS:
                workdir = str(FULL_JOBS[job]["workdir"])
            if not workdir or not name:
                self._send(400, b"workdir and name required", "text/plain; charset=utf-8")
                return
            target = _full_file_path(workdir, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".pdf": "application/pdf",
                ".csv": "text/csv; charset=utf-8",
                ".json": "application/json",
                ".html": "text/html; charset=utf-8",
                ".css": "text/css; charset=utf-8",
                ".md": "text/markdown; charset=utf-8",
                ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/dock/results":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(
                dock_results(info),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/dock/file":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            name = Path(query.get("name", [""])[0])
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            target = _dock_file_path(info, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".pdf": "application/pdf",
                ".csv": "text/csv; charset=utf-8",
                ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                ".json": "application/json",
                ".pdbqt": "chemical/x-pdbqt",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/dock/history":
            body = json.dumps(
                load_dock_history(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path in ("/dock/network/file", "/dock/faers/file"):
            query = parse_qs(parsed.query)
            workdir = query.get("workdir", [""])[0]
            name = query.get("name", [""])[0]
            kind = (
                "network"
                if parsed.path == "/dock/network/file"
                else "faers"
            )
            if not workdir or not name:
                self._send(
                    400,
                    b"workdir and name required",
                    "text/plain; charset=utf-8",
                )
                return
            target = _analysis_file_path(workdir, name, kind)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".pdf": "application/pdf",
                ".csv": "text/csv; charset=utf-8",
                ".html": "text/html; charset=utf-8",
                ".json": "application/json",
                ".md": "text/markdown; charset=utf-8",
                ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/dock/knockout/file":
            query = parse_qs(parsed.query)
            workdir = query.get("workdir", [""])[0]
            name = query.get("name", [""])[0]
            target = _ko_file_path(workdir, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".jpeg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".pdf": "application/pdf",
                ".csv": "text/csv; charset=utf-8",
                ".md": "text/markdown; charset=utf-8",
                ".json": "application/json",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        if parsed.path == "/dock/validation-report":
            body = json.dumps(
                {
                    "report": validation_report_text(),
                    "exists": VALIDATION_REPORT_PATH.exists(),
                },
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/dock/validation-status":
            body = json.dumps(
                validation_job_status(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/history":
            body = json.dumps(
                load_history(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/environment":
            try:
                result = subprocess.run(
                    [
                        sys.executable,
                        str(APP_ROOT / "launchers" / "check_pipeline_environment.py"),
                    ],
                    cwd=APP_ROOT,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=180,
                )
                body = json.dumps({
                    "ok": result.returncode == 0,
                    "output": result.stdout + result.stderr,
                }, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                self._send(200, json.dumps({
                    "ok": False,
                    "output": str(exc),
                }, ensure_ascii=False).encode("utf-8"), "application/json")
            return
        if parsed.path == "/dock-environment":
            try:
                result = subprocess.run(
                    [
                        sys.executable,
                        str(APP_ROOT / "launchers" / "check_dock_environment.py"),
                    ],
                    cwd=APP_ROOT,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=180,
                )
                body = json.dumps({
                    "ok": result.returncode == 0,
                    "output": result.stdout + result.stderr,
                }, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                self._send(200, json.dumps({
                    "ok": False,
                    "output": str(exc),
                }, ensure_ascii=False).encode("utf-8"), "application/json")
            return
        if parsed.path == "/install-status":
            proc = INSTALL_JOB.get("proc")
            if not proc:
                body = json.dumps({
                    "running": False,
                    "ok": False,
                    "log": "",
                }).encode("utf-8")
                self._send(200, body, "application/json")
                return
            running = proc.poll() is None
            log_text = ""
            if INSTALL_LOG.exists():
                log_text = INSTALL_LOG.read_text(
                    encoding="utf-8",
                    errors="replace",
                )[-6000:]
            body = json.dumps({
                "running": running,
                "ok": not running and proc.returncode == 0,
                "log": log_text,
            }, ensure_ascii=False).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/versions":
            body = json.dumps(get_versions(), ensure_ascii=False).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path.startswith("/static/"):
            name = Path(parsed.path.removeprefix("/static/")).name
            target = (STATIC_DIR / name).resolve()
            if not target.is_file() or target.parent != STATIC_DIR.resolve():
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            content_type = {
                ".png": "image/png",
                ".jpg": "image/jpeg",
                ".svg": "image/svg+xml",
                ".css": "text/css; charset=utf-8",
                ".json": "application/json; charset=utf-8",
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self):
        parsed = urlparse(self.path)
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8", "replace")
        data = parse_qs(body)

        if parsed.path == "/heartbeat":
            params = parse_qs(parsed.query)
            params.update(data)
            self._handle_heartbeat(params)
            return

        if parsed.path == "/pause":
            job = data.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            info["out"].mkdir(parents=True, exist_ok=True)
            (info["out"] / "pause_request.flag").write_text(
                "pause", encoding="utf-8"
            )
            self._send(200, json.dumps({"paused_requested": True}).encode("utf-8"), "application/json")
            return

        if parsed.path == "/resume":
            job = data.get("job", [""])[0]
            try:
                result = resume_job(job)
                self._send(200, json.dumps(result).encode("utf-8"), "application/json")
            except Exception as exc:
                self._send(400, json.dumps({"error": str(exc)}).encode("utf-8"), "application/json")
            return

        if parsed.path == "/dock/pause":
            job = data.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            info["workdir"].mkdir(parents=True, exist_ok=True)
            (info["workdir"] / "pause_request.flag").write_text(
                "pause",
                encoding="utf-8",
            )
            proc = info.get("proc")
            if proc and proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
            info["paused"] = True
            self._send(
                200,
                json.dumps({"paused": True}).encode("utf-8"),
                "application/json",
            )
            return

        if parsed.path == "/dock/resume":
            job = data.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            flag = info["workdir"] / "pause_request.flag"
            if flag.exists():
                flag.unlink()
            info["paused"] = False
            info["notified"] = False
            log_handle = info["log"].open("a", encoding="utf-8", errors="replace")
            log_handle.write("\n[pipeline] resume requested\n")
            log_handle.flush()
            proc = subprocess.Popen(
                info["cmd"],
                cwd=APP_ROOT,
                env=info.get("env", os.environ.copy()),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            info["proc"] = proc
            info["started"] = time.time()
            self._send(
                200,
                json.dumps({"resumed": True}).encode("utf-8"),
                "application/json",
            )
            return

        if parsed.path == "/full/pause":
            job = data.get("job", [""])[0]
            info = FULL_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            proc = info.get("proc")
            if proc and proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
            info["paused"] = True
            self._send(
                200,
                json.dumps({"paused": True}).encode("utf-8"),
                "application/json",
            )
            return

        if parsed.path == "/full/resume":
            job = data.get("job", [""])[0]
            info = FULL_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            info["paused"] = False
            info["notified"] = False
            log_handle = info["log"].open("a", encoding="utf-8", errors="replace")
            log_handle.write("\n[full pipeline] resume requested\n")
            log_handle.flush()
            proc = subprocess.Popen(
                info["cmd"],
                cwd=APP_ROOT,
                env=info.get("env", os.environ.copy()),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            info["proc"] = proc
            info["started"] = time.time()
            self._send(
                200,
                json.dumps({"resumed": True}).encode("utf-8"),
                "application/json",
            )
            return

        if parsed.path == "/tasks/history/clear":
            body = json.dumps(
                clear_task_history(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return

        if parsed.path == "/full/start":
            try:
                result = start_full_job(data)
                self._send(
                    200,
                    json.dumps(result).encode("utf-8"),
                    "application/json",
                )
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}).encode("utf-8"),
                    "application/json",
                )
            return

        if parsed.path == "/datasets/search":
            try:
                result = dataset_search_request(data)
                self._send(
                    200,
                    json.dumps(result, ensure_ascii=False).encode("utf-8"),
                    "application/json",
                )
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}, ensure_ascii=False).encode(
                        "utf-8"
                    ),
                    "application/json",
                )
            return

        if parsed.path == "/datasets/download":
            try:
                result = start_dataset_download(data)
                self._send(
                    200,
                    json.dumps(result, ensure_ascii=False).encode("utf-8"),
                    "application/json",
                )
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}, ensure_ascii=False).encode(
                        "utf-8"
                    ),
                    "application/json",
                )
            return

        if parsed.path == "/dock/check-env":
            from docking.environment import check_environment

            checks = check_environment()
            self._send(
                200,
                json.dumps({"checks": checks}, ensure_ascii=False).encode("utf-8"),
                "application/json",
            )
            return

        if parsed.path == "/dock/install":
            try:
                proc = subprocess.run(
                    [
                        sys.executable,
                        str(APP_ROOT / "launchers" / "install_dock_dependencies.py"),
                    ],
                    cwd=APP_ROOT,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=600,
                )
                output = (proc.stdout or "") + (proc.stderr or "")
                body = json.dumps(
                    {
                        "ok": proc.returncode == 0,
                        "output": output[-8000:],
                    },
                    ensure_ascii=False,
                ).encode("utf-8")
                self._send(200, body, "application/json")
            except subprocess.TimeoutExpired:
                self._send(
                    200,
                    json.dumps(
                        {"ok": False, "output": "install timed out after 600s"}
                    ).encode("utf-8"),
                    "application/json",
                )
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}).encode("utf-8"),
                    "application/json",
                )
            return

        if parsed.path == "/dock/detect-box":
            try:
                workdir = _first(data, "workdir", "")
                receptor = _first(data, "receptor", "")
                if not workdir or not receptor:
                    raise ValueError("workdir and receptor are required")
                result = _detect_and_save_box(workdir, receptor)
                self._send(
                    200,
                    json.dumps(result, ensure_ascii=False).encode("utf-8"),
                    "application/json",
                )
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}).encode("utf-8"),
                    "application/json",
                )
            return

        if parsed.path == "/dock/network":
            try:
                result = run_network_request(data)
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/dock/faers":
            try:
                result = run_faers_request(data)
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/dock/knockout":
            try:
                result = run_knockout_request(data)
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/dock/knockout/validate":
            try:
                result = run_validation_request(data)
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/dock/validation/run":
            try:
                result = start_validation_job()
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/dock/start":
            try:
                result = start_dock_job(data)
                self._send(200, json.dumps(result).encode("utf-8"), "application/json")
            except Exception as exc:
                self._send(
                    400,
                    json.dumps({"error": str(exc)}).encode("utf-8"),
                    "application/json",
                )
            return

        if parsed.path == "/install":
            if INSTALL_JOB.get("proc") and INSTALL_JOB["proc"].poll() is None:
                self._send(200, json.dumps({
                    "running": True,
                    "message": "环境补全已在运行",
                }).encode("utf-8"), "application/json")
                return
            project = _first(data, "project", "single")
            target = _first(data, "target", "")
            script_name = (
                "install_dock_dependencies.py"
                if project == "dock"
                else "install_pipeline_dependencies.py"
            )
            cmd = [
                sys.executable,
                str(APP_ROOT / "launchers" / script_name),
            ]
            if target:
                cmd += ["--target", target]
            log_handle = INSTALL_LOG.open("w", encoding="utf-8", errors="replace")
            proc = subprocess.Popen(
                cmd,
                cwd=APP_ROOT,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            INSTALL_JOB["proc"] = proc
            INSTALL_JOB["project"] = project
            self._send(200, json.dumps({
                "running": True,
                "message": f"{project} 环境自动补全已启动",
            }).encode("utf-8"), "application/json")
            return

        if parsed.path != "/start":
            self._send(404, b"not found", "text/plain; charset=utf-8")
            return

        accession = data.get("accession", [""])[0]
        output = data.get("output", [""])[0]
        species = data.get("species", ["hs"])[0]
        skip_figs = [
            name for name in FIGURE_NAMES
            if data.get(f"fig_{name}", ["yes"])[0].strip().lower() != "yes"
        ]
        figure_styles = {}
        for item in FIGURES:
            if "styles" not in item:
                continue
            style = data.get(
                f"style_{item['file']}", [item["styles"][0]]
            )[0]
            if style in item["styles"]:
                figure_styles[item["file"]] = style
        params = {}
        for key in [
            "LIVER_QC_MIN_FEATURES",
            "LIVER_QC_MAX_FEATURES",
            "LIVER_QC_MIN_COUNTS",
            "LIVER_QC_MAX_COUNTS",
            "LIVER_QC_MAX_MT",
            "LIVER_CLUSTER_RESOLUTION",
            "LIVER_CLUSTER_ALGORITHM",
            "LIVER_DE_LOGFc",
            "LIVER_DE_PADJ",
            "LIVER_DE_VIOLIN_TOP_N",
            "LIVER_DE_VIOLIN_MAX_CELLS",
            "LIVER_ML_MODEL",
        ]:
            if data.get(key):
                params[key] = data[key][0]

        try:
            result = start_job(
                accession,
                output,
                species,
                skip_figs,
                figure_styles,
                params,
            )
            self._send(200, json.dumps(result).encode("utf-8"), "application/json")
        except Exception as exc:
            self._send(
                400,
                json.dumps({"error": str(exc)}).encode("utf-8"),
                "application/json",
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Local web UI for pipeline")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--page",
        choices=["single", "dock", "full", "results", "tasks", "datasets"],
        default="full",
        help="page to open in the browser",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="start the server without opening a browser window",
    )
    args = parser.parse_args()

    if not _cleanup_stale_web_ui(args.host, args.port):
        print(f"ERROR: port {args.port} is still in use by another process.")
        return 1

    INDEX_PATH.write_text(render_page(), encoding="utf-8")
    # On Windows, SO_REUSEADDR allows a second instance to bind the same port
    # and steal incoming connections, which surfaces as "connection refused".
    ThreadingHTTPServer.allow_reuse_address = False
    try:
        server = ThreadingHTTPServer((args.host, args.port), Handler)
    except OSError as exc:
        print(f"ERROR: cannot bind {args.host}:{args.port}: {exc}")
        return 1
    url = f"http://{args.host}:{args.port}"
    print(f"Web UI started: {url}")
    threading.Thread(
        target=_run_idle_shutdown_monitor,
        args=(server, time.monotonic()),
        daemon=True,
        name="web-ui-idle-shutdown",
    ).start()
    if args.page == "dock":
        open_url = url + "/dock"
    elif args.page == "full":
        open_url = url + "/full"
    elif args.page == "results":
        open_url = url + "/results"
    elif args.page == "tasks":
        open_url = url + "/tasks"
    elif args.page == "datasets":
        open_url = url + "/datasets"
    else:
        open_url = url
    if not args.no_browser:
        threading.Timer(1.0, lambda: webbrowser.open(open_url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping web UI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
