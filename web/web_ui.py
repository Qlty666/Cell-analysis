#!/usr/bin/env python3
"""Local HTML web UI for the single-cell pipeline."""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

WEB_DIR = Path(__file__).resolve().parent
APP_ROOT = WEB_DIR.parent
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

SRC_DIR = APP_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

DOCK_TEMPLATE_PATH = TEMPLATE_DIR / "dock_page_template.html"
DOCK_JOBS = {}
DOCK_QUEUE = []
DOCK_QUEUE_LOCK = threading.Lock()
DOCK_HISTORY_PATH = WEB_DIR / "dock_history.json"
NAV_HTML = (
    '<div class="topnav">'
    '<a href="/">单细胞分析</a>'
    '<a href="/dock">虚拟筛选</a>'
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
)

FIGURES = [
    {"file": "fig_01_qc_raw_violin.png", "label": "QC 小提琴图（原始）"},
    {"file": "fig_01_qc_filtered_violin.png", "label": "QC 小提琴图（过滤后）"},
    {"file": "fig_02_doublet_scores.png", "label": "双细胞得分图"},
    {"file": "fig_03_umap_clusters.png", "label": "UMAP 聚类图"},
    {"file": "fig_04_umap_condition.png", "label": "UMAP 分组图"},
    {"file": "fig_05_umap_annotation.png", "label": "UMAP 注释图"},
    {"file": "fig_06_dotplot_markers.png", "label": "Marker 基因 DotPlot"},
    {"file": "fig_07_annotation_confusion_heatmap.png", "label": "注释混淆矩阵热图"},
    {"file": "fig_08_volcano.png", "label": "差异表达图", "styles": ["volcano", "maplot"]},
    {"file": "fig_09_deg_heatmap.png", "label": "Top DEG 热图"},
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
    {"file": "fig_22_go_network.png", "label": "GO BP 通路网络图", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_23_kegg_network.png", "label": "KEGG 通路网络图", "styles": ["cnetplot", "emapplot"]},
    {"file": "fig_24_ml_feature_importance.png", "label": "ML 特征重要性图"},
    {"file": "fig_25_ml_shap.png", "label": "SHAP 可解释性图"},
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
    if not re.fullmatch(r"GSE\d+", acc):
        raise ValueError("GSE accession must look like GSE125449")
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
        str(APP_ROOT / "run_pipeline.py"),
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
        figure_count = len(list(fig_dir.glob("*")))
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
    info["started"] = time.time()
    return {"job": job_id, "resumed": True}


PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>单细胞分析流水线</title>
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
</style>
</head>
<body>
<div class="wrap">
<h1>单细胞分析流水线</h1>
<div class="card">
  <form id="form">
    <label for="acc">GSE 数据集编号</label>
    <input id="acc" name="accession" placeholder="GSE125449" required>

    <label for="out">结果保存路径</label>
    <input id="out" name="output" placeholder="C:\\results\\out" required>

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
let currentJob = null;
let pollTimer = null;

async function startRun() {
  const form = document.getElementById('form');
  const data = new URLSearchParams(new FormData(form));
  const btn = document.getElementById('startBtn');
  const msg = document.getElementById('message');
  btn.disabled = true;
  msg.className = '';
  msg.textContent = '正在启动...';

  try {
    const resp = await fetch('/start', {method: 'POST', body: data});
    const result = await resp.json();
    if (!resp.ok) throw new Error(result.error || '启动失败');
    currentJob = result.job;
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
    if (!status.running) {
      clearInterval(pollTimer);
      const msg = document.getElementById('message');
      msg.className = status.ok ? 'ok' : 'error';
      msg.textContent = status.ok ? '流水线已完成' : '流水线运行失败，请查看日志';
      if (status.ok) await loadResults(currentJob);
      currentJob = null;
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
    gallery.innerHTML = '<p class="error">结果图加载失败：' + String(e.message || e) + '</p>';
  }
}
</script>
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
    }
    cfg = load_config(APP_ROOT / "config" / "docking_config.json", overrides)
    save_config(cfg, cfg_path)

    log_path = workdir / "logs" / f"web_dock_{job_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(APP_ROOT / "run_docking.py"),
        stage,
        "--config",
        str(cfg_path),
    ]
    if _first(data, "force", "") in ("1", "true", "on", "yes"):
        cmd.append("--force")

    env = os.environ.copy()
    env["DOCK_WORKDIR"] = str(workdir)
    DOCK_JOBS[job_id] = {
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
        return {"running": False, "ok": False, "queued": True, "paused": False}
    running = info["proc"].poll() is None
    paused = bool(info.get("paused"))
    ok = False
    if not running and not paused:
        ok = info["proc"].returncode == 0
        if not info.get("recorded"):
            record_dock_job(info, ok)
            info["recorded"] = True
        _drain_dock_queue()
    return {"running": running, "ok": ok, "queued": False, "paused": paused}


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

    reports = info["output_dir"] / "reports"
    summary = {}
    summary_path = reports / "summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            summary = {}
    rows: list[dict] = []
    ranked = reports / "ranked_results.csv"
    if ranked.exists():
        with ranked.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= 200:
                    break
                rows.append(row)
    figures = (
        sorted(p.name for p in reports.glob("*.png"))
        if reports.exists()
        else []
    )
    files = (
        sorted(p.name for p in reports.iterdir() if p.is_file())
        if reports.exists()
        else []
    )
    return {
        "summary": summary,
        "rows": rows,
        "figures": figures,
        "files": files,
        "output_dir": str(info["output_dir"]),
        "stage": info.get("stage", ""),
    }


def _dock_file_path(info: dict, name: str):
    out = info["output_dir"].resolve()
    base = Path(name).name
    for folder in [(out / "reports").resolve(), (out / "docked").resolve()]:
        target = (folder / base).resolve()
        if target.is_file() and target.parent == folder:
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
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self._send(200, get_page().encode("utf-8"), "text/html; charset=utf-8")
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
            queued = info.get("proc") is None
            running = False if queued else info["proc"].poll() is None
            ok = False
            paused = False
            if queued:
                body = json.dumps(
                    {"running": False, "ok": False, "paused": False, "queued": True}
                ).encode("utf-8")
                self._send(200, body, "application/json")
                return
            if not running:
                ok = info["proc"].returncode == 0
                paused = info["proc"].returncode == 98
                if not info.get("recorded"):
                    record_job(info, ok)
                    info["recorded"] = True
                _drain_queue()
            body = json.dumps(
                {"running": running, "ok": ok, "paused": paused}
            ).encode("utf-8")
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
                p.name for p in fig_dir.glob("*") if p.is_file()
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
            target = (fig_dir / name).resolve()
            if not target.is_file() or target.parent != fig_dir:
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
            name = Path(query.get("name", [""])[0]).name
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
            }.get(target.suffix.lower(), "application/octet-stream")
            self._send(200, target.read_bytes(), content_type)
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self):
        parsed = urlparse(self.path)
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8", "replace")
        data = parse_qs(body)

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
        choices=["single", "dock"],
        default="single",
        help="page to open in the browser",
    )
    args = parser.parse_args()

    INDEX_PATH.write_text(render_page(), encoding="utf-8")
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{args.host}:{args.port}"
    print(f"Web UI started: {url}")
    open_url = url + ("/dock" if args.page == "dock" else "")
    threading.Timer(1.0, lambda: webbrowser.open(open_url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping web UI")
    return 0


if __name__ == "__main__":
    sys.exit(main())
