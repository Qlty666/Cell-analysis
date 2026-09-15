#!/usr/bin/env python3
"""Local HTML web UI for the single-cell pipeline."""

import argparse
import json
import logging
import os
import re
import secrets
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
TEMPLATE_DIR = WEB_DIR / "templates"
STATIC_DIR = WEB_DIR / "static"
PAGE_TEMPLATE_PATH = TEMPLATE_DIR / "web_page_template.html"
GUIDE_TEMPLATE_PATH = TEMPLATE_DIR / "guide_page_template.html"
ENVIRONMENT_TEMPLATE_PATH = TEMPLATE_DIR / "environment_page_template.html"
INSTALL_LOG = WEB_DIR / "install_log.txt"

log = logging.getLogger("web_ui")

HEARTBEAT_LAST_SEEN_AT: float | None = None

SRC_DIR = APP_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from web_state import (  # noqa: E402
    JOB_STORE_MAX_RECORDS,
    JOBS,
    QUEUE,
    QUEUE_LOCK,
    HISTORY_PATH,
    HISTORY_LOCK,
    INSTALL_JOB,
    FINISHED_NOTIFICATIONS,
    NOTIFY_LOCK,
    TASK_HISTORY_PATH,
    TASK_HISTORY_LOCK,
    JOB_RECORD_LOCK,
    HEARTBEAT_CLIENTS,
    HEARTBEAT_LOCK,
    HEARTBEAT_INTERVAL_SECONDS,
    HEARTBEAT_IDLE_TIMEOUT_SECONDS,
    HEARTBEAT_START_GRACE_SECONDS,
    HEARTBEAT_SHUTDOWN_GRACE_SECONDS,
    DOCK_JOBS,
    DOCK_QUEUE,
    DOCK_QUEUE_LOCK,
    DOCK_HISTORY_PATH,
    DOCK_HISTORY_LOCK,
    MOLECULAR_DOCK_JOBS,
    MOLECULAR_DOCK_QUEUE,
    MOLECULAR_DOCK_QUEUE_LOCK,
    MOLECULAR_DOCK_HISTORY_PATH,
    MOLECULAR_DOCK_HISTORY_LOCK,
    FULL_JOBS,
    FULL_QUEUE,
    FULL_QUEUE_LOCK,
    ANALYSIS_JOBS,
    ANALYSIS_QUEUE,
    DATASET_DOWNLOAD_JOBS,
    DATASET_DOWNLOAD_LOCK,
    VALIDATION_JOB,
    _has_active_jobs,
    _prune_job_stores,
    _read_history_file,
    _write_json_atomic,
    _write_history_file,
    load_history,
    save_history,
    load_dock_history,
    save_dock_history,
    load_molecular_docking_history,
    save_molecular_docking_history,
    _load_task_history,
    _save_task_history,
    _spawn_process,
    _drain_store,
)
from web_results import (  # noqa: E402
    _marker_progress,
    _current_stage,
    _log_tail,
    _stage_from_log,
    _full_stage_from_log,
    _full_stage_label,
    _extract_error,
    _finished_info,
    _append_task_history,
    task_history_data,
    clear_task_history,
    _notify_finished,
    record_job,
    record_dock_job,
    record_molecular_docking_job,
    record_analysis_job,
    _full_log_paths,
    _single_cell_root_from_workdir,
    _read_json,
)
from common.env import require_rscript  # noqa: E402

DOCK_TEMPLATE_PATH = TEMPLATE_DIR / "dock_page_template.html"
MD_TEMPLATE_PATH = TEMPLATE_DIR / "md_simulation_page_template.html"
KNOCKOUT_TEMPLATE_PATH = TEMPLATE_DIR / "knockout_page_template.html"
NETWORK_TEMPLATE_PATH = TEMPLATE_DIR / "network_page_template.html"
FAERS_TEMPLATE_PATH = TEMPLATE_DIR / "faers_page_template.html"
VALIDATION_TEMPLATE_PATH = TEMPLATE_DIR / "validation_page_template.html"
MOLECULAR_DOCK_TEMPLATE_PATH = TEMPLATE_DIR / "molecular_docking_template.html"
FULL_TEMPLATE_PATH = TEMPLATE_DIR / "full_page_template.html"
RESULTS_TEMPLATE_PATH = TEMPLATE_DIR / "results_manifest_optimized.html"
RESULT_GUIDE_PATH = APP_ROOT / "docs" / "result_figure_guide.md"
RESULT_DETAILS_PATH = STATIC_DIR / "result_details.json"
TASKS_TEMPLATE_PATH = TEMPLATE_DIR / "tasks_template.html"
DATASET_TEMPLATE_PATH = TEMPLATE_DIR / "datasets_template.html"
ANALYSIS_TEMPLATE_PATH = TEMPLATE_DIR / "analysis_tools_page_template.html"
DATASET_SEARCH_DIR = APP_ROOT / "data_cache" / "dataset_search"
DATASET_DATABASES = ("geo", "biostudies", "atlas")
MAX_POST_BODY_BYTES = 1_000_000
MAX_SERVED_FILE_BYTES = 64 * 1024 * 1024
LOCAL_ORIGIN_HOSTS = {"127.0.0.1", "localhost", "[::1]", "::1"}
LOOPBACK_HOSTS = {"127.0.0.1", "localhost", "::1", "[::1]"}

# Runtime security state. ``AUTH_TOKEN`` is only created when the server binds a
# non-loopback host; loopback usage keeps the zero-friction local workflow.
SERVER_PORT: int | None = None
AUTH_TOKEN: str | None = None

# Directories the web console may read results from. Anything else must be a
# workdir already registered by a job started through this console, or an
# explicit ``--allow-path`` given on the command line.
DEFAULT_WORKDIR_ROOTS: tuple[Path, ...] = tuple(
    Path(item).expanduser().resolve()
    for item in (
        os.environ.get("LIVER_OUTPUT_ROOT")
        or str(APP_ROOT.parent / "liver_cancer"),
        os.environ.get("LIVER_VALIDATION_ROOT")
        or str(APP_ROOT / "data_cache"),
        str(APP_ROOT / "dock"),
        str(APP_ROOT / "molecular_docking"),
        str(APP_ROOT / "data_cache"),
    )
)
EXTRA_WORKDIR_ROOTS: list[Path] = []

# One table for every result/static file the console serves (previously copied
# inline in eight request handlers).
CONTENT_TYPES = {
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".svg": "image/svg+xml",
    ".pdf": "application/pdf",
    ".csv": "text/csv; charset=utf-8",
    ".tsv": "text/tab-separated-values; charset=utf-8",
    ".json": "application/json",
    ".html": "text/html; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".md": "text/markdown; charset=utf-8",
    ".mtx": "text/plain; charset=utf-8",
    ".xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    ".pdbqt": "chemical/x-pdbqt",
}


def _content_type(suffix: str) -> str:
    return CONTENT_TYPES.get((suffix or "").lower(), "application/octet-stream")


from web_data import (  # noqa: E402
    NAV_HTML,
    NAV_CSS,
    ENV_MODULES,
    FIGURES,
    FIGURE_NAMES,
    STYLE_LABELS,
    SOFTWARE,
    SINGLE_STAGE_LABELS,
    DOCK_STAGE_LABELS,
    MOLECULAR_DOCK_STAGE_LABELS,
    FULL_STAGE_LABELS,
    RESULT_IMAGE_SUFFIXES,
)
from web_analysis import (  # noqa: E402
    _advanced_analysis_file_path,
    _analysis_status,
    _drain_analysis_queue,
    analysis_results,
    start_analysis_job,
)
from web_files import (  # noqa: E402
    analysis_files as _analysis_files,
    is_result_file as _is_result_file,
    list_result_files as _list_result_files,
    list_result_images as _list_result_images,
)
from web_utils import (  # noqa: E402
    cli_path,
    first_value as _first,
    float3 as _float3,
    float_field as _float_field,
    integer_field as _int_field,
    integer_list_field as _int_list_field,
    raw_count_flag as _raw_count_flag,
)
from web_validation import (  # noqa: E402
    _validation_report_path,
    start_validation_job,
    validation_job_status,
    validation_report_text,
)


def _cli_path(value: str) -> str:
    return cli_path(value, APP_ROOT)


def environment_module_cards() -> str:
    """Return the per-module environment cards for the web board."""
    cards = []
    for name, meta in ENV_MODULES.items():
        badges = []
        if meta["r_deps"]:
            badges.append('<span class="badge env-r">R 包</span>')
        if meta["dock_tools"]:
            badges.append('<span class="badge env-dock">对接工具</span>')
        if meta["skills"]:
            badges.append('<span class="badge env-skills">Codex Skills</span>')
        badge_html = (
            '<div class="env-badges">' + "".join(badges) + "</div>"
            if badges
            else ""
        )
        install_bat = meta.get("install_bat")
        check_bat = meta.get("check_bat")
        if install_bat and check_bat:
            install_cmd = install_bat.replace("/", "\\")
            check_cmd = meta["check_bat"].replace("/", "\\")
            launcher_html = (
                '<div class="env-cmd">一键补全：<code>'
                + install_cmd
                + "</code>　检查：<code>"
                + check_cmd
                + "</code></div>"
            )
        elif install_bat:
            install_cmd = install_bat.replace("/", "\\")
            launcher_html = (
                '<div class="env-cmd">一键补全：<code>'
                + install_cmd
                + "</code>　检查：统一命令</div>"
            )
        else:
            launcher_html = (
                '<div class="env-cmd">此模块不自动安装外部软件，'
                "请在环境检查结果和版本文档中确认。</div>"
            )
        install_button = (
            '<button type="button" class="secondary" data-action="install" data-name="'
            + name
            + '">一键补全</button>'
            if install_bat
            else ""
        )
        cards.append(
            '<section class="card env-module" data-module="'
            + name
            + '">'
            + '<div class="env-head">'
            + "<h2>"
            + meta["title"]
            + "</h2>"
            + badge_html
            + "</div>"
            + '<p class="muted">'
            + meta["summary"]
            + "</p>"
            + launcher_html
            + '<p class="env-note">'
            + meta["note"]
            + "</p>"
            + '<div class="actions env-actions">'
            + '<button type="button" data-action="check" data-name="'
            + name
            + '">检查环境</button>'
            + install_button
            + "</div>"
            + '<pre id="envLog-'
            + name
            + '" class="env-log">尚未检查此板块。</pre>'
            + "</section>"
        )
    return "\n".join(cards)

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












def _is_loopback_host(host: str) -> bool:
    return (host or "").strip().lower() in LOOPBACK_HOSTS


def _configure_runtime(host: str, port: int) -> str | None:
    """Record the bound port and mint an auth token for non-loopback binds."""
    global SERVER_PORT, AUTH_TOKEN
    SERVER_PORT = int(port)
    AUTH_TOKEN = None if _is_loopback_host(host) else secrets.token_urlsafe(24)
    return AUTH_TOKEN


def _known_workdirs() -> set[str]:
    """Workdirs registered by jobs that this console actually started."""
    known: set[str] = set()
    for store in (JOBS, DOCK_JOBS, MOLECULAR_DOCK_JOBS, FULL_JOBS):
        for info in list(store.values()):
            for key in ("out", "workdir"):
                value = info.get(key)
                if value:
                    try:
                        known.add(str(Path(value).resolve()))
                    except (OSError, ValueError):
                        continue
    return known


def _workdir_allowed(value: str | Path) -> Path | None:
    """Return the resolved workdir when it is inside an allowed root, else None."""
    try:
        candidate = Path(value).expanduser().resolve()
    except (OSError, RuntimeError, ValueError):
        return None
    if str(candidate) in _known_workdirs():
        return candidate
    for root in (*DEFAULT_WORKDIR_ROOTS, *EXTRA_WORKDIR_ROOTS):
        try:
            if candidate == root or candidate.is_relative_to(root):
                return candidate
        except (OSError, ValueError):
            continue
    return None


def _request_token(handler) -> str:
    headers = getattr(handler, "headers", None)
    if headers is not None:
        header = (headers.get("X-Auth-Token") or "").strip()
        if header:
            return header
        cookie = headers.get("Cookie") or ""
        for chunk in cookie.split(";"):
            name, _, value = chunk.strip().partition("=")
            if name == "liverbio_token" and value:
                return value.strip()
    try:
        query = parse_qs(urlparse(handler.path).query)
    except (AttributeError, ValueError):
        return ""
    return (query.get("token", [""])[0] or "").strip()


def _request_authorized(handler) -> bool:
    """Loopback requests are unauthenticated; remote binds require the token."""
    if AUTH_TOKEN is None:
        return True
    return secrets.compare_digest(_request_token(handler), AUTH_TOKEN)


def _origin_allowed(origin: str) -> bool:
    origin = (origin or "").strip()
    if not origin:
        return True
    try:
        parts = urlparse(origin)
    except Exception:
        return False
    if parts.scheme.lower() not in ("http", "https"):
        return False
    if (parts.hostname or "").lower() not in LOCAL_ORIGIN_HOSTS:
        return False
    if SERVER_PORT is not None and parts.port is not None and parts.port != SERVER_PORT:
        return False
    return True


def _fetch_site_allowed(site: str) -> bool:
    site = (site or "").strip().lower()
    return not site or site in ("same-origin", "same-site", "none")


def _inject_heartbeat_script(body: bytes) -> bytes:
    if b"</body>" not in body:
        return body
    return body.replace(b"</body>", HEARTBEAT_SCRIPT + b"</body>", 1)


def _run_idle_shutdown_monitor(
    server: ThreadingHTTPServer,
    started_at: float,
) -> None:
    idle_since: float | None = None
    last_prune = time.monotonic()
    while True:
        now = time.monotonic()
        _purge_stale_heartbeats(now)
        if now - last_prune >= 60.0:
            _prune_job_stores()
            last_prune = now
        if _heartbeat_client_ids() or _has_active_jobs():
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
    return require_rscript()


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
    try:
        from docking import __version__ as project_version
    except Exception:
        project_version = "unknown"
    out.insert(
        0,
        {
            "name": "项目版本",
            "version": project_version,
            "kind": "software",
            "url": "https://github.com/Qlty666/Cell-analysis",
            "install": "",
        },
    )
    try:
        from experiment_plan_one.coverage import audit_plan_environment

        environment = audit_plan_environment()
        vina = environment.get("vina") or {}
        gromacs = environment.get("gromacs") or {}
        mmpbsa = environment.get("gmx_mmpbsa") or {}
        cellchat = environment.get("cellchat") or {}
        packages = environment.get("python_packages") or {}
        out.extend(
            [
                {
                    "name": "Plan One Vina",
                    "version": vina.get("version") or "not installed",
                    "kind": "software",
                    "url": "https://github.com/ccsb-scripps/AutoDock-Vina",
                    "install": f"expected {vina.get('expected', '1.2.3')}",
                },
                {
                    "name": "Plan One GROMACS",
                    "version": gromacs.get("version") or "not installed",
                    "kind": "software",
                    "url": "https://www.gromacs.org/",
                    "install": f"expected {gromacs.get('expected', '2022')}",
                },
                {
                    "name": "gmx_MMPBSA",
                    "version": (
                        mmpbsa.get("version") or "not installed"
                        if mmpbsa.get("available")
                        else "not installed"
                    ),
                    "kind": "software",
                    "url": "https://valdes-tresanco-ms.github.io/gmx_MMPBSA/",
                    "install": "required for Figure 5h",
                },
                {
                    "name": "R CellChat",
                    "version": (
                        "available"
                        if cellchat.get("available")
                        else "not installed"
                    ),
                    "kind": "package",
                    "url": "https://github.com/jinworks/CellChat",
                    "install": "install.packages('CellChat')",
                },
                {
                    "name": "PoseBusters",
                    "version": (
                        "available"
                        if packages.get("posebusters")
                        else "not installed"
                    ),
                    "kind": "package",
                    "url": "https://github.com/maabuu/posebusters",
                    "install": "pip install posebusters",
                },
                {
                    "name": "Meeko",
                    "version": (
                        "available" if packages.get("meeko") else "not installed"
                    ),
                    "kind": "package",
                    "url": "https://github.com/forlilab/Meeko",
                    "install": "pip install meeko",
                },
            ]
        )
    except Exception as exc:
        log.debug("plan one environment versions unavailable: %s", exc)
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
    with QUEUE_LOCK:
        QUEUE.append(JOBS[job_id])
    _drain_queue()
    return {
        "job": job_id,
        "log_url": f"/log?job={job_id}",
        "status_url": f"/status?job={job_id}",
    }






def _start_process(info: dict) -> None:
    _spawn_process(info)


def _drain_queue() -> None:
    _drain_store(QUEUE, QUEUE_LOCK)








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
<script src="/static/nav.js" defer></script>
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


def _with_shared_nav(html: str) -> str:
    """Replace a page's local nav block with the shared navigation."""
    return re.sub(
        r'<div class="topnav">.*?</div>',
        lambda _: NAV_HTML,
        html,
        count=1,
        flags=re.S,
    )


def render_dock_page() -> str:
    if DOCK_TEMPLATE_PATH.exists():
        return _with_shared_nav(DOCK_TEMPLATE_PATH.read_text(encoding="utf-8"))
    return (
        "<html><body><h1>dock template missing</h1>"
        "<p>web/templates/dock_page_template.html not found</p></body></html>"
    )


def _render_template(path: Path, missing: str) -> str:
    if path.exists():
        return _with_shared_nav(path.read_text(encoding="utf-8"))
    return f"<html><body><h1>{missing}</h1></body></html>"


def render_md_simulation_page() -> str:
    return _render_template(
        MD_TEMPLATE_PATH,
        "md simulation template missing",
    )


def render_knockout_page() -> str:
    return _render_template(
        KNOCKOUT_TEMPLATE_PATH,
        "knockout template missing",
    )


def render_network_page() -> str:
    return _render_template(
        NETWORK_TEMPLATE_PATH,
        "network toxicology template missing",
    )


def render_faers_page() -> str:
    return _render_template(
        FAERS_TEMPLATE_PATH,
        "faers template missing",
    )


def render_validation_page() -> str:
    return _render_template(
        VALIDATION_TEMPLATE_PATH,
        "validation template missing",
    )


def render_molecular_docking_page() -> str:
    return _render_template(
        MOLECULAR_DOCK_TEMPLATE_PATH,
        "molecular docking template missing",
    )


def render_full_page() -> str:
    return _render_template(
        FULL_TEMPLATE_PATH,
        "full pipeline template missing",
    )


def render_results_page() -> str:
    return _render_template(
        RESULTS_TEMPLATE_PATH,
        "results template missing",
    )


def render_guide_page() -> str:
    return _render_template(
        GUIDE_TEMPLATE_PATH,
        "guide template missing",
    )


def render_environment_page() -> str:
    if not ENVIRONMENT_TEMPLATE_PATH.exists():
        return (
            "<html><body><h1>environment template missing</h1>"
            "<p>web/templates/environment_page_template.html not found</p>"
            "</body></html>"
        )
    html = ENVIRONMENT_TEMPLATE_PATH.read_text(encoding="utf-8")
    html = html.replace(
        "<!--ENV_MODULE_CARDS-->",
        environment_module_cards(),
    )
    return _with_shared_nav(html)


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
    return _render_template(
        TASKS_TEMPLATE_PATH,
        "tasks template missing",
    )


def render_datasets_page() -> str:
    return _render_template(
        DATASET_TEMPLATE_PATH,
        "datasets template missing",
    )


def render_analysis_page() -> str:
    return _render_template(
        ANALYSIS_TEMPLATE_PATH,
        "advanced analysis template missing",
    )


def run_environment_check(module: str, with_ml: bool = False) -> dict:
    if module == "experiment-plan-one":
        try:
            from experiment_plan_one.coverage import audit_plan_environment

            environment = audit_plan_environment()
            packages = environment.get("python_packages") or {}
            vina = environment.get("vina") or {}
            gromacs = environment.get("gromacs") or {}
            mmpbsa = environment.get("gmx_mmpbsa") or {}
            cellchat = environment.get("cellchat") or {}
            required_ok = bool(
                vina.get("path")
                and gromacs.get("path")
                and packages.get("rdkit")
                and packages.get("scanpy")
                and packages.get("shap")
            )
            lines = [
                "实验方案一环境检查",
                f"Vina: {vina.get('version') or 'not found'} "
                f"(expected {vina.get('expected', '1.2.3')})",
                f"GROMACS: {gromacs.get('version') or 'not found'} "
                f"(expected {gromacs.get('expected', '2022')})",
                "gmx_MMPBSA: "
                + ("available" if mmpbsa.get("available") else "not installed"),
                "R CellChat: "
                + ("available" if cellchat.get("available") else "not installed"),
                "Python packages: "
                + ", ".join(
                    f"{name}={'yes' if available else 'no'}"
                    for name, available in packages.items()
                ),
            ]
            return {
                "module": module,
                "ok": required_ok,
                "output": "\n".join(lines),
                "optional": {
                    "gmx_MMPBSA": bool(mmpbsa.get("available")),
                    "cellchat": bool(cellchat.get("available")),
                    "vina_version_match": bool(vina.get("match")),
                    "gromacs_version_match": bool(gromacs.get("match")),
                },
            }
        except Exception as exc:  # noqa: BLE001
            return {
                "module": module,
                "ok": False,
                "output": f"实验方案一环境检查失败：{exc}",
            }
    cmd = [
        sys.executable,
        str(APP_ROOT / "launchers" / "install_environment.py"),
        "check",
        module,
    ]
    if not with_ml:
        cmd += ["--skip-ml"]
    try:
        result = subprocess.run(
            cmd,
            cwd=APP_ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=420,
        )
        output = (result.stdout or "") + (result.stderr or "")
        optional = {}
        if module == "full":
            import importlib.util

            posebusters = importlib.util.find_spec("posebusters") is not None
            optional["posebusters"] = posebusters
            output += (
                "\n[optional] PoseBusters: "
                + ("available" if posebusters else "not installed")
            )
        return {
            "module": module,
            "ok": result.returncode == 0,
            "output": output[-8000:] or "环境检查完成。",
            "optional": optional,
        }
    except subprocess.TimeoutExpired:
        return {
            "module": module,
            "ok": False,
            "output": "环境检查超时，请确认本机网络或稍后重试。",
        }
    except Exception as exc:
        return {
            "module": module,
            "ok": False,
            "output": f"环境检查失败：{exc}",
        }


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
        log.warning("invalid max_results %r; using 20", _first(data, "max_results", ""))
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
        # Only load models produced by this console. joblib/pickle deserializes
        # arbitrary code, so an unconstrained path is remote code execution.
        allowed_model = _model_path_allowed(model_file)
        if allowed_model is None:
            raise ValueError(
                "模型文件必须位于数据检索缓存目录内：" + str(DATASET_SEARCH_DIR)
            )
        from dataset_search_ml import load_model, rerank

        try:
            model = load_model(allowed_model, allow_root=DATASET_SEARCH_DIR)
        except TypeError:  # older signature without the allowlist argument
            model = load_model(allowed_model)
        rows = rerank(
            rows,
            disease,
            research_direction,
            model=model,
        )
        model_applied = True
        model_path = str(allowed_model)

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
                _write_json_atomic(DATASET_SEARCH_DIR / "download_results.json", results)
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


def _model_path_allowed(path: Path) -> Path | None:
    """Only accept ML model files inside the dataset-search cache directory."""
    try:
        candidate = path.expanduser().resolve()
    except (OSError, ValueError):
        return None
    root = DATASET_SEARCH_DIR.resolve()
    if candidate.is_file() and candidate.is_relative_to(root):
        return candidate
    return None


def dataset_full_pipeline_url(row: dict) -> str:
    """Build the full-pipeline URL; run_supported is a candidate, verified on download."""
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
        "seed": _int_field(data, "docking_seed"),
        "seeds": _int_list_field(data, "docking_seeds"),
        "replicates": _int_field(data, "docking_replicates"),
        "positive_control_pdbqt": (
            _first(data, "docking_positive_control", "") or None
        ),
        "positive_control_id": (
            _first(data, "docking_positive_control_id", "") or None
        ),
        "positive_control_max_affinity": _float_field(
            data,
            "docking_positive_control_max_affinity",
        ),
        "cutoff": _float_field(data, "cutoff"),
        "moderate_cutoff": _float_field(data, "moderate_cutoff"),
        "top_n": _int_field(data, "top_n"),
        "model": _first(data, "model", "") or None,
        "training_csv": _first(data, "training_csv", "") or None,
        "label_column": _first(data, "label_column", "") or None,
        "md_mode": _first(data, "md_mode", "") or None,
        "md_receptor_pdb": _first(data, "md_receptor_pdb", "") or None,
        "md_gmx_data": _first(data, "md_gmx_data", "") or None,
        "md_top_n": _int_field(data, "md_top_n"),
        "md_protein_ff": _first(data, "md_protein_ff", "") or None,
        "md_water": _first(data, "md_water", "") or None,
        "md_em_steps": _int_field(data, "md_em_steps"),
        "md_equil_steps": _int_field(data, "md_equil_steps"),
        "md_prod_steps": _int_field(data, "md_prod_steps"),
        "md_gen_seed": _int_field(data, "md_gen_seed"),
        "md_temperature": _float_field(data, "md_temperature"),
        "md_ligand_charge": _int_field(data, "md_ligand_charge"),
        "md_cpu": _int_field(data, "md_cpu"),
        "md_gpu": _first(data, "md_gpu", "0") in ("1", "true", "on", "yes"),
        "md_topology_dir": _first(data, "md_topology_dir", "") or None,
        "md_mmpbsa_command": _first(data, "md_mmpbsa_command", "") or None,
        "md_mmpbsa_timeout": _int_field(data, "md_mmpbsa_timeout"),
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
    with DOCK_QUEUE_LOCK:
        DOCK_QUEUE.append(DOCK_JOBS[job_id])
    _drain_dock_queue()
    return {
        "job": job_id,
        "log_url": f"/dock/log?job={job_id}",
        "status_url": f"/dock/status?job={job_id}",
    }


def _start_dock_process(info: dict) -> None:
    _spawn_process(info)


def _drain_dock_queue() -> None:
    _drain_store(DOCK_QUEUE, DOCK_QUEUE_LOCK)


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
            if stage in ("", "未知") and info.get("stage") == "md-simulation":
                stage = "MD 模拟"
            return {
                "running": False,
                "ok": False,
                "queued": False,
                "paused": True,
                "stage": stage,
                "error": error,
            }
        ok = info["proc"].returncode == 0
        with JOB_RECORD_LOCK:
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
        if stage in ("", "未知") and info.get("stage") == "md-simulation":
            stage = "MD 模拟"
        return {
            "running": False,
            "ok": ok,
            "queued": False,
            "paused": False,
            "stage": stage,
            "error": error,
        }
    return {"running": True, "ok": False, "queued": False, "paused": False, "stage": "", "error": ""}


def start_molecular_docking_job(data: dict) -> dict:
    from molecular_docking.config import load_config, save_config

    workdir = Path(
        _first(data, "workdir", str(APP_ROOT / "molecular_docking"))
    ).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    stage = _first(data, "stage", "pipeline")
    job_id = uuid.uuid4().hex[:8]
    cfg_path = workdir / "config" / f"molecular_docking_web_{job_id}.json"

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
        "moderate_cutoff": _float_field(data, "moderate_cutoff"),
        "top_n": _int_field(data, "top_n"),
        "executable": _first(data, "executable", "") or None,
    }
    cfg = load_config(
        APP_ROOT / "config" / "molecular_docking_config.json",
        overrides,
    )
    save_config(cfg, cfg_path)

    log_path = workdir / "logs" / f"molecular_docking_{job_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPTS_DIR / "run_molecular_docking.py"),
        stage,
        "--config",
        str(cfg_path),
    ]
    if _first(data, "force", "") in ("1", "true", "on", "yes"):
        cmd.append("--force")

    env = os.environ.copy()
    env["DOCK_WORKDIR"] = str(workdir)
    MOLECULAR_DOCK_JOBS[job_id] = {
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
    with MOLECULAR_DOCK_QUEUE_LOCK:
        MOLECULAR_DOCK_QUEUE.append(MOLECULAR_DOCK_JOBS[job_id])
    _drain_molecular_docking_queue()
    return {
        "job": job_id,
        "log_url": f"/molecular-docking/log?job={job_id}",
        "status_url": f"/molecular-docking/status?job={job_id}",
    }


def _start_molecular_docking_process(info: dict) -> None:
    _spawn_process(info)


def _drain_molecular_docking_queue() -> None:
    _drain_store(MOLECULAR_DOCK_QUEUE, MOLECULAR_DOCK_QUEUE_LOCK)


def _molecular_docking_status(info: dict) -> dict:
    queued = info.get("proc") is None
    if queued:
        return {
            "running": False,
            "ok": False,
            "queued": True,
            "paused": False,
            "stage": "",
            "error": "",
        }
    running = info["proc"].poll() is None
    paused = bool(info.get("paused"))
    ok = False
    if not running:
        marker_dir = info["output_dir"] / ".stages"
        log_paths = [info["log"]]
        if paused:
            stage, error = _finished_info(
                marker_dir,
                MOLECULAR_DOCK_STAGE_LABELS,
                log_paths,
                True,
            )
            _notify_finished(
                info,
                "molecular-docking",
                "分子对接",
                f"{info.get('stage', 'pipeline')} 分子对接",
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
        with JOB_RECORD_LOCK:
            if not info.get("recorded"):
                record_molecular_docking_job(info, ok)
                info["recorded"] = True
        _drain_molecular_docking_queue()
        stage, error = _finished_info(
            marker_dir,
            MOLECULAR_DOCK_STAGE_LABELS,
            log_paths,
            False,
        )
        status = "completed" if ok else "interrupted"
        if status == "completed":
            error = ""
        _notify_finished(
            info,
            "molecular-docking",
            "分子对接",
            f"{info.get('stage', 'pipeline')} 分子对接",
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
    return {
        "running": True,
        "ok": False,
        "queued": False,
        "paused": False,
        "stage": "",
        "error": "",
    }


def _start_plan_one_job(data: dict) -> dict:
    output_value = _first(data, "plan_output_root", "").strip()
    if not output_value:
        raise ValueError("实验方案一结果目录不能为空")
    output_root = Path(output_value).expanduser()
    if not output_root.is_absolute():
        output_root = APP_ROOT / output_root
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    base_config_path = APP_ROOT / "config" / "experiment_plan_one.json"
    config = json.loads(base_config_path.read_text(encoding="utf-8"))
    config["output_root"] = str(output_root)
    compound = config.setdefault("compound", {})
    disease = config.setdefault("disease", {})
    md = config.setdefault("md", {})
    docking = config.setdefault("docking", {})
    single_cell = config.setdefault("single_cell", {}).setdefault("mouse", {})

    for field, key in (
        ("plan_compound_name", "name"),
        ("plan_pubchem_cid", "pubchem_cid"),
        ("plan_target_prediction_file", "target_prediction_file"),
        ("plan_target_prediction_gene_column", "target_prediction_gene_column"),
        ("plan_target_prediction_score_column", "target_prediction_score_column"),
    ):
        value = _first(data, field, "").strip()
        if value:
            compound[key] = value
    for field, key in (
        ("plan_disease_name", "name"),
        ("plan_gene_cards_file", "gene_cards_file"),
        ("plan_omim_file", "omim_file"),
        ("plan_ttd_file", "ttd_file"),
    ):
        value = _first(data, field, "").strip()
        if value:
            disease[key] = value
    disease["use_open_evidence_sources"] = _first(
        data,
        "plan_use_open_evidence",
        "",
    ) in ("1", "true", "on", "yes")

    for field, key, default in (
        ("plan_docking_targets", "targets", 5),
        ("plan_docking_exhaustiveness", "exhaustiveness", 16),
        ("plan_docking_cpu", "cpu", 4),
    ):
        value = _int_field(data, field)
        if value is not None:
            docking[key] = value
        elif key not in docking:
            docking[key] = default
    md["run"] = _first(data, "plan_md_run", "") in (
        "1",
        "true",
        "on",
        "yes",
    )
    md["gpu"] = _first(data, "plan_md_gpu", "") in (
        "1",
        "true",
        "on",
        "yes",
    )
    md_cpu = _int_field(data, "plan_md_cpu")
    if md_cpu is not None:
        md["cpu"] = md_cpu
    mouse_cells = _int_field(data, "plan_mouse_max_cells")
    if mouse_cells is not None:
        single_cell["max_cells"] = mouse_cells
    permutations = _int_field(data, "plan_cellchat_permutations")
    if permutations is not None:
        single_cell["cellchat_permutations"] = permutations

    config_path = output_root / "web_experiment_plan_one_config.json"
    config_path.write_text(
        json.dumps(config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    job_id = uuid.uuid4().hex[:8]
    log_path = output_root / "logs" / f"web_plan_one_{job_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPTS_DIR / "run_experiment_plan_one.py"),
        "--config",
        str(config_path),
        "--output-root",
        str(output_root),
    ]
    start_stage = _first(data, "plan_start_stage", "").strip()
    if start_stage:
        cmd += ["--start-stage", start_stage]
    for value in _first(data, "plan_stages", "").replace("\n", ",").split(","):
        if value.strip():
            cmd += ["--stage", value.strip()]
    for value in _first(data, "plan_skip_stages", "").replace(
        "\n",
        ",",
    ).split(","):
        if value.strip():
            cmd += ["--skip-stage", value.strip()]
    if _first(data, "plan_force", "") in ("1", "true", "on", "yes"):
        cmd.append("--force")
    if _first(data, "plan_verbose", "") in ("1", "true", "on", "yes"):
        cmd.append("--verbose")
    if _first(data, "plan_dry_run", "") in ("1", "true", "on", "yes"):
        cmd.append("--dry-run")
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    FULL_JOBS[job_id] = {
        "job_id": job_id,
        "mode": "experiment_plan_one",
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "workdir": output_root,
        "output": str(output_root),
        "cmd": cmd,
        "env": env,
        "queued": True,
        "notified": False,
    }
    with FULL_QUEUE_LOCK:
        FULL_QUEUE.append(FULL_JOBS[job_id])
    _drain_full_queue()
    return {
        "job": job_id,
        "log_url": f"/full/log?job={job_id}",
        "status_url": f"/full/status?job={job_id}",
    }


def start_full_job(data: dict) -> dict:
    if _first(data, "run_mode", "").strip() == "experiment_plan_one":
        return _start_plan_one_job(data)
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
    output = Path(output_value).expanduser()
    if not output.is_absolute():
        output = APP_ROOT / output
    output = output.resolve()
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
    cmd += ["--output", str(output)]
    species = _first(data, "species", "").strip()
    if species:
        cmd += ["--species", species]
    ml_model = _first(data, "ml_model", "").strip()
    if ml_model:
        cmd += ["--ml-model", ml_model]
    top_genes = _int_field(data, "top_genes")
    if top_genes:
        cmd += ["--top-genes", str(top_genes)]
    candidate_universe_size = _int_field(data, "candidate_universe_size")
    if candidate_universe_size is not None:
        cmd += [
            "--candidate-universe-size",
            str(candidate_universe_size),
        ]
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
    advanced_priority_csv = _first(data, "advanced_priority_csv", "").strip()
    if advanced_priority_csv:
        cmd += [
            "--advanced-priority-csv",
            _cli_path(advanced_priority_csv),
        ]
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
        "skip_evidence_hub",
        "evidence_hub_offline",
        "evidence_hub_strict",
        "skip_pseudobulk",
        "skip_knockout",
        "skip_docking",
        "allow_review_docking",
        "skip_md",
        "skip_handoff",
        "skip_docking_ml",
        "skip_network",
        "skip_faers",
        "skip_cell_feedback",
        "skip_qc_gate",
        "skip_differential_abundance",
        "keep_all_genes",
        "force",
        "dry_run",
    ]:
        if _first(data, flag, "") in ("1", "true", "on", "yes"):
            cmd.append("--" + flag.replace("_", "-"))

    md_mode = _first(data, "md_mode", "").strip()
    if md_mode:
        cmd += ["--md-mode", md_mode]
    md_top_n = _int_field(data, "md_top_n")
    if md_top_n:
        cmd += ["--md-top-n", str(md_top_n)]
    docking_ml_model = _first(data, "docking_ml_model", "").strip()
    if docking_ml_model:
        cmd += ["--docking-ml-model", docking_ml_model]
    docking_ml_training_csv = _first(
        data, "docking_ml_training_csv", ""
    ).strip()
    if docking_ml_training_csv:
        cmd += ["--docking-ml-training-csv", docking_ml_training_csv]
    docking_ml_label_column = _first(
        data, "docking_ml_label_column", ""
    ).strip()
    if docking_ml_label_column:
        cmd += ["--docking-ml-label-column", docking_ml_label_column]
    for attr in (
        "evidence_hub_config",
        "evidence_disease",
        "evidence_disease_id",
    ):
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_", "-"), value]
    for attr in (
        "evidence_max_targets",
        "evidence_max_records",
        "evidence_hub_timeout",
        "evidence_legacy_pool_size",
        "candidate_expansion_max_targets",
        "benchmark_top_n",
    ):
        value = _int_field(data, attr)
        if value is not None:
            cmd += ["--" + attr.replace("_", "-"), str(value)]
    for attr in (
        "benchmark_positive_targets",
        "benchmark_negative_targets",
    ):
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_targets", "").replace("_", "-"), value]
    for attr in (
        "benchmark_source",
        "benchmark_version",
        "benchmark_exclude_sources",
    ):
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_", "-"), value]
    if _first(data, "benchmark_independent", "") in (
        "1",
        "true",
        "on",
        "yes",
    ):
        cmd.append("--benchmark-independent")
    external_validation_path = _first(
        data,
        "external_validation_path",
        "",
    ).strip()
    if external_validation_path:
        cmd += ["--external-validation-path", external_validation_path]
    for attr in (
        "external_validation_target_column",
        "external_validation_score_column",
        "external_validation_label_column",
    ):
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_", "-"), value]
    external_validation_threshold = _first(
        data,
        "external_validation_threshold",
        "",
    ).strip()
    if external_validation_threshold:
        cmd += [
            "--external-validation-threshold",
            external_validation_threshold,
        ]
    external_validation_bootstrap = _int_field(
        data,
        "external_validation_bootstrap",
    )
    if external_validation_bootstrap is not None:
        cmd += [
            "--external-validation-bootstrap",
            str(external_validation_bootstrap),
        ]
    if _first(
        data,
        "external_validation_allow_external_score",
        "",
    ) in ("1", "true", "on", "yes"):
        cmd.append("--allow-external-validation-score")
    for attr in [
        "network_compound_targets_csv",
        "network_disease_genes_csv",
        "network_ppi_network_csv",
        "network_target_sources_dir",
    ]:
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_", "-"), value]
    network_disease_gene_column = _first(
        data, "network_disease_gene_column", ""
    ).strip()
    if network_disease_gene_column:
        cmd += ["--network-disease-gene-column", network_disease_gene_column]
    network_cytoscape = _first(data, "network_cytoscape", "").strip()
    if network_cytoscape in ("auto", "on", "off"):
        cmd += ["--network-cytoscape", network_cytoscape]
    network_cytoscape_url = _first(
        data,
        "network_cytoscape_url",
        "",
    ).strip()
    if network_cytoscape_url:
        cmd += ["--network-cytoscape-url", network_cytoscape_url]
    network_cytoscape_layout = _first(
        data,
        "network_cytoscape_layout",
        "",
    ).strip()
    if network_cytoscape_layout:
        cmd += ["--network-cytoscape-layout", network_cytoscape_layout]
    if _first(data, "network_cytoscape_session", "") in (
        "1",
        "true",
        "on",
        "yes",
    ):
        cmd += ["--network-cytoscape-session"]
    network_max_ppi_edges = _int_field(data, "network_max_ppi_edges")
    if network_max_ppi_edges is not None:
        cmd += ["--network-max-ppi-edges", str(network_max_ppi_edges)]
    if _first(data, "network_run_enrichment", "") in (
        "1",
        "true",
        "on",
        "yes",
    ):
        cmd.append("--network-run-enrichment")
    network_enrichment_timeout = _int_field(
        data,
        "network_enrichment_timeout",
    )
    if network_enrichment_timeout:
        cmd += [
            "--network-enrichment-timeout",
            str(network_enrichment_timeout),
        ]
    faers_input = _first(data, "faers_input", "").strip()
    if faers_input:
        cmd += ["--faers-input", faers_input]
    for attr in ["faers_drug_column", "faers_event_column"]:
        value = _first(data, attr, "").strip()
        if value:
            cmd += ["--" + attr.replace("_", "-"), value]
    faers_count_column = _first(data, "faers_count_column", "").strip()
    if faers_count_column:
        cmd += ["--faers-count-column", faers_count_column]
    faers_min_count = _int_field(data, "faers_min_count")
    if faers_min_count:
        cmd += ["--faers-min-count", str(faers_min_count)]

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    advanced_config = _cli_path(_first(data, "advanced_config", ""))
    if advanced_config:
        env["LIVER_ADVANCED_CONFIG"] = advanced_config
    mr_config = _cli_path(_first(data, "mr_config", ""))
    if mr_config:
        env["LIVER_MR_CONFIG"] = mr_config
    FULL_JOBS[job_id] = {
        "job_id": job_id,
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "workdir": workdir,
        "accession": accession,
        "output": str(output),
        "cmd": cmd,
        "env": env,
        "queued": True,
        "notified": False,
    }
    with FULL_QUEUE_LOCK:
        FULL_QUEUE.append(FULL_JOBS[job_id])
    _drain_full_queue()
    return {
        "job": job_id,
        "log_url": f"/full/log?job={job_id}",
        "status_url": f"/full/status?job={job_id}",
    }


def _start_full_process(info: dict) -> None:
    _spawn_process(info)


def _drain_full_queue() -> None:
    _drain_store(FULL_QUEUE, FULL_QUEUE_LOCK)


def _plan_one_status(info: dict) -> dict:
    queued = info.get("proc") is None
    if queued:
        return {
            "running": False,
            "ok": False,
            "queued": True,
            "paused": False,
            "stage": "",
            "error": "",
        }
    proc = info["proc"]
    running = proc.poll() is None
    log_text = (
        info["log"].read_text(encoding="utf-8", errors="replace")[-12000:]
        if info["log"].exists()
        else ""
    )
    stage = ""
    for line in reversed(log_text.splitlines()):
        match = re.search(r"(?:starting stage|completed stage)\s+([A-Za-z_]+)", line)
        if match:
            stage = match.group(1)
            break
    error = ""
    if not running and proc.returncode != 0:
        error_lines = [
            line
            for line in log_text.splitlines()
            if line.strip().startswith("ERROR:")
        ]
        error = error_lines[-1] if error_lines else "experiment plan one failed"
    if not running:
        _notify_finished(
            info,
            "full",
            "实验方案一",
            "6PPD-Q / NAFLD 实验方案一",
            "finished" if proc.returncode == 0 else "failed",
            stage,
            error,
            exit_code=proc.returncode,
        )
    return {
        "running": running,
        "ok": not running and proc.returncode == 0,
        "queued": False,
        "paused": False,
        "stage": stage,
        "error": error,
    }


def _full_status(info: dict) -> dict:
    if info.get("mode") == "experiment_plan_one":
        return _plan_one_status(info)
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
        with JOB_RECORD_LOCK:
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
































def running_tasks_data(include_logs: bool = True) -> dict:
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
        log_text = (
            "\n".join(
                [
                    _log_tail(info["log"]),
                    _log_tail(info["out"] / "logs" / "pipeline_r.log"),
                ]
            ).strip()
            if include_logs
            else ""
        )
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
            else (
                "MD 模拟"
                if info.get("stage") == "md-simulation"
                else _current_stage(marker_dir, DOCK_STAGE_LABELS) or "准备中"
            )
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
                "log_tail": _log_tail(info["log"]) if include_logs else "",
            }
        )

    for job_id, info in list(MOLECULAR_DOCK_JOBS.items()):
        status = _molecular_docking_status(info)
        if not (
            status.get("running")
            or status.get("queued")
            or status.get("paused")
        ):
            continue
        state = (
            "queued"
            if status.get("queued")
            else "paused"
            if status.get("paused")
            else "running"
        )
        marker_dir = info["output_dir"] / ".stages"
        progress = (
            0
            if state == "queued"
            else _marker_progress(marker_dir, len(MOLECULAR_DOCK_STAGE_LABELS))
        )
        stage_label = (
            "排队中"
            if state == "queued"
            else "已暂停"
            if state == "paused"
            else _current_stage(marker_dir, MOLECULAR_DOCK_STAGE_LABELS)
            or "准备中"
        )
        started = float(info.get("started") or now)
        tasks.append(
            {
                "page": "molecular-docking",
                "page_label": "分子对接",
                "job": job_id,
                "url": f"/molecular-docking?job={job_id}",
                "title": f"{info.get('stage', 'pipeline')} · 分子对接",
                "detail": str(info["workdir"]),
                "status": state,
                "started": started,
                "elapsed": int(now - started),
                "progress": progress,
                "stage_label": stage_label,
                "log_tail": _log_tail(info["log"]) if include_logs else "",
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
        progress = 0 if state == "queued" else _marker_progress(marker_dir, 11)
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
                "log_tail": (
                    "\n".join(
                        _log_tail(path, 4000) for path in _full_log_paths(info)
                    ).strip()
                    if include_logs
                    else ""
                ),
            }
        )

    for job_id, info in list(ANALYSIS_JOBS.items()):
        status = _analysis_status(info)
        if not (status.get("running") or status.get("queued")):
            continue
        state = "queued" if status.get("queued") else "running"
        started = float(info.get("started") or now)
        tasks.append(
            {
                "page": "analysis",
                "page_label": "高级分析",
                "job": job_id,
                "url": (
                    f"/analysis?kind={info.get('kind', '')}&job={job_id}"
                ),
                "title": info.get("title", "高级分析"),
                "detail": str(info.get("output", "")),
                "status": state,
                "started": started,
                "elapsed": int(now - started),
                "progress": 0,
                "stage_label": "排队中" if state == "queued" else "运行中",
                "log_tail": _log_tail(info["log"]) if include_logs else "",
            }
        )

    tasks.sort(key=lambda item: float(item.get("started") or 0))
    return {"tasks": tasks}


def running_task_counts() -> dict:
    tasks = running_tasks_data(include_logs=False)["tasks"]
    return {
        "count": len(tasks),
        "running": sum(1 for task in tasks if task.get("status") == "running"),
        "queued": sum(1 for task in tasks if task.get("status") == "queued"),
        "paused": sum(1 for task in tasks if task.get("status") == "paused"),
    }






def _full_result_files(workdir: Path) -> list[str]:
    files: list[str] = []
    roots: list[tuple[Path, str]] = [
        (workdir / "outputs" / "integration", ""),
        (workdir / "outputs" / "run_001" / "results", "outputs/run_001/results"),
        (workdir / "outputs" / "run_001" / "docked", "outputs/run_001/docked"),
        (
            workdir / "outputs" / "run_001" / "external",
            "outputs/run_001/external",
        ),
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
            roots.append(
                (
                    gene_dir / "outputs" / "md",
                    f"work/{gene_dir.name}/outputs/md",
                )
            )
    for root, prefix in roots:
        for rel in _list_result_files(root):
            files.append(f"{prefix}/{rel}" if prefix else rel)
    return sorted(set(files))


PLAN_ONE_RESULT_SUFFIXES = {
    ".png",
    ".jpg",
    ".jpeg",
    ".svg",
    ".pdf",
    ".csv",
    ".json",
    ".md",
    ".txt",
    ".html",
    ".tsv",
    ".xlsx",
    ".dat",
    ".xvg",
    ".gro",
    ".tpr",
    ".xtc",
    ".pdb",
    ".pdbqt",
    ".itp",
    ".top",
    ".log",
}


def _plan_one_result_files(workdir: Path) -> list[str]:
    files: list[str] = []
    roots = [
        "01_compound_characterization",
        "02_disease_targets",
        "02b_evidence",
        "03_intersection_ppi",
        "04_bulk_training",
        "05_machine_learning",
        "06_single_cell_mouse",
        "07_single_cell_human",
        "08_docking",
        "09_md_mmpbsa",
        "10_reports",
        "按方案分类",
    ]
    summary = workdir / "RESULTS_SUMMARY.md"
    if summary.exists():
        files.append(summary.name)
    for name in roots:
        root = workdir / name
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if (
                path.is_file()
                and path.suffix.lower() in PLAN_ONE_RESULT_SUFFIXES
            ):
                files.append(
                    f"{name}/{path.relative_to(root).as_posix()}"
                )
    return sorted(set(files))


def _plan_one_results(workdir: Path) -> dict:
    workdir = Path(workdir).expanduser().resolve()
    coverage = _read_json(
        workdir / "10_reports" / "plan_coverage" / "plan_coverage.json"
    )
    figure_audit = _read_json(
        workdir
        / "10_reports"
        / "figure_quality_audit"
        / "figure_quality_audit.json"
    )
    status = _read_json(workdir / "10_reports" / "analysis_status.json")
    summary_path = workdir / "RESULTS_SUMMARY.md"
    summary_text = (
        summary_path.read_text(encoding="utf-8", errors="replace")
        if summary_path.exists()
        else ""
    )
    return {
        "workdir": str(workdir),
        "exists": bool(coverage or figure_audit or status),
        "mode": "experiment_plan_one",
        "files": _plan_one_result_files(workdir),
        "plan_one": {
            "coverage": coverage,
            "figure_audit": figure_audit,
            "status": status,
            "summary_markdown": summary_text,
        },
        "summary": {},
    }


def full_results(workdir: Path) -> dict:
    workdir = Path(workdir).expanduser().resolve()
    out = workdir / "outputs" / "integration"
    if (
        (workdir / "10_reports" / "plan_coverage" / "plan_coverage.json").exists()
        or (not out.exists() and (workdir / "10_reports").exists())
    ):
        return _plan_one_results(workdir)
    result = {
        "workdir": str(workdir),
        "exists": out.exists(),
        "files": [],
        "summary": {},
        "qc_metrics": {},
        "differential_abundance": [],
        "key_genes": [],
        "target_priority": [],
        "target_priority_summary": {},
        "target_validation": [],
        "external_validation": {},
        "target_decisions": [],
        "target_decision_summary": {},
        "omics_qc": {},
        "pose_qc": {},
        "structural_quality": {},
        "publication_readiness": {},
        "benchmark": {},
        "reproducibility": {},
        "knockout": [],
        "docking": [],
        "cadd_downstream_summary": {},
        "cadd_downstream": [],
        "network_summary": {},
        "network": [],
        "faers_summary": {},
        "faers": [],
        "cell_feedback": [],
        "cell_feedback_deg": [],
        "cell_feedback_enrichment_go": [],
        "cell_feedback_enrichment_kegg": [],
        "evidence": [],
        "evidence_coverage": [],
        "source_ablation": [],
    }
    if not out.exists():
        return result
    result["summary"] = _read_json(out / "integration_summary.json")
    result["qc_metrics"] = _read_json(out / "qc_metrics.json")
    result["target_priority_summary"] = (
        result["summary"].get("target_priority") or {}
    )
    result["publication_readiness"] = (
        result["summary"].get("publication_readiness") or {}
    )
    result["benchmark"] = (
        result["target_priority_summary"].get("benchmark")
        or (result["summary"].get("evidence_hub") or {}).get("benchmark")
        or {}
    )
    result["target_decision_summary"] = (
        result["summary"].get("target_decisions") or {}
    )
    result["omics_qc"] = _read_json(out / "omics_qc_summary.json")
    result["pose_qc"] = _read_json(out / "pose_qc_summary.json")
    result["structural_quality"] = _read_json(
        out / "structural_quality_summary.json"
    )
    result["reproducibility"] = _read_json(
        out / "reproducibility_manifest.json"
    )
    result["files"] = _full_result_files(workdir)
    try:
        result["key_genes"] = json.loads(
            pd_read_csv(out / "key_genes.csv").to_json(orient="records")
        )
    except Exception:
        result["key_genes"] = []
    priority_path = out / "integrated_target_priority.csv"
    if not priority_path.exists():
        priority_path = out / "target_priority.csv"
    try:
        result["target_priority"] = json.loads(
            pd_read_csv(priority_path).head(200).to_json(orient="records")
        )
    except Exception:
        result["target_priority"] = []
    try:
        result["target_validation"] = json.loads(
            pd_read_csv(out / "target_validation_scores.csv")
            .head(200)
            .to_json(orient="records")
        )
    except Exception:
        result["target_validation"] = []
    result["external_validation"] = _read_json(
        out / "external_validation_summary.json"
    )
    try:
        result["target_decisions"] = json.loads(
            pd_read_csv(out / "target_decision_report.csv")
            .head(200)
            .to_json(orient="records")
        )
    except Exception:
        result["target_decisions"] = []
    try:
        # The page renders only the top knockout rows; cap the JSON payload.
        result["knockout"] = json.loads(
            pd_read_csv(
                workdir
                / "outputs"
                / "run_001"
                / "results"
                / "04_knockout"
                / "data"
                / "fig_52_53_ranked_knockout.csv"
            )
            .head(200)
            .to_json(orient="records")
        )
    except Exception:
        result["knockout"] = []
    try:
        result["docking"] = json.loads(
            pd_read_csv(out / "docking_targets.csv").to_json(orient="records")
        )
    except Exception:
        result["docking"] = []
    result["cadd_downstream_summary"] = _read_json(
        out / "cadd_downstream_summary.json"
    )
    try:
        result["cadd_downstream"] = json.loads(
            pd_read_csv(out / "cadd_targets.csv").to_json(orient="records")
        )
    except Exception:
        result["cadd_downstream"] = []
    result["network_summary"] = _read_json(out / "network_summary.json")
    network_csv = (
        (result["network_summary"].get("outputs") or {}).get("overlap_csv")
        or result["network_summary"].get("overlap_csv")
        or ""
    )
    if network_csv:
        try:
            result["network"] = json.loads(
                pd_read_csv(Path(str(network_csv)))
                .head(200)
                .to_json(orient="records")
            )
        except Exception:
            result["network"] = []
    result["faers_summary"] = _read_json(out / "faers_summary.json")
    faers_csv = result["faers_summary"].get("output_csv") or ""
    if faers_csv:
        try:
            result["faers"] = json.loads(
                pd_read_csv(Path(str(faers_csv)))
                .head(200)
                .to_json(orient="records")
            )
        except Exception:
            result["faers"] = []
    try:
        result["cell_feedback"] = json.loads(
            pd_read_csv(
                out / "cell_feedback" / "data" / "feedback_targets.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["cell_feedback"] = []
    try:
        result["cell_feedback_deg"] = json.loads(
            pd_read_csv(
                out / "cell_feedback" / "data" / "feedback_deg.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["cell_feedback_deg"] = []
    try:
        result["cell_feedback_enrichment_go"] = json.loads(
            pd_read_csv(
                out / "cell_feedback" / "data" / "feedback_enrichment_go.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["cell_feedback_enrichment_go"] = []
    try:
        result["cell_feedback_enrichment_kegg"] = json.loads(
            pd_read_csv(
                out / "cell_feedback" / "data" / "feedback_enrichment_kegg.csv"
            ).to_json(orient="records")
        )
    except Exception:
        result["cell_feedback_enrichment_kegg"] = []
    try:
        result["evidence"] = json.loads(
            pd_read_csv(out / "gene_evidence.csv").to_json(orient="records")
        )
    except Exception:
        result["evidence"] = []
    try:
        result["evidence_coverage"] = json.loads(
            pd_read_csv(out / "evidence_hub" / "evidence_coverage.csv")
            .head(200)
            .to_json(orient="records")
        )
    except Exception:
        result["evidence_coverage"] = []
    try:
        result["source_ablation"] = json.loads(
            pd_read_csv(out / "evidence_hub" / "source_ablation.csv")
            .head(200)
            .to_json(orient="records")
        )
    except Exception:
        result["source_ablation"] = []
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


def _full_workdir_for_query(job: str, workdir: str) -> str:
    """Resolve a full-pipeline workdir from an explicit query or a live job."""
    if workdir:
        return workdir
    info = FULL_JOBS.get(job or "")
    if info:
        return str(info["workdir"])
    return ""


def _full_file_path(workdir: Path, name: str) -> Path | None:
    workdir = Path(workdir).expanduser().resolve()
    name_path = Path(name)
    if name_path.is_absolute():
        return None
    if (
        (workdir / "10_reports" / "plan_coverage" / "plan_coverage.json").exists()
        or (
            not (workdir / "outputs" / "integration").exists()
            and (workdir / "10_reports").exists()
        )
    ):
        target = (workdir / name_path).resolve()
        if (
            target.is_relative_to(workdir)
            and target.is_file()
            and target.suffix.lower() in PLAN_ONE_RESULT_SUFFIXES
        ):
            return target
        return None
    single_cell_root = _single_cell_root_from_workdir(workdir)
    allowed_roots = [
        (workdir / "outputs" / "integration").resolve(),
        (workdir / "outputs" / "run_001" / "results").resolve(),
        (workdir / "outputs" / "run_001" / "docked").resolve(),
        (workdir / "outputs" / "run_001" / "external").resolve(),
        (workdir / "outputs" / "run_001" / "network_toxicology").resolve(),
        (workdir / "outputs" / "run_001" / "faers").resolve(),
        (workdir / "work").resolve(),
        (workdir / "data" / "knockout").resolve(),
    ]
    if single_cell_root:
        allowed_roots.append(single_cell_root.resolve())
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
    md = _md_results(info)
    return {
        "summary": summary,
        "rows": rows,
        "figures": figures,
        "files": files,
        "md": md,
        "output_dir": str(reports),
        "stage": info.get("stage", ""),
    }


def _md_results(info: dict) -> dict:
    md_dir = info["output_dir"] / "results" / "06_md"
    if not md_dir.exists():
        return {}
    summary = {}
    summary_path = md_dir / "md_simulation_summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            summary = {}
    rows: list[dict] = []
    results_csv = md_dir / "md_simulation_results.csv"
    if results_csv.exists():
        import csv as csv_module

        with results_csv.open("r", newline="", encoding="utf-8") as fh:
            for i, row in enumerate(csv_module.DictReader(fh)):
                if i >= 200:
                    break
                rows.append(row)
    return {
        "summary": summary,
        "rows": rows,
        "figures": _list_result_images(md_dir),
        "files": _list_result_files(md_dir),
        "output_dir": str(md_dir),
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








def molecular_docking_results(info: dict) -> dict:
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
    docked_root = info["output_dir"] / "docked"
    if docked_root.exists():
        for p in sorted(docked_root.glob("*.pdbqt")):
            if p.is_file():
                files.append(f"docked/{p.name}")
    report_rel = "molecular_docking_report.html"
    html_report = report_rel if (reports / report_rel).exists() else ""
    return {
        "summary": summary,
        "rows": rows,
        "figures": figures,
        "files": sorted(set(files)),
        "html_report": html_report,
        "output_dir": str(reports),
        "stage": info.get("stage", ""),
    }


def _molecular_docking_file_path(info: dict, name: str):
    out = info["output_dir"].resolve()
    results_root = (out / "results").resolve()
    docked_root = (out / "docked").resolve()
    clean = Path(name)
    target = (results_root / clean).resolve()
    if (
        (target.is_file() or target.suffix.lower() == ".html")
        and results_root in target.parents
    ):
        return target
    base = Path(name).name
    target = (docked_root / base).resolve()
    if target.is_file() and target.parent == docked_root:
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
        "group_column": _first(data, "ko_group_column", "") or None,
        "cell_type_column": _first(data, "ko_cell_type_column", "") or None,
        "ko_top_n": _int_field(data, "ko_top_n"),
        "insilico_gene": _first(data, "ko_insilico_gene", "") or None,
        "insilico_engine": _first(data, "ko_insilico_engine", "") or None,
        "insilico_raw_count_input": _raw_count_flag(
            _first(data, "ko_insilico_raw_count_input", "")
        ),
        "insilico_species": _first(data, "ko_insilico_species", "") or None,
        "insilico_embedding_csv": _first(
            data, "ko_insilico_embedding", ""
        ) or None,
        "insilico_regulators_csv": _first(
            data, "ko_insilico_regulators", ""
        ) or None,
        "insilico_photo_dir": _first(
            data, "ko_insilico_photo", ""
        ) or None,
    }
    cfg = load_config(APP_ROOT / "config" / "docking_config.json", overrides)
    insilico = cfg.data.setdefault("insilico_knockout", {})
    if insilico.get("ko_gene"):
        insilico["enabled"] = True
    for field, config_key, cast in (
        ("ko_insilico_max_genes", "max_genes", int),
        ("ko_insilico_max_cells", "max_cells", int),
        ("ko_insilico_n_propagation", "n_propagation", int),
        ("ko_drugreflector_top_n", "drugreflector_top_n", int),
    ):
        value = _first(data, field, "")
        if str(value).strip():
            try:
                insilico[config_key] = cast(value)
            except (TypeError, ValueError):
                log.warning(
                    "ignoring invalid %s value %r (expected %s)",
                    field,
                    value,
                    cast.__name__,
                )
    checkpoint = _first(data, "ko_drugreflector_checkpoint_dir", "").strip()
    if checkpoint:
        insilico["drugreflector_checkpoint_dir"] = checkpoint
    enrichment_value = _first(data, "ko_run_enrichment", "0").strip().lower()
    insilico["run_enrichment"] = enrichment_value in ("1", "true", "on", "yes")
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
    report_rel = "in_silico/in_silico_knockout_report.html"
    if (ko_dir / report_rel).exists():
        files.append(report_rel)
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
        "target_sources_dir": (
            _first(data, "net_target_sources_dir", "") or None
        ),
        "disease_genes_csv": _first(data, "net_disease_genes", "") or None,
        "disease_gene_column": (
            _first(data, "net_disease_gene_column", "") or None
        ),
        "ppi_network_csv": _first(data, "net_ppi", "") or None,
        "venn": _first(data, "net_venn", "1")
        in ("1", "true", "on", "yes"),
        "network_cytoscape": (
            _first(data, "net_cytoscape", "auto") or "auto"
        ),
        "network_cytoscape_url": (
            _first(data, "net_cytoscape_url", "http://127.0.0.1:1234")
            or "http://127.0.0.1:1234"
        ),
        "network_cytoscape_layout": _first(data, "net_cytoscape_layout", "") or None,
        "network_cytoscape_session": _first(
            data,
            "net_cytoscape_session",
            "",
        )
        in ("1", "true", "on", "yes"),
        "network_max_ppi_edges": _int_field(data, "net_max_ppi_edges"),
        "run_enrichment": _first(data, "net_run_enrichment", "")
        in ("1", "true", "on", "yes"),
        "enrichment_timeout": _int_field(
            data,
            "net_enrichment_timeout",
        ),
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
        ".xgmml",
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
        if (
            (_is_result_file(target) or target.suffix.lower() == ".html")
            and folder in target.parents
        ):
            return target
    return None


def _detect_and_save_molecular_box(workdir: str, receptor: str) -> dict:
    from docking.box import detect_box_data
    from molecular_docking.config import load_config, save_config

    workdir = Path(workdir).expanduser().resolve()
    receptor_path = Path(receptor).expanduser()
    if not receptor_path.is_absolute():
        receptor_path = workdir / receptor_path
    receptor_path = receptor_path.resolve()
    if not receptor_path.is_file():
        raise ValueError(f"receptor file not found: {receptor_path}")
    center, size, mode = detect_box_data(receptor_path)
    config_path = workdir / "config" / "molecular_docking_config.json"
    if config_path.exists():
        cfg = load_config(config_path)
    else:
        cfg = load_config(
            APP_ROOT / "config" / "molecular_docking_config.json",
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




def main() -> int:
    parser = argparse.ArgumentParser(description="Local web UI for pipeline")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument(
        "--page",
        choices=[
            "single",
            "dock",
            "md-simulation",
            "molecular-docking",
            "knockout",
            "network",
            "faers",
            "validation",
            "analysis",
            "full",
            "results",
            "tasks",
            "datasets",
            "guide",
            "environment",
        ],
        default="full",
        help="page to open in the browser",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="start the server without opening a browser window",
    )
    parser.add_argument(
        "--allow-path",
        action="append",
        default=[],
        metavar="DIR",
        help=(
            "additional directory the results browser may read from "
            "(repeatable); by default only project output roots and workdirs "
            "registered by jobs started in this console are allowed"
        ),
    )
    args = parser.parse_args()

    for extra in args.allow_path:
        EXTRA_WORKDIR_ROOTS.append(Path(extra).expanduser().resolve())

    if not _cleanup_stale_web_ui(args.host, args.port):
        print(f"ERROR: port {args.port} is still in use by another process.")
        return 1

    token = _configure_runtime(args.host, args.port)
    from web_handler import Handler  # noqa: E402 - lazy, breaks import cycle
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
    if token:
        print(
            "Non-loopback bind: a session token is required. Open this URL:\n"
            f"  {url}/?token={token}\n"
            "or send the header 'X-Auth-Token: <token>' with every request."
        )
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print(
            "WARNING: Web UI exposes pipeline command and file endpoints "
            "on a non-loopback host; only use this on a trusted network."
        )
    threading.Thread(
        target=_run_idle_shutdown_monitor,
        args=(server, time.monotonic()),
        daemon=True,
        name="web-ui-idle-shutdown",
    ).start()
    if args.page == "dock":
        open_url = url + "/dock"
    elif args.page == "md-simulation":
        open_url = url + "/md-simulation"
    elif args.page == "molecular-docking":
        open_url = url + "/molecular-docking"
    elif args.page == "knockout":
        open_url = url + "/knockout"
    elif args.page == "network":
        open_url = url + "/network"
    elif args.page == "faers":
        open_url = url + "/faers"
    elif args.page == "validation":
        open_url = url + "/validation"
    elif args.page == "analysis":
        open_url = url + "/analysis"
    elif args.page == "full":
        open_url = url + "/full"
    elif args.page == "results":
        open_url = url + "/results"
    elif args.page == "tasks":
        open_url = url + "/tasks"
    elif args.page == "datasets":
        open_url = url + "/datasets"
    elif args.page == "guide":
        open_url = url + "/guide"
    elif args.page == "environment":
        open_url = url + "/environment"
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
    # Running this file as a script makes it ``__main__``; register it under its
    # canonical name first so ``web_handler``'s ``import web_ui`` binds to THIS
    # module instead of importing a second, empty copy (which would split the
    # job registry and silently disable the auth token / allow-path checks).
    sys.modules.setdefault("web_ui", sys.modules["__main__"])
    sys.exit(main())
