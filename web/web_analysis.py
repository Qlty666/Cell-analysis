"""Advanced-analysis, MR and export jobs for the web console."""

from __future__ import annotations

import os
import sys
import time
import uuid
from pathlib import Path

from web_data import RESULT_FILE_SUFFIXES
from web_files import analysis_files
from web_results import (
    _extract_error,
    _log_tail,
    _notify_finished,
    record_analysis_job,
)
from web_state import (
    ANALYSIS_JOBS,
    ANALYSIS_QUEUE,
    ANALYSIS_QUEUE_LOCK,
    JOB_RECORD_LOCK,
    _drain_store,
)
from web_utils import (
    analysis_flags,
    cli_path,
    first_value,
    integer_field,
    require_existing_file,
    truthy,
)

WEB_DIR = Path(__file__).resolve().parent
APP_ROOT = WEB_DIR.parent
SCRIPTS_DIR = APP_ROOT / "scripts"
ANALYSIS_LOG_DIR = APP_ROOT / "logs" / "analysis"

ANALYSIS_KIND_LABELS = {
    "advanced": "多队列分析与靶点排序",
    "mr": "孟德尔随机化与共定位",
    "export": "分析工作区导出",
}


def _cli_path(value: str) -> str:
    return cli_path(value, APP_ROOT)


def _analysis_output_path(data: dict, kind: str) -> Path:
    if kind == "export":
        value = first_value(data, "analysis_root", "").strip()
    else:
        value = first_value(data, "output", "").strip()
    if not value:
        field = "分析工作区根目录" if kind == "export" else "输出目录"
        raise ValueError(f"{field}不能为空")
    output = Path(value).expanduser()
    if not output.is_absolute():
        output = APP_ROOT / output
    return output.resolve()


def start_analysis_job(data: dict) -> dict:
    """Start an Advanced, MR/coloc or analysis-export job."""
    kind = first_value(data, "kind", "").strip().lower()
    if kind not in ANALYSIS_KIND_LABELS:
        raise ValueError("分析类型必须是 advanced、mr 或 export")
    output = _analysis_output_path(data, kind)
    dry_run = truthy(first_value(data, "dry_run", ""))

    if kind == "advanced":
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "run_advanced_analysis.py"),
            "--output",
            str(output),
        ]
        config = _cli_path(first_value(data, "config", ""))
        expression = first_value(data, "expression", "").strip()
        if not config and not expression:
            raise ValueError("请填写高级分析配置文件或发现队列表达矩阵")
        if config:
            config = require_existing_file(
                config,
                "高级分析配置文件",
                APP_ROOT,
            )
            cmd += ["--config", config]
        for key in (
            "expression",
            "metadata",
            "cohort_name",
            "condition_column",
            "gene_column",
            "case_label",
            "control_label",
        ):
            value = first_value(data, key, "").strip()
            if value:
                if key in {"expression", "metadata"}:
                    value = require_existing_file(
                        value,
                        "发现队列表",
                        APP_ROOT,
                    )
                cmd += ["--" + key.replace("_", "-"), value]
        for key in ("feature_cap", "cv_folds", "cv_repeats", "seed"):
            value = integer_field(data, key)
            if value is not None:
                cmd += ["--" + key.replace("_", "-"), str(value)]
        cmd += analysis_flags(
            data,
            (
                "skip_r",
                "skip_ml",
                "skip_wgcna",
                "skip_immune",
                "skip_survival",
                "verbose",
            ),
        )
    elif kind == "mr":
        config = _cli_path(first_value(data, "config", ""))
        if not config:
            raise ValueError("MR/共定位配置文件不能为空")
        config = require_existing_file(
            config,
            "MR/共定位配置文件",
            APP_ROOT,
        )
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "run_mr_coloc.py"),
            "--config",
            config,
            "--output",
            str(output),
        ]
        cmd += analysis_flags(data, ("skip_r", "verbose"))
    else:
        cmd = [sys.executable, str(SCRIPTS_DIR / "export_to_analysis.py")]
        source = first_value(data, "source", "").strip()
        run = first_value(data, "run", "").strip()
        if source:
            source_path = Path(_cli_path(source))
            if not source_path.exists():
                raise ValueError(f"运行结果目录不存在：{source}")
            cmd += ["--source", str(source_path)]
        elif run:
            cmd.append(run)
        else:
            raise ValueError("请填写数据集编号、运行目录或来源目录")
        cmd += ["--analysis-root", str(output)]
        name = first_value(data, "name", "").strip()
        if name:
            cmd += ["--name", name]
        inventory_script = first_value(data, "inventory_script", "").strip()
        if inventory_script:
            cmd += [
                "--inventory-script",
                require_existing_file(
                    inventory_script,
                    "清单脚本",
                    APP_ROOT,
                ),
            ]
        cmd += analysis_flags(
            data,
            ("remember_analysis_root", "no_inventory", "dry_run"),
        )

    if kind != "export" or not dry_run:
        output.mkdir(parents=True, exist_ok=True)
    ANALYSIS_LOG_DIR.mkdir(parents=True, exist_ok=True)
    job_id = uuid.uuid4().hex[:8]
    log_path = ANALYSIS_LOG_DIR / f"{kind}_{job_id}.log"
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    info = {
        "job_id": job_id,
        "kind": kind,
        "title": ANALYSIS_KIND_LABELS[kind],
        "log": log_path,
        "proc": None,
        "started": time.time(),
        "output": output,
        "cmd": cmd,
        "env": env,
        "queued": True,
        "recorded": False,
        "notified": False,
    }
    ANALYSIS_JOBS[job_id] = info
    with ANALYSIS_QUEUE_LOCK:
        ANALYSIS_QUEUE.append(info)
    try:
        _drain_analysis_queue()
    except Exception:
        with ANALYSIS_QUEUE_LOCK:
            ANALYSIS_QUEUE[:] = [
                item for item in ANALYSIS_QUEUE if item is not info
            ]
        ANALYSIS_JOBS.pop(job_id, None)
        raise
    return {
        "job": job_id,
        "log_url": f"/analysis/log?job={job_id}",
        "status_url": f"/analysis/status?job={job_id}",
    }


def _drain_analysis_queue() -> None:
    _drain_store(ANALYSIS_QUEUE, ANALYSIS_QUEUE_LOCK)


def _analysis_status(info: dict) -> dict:
    if info.get("proc") is None:
        return {
            "running": False,
            "ok": False,
            "queued": True,
            "paused": False,
            "stage": "排队中",
            "error": "",
            "kind": info.get("kind", ""),
            "output": str(info.get("output", "")),
        }
    running = info["proc"].poll() is None
    if running:
        return {
            "running": True,
            "ok": False,
            "queued": False,
            "paused": False,
            "stage": ANALYSIS_KIND_LABELS.get(
                info.get("kind", ""),
                "分析中",
            ),
            "error": "",
            "kind": info.get("kind", ""),
            "output": str(info.get("output", "")),
        }
    ok = info["proc"].returncode == 0
    with JOB_RECORD_LOCK:
        if not info.get("recorded"):
            record_analysis_job(info, ok)
            info["recorded"] = True
    _drain_analysis_queue()
    log_text = _log_tail(info["log"], 30000)
    error = "" if ok else (
        _extract_error(log_text) or "进程已退出，请查看日志"
    )
    _notify_finished(
        info,
        "analysis",
        "高级分析",
        info.get("title", "高级分析"),
        "completed" if ok else "interrupted",
        info.get("title", ""),
        error,
        exit_code=info["proc"].returncode,
    )
    return {
        "running": False,
        "ok": ok,
        "queued": False,
        "paused": False,
        "stage": info.get("title", ""),
        "error": error,
        "kind": info.get("kind", ""),
        "output": str(info.get("output", "")),
    }


def _analysis_result_root(info: dict) -> Path:
    """Return the directory whose files are exposed by the analysis page."""
    root = Path(info["output"]).resolve()
    if info.get("kind") == "export":
        imported = root / "data" / "imported_results"
        if imported.is_dir():
            return imported.resolve()
    return root


def analysis_results(info: dict) -> dict:
    root = _analysis_result_root(info)
    files = analysis_files(root)
    if len(files) > 500:
        files = files[:500]
    return {
        "job": info.get("job_id", ""),
        "kind": info.get("kind", ""),
        "output": str(root),
        "files": files,
        "count": len(files),
    }


def _advanced_analysis_file_path(info: dict, name: str) -> Path | None:
    root = _analysis_result_root(info)
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        return None
    target = (root / relative).resolve()
    allowed_suffixes = RESULT_FILE_SUFFIXES | {
        ".json",
        ".md",
        ".html",
        ".txt",
    }
    if (
        not target.is_relative_to(root)
        or not target.is_file()
        or target.suffix.lower() not in allowed_suffixes
    ):
        return None
    return target
