"""Status, task-history and job-record helpers for the web console.

These functions only read the shared state in :mod:`web_state` and the static
stage labels in :mod:`web_data`, so they can live outside ``web_ui.py``.
"""

from __future__ import annotations

import json
import logging
import re
import time
from pathlib import Path

from web_data import FULL_STAGE_LABELS, SINGLE_STAGE_LABELS
from web_state import (
    DOCK_HISTORY_LOCK,
    FINISHED_NOTIFICATIONS,
    HISTORY_LOCK,
    MOLECULAR_DOCK_HISTORY_LOCK,
    NOTIFY_LOCK,
    TASK_HISTORY_LOCK,
    _load_task_history,
    _save_task_history,
    load_dock_history,
    load_history,
    load_molecular_docking_history,
    save_dock_history,
    save_history,
    save_molecular_docking_history,
)


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
        size = path.stat().st_size
        with path.open("rb") as handle:
            if size > limit:
                handle.seek(size - limit)
            data = handle.read(limit)
        return data.decode("utf-8", errors="replace")
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
    except Exception as exc:  # noqa: BLE001 - notification is best effort
        logging.getLogger("web_ui").warning(
            "could not append task history for %s: %s", item.get("job", "?"), exc
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
    with HISTORY_LOCK:
        records = load_history()
        records.insert(0, record)
        save_history(records)


def record_dock_job(info: dict, ok: bool) -> None:
    with DOCK_HISTORY_LOCK:
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


def record_molecular_docking_job(info: dict, ok: bool) -> None:
    with MOLECULAR_DOCK_HISTORY_LOCK:
        records = load_molecular_docking_history()
        records.insert(
            0,
            {
                "job": info["log"].stem.replace("molecular_docking_", ""),
                "stage": info.get("stage", ""),
                "workdir": str(info["workdir"]),
                "output": str(info["output_dir"]),
                "status": "success" if ok else "failed",
                "started": info["started"],
                "finished": time.time(),
            },
        )
        save_molecular_docking_history(records)


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


def _single_cell_root_from_workdir(workdir: Path) -> Path | None:
    context = _read_json(
        workdir / "outputs" / "integration" / ".stages" / "run_context.json"
    )
    value = context.get("single_cell_root")
    if not value:
        return None
    return Path(str(value)).expanduser().resolve()


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}
