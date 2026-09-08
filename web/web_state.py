"""Shared runtime state for the web console.

Job registries, queues, locks and JSON history helpers live here so
``web_ui.py`` can focus on request handling. Everything in this module is
process-local state; it must not import ``web_ui`` (that would be circular).
"""

from __future__ import annotations

import json
import os
import subprocess
import threading
import time
from pathlib import Path

WEB_DIR = Path(__file__).resolve().parent
APP_ROOT = WEB_DIR.parent


JOB_STORE_MAX_RECORDS = 200
JOBS = {}
QUEUE = []
QUEUE_LOCK = threading.RLock()
HISTORY_PATH = WEB_DIR / "history.json"
HISTORY_LOCK = threading.RLock()
INSTALL_JOB = {}
FINISHED_NOTIFICATIONS: list[dict] = []
NOTIFY_LOCK = threading.Lock()
TASK_HISTORY_PATH = WEB_DIR / "task_history.json"
TASK_HISTORY_LOCK = threading.RLock()
JOB_RECORD_LOCK = threading.Lock()
HEARTBEAT_CLIENTS: dict[str, float] = {}
HEARTBEAT_LOCK = threading.Lock()
HEARTBEAT_INTERVAL_SECONDS = 5
HEARTBEAT_IDLE_TIMEOUT_SECONDS = 120
HEARTBEAT_START_GRACE_SECONDS = 20
HEARTBEAT_SHUTDOWN_GRACE_SECONDS = 5
DOCK_JOBS = {}
DOCK_QUEUE = []
DOCK_QUEUE_LOCK = threading.RLock()
DOCK_HISTORY_PATH = WEB_DIR / "dock_history.json"
DOCK_HISTORY_LOCK = threading.RLock()
MOLECULAR_DOCK_JOBS = {}
MOLECULAR_DOCK_QUEUE = []
MOLECULAR_DOCK_QUEUE_LOCK = threading.RLock()
MOLECULAR_DOCK_HISTORY_PATH = WEB_DIR / "molecular_docking_history.json"
MOLECULAR_DOCK_HISTORY_LOCK = threading.RLock()
FULL_JOBS = {}
FULL_QUEUE = []
FULL_QUEUE_LOCK = threading.RLock()
DATASET_DOWNLOAD_JOBS = {}
DATASET_DOWNLOAD_LOCK = threading.Lock()
VALIDATION_JOB = {"proc": None, "log": None, "handle": None, "started": None}


def _has_active_jobs() -> bool:
    for store in (JOBS, DOCK_JOBS, MOLECULAR_DOCK_JOBS, FULL_JOBS):
        # Iterate over a snapshot: handlers insert into these dicts from other
        # threads without holding a lock.
        for info in list(store.values()):
            proc = info.get("proc")
            if proc is None:
                if info.get("queued"):
                    return True
                continue
            exit_code = proc.poll()
            if exit_code is None or info.get("paused") or exit_code == 98:
                return True
    with DATASET_DOWNLOAD_LOCK:
        if any(info.get("running") for info in DATASET_DOWNLOAD_JOBS.values()):
            return True
    for proc in (VALIDATION_JOB.get("proc"), INSTALL_JOB.get("proc")):
        if proc is not None and proc.poll() is None:
            return True
    return False


def _prune_job_stores(max_records: int = JOB_STORE_MAX_RECORDS) -> None:
    """Drop finished job records beyond ``max_records`` so memory stays bounded.

    Running/queued/paused jobs are never removed.
    """
    for store in (JOBS, DOCK_JOBS, MOLECULAR_DOCK_JOBS, FULL_JOBS):
        if len(store) <= max_records:
            continue
        removable: list[str] = []
        for job_id, info in list(store.items()):
            if len(store) - len(removable) <= max_records:
                break
            proc = info.get("proc")
            if info.get("queued") or info.get("paused"):
                continue
            if proc is not None and proc.poll() is None:
                continue
            removable.append(job_id)
        for job_id in removable:
            store.pop(job_id, None)
    with DATASET_DOWNLOAD_LOCK:
        if len(DATASET_DOWNLOAD_JOBS) > max_records:
            for job_id, info in list(DATASET_DOWNLOAD_JOBS.items()):
                if len(DATASET_DOWNLOAD_JOBS) <= max_records:
                    break
                if info.get("running"):
                    continue
                DATASET_DOWNLOAD_JOBS.pop(job_id, None)


def _read_history_file(path: Path, lock) -> list[dict]:
    with lock:
        if not path.exists():
            return []
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return []
        return data if isinstance(data, list) else []


def _write_json_atomic(path: Path, data) -> None:
    """Write JSON via a temp file + replace so readers never see a partial file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    os.replace(tmp, path)


def _write_history_file(path: Path, records: list[dict], lock) -> None:
    """Write JSON history atomically so a crash cannot truncate the file."""
    with lock:
        _write_json_atomic(path, records)


def load_history() -> list[dict]:
    return _read_history_file(HISTORY_PATH, HISTORY_LOCK)


def save_history(records: list[dict]) -> None:
    _write_history_file(HISTORY_PATH, records, HISTORY_LOCK)


def load_dock_history() -> list[dict]:
    return _read_history_file(DOCK_HISTORY_PATH, DOCK_HISTORY_LOCK)


def save_dock_history(records: list[dict]) -> None:
    _write_history_file(DOCK_HISTORY_PATH, records, DOCK_HISTORY_LOCK)


def load_molecular_docking_history() -> list[dict]:
    return _read_history_file(MOLECULAR_DOCK_HISTORY_PATH, MOLECULAR_DOCK_HISTORY_LOCK)


def save_molecular_docking_history(records: list[dict]) -> None:
    _write_history_file(
        MOLECULAR_DOCK_HISTORY_PATH, records, MOLECULAR_DOCK_HISTORY_LOCK
    )


def _load_task_history() -> list[dict]:
    return _read_history_file(TASK_HISTORY_PATH, TASK_HISTORY_LOCK)


def _save_task_history(records: list[dict]) -> None:
    _write_history_file(TASK_HISTORY_PATH, records, TASK_HISTORY_LOCK)


def _spawn_process(info: dict) -> None:
    """Launch a queued job process; identical for every job domain."""
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


def _drain_store(queue: list, lock) -> None:
    """Start at most one queued job per domain (each domain is serial)."""
    with lock:
        for info in queue:
            if info.get("proc") is None:
                _spawn_process(info)
                break
