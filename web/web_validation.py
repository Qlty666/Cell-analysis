"""Validation jobs and report rendering for the web console."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

from web_state import VALIDATION_JOB
from web_utils import first_value, integer_field, truthy

WEB_DIR = Path(__file__).resolve().parent
APP_ROOT = WEB_DIR.parent
SCRIPTS_DIR = APP_ROOT / "scripts"

VALIDATION_ROOT = Path(
    os.environ.get(
        "LIVER_VALIDATION_ROOT",
        str(APP_ROOT / "data_cache" / "validation_runs"),
    )
).resolve()
VALIDATION_REPORT_DIR = VALIDATION_ROOT / "validation"
VALIDATION_REPORT_PATH = VALIDATION_REPORT_DIR / "validation_summary.json"
VALIDATION_LOG = WEB_DIR / "validation_run.log"
VALIDATION_FEATURES_ROOT = (
    APP_ROOT / "dock" / "validation_real" / "pan_cancer_20"
)
VALIDATION_EVIDENCE_ROOT = APP_ROOT / "dock" / "validation_real"
VALIDATION_RANDOM_EVIDENCE_ROOT = (
    APP_ROOT / "dock" / "validation_real_random"
)

VALIDATION_KIND_LABELS = {
    "random": "随机真实 GSE 全流程验证",
    "features": "多队列靶点功能验证",
    "evidence": "真实 PDB 证据验证",
    "random-evidence": "随机真实证据与对接盒验证",
}


def _validation_mode(value: object = "random") -> str:
    mode = str(value or "random").strip().lower()
    return mode if mode in VALIDATION_KIND_LABELS else "random"


def _validation_report_path(mode: str) -> Path:
    mode = _validation_mode(mode)
    if mode == "features":
        return VALIDATION_FEATURES_ROOT / "validation_report.md"
    if mode == "evidence":
        return VALIDATION_EVIDENCE_ROOT / "summary.json"
    if mode == "random-evidence":
        return VALIDATION_RANDOM_EVIDENCE_ROOT / "summary.json"
    return VALIDATION_REPORT_PATH


def validation_report_text(mode: str = "random") -> str:
    mode = _validation_mode(mode)
    report_path = _validation_report_path(mode)
    if mode in {"evidence", "random-evidence"} and report_path.exists():
        try:
            data = json.loads(
                report_path.read_text(encoding="utf-8", errors="replace")
            )
        except Exception:
            data = {}
        lines = [
            f"# {VALIDATION_KIND_LABELS[mode]}",
            "",
            f"- 状态：{data.get('status', '')}",
            f"- 通过：{data.get('passed', '')}",
            f"- 最少成功靶点：{data.get('min_ok_targets', '')}",
            f"- 最少配体记录：{data.get('min_ligands', '')}",
            f"- 完成时间：{data.get('finished_at', '')}",
            "",
            "## 原始汇总",
            "",
            "```json",
            json.dumps(data, ensure_ascii=False, indent=2)[-30000:],
            "```",
        ]
        return "\n".join(lines)
    if mode == "features" and report_path.exists():
        return report_path.read_text(encoding="utf-8", errors="replace")
    if mode == "random" and VALIDATION_REPORT_PATH.exists():
        try:
            data = json.loads(
                VALIDATION_REPORT_PATH.read_text(
                    encoding="utf-8",
                    errors="replace",
                )
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
        f"请运行“{VALIDATION_KIND_LABELS[mode]}”，或在命令行运行对应验证脚本。"
    )


def start_validation_job(data: dict | None = None) -> dict:
    proc = VALIDATION_JOB.get("proc")
    if proc is not None and proc.poll() is None:
        return {"running": True, "message": "验证任务已在运行"}
    old_handle = VALIDATION_JOB.get("handle")
    if old_handle is not None and not old_handle.closed:
        old_handle.close()

    payload = data or {}
    mode = _validation_mode(first_value(payload, "mode", "random"))
    count = integer_field(payload, "count") or 10
    count = max(1, min(20 if mode == "features" else 30, count))
    raw_seed = first_value(payload, "seed", "20260813") or "20260813"
    try:
        seed = int(raw_seed)
    except (TypeError, ValueError):
        seed = 20260813

    VALIDATION_LOG.parent.mkdir(parents=True, exist_ok=True)
    handle = VALIDATION_LOG.open("w", encoding="utf-8", errors="replace")
    if mode == "random":
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "validate_random_real_full_pipeline.py"),
            "--result-root",
            str(VALIDATION_ROOT),
            "--count",
            str(count),
            "--seed",
            str(seed),
        ]
    elif mode == "features":
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "validate_new_features.py"),
            "--max-studies",
            str(count),
        ]
        if truthy(first_value(payload, "skip_gse", "")):
            cmd.append("--skip-gse")
        if truthy(first_value(payload, "skip_build", "")):
            cmd.append("--skip-build")
    elif mode == "evidence":
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "validate_real_evidence.py"),
        ]
    else:
        cmd = [
            sys.executable,
            str(SCRIPTS_DIR / "validate_real_random.py"),
        ]
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=APP_ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except Exception:
        handle.close()
        raise
    VALIDATION_JOB.update(
        {
            "proc": proc,
            "log": VALIDATION_LOG,
            "handle": handle,
            "started": time.time(),
            "mode": mode,
            "label": VALIDATION_KIND_LABELS[mode],
        }
    )
    return {
        "running": True,
        "job": "validation",
        "mode": mode,
        "label": VALIDATION_KIND_LABELS[mode],
    }


def validation_job_status() -> dict:
    proc = VALIDATION_JOB.get("proc")
    if proc is None:
        return {
            "running": False,
            "ok": False,
            "started": False,
            "log": "",
            "mode": "",
            "label": "",
        }
    running = proc.poll() is None
    handle = VALIDATION_JOB.get("handle")
    if not running and handle is not None and not handle.closed:
        handle.close()
    log_text = ""
    if VALIDATION_LOG.exists():
        log_text = VALIDATION_LOG.read_text(
            encoding="utf-8",
            errors="replace",
        )[-6000:]
    return {
        "running": running,
        "ok": not running and proc.returncode == 0,
        "started": True,
        "log": log_text,
        "mode": VALIDATION_JOB.get("mode", "random"),
        "label": VALIDATION_JOB.get(
            "label",
            VALIDATION_KIND_LABELS["random"],
        ),
    }
