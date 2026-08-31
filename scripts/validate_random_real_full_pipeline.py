#!/usr/bin/env python3
"""Random real-GSE validation for the fully automated pipeline.

The runner samples real GEO accessions from a curated pool, runs the full
pipeline for each one, repairs common environment/download/pseudobulk issues
automatically, and writes all outputs plus a validation summary to the
requested result root.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from common.env import find_rscript as _common_find_rscript  # noqa: E402

DEFAULT_RESULT_ROOT = ROOT.parent / "y3"

# Real GEO expression datasets with public count matrices small enough for
# automated validation. The pool is deliberately conservative: only datasets
# whose supplementary files expose 10x MTX, UMI matrices, or RAW tar counts.
POOL = [
    "GSE125449",
    "GSE146409",
    "GSE165816",
    "GSE181919",
    "GSE178318",
    "GSE171306",
    "GSE144735",
    "GSE156728",
    "GSE140228",
    "GSE149614",
]


def series_prefix(accession: str) -> str:
    digits = accession[3:]
    if len(digits) <= 3:
        return "GSE" + digits
    return "GSE" + digits[: len(digits) - 3] + "nnn"


def normalized_only_accessions(values: list[str]) -> list[str]:
    order = [acc.upper() for acc in values]
    invalid = [acc for acc in order if not re.fullmatch(r"GSE\d+", acc)]
    if invalid:
        raise ValueError(f"invalid GEO accessions: {', '.join(invalid)}")
    return order


def log(message: str) -> None:
    print(f"[validation] {message}", flush=True)


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def read_json(path: Path, default=None):
    if not path.exists():
        return default or {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else default or {}
    except Exception:
        return default or {}


def safe_delete(path: Path) -> None:
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path, ignore_errors=True)
    else:
        path.unlink(missing_ok=True)


def find_rscript() -> str:
    found = _common_find_rscript()
    if found is None:
        raise RuntimeError("Rscript not found")
    return found


def install_r_deps() -> bool:
    try:
        result = subprocess.run(
            [
                find_rscript(),
                str(ROOT / "src" / "analysis" / "install_deps.R"),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=1800,
        )
        return result.returncode == 0
    except Exception:
        return False


def clean_download(accession: str, single_cell_root: Path) -> None:
    safe_delete(single_cell_root / "data" / "raw" / accession)
    safe_delete(single_cell_root / "data" / f"{accession}_manifest.json")
    safe_delete(ROOT / "data_cache" / accession)


def clean_integration_markers(workdir: Path) -> None:
    stage_dir = workdir / "outputs" / "integration" / ".stages"
    safe_delete(stage_dir)


def rebuild_pseudobulk(workdir: Path) -> None:
    safe_delete(workdir / "data" / "knockout")
    stage_dir = workdir / "outputs" / "integration" / ".stages"
    for marker in [
        "04_knockout_inputs.done",
        "05_knockout.done",
        "06_docking.done",
        "07_cell_feedback.done",
        "08_report.done",
    ]:
        safe_delete(stage_dir / marker)


def diagnose(log_text: str, accession: str, single_cell_root: Path, workdir: Path) -> list[str]:
    """Return ordered repair actions based on the latest pipeline log."""
    actions: list[str] = []
    if re.search(
        r"there is no package called|package or namespace load failed for",
        log_text,
        re.IGNORECASE,
    ):
        actions.append("install-deps")
    if re.search(
        r"No count matrix files found|Downloaded file too small|"
        r"cannot open file|No such file or directory|HTTP Error|"
        r"curl.*failed|Connection timed out|Service Unavailable",
        log_text,
        re.IGNORECASE,
    ):
        actions.append("clean-download")
    if re.search(
        r"expression matrix has only|pseudobulk.*missing|"
        r"pseudobulk.*empty|exporting pseudobulk.*failed|"
        r"must contain at least two conditions",
        log_text,
        re.IGNORECASE,
    ):
        actions.append("rebuild-pseudobulk")
    if re.search(
        r"KeyError|Traceback \(most recent call last\)|"
        r"missing; run stage|stage .* crashed",
        log_text,
        re.IGNORECASE,
    ):
        actions.append("reset-integration")
    if not actions:
        actions.append("reset-integration")
    return actions


def apply_action(
    action: str,
    accession: str,
    single_cell_root: Path,
    workdir: Path,
) -> str:
    if action == "install-deps":
        if install_r_deps():
            return "installed R dependencies"
        return "R dependency install failed"
    if action == "clean-download":
        clean_download(accession, single_cell_root)
        rebuild_pseudobulk(workdir)
        return "cleaned GEO cache and pseudobulk"
    if action == "rebuild-pseudobulk":
        rebuild_pseudobulk(workdir)
        return "cleared stale pseudobulk outputs"
    if action == "reset-integration":
        clean_integration_markers(workdir)
        return "reset integration stage markers"
    return "no repair action"


def verify_passed(single_cell_root: Path, workdir: Path, skip_docking: bool) -> tuple[bool, list[str]]:
    integration = workdir / "outputs" / "integration"
    missing: list[str] = []
    required = [
        integration / "integration_summary.json",
        integration / "integration_report.html",
        integration / "key_genes.csv",
        integration / "gene_evidence.csv",
        integration / "knockout_summary.json",
        integration / "docking_summary.json",
    ]
    if not skip_docking:
        required.append(integration / "docking_targets.csv")
    for path in required:
        if not path.exists():
            missing.append(str(path.relative_to(workdir)))
    if not (single_cell_root / "results" / "pipeline_complete.json").exists():
        missing.append("single-cell pipeline_complete.json")
    return not missing, missing


def run_one(
    accession: str,
    single_cell_root: Path,
    workdir: Path,
    log_path: Path,
    args,
) -> dict:
    if verify_passed(single_cell_root, workdir, args.skip_docking)[0]:
        log(f"{accession} already passed, skipping")
        return {
            "accession": accession,
            "status": "passed",
            "attempts": [{"attempt": 0, "note": "already complete"}],
            "repair_actions": [],
            "single_cell_output": str(single_cell_root),
            "workdir": str(workdir),
            "log": str(log_path),
            "elapsed_seconds": 0.0,
            "summary": read_json(
                workdir / "outputs" / "integration" / "integration_summary.json"
            ),
        }

    started = time.time()
    attempts: list[dict] = []
    status = "failed"
    summary: dict = {}
    actions_used: list[str] = []

    for attempt in range(1, args.max_attempts + 1):
        log(f"{accession} attempt {attempt}/{args.max_attempts}")
        cmd = [
            sys.executable,
            str(ROOT / "scripts" / "run_full_pipeline.py"),
            "--accession",
            accession,
            "--output",
            str(single_cell_root),
            "--workdir",
            str(workdir),
            "--species",
            args.species,
            "--skip-deps",
            "--top-genes",
            str(args.top_genes),
            "--docking-targets",
            str(args.docking_targets),
            "--evidence-workers",
            "4",
            "--evidence-timeout",
            "180",
        ]
        if args.skip_evidence_fetch:
            cmd.append("--skip-evidence-fetch")
        if args.skip_docking:
            cmd.append("--skip-docking")

        log_path.parent.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env["LIVER_RUN_CELLCYCLE"] = "no"
        env["LIVER_RUN_CLUSTER_MARKERS"] = "no"
        env["LIVER_RUN_SIGNATURES"] = "no"
        env["LIVER_RUN_CNV"] = "no"
        env["LIVER_RUN_SINGLER"] = "no"
        try:
            with log_path.open("w", encoding="utf-8") as fh:
                proc = subprocess.run(
                    cmd,
                    cwd=ROOT,
                    stdout=fh,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=args.timeout,
                    env=env,
                )
        except subprocess.TimeoutExpired:
            attempts.append(
                {
                    "attempt": attempt,
                    "returncode": None,
                    "passed_verification": False,
                    "missing": ["timeout"],
                    "log": str(log_path),
                }
            )
            log_text = ""
            try:
                log_text = log_path.read_text(encoding="utf-8", errors="replace")
            except Exception:
                pass
            if attempt < args.max_attempts:
                actions = diagnose(
                    log_text[-40000:],
                    accession,
                    single_cell_root,
                    workdir,
                )
                for action in actions:
                    result = apply_action(
                        action,
                        accession,
                        single_cell_root,
                        workdir,
                    )
                    actions_used.append(f"attempt{attempt}:{result}")
                    log(f"{accession} auto-repair: {result}")
            continue
        ok, missing = verify_passed(single_cell_root, workdir, args.skip_docking)
        attempt_data = {
            "attempt": attempt,
            "returncode": proc.returncode,
            "passed_verification": ok,
            "missing": missing,
            "log": str(log_path),
        }
        attempts.append(attempt_data)
        if ok:
            status = "passed"
            summary = read_json(
                workdir / "outputs" / "integration" / "integration_summary.json"
            )
            break

        log_text = ""
        try:
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            pass
        if attempt >= args.max_attempts:
            break
        actions = diagnose(log_text[-40000:], accession, single_cell_root, workdir)
        for action in actions:
            result = apply_action(action, accession, single_cell_root, workdir)
            actions_used.append(f"attempt{attempt}:{result}")
            log(f"{accession} auto-repair: {result}")

    elapsed = round(time.time() - started, 1)
    return {
        "accession": accession,
        "status": status,
        "attempts": attempts,
        "repair_actions": actions_used,
        "single_cell_output": str(single_cell_root),
        "workdir": str(workdir),
        "log": str(log_path),
        "elapsed_seconds": elapsed,
        "summary": summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate the full pipeline with random real GSE datasets."
    )
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--result-root", default=str(DEFAULT_RESULT_ROOT))
    parser.add_argument("--species", default="auto", choices=["hs", "mm", "auto"])
    parser.add_argument("--top-genes", type=int, default=50)
    parser.add_argument("--docking-targets", type=int, default=1)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=10800)
    parser.add_argument("--skip-docking", action="store_true")
    parser.add_argument("--skip-evidence-fetch", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--only", nargs="*", help="explicit GSE list (overrides random)")
    args = parser.parse_args()

    result_root = Path(args.result_root).resolve()
    single_cell_root = result_root / "single_cell"
    work_root = result_root / "work"
    log_root = result_root / "logs"
    summary_root = result_root / "validation"
    for path in [single_cell_root, work_root, log_root, summary_root]:
        path.mkdir(parents=True, exist_ok=True)

    if args.only:
        try:
            order = normalized_only_accessions(args.only)
        except ValueError as exc:
            parser.error(str(exc))
    else:
        random.seed(args.seed)
        order = list(random.sample(POOL, min(args.count, len(POOL))))
        backups = [acc for acc in POOL if acc not in order]
        order += backups

    log(f"validation root: {result_root}")
    log(f"run order: {', '.join(order[:args.count])}"
        + (f" (+{len(order) - args.count} replacements)" if len(order) > args.count else ""))
    if args.dry_run:
        return 0

    results = []
    passed = 0
    for accession in order:
        if passed >= args.count:
            break
        acc_root = single_cell_root / accession
        acc_work = work_root / accession
        acc_log = log_root / f"{accession}.log"
        record = run_one(accession, acc_root, acc_work, acc_log, args)
        results.append(record)
        if record["status"] == "passed":
            passed += 1
            log(f"{accession} PASS ({record['elapsed_seconds']}s)")
        else:
            log(f"{accession} FAIL; trying next replacement if available")

    summary_path = summary_root / "validation_summary.json"
    write_json(
        summary_path,
        {
            "requested": args.count,
            "passed": passed,
            "seed": args.seed,
            "pool": POOL,
            "result_root": str(result_root),
            "finished_at": datetime.now().isoformat(timespec="seconds"),
            "results": results,
        },
    )
    csv_path = summary_root / "results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "accession",
                "status",
                "attempts",
                "elapsed_seconds",
                "single_cell_output",
                "workdir",
                "log",
                "repair_actions",
            ]
        )
        for record in results:
            writer.writerow(
                [
                    record["accession"],
                    record["status"],
                    len(record["attempts"]),
                    record["elapsed_seconds"],
                    record["single_cell_output"],
                    record["workdir"],
                    record["log"],
                    "; ".join(record["repair_actions"]),
                ]
            )

    log(f"PASSED {passed}/{args.count}")
    log(f"summary: {summary_path}")
    return 0 if passed >= args.count else 1


if __name__ == "__main__":
    sys.exit(main())
