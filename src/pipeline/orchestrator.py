#!/usr/bin/env python3
"""Orchestrator for the liver cancer single-cell pipeline.

Auto-downloads data, checks R dependencies, runs the Seurat analysis with
stage markers, watches for stalls every 10 minutes, diagnoses failures, and
regenerates the HTML report.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

try:
    import psutil
except ImportError:  # pragma: no cover - fallback keeps stall detection working
    psutil = None

ROOT = Path(__file__).resolve().parent.parent.parent
OUTPUT_ROOT = Path(
    os.environ.get("LIVER_OUTPUT_ROOT", str(ROOT.parent / "liver_cancer"))
).resolve()
CONFIG_PATH = ROOT / "config" / "project_config.json"
STALL_SECONDS = int(os.environ.get("LIVER_STALL_SECONDS", "1800"))
POLL_INTERVAL = 30
MAX_ATTEMPTS = 4
STALL_LOG_INTERVAL = 600

DETERMINISTIC_FAILURE_MARKERS = (
    "Error: unexpected symbol",
    "Error: unexpected ')'",
    "Error: unexpected input",
    "Error: unexpected string constant",
    "argument 1 is not a vector",
)

SINGLE_CELL_STAGES = [
    "01_load_data",
    "02_qc_filter",
    "03_doublets",
    "04_cluster",
    "05_annotation",
    "06_differential_expression",
    "07_enrichment",
    "08_publication_analyses",
    "09_summary_outputs",
]

try:
    from data.download_data import ensure_data_for_accession
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from data.download_data import ensure_data_for_accession

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


def log(message: str) -> None:
    print(f"[pipeline] {message}", flush=True)


def load_config() -> dict:
    with CONFIG_PATH.open("r", encoding="utf-8") as fh:
        return json.load(fh)


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


def install_deps() -> None:
    log("checking/installing R dependencies")
    subprocess.run(
        [find_rscript(), str(ROOT / "src" / "analysis" / "install_deps.R")],
        cwd=ROOT,
        check=True,
    )
    log("R dependencies ready")


def read_tail(path: Path, limit: int = 4000) -> str:
    if not path.exists():
        return ""
    data = path.read_text(encoding="utf-8", errors="replace")
    return data[-limit:]


def _snapshot_r_script(log_dir: Path) -> Path:
    """Copy the R script so live edits cannot corrupt a running parse."""
    log_dir.mkdir(parents=True, exist_ok=True)
    source = ROOT / "src" / "analysis" / "analysis_pipeline.R"
    snapshot = log_dir / "pipeline_analysis.R"
    shutil.copy2(source, snapshot)
    log(f"R script snapshot: {snapshot}")
    return snapshot


def _check_r_script_syntax(script: Path) -> None:
    proc = subprocess.run(
        [
            find_rscript(),
            "-e",
            "invisible(parse(file=commandArgs(TRUE)[1]))",
            str(script),
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=120,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "R script syntax check failed:\n"
            + (proc.stderr or proc.stdout)[-2000:]
        )


def _cpu_seconds(pid: int) -> float:
    if psutil is None:
        return 0.0
    try:
        proc = psutil.Process(pid)
        total = sum(proc.cpu_times()[:2])
        for child in proc.children(recursive=True):
            try:
                times = child.cpu_times()
                total += times.user + times.system
            except (psutil.Error, OSError):
                continue
        return total
    except (psutil.Error, OSError):
        return 0.0


def _process_using_cpu(pid: int, sample_seconds: float = 2.0) -> bool:
    if psutil is None:
        return False
    start = _cpu_seconds(pid)
    time.sleep(sample_seconds)
    end = _cpu_seconds(pid)
    return end - start >= 0.2


def _terminate_process_tree(proc: subprocess.Popen, timeout: float = 20) -> None:
    if psutil is None:
        proc.terminate()
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        return
    try:
        root = psutil.Process(proc.pid)
        children = root.children(recursive=True)
        for child in children:
            try:
                child.terminate()
            except (psutil.Error, OSError):
                continue
        try:
            root.terminate()
        except (psutil.Error, OSError):
            pass
        try:
            proc.wait(timeout=timeout)
            return
        except subprocess.TimeoutExpired:
            for process in [*children, root]:
                try:
                    process.kill()
                except (psutil.Error, OSError):
                    continue
            proc.wait(timeout=5)
            return
    except (psutil.Error, OSError):
        proc.terminate()
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()


def run_r_pipeline(
    accession: str,
    force: bool,
    start_stage: str = "01",
    species: str = "hs",
):
    log_dir = OUTPUT_ROOT / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "pipeline_r.log"

    script = _snapshot_r_script(log_dir)
    _check_r_script_syntax(script)
    cmd = [find_rscript(), str(script)]
    if force:
        cmd.append("--force")
    cmd.append(f"--start-stage={start_stage}")

    env = os.environ.copy()
    env["LIVER_ROOT"] = str(OUTPUT_ROOT)
    env["LIVER_ACCESSION"] = accession
    env["LIVER_SPECIES"] = species
    pause_path = OUTPUT_ROOT / "pause_request.flag"
    log(f"running R pipeline (log: {log_path})")

    with log_path.open("w", encoding="utf-8") as fh:
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )

        stalled = {"flag": False}
        paused = {"flag": False}

        def watch_stall() -> None:
            last_stall_log = float("-inf")
            while proc.poll() is None:
                if pause_path.exists():
                    log("pause requested; stopping current run")
                    paused["flag"] = True
                    _terminate_process_tree(proc)
                    return
                time.sleep(POLL_INTERVAL)
                try:
                    age = time.time() - log_path.stat().st_mtime
                except OSError:
                    continue
                if age > STALL_SECONDS:
                    if _process_using_cpu(proc.pid):
                        if age - last_stall_log >= STALL_LOG_INTERVAL:
                            log(
                                f"no log output for {int(age)}s but R is still "
                                "computing; extending stall grace"
                            )
                            last_stall_log = age
                    else:
                        log(
                            f"progress stalled for {int(age)}s; "
                            "killing and restarting"
                        )
                        stalled["flag"] = True
                        _terminate_process_tree(proc)
                        return

        watcher = threading.Thread(target=watch_stall, daemon=True)
        watcher.start()
        proc.wait()
        watcher.join(timeout=5)

    if stalled["flag"]:
        return 99, log_path
    if paused["flag"]:
        return 98, log_path
    return proc.returncode, log_path


def diagnose_failure(log_text: str) -> str | None:
    for marker in DETERMINISTIC_FAILURE_MARKERS:
        if marker in log_text:
            return f"deterministic R error: {marker}"

    missing = re.search(
        r"there is no package called ['\"]([^'\"]+)", log_text
    ) or re.search(
        r"package or namespace load failed for ['\"]([^'\"]+)", log_text
    )
    if missing:
        return f"missing R package: {missing.group(1)}"

    if "cannot open file" in log_text or "No such file or directory" in log_text:
        return "missing input file"

    return None


def verify_outputs() -> None:
    required = [
        OUTPUT_ROOT / "results" / "summary.json",
        OUTPUT_ROOT / "results" / "pipeline_complete.json",
        OUTPUT_ROOT / "results" / "data" / "05_deg" / "fig_08_deg_all.csv",
        OUTPUT_ROOT / "results" / "data" / "05_deg" / "fig_09_deg_significant.csv",
        OUTPUT_ROOT / "results" / "data" / "05_deg" / "fig_09_deg_horizontal_violin.csv",
        OUTPUT_ROOT / "results" / "data" / "01_qc" / "fig_01_qc_metrics.csv",
        OUTPUT_ROOT
        / "results"
        / "data"
        / "04_annotation"
        / "fig_05_16_17_cell_annotations.csv",
        OUTPUT_ROOT / "results" / "figures" / "05_deg" / "fig_08_volcano.png",
        OUTPUT_ROOT / "results" / "figures" / "05_deg" / "fig_09_deg_horizontal_violin.png",
        OUTPUT_ROOT / "results" / "figures" / "06_enrichment" / "fig_10_go_up.png",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise RuntimeError("Missing required outputs: " + ", ".join(missing))


def next_missing_stage(force: bool) -> str:
    if force:
        return "01"
    stage_dir = OUTPUT_ROOT / "results" / ".stages"
    if not stage_dir.exists():
        return "01"
    done = {path.stem for path in stage_dir.glob("*.done")}
    for name in SINGLE_CELL_STAGES:
        if name not in done:
            return name.split("_", 1)[0]
    return "01"


def generate_report() -> None:
    log("generating HTML report")
    env = os.environ.copy()
    env["LIVER_OUTPUT_ROOT"] = str(OUTPUT_ROOT)
    subprocess.run(
        [sys.executable, str(ROOT / "src" / "report" / "generate_report.py")],
        cwd=ROOT,
        env=env,
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "src" / "report" / "export_report.py"),
            str(OUTPUT_ROOT),
        ],
        cwd=ROOT,
        env=env,
        check=False,
    )


def run_cellchat(species: str) -> None:
    log("running optional CellChat analysis")
    subprocess.run(
        [
            find_rscript(),
            str(ROOT / "src" / "analysis" / "cellchat_analysis.R"),
            str(OUTPUT_ROOT),
            species,
        ],
        cwd=ROOT,
        check=False,
    )


def run_ml_analysis() -> None:
    log("running optional ML analysis")
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "src" / "analysis" / "ml_analysis.py"),
            str(OUTPUT_ROOT),
        ],
        cwd=ROOT,
        check=False,
    )


def run_pipeline(
    force: bool,
    skip_download: bool,
    skip_deps: bool,
    accession: str = "GSE125449",
    output_root: str | None = None,
    species: str = "hs",
) -> int:
    global OUTPUT_ROOT
    accession = accession.strip().upper()
    if not re.fullmatch(r"GSE\d+", accession):
        raise RuntimeError("GSE accession must look like GSE125449")
    if output_root:
        OUTPUT_ROOT = Path(output_root).resolve()
        OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    cfg = load_config()

    if not skip_download:
        ensure_data_for_accession(accession, cfg, OUTPUT_ROOT, log)

    if species == "auto":
        manifest_path = OUTPUT_ROOT / "data" / f"{accession.upper()}_manifest.json"
        if manifest_path.exists():
            try:
                manifest = json.loads(
                    manifest_path.read_text(encoding="utf-8")
                )
            except (OSError, ValueError):
                manifest = {}
            species = manifest.get("organism", "hs")
            if manifest.get("mode") == "bulk":
                raise RuntimeError(
                    f"{accession} is a bulk RNA-seq dataset, not "
                    "single-cell; the single-cell pipeline and the full "
                    "pipeline only support single-cell datasets. Use a "
                    "single-cell GSE accession instead."
                )
            if species not in ("hs", "mm"):
                species = "hs"
        else:
            species = "hs"

    if species not in ("hs", "mm"):
        raise RuntimeError("--species must be hs, mm, or auto")

    if not skip_deps:
        install_deps()

    complete_file = OUTPUT_ROOT / "results" / "pipeline_complete.json"
    if complete_file.exists() and not force:
        log("pipeline already complete; use --force to rerun")
        generate_report()
        return 0

    success = False
    for attempt in range(1, MAX_ATTEMPTS + 1):
        log(f"pipeline attempt {attempt}/{MAX_ATTEMPTS}")
        start_stage = next_missing_stage(force)
        log(f"resuming from stage {start_stage}")
        code, log_path = run_r_pipeline(accession, force, start_stage, species)

        if code == 98:
            log("pipeline paused; run again to resume from checkpoint")
            return 98

        if code == 99:
            log("stall detected; restarting immediately")
            continue

        if code == 0:
            success = True
            break

        tail = read_tail(log_path)
        issue = diagnose_failure(tail)
        if issue:
            log(f"issue detected: {issue}")
            if issue.startswith("deterministic"):
                print(tail[-2000:])
                raise RuntimeError(
                    "Pipeline failed with a deterministic error; "
                    "retrying cannot fix it"
                )
            if "package" in issue:
                install_deps()
        else:
            log("pipeline failed without a recognized fix")
            print(tail[-2000:])

        if attempt < MAX_ATTEMPTS:
            time.sleep(3)

    if not success:
        raise RuntimeError("Pipeline failed after repeated attempts")

    verify_outputs()
    run_ml_analysis()
    if os.environ.get("LIVER_RUN_CELLCHAT", "").lower() == "yes":
        run_cellchat(species)
    generate_report()
    log("pipeline finished successfully")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Liver cancer scRNA-seq pipeline")
    parser.add_argument(
        "accession",
        help="GEO dataset accession (e.g. GSE125449)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="output folder for results",
    )
    parser.add_argument(
        "--species",
        default="auto",
        choices=["hs", "mm", "auto"],
        help="organism: hs, mm, or auto",
    )
    parser.add_argument("--force", action="store_true", help="rerun all R stages")
    parser.add_argument("--skip-download", action="store_true", help="skip data download")
    parser.add_argument("--skip-deps", action="store_true", help="skip R dependency check")
    args = parser.parse_args()

    try:
        return run_pipeline(
            args.force,
            args.skip_download,
            args.skip_deps,
            args.accession,
            args.output,
            args.species,
        )
    except Exception as exc:  # noqa: BLE001
        log(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
