#!/usr/bin/env python3
"""Drive the web UI full pipeline for three new real GEO datasets.

Starts web/web_ui.py in a child process, keeps it alive with page heartbeats,
submits one full-pipeline job per dataset through /full/start, and then polls
all jobs until they finish so the helper owns the server lifecycle.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
WEB_UI = APP_ROOT / "web" / "web_ui.py"
OUTPUT_ROOT = Path(r"D:\AAA Liver cancer\y2")

DATASETS = [
    {
        "accession": "GSE235863",
        "title": "anti-PD-1 + lenvatinib HCC (bulk + scRNA subset)",
        "species": "hs",
    },
    {
        "accession": "GSE307831",
        "title": "HCC Foxp3-high Treg CITE-seq",
        "species": "hs",
    },
    {
        "accession": "GSE249843",
        "title": "liver cancer lymphatic metastasis scRNA-seq",
        "species": "hs",
    },
]

HEARTBEAT_INTERVAL = 4.0
POLL_INTERVAL = 20.0


def log(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    (APP_ROOT / "logs").mkdir(exist_ok=True)
    with (APP_ROOT / "logs" / "web_full_driver.log").open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def post(url: str, data: dict) -> dict:
    body = urllib.parse.urlencode(data).encode("utf-8")
    req = urllib.request.Request(url, data=body, headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def get(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def heartbeat_loop(base: str, stop: list[bool]) -> None:
    client_id = "codex-full-driver"
    while not stop[0]:
        try:
            urllib.request.urlopen(
                f"{base}/heartbeat?client={urllib.parse.quote(client_id)}",
                timeout=10,
            ).read()
        except Exception:
            pass
        time.sleep(HEARTBEAT_INTERVAL)


def submit(base: str, accession: str, species: str) -> str:
    workdir = OUTPUT_ROOT / accession
    output = workdir / "outputs"
    data = {
        "accession": accession,
        "output": str(output),
        "workdir": str(workdir),
        "species": species,
        "top_genes": "50",
        "docking_targets": "3",
        "ko_top_n": "50",
        "feedback_top_n": "12",
        "feedback_max_features": "8",
    }
    result = post(f"{base}/full/start", data)
    return str(result["job"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--poll-interval", type=float, default=POLL_INTERVAL)
    parser.add_argument("--max-hours", type=float, default=20)
    args = parser.parse_args()

    log_dir = APP_ROOT / "logs"
    log_dir.mkdir(exist_ok=True)
    server_log = (log_dir / "web_ui_driver_server.log").open("w", encoding="utf-8", errors="replace")
    server = subprocess.Popen(
        [sys.executable, str(WEB_UI), "--host", args.host, "--port", str(args.port), "--no-browser"],
        cwd=str(APP_ROOT),
        stdout=server_log,
        stderr=subprocess.STDOUT,
    )

    base = f"http://{args.host}:{args.port}"
    stop = [False]
    import threading

    threading.Thread(target=heartbeat_loop, args=(base, stop), daemon=True).start()

    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        if server.poll() is not None:
            log(f"web UI exited early with code {server.returncode}")
            return 1
        try:
            urllib.request.urlopen(base + "/full", timeout=5).read()
            break
        except Exception:
            time.sleep(1)
    else:
        log("web UI did not become ready in time")
        server.terminate()
        return 1

    jobs: dict[str, str] = {}
    attempts: dict[str, int] = {}
    for ds in DATASETS:
        job = submit(base, ds["accession"], ds["species"])
        jobs[ds["accession"]] = job
        attempts[ds["accession"]] = 0
        log(f"submitted {ds['accession']} -> job {job} ({ds['title']})")

    last_stage: dict[str, str] = {}
    finished: dict[str, bool] = {}
    started_at = time.monotonic()
    while time.monotonic() - started_at < args.max_hours * 3600:
        if server.poll() is not None:
            log(f"web UI exited with code {server.returncode} while jobs are running")
            return 1
        running_count = 0
        for accession, job in jobs.items():
            if finished.get(accession):
                continue
            try:
                status = get(f"{base}/full/status?job={job}")
            except Exception as exc:
                log(f"status query failed for {accession}: {exc!r}")
                continue
            stage = status.get("stage") or ""
            running = status.get("running", False)
            queued = status.get("queued", False)
            paused = status.get("paused", False)
            ok = status.get("ok", False)
            error = status.get("error") or ""
            if stage != last_stage.get(accession):
                last_stage[accession] = stage
                log(
                    f"{accession} {job}: stage={stage!r} running={running} "
                    f"queued={queued} paused={paused} ok={ok}"
                    + (f" error={error[:300]}" if error else "")
                )
            if running or queued:
                running_count += 1
                continue
            if not status.get("ok") and attempts[accession] < 3:
                attempts[accession] += 1
                log(f"resubmitting {accession} after failed job {job} (attempt {attempts[accession]})")
                jobs[accession] = submit(base, accession, next(
                    ds["species"] for ds in DATASETS if ds["accession"] == accession
                ))
                running_count += 1
                continue
            finished[accession] = True
        if running_count == 0:
            break
        time.sleep(args.poll_interval)

    stop[0] = True
    summary = {}
    for accession, job in jobs.items():
        try:
            summary[accession] = get(f"{base}/full/status?job={job}")
        except Exception as exc:
            summary[accession] = {"error": str(exc)}
    log("driver summary: " + json.dumps(summary, ensure_ascii=False))

    # Give the web UI a few seconds to write task history, then shut it down.
    try:
        server.terminate()
        server.wait(timeout=10)
    except Exception:
        server.kill()
    return 0


if __name__ == "__main__":
    sys.exit(main())
