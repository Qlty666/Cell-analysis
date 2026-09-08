"""HTTP request handler for the web console.

Moved out of ``web_ui.py``. This module reads the console's module-level
helpers through ``import web_ui as _ui``; ``web_ui`` imports ``Handler``
lazily inside ``main()`` so the import graph stays acyclic.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import web_ui as _ui
from web_data import ENV_MODULES, FIGURES, FIGURE_NAMES
from web_results import (
    _log_tail,
    clear_task_history,
    load_dock_history,
    load_history,
    load_molecular_docking_history,
    task_history_data,
)
from web_state import (
    DOCK_JOBS,
    FINISHED_NOTIFICATIONS,
    FULL_JOBS,
    INSTALL_JOB,
    JOBS,
    MOLECULAR_DOCK_JOBS,
    NOTIFY_LOCK,
)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass

    def _send(self, code: int, body: bytes, content_type: str) -> None:
        if content_type.startswith("text/html"):
            body = _ui._inject_heartbeat_script(body)
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store, no-cache, must-revalidate")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        if _ui.AUTH_TOKEN is not None and code == 200 and _ui._request_authorized(self):
            self.send_header(
                "Set-Cookie",
                f"liverbio_token={_ui.AUTH_TOKEN}; Path=/; SameSite=Strict; HttpOnly",
            )
        self.end_headers()
        self.wfile.write(body)

    def _send_file(self, target: Path, content_type: str) -> None:
        """Serve a result file, refusing oversized payloads."""
        try:
            if target.stat().st_size > _ui.MAX_SERVED_FILE_BYTES:
                self._send(413, b"file too large to serve", "text/plain; charset=utf-8")
                return
            body = target.read_bytes()
        except OSError as exc:
            self._send(
                404,
                f"file not readable: {exc}".encode("utf-8"),
                "text/plain; charset=utf-8",
            )
            return
        self._send(200, body, content_type)

    def _handle_heartbeat(self, params: dict) -> None:
        client_id = _ui._first(params, "client")
        if _ui._first(params, "left", "0") == "1":
            _ui.register_heartbeat(client_id)
            _ui.unregister_heartbeat(client_id)
        else:
            _ui.register_heartbeat(client_id)
        self._send(200, b'{"ok": true}', "application/json")

    def do_GET(self):
        parsed = urlparse(self.path)
        if not _ui._request_authorized(self):
            self._send(403, b"unauthorized: invalid or missing token", "text/plain; charset=utf-8")
            return
        if not _ui._fetch_site_allowed(self.headers.get("Sec-Fetch-Site", "")):
            self._send(403, b"cross-site request blocked", "text/plain; charset=utf-8")
            return
        if parsed.path == "/heartbeat":
            self._handle_heartbeat(parse_qs(parsed.query))
            return
        if parsed.path == "/":
            self._send(200, _ui.get_page().encode("utf-8"), "text/html; charset=utf-8")
            return
        if parsed.path == "/guide":
            self._send(
                200,
                _ui.render_guide_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/environment":
            self._send(
                200,
                _ui.render_environment_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/datasets":
            self._send(
                200,
                _ui.render_datasets_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/dock":
            self._send(200, _ui.render_dock_page().encode("utf-8"), "text/html; charset=utf-8")
            return
        if parsed.path == "/md-simulation":
            self._send(
                200,
                _ui.render_md_simulation_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/knockout":
            self._send(
                200,
                _ui.render_knockout_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/network":
            self._send(
                200,
                _ui.render_network_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/faers":
            self._send(
                200,
                _ui.render_faers_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/validation":
            self._send(
                200,
                _ui.render_validation_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/molecular-docking":
            self._send(
                200,
                _ui.render_molecular_docking_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/log":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            parts = []
            parts.append(_log_tail(info["log"], limit=20000))
            r_log = info["out"] / "logs" / "pipeline_r.log"
            parts.append(_log_tail(r_log, limit=20000))
            text = "\n".join(part for part in parts if part)
            self._send(200, text.encode("utf-8"), "text/plain; charset=utf-8")
            return
        if parsed.path == "/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(_ui._single_status(info)).encode("utf-8")
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
            candidate = (fig_dir / name).resolve()
            if not candidate.is_relative_to(fig_dir) or not candidate.is_file():
                self._send(404, b"figure not found", "text/plain; charset=utf-8")
                return
            target = candidate
            suffix = target.suffix.lower()
            self._send_file(target, _ui._content_type(suffix))
            return
        if parsed.path == "/report":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            output = query.get("output", [""])[0]
            info = JOBS.get(job) if job else None
            if info is None and output and not _ui._workdir_allowed(output):
                self._send(403, b"output directory not allowed", "text/plain; charset=utf-8")
                return
            target = (
                _ui._single_report_path(info)
                if info is not None
                else _ui._single_report_path(None, output)
            )
            if not target:
                self._send(404, b"report not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, "text/html; charset=utf-8")
            return
        if parsed.path == "/datasets/file":
            query = parse_qs(parsed.query)
            name = query.get("name", [""])[0]
            target = _ui.dataset_file_path(name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        if parsed.path == "/datasets/download/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            try:
                body = json.dumps(
                    _ui.dataset_download_status(job),
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
            body = json.dumps(_ui._dock_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/molecular-docking/log":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = MOLECULAR_DOCK_JOBS.get(job)
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
        if parsed.path == "/molecular-docking/status":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = MOLECULAR_DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(_ui._molecular_docking_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/molecular-docking/results":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = MOLECULAR_DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(
                _ui.molecular_docking_results(info),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/molecular-docking/file":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            name = Path(query.get("name", [""])[0])
            info = MOLECULAR_DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "text/plain; charset=utf-8")
                return
            target = _ui._molecular_docking_file_path(info, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        if parsed.path == "/molecular-docking/history":
            body = json.dumps(
                load_molecular_docking_history(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/molecular-docking/check-env":
            from docking.environment import check_environment

            self._send(
                200,
                json.dumps(
                    {"checks": check_environment()},
                    ensure_ascii=False,
                ).encode("utf-8"),
                "application/json",
            )
            return
        if parsed.path == "/full":
            self._send(
                200,
                _ui.render_full_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/results":
            self._send(
                200,
                _ui.render_results_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/result-guide":
            self._send(
                200,
                json.dumps(_ui.result_guide_data(), ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
            )
            return
        if parsed.path == "/result-details":
            self._send(
                200,
                json.dumps(_ui.result_details_data(), ensure_ascii=False).encode("utf-8"),
                "application/json; charset=utf-8",
            )
            return
        if parsed.path == "/tasks":
            self._send(
                200,
                _ui.render_tasks_page().encode("utf-8"),
                "text/html; charset=utf-8",
            )
            return
        if parsed.path == "/tasks/count":
            body = json.dumps(
                _ui.running_task_counts(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/tasks/data":
            body = json.dumps(
                _ui.running_tasks_data(),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/tasks/notifications":
            # Read-only peek: draining is a POST so a prefetch/retry cannot
            # silently swallow notifications.
            with NOTIFY_LOCK:
                items = FINISHED_NOTIFICATIONS[:]
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
            body = json.dumps(_ui._full_status(info)).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full/results":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            workdir = query.get("workdir", [""])[0]
            workdir = _ui._full_workdir_for_query(job, workdir)
            if not workdir:
                self._send(400, b"job or workdir required", "application/json")
                return
            if job not in FULL_JOBS and not _ui._workdir_allowed(workdir):
                self._send(403, b"workdir not allowed", "text/plain; charset=utf-8")
                return
            body = json.dumps(
                _ui.full_results(workdir),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/full/file":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            name = query.get("name", [""])[0]
            workdir = query.get("workdir", [""])[0]
            workdir = _ui._full_workdir_for_query(job, workdir)
            if not workdir or not name:
                self._send(400, b"workdir and name required", "text/plain; charset=utf-8")
                return
            if job not in FULL_JOBS and not _ui._workdir_allowed(workdir):
                self._send(403, b"workdir not allowed", "text/plain; charset=utf-8")
                return
            target = _ui._full_file_path(workdir, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        if parsed.path == "/dock/results":
            query = parse_qs(parsed.query)
            job = query.get("job", [""])[0]
            info = DOCK_JOBS.get(job)
            if not info:
                self._send(404, b"job not found", "application/json")
                return
            body = json.dumps(
                _ui.dock_results(info),
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
            target = _ui._dock_file_path(info, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
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
            if not _ui._workdir_allowed(workdir):
                self._send(403, b"workdir not allowed", "text/plain; charset=utf-8")
                return
            target = _ui._analysis_file_path(workdir, name, kind)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        if parsed.path == "/dock/knockout/file":
            query = parse_qs(parsed.query)
            workdir = query.get("workdir", [""])[0]
            name = query.get("name", [""])[0]
            if not _ui._workdir_allowed(workdir):
                self._send(403, b"workdir not allowed", "text/plain; charset=utf-8")
                return
            target = _ui._ko_file_path(workdir, name)
            if not target:
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        if parsed.path == "/dock/validation-report":
            body = json.dumps(
                {
                    "report": _ui.validation_report_text(),
                    "exists": _ui.VALIDATION_REPORT_PATH.exists(),
                },
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/dock/validation-status":
            body = json.dumps(
                _ui.validation_job_status(),
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
        if parsed.path == "/environment/check":
            query = parse_qs(parsed.query)
            module = _ui._first(query, "module", "").strip().lower()
            if module not in ENV_MODULES:
                self._send(
                    400,
                    json.dumps(
                        {"error": f"unknown module: {module}"},
                        ensure_ascii=False,
                    ).encode("utf-8"),
                    "application/json",
                )
                return
            with_ml = _ui._first(query, "with_ml", "").strip().lower() in (
                "1",
                "yes",
                "true",
                "on",
            )
            body = json.dumps(
                _ui.run_environment_check(module, with_ml=with_ml),
                ensure_ascii=False,
            ).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/dock-environment":
            try:
                result = subprocess.run(
                    [
                        sys.executable,
                        str(_ui.APP_ROOT / "launchers" / "check_dock_environment.py"),
                    ],
                    cwd=_ui.APP_ROOT,
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
                    "module": "",
                }).encode("utf-8")
                self._send(200, body, "application/json")
                return
            running = proc.poll() is None
            log_text = ""
            if _ui.INSTALL_LOG.exists():
                log_text = _ui.INSTALL_LOG.read_text(
                    encoding="utf-8",
                    errors="replace",
                )[-6000:]
            body = json.dumps({
                "running": running,
                "ok": not running and proc.returncode == 0,
                "log": log_text,
                "module": INSTALL_JOB.get("module", ""),
            }, ensure_ascii=False).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path == "/versions":
            body = json.dumps(_ui.get_versions(), ensure_ascii=False).encode("utf-8")
            self._send(200, body, "application/json")
            return
        if parsed.path.startswith("/static/"):
            name = Path(parsed.path.removeprefix("/static/")).name
            target = (_ui.STATIC_DIR / name).resolve()
            if not target.is_file() or target.parent != _ui.STATIC_DIR.resolve():
                self._send(404, b"file not found", "text/plain; charset=utf-8")
                return
            self._send_file(target, _ui._content_type(target.suffix))
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self):
        parsed = urlparse(self.path)
        if not _ui._request_authorized(self):
            self._send(403, b"unauthorized: invalid or missing token", "text/plain; charset=utf-8")
            return
        if not _ui._fetch_site_allowed(self.headers.get("Sec-Fetch-Site", "")):
            self._send(403, b"cross-site request blocked", "text/plain; charset=utf-8")
            return
        if not _ui._origin_allowed(self.headers.get("Origin", "")):
            self._send(403, b"cross-origin request blocked", "text/plain; charset=utf-8")
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            self._send(400, b"invalid Content-Length", "text/plain; charset=utf-8")
            return
        if length < 0 or length > _ui.MAX_POST_BODY_BYTES:
            self._send(413, b"request body too large", "text/plain; charset=utf-8")
            return
        body = self.rfile.read(length).decode("utf-8", "replace")
        data = parse_qs(body)

        if parsed.path == "/heartbeat":
            params = parse_qs(parsed.query)
            params.update(data)
            self._handle_heartbeat(params)
            return

        if parsed.path == "/tasks/notifications/clear":
            with NOTIFY_LOCK:
                items = FINISHED_NOTIFICATIONS[:]
                FINISHED_NOTIFICATIONS.clear()
            self._send(
                200,
                json.dumps({"cleared": len(items)}, ensure_ascii=False).encode("utf-8"),
                "application/json",
            )
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
                result = _ui.resume_job(job)
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
                cwd=_ui.APP_ROOT,
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

        if parsed.path == "/molecular-docking/pause":
            job = data.get("job", [""])[0]
            info = MOLECULAR_DOCK_JOBS.get(job)
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

        if parsed.path == "/molecular-docking/resume":
            job = data.get("job", [""])[0]
            info = MOLECULAR_DOCK_JOBS.get(job)
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
                cwd=_ui.APP_ROOT,
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
                cwd=_ui.APP_ROOT,
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
                result = _ui.start_full_job(data)
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
                result = _ui.dataset_search_request(data)
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
                result = _ui.start_dataset_download(data)
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

        if parsed.path == "/molecular-docking/check-env":
            from docking.environment import check_environment

            checks = check_environment()
            self._send(
                200,
                json.dumps({"checks": checks}, ensure_ascii=False).encode("utf-8"),
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
                target = _ui._first(data, "target", "")
                cmd = [
                    sys.executable,
                    str(_ui.APP_ROOT / "launchers" / "install_environment.py"),
                    "install",
                    "docking",
                ]
                if target:
                    cmd += ["--target", target]
                proc = subprocess.run(
                    cmd,
                    cwd=_ui.APP_ROOT,
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

        if parsed.path == "/molecular-docking/install":
            try:
                target = _ui._first(data, "target", "")
                cmd = [
                    sys.executable,
                    str(_ui.APP_ROOT / "launchers" / "install_environment.py"),
                    "install",
                    "molecular-docking",
                ]
                if target:
                    cmd += ["--target", target]
                proc = subprocess.run(
                    cmd,
                    cwd=_ui.APP_ROOT,
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

        if parsed.path == "/molecular-docking/detect-box":
            try:
                workdir = _ui._first(data, "workdir", "")
                receptor = _ui._first(data, "receptor", "")
                if not workdir or not receptor:
                    raise ValueError("workdir and receptor are required")
                result = _ui._detect_and_save_molecular_box(workdir, receptor)
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

        if parsed.path == "/dock/detect-box":
            try:
                workdir = _ui._first(data, "workdir", "")
                receptor = _ui._first(data, "receptor", "")
                if not workdir or not receptor:
                    raise ValueError("workdir and receptor are required")
                result = _ui._detect_and_save_box(workdir, receptor)
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
                result = _ui.run_network_request(data)
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
                result = _ui.run_faers_request(data)
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
                result = _ui.run_knockout_request(data)
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
                result = _ui.run_validation_request(data)
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
                result = _ui.start_validation_job(data)
                body = json.dumps(result, ensure_ascii=False).encode("utf-8")
                self._send(200, body, "application/json")
            except Exception as exc:
                body = json.dumps(
                    {"error": str(exc)}, ensure_ascii=False
                ).encode("utf-8")
                self._send(400, body, "application/json")
            return

        if parsed.path == "/molecular-docking/start":
            try:
                result = _ui.start_molecular_docking_job(data)
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

        if parsed.path == "/dock/start":
            try:
                result = _ui.start_dock_job(data)
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
                    "module": INSTALL_JOB.get("module", ""),
                }).encode("utf-8"), "application/json")
                return
            project = _ui._first(data, "project", "single")
            target = _ui._first(data, "target", "")
            module = {
                "single": "expression",
                "expression": "expression",
                "dock": "docking",
                "docking": "docking",
                "molecular-docking": "molecular-docking",
                "md": "md",
                "full": "full",
                "datasets": "datasets",
                "web": "web",
                "skills": "skills",
            }.get(project, project)
            if module not in ENV_MODULES:
                self._send(
                    400,
                    json.dumps(
                        {"error": f"unknown module: {module}"},
                        ensure_ascii=False,
                    ).encode("utf-8"),
                    "application/json",
                )
                return
            cmd = [
                sys.executable,
                str(_ui.APP_ROOT / "launchers" / "install_environment.py"),
                "install",
                module,
            ]
            if target:
                cmd += ["--target", target]
            if _ui._first(data, "with_ml", "").strip().lower() in (
                "1",
                "yes",
                "true",
                "on",
            ):
                cmd += ["--with-ml"]
            log_handle = _ui.INSTALL_LOG.open("w", encoding="utf-8", errors="replace")
            proc = subprocess.Popen(
                cmd,
                cwd=_ui.APP_ROOT,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            INSTALL_JOB["proc"] = proc
            INSTALL_JOB["project"] = project
            INSTALL_JOB["module"] = module
            self._send(200, json.dumps({
                "running": True,
                "message": f"{project} 环境自动补全已启动",
                "module": module,
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
            "LIVER_QC_MAX_RIBO",
            "LIVER_QC_MAX_HB",
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
            result = _ui.start_job(
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
