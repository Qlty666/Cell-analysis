#!/usr/bin/env python3
"""Regression test: running web_ui.py as a script must not fork its state.

When ``python web/web_ui.py`` runs, the file is ``__main__``. ``web_handler``
imports it as ``web_ui``; without the alias registered in the ``__main__``
block that import creates a second, empty module, which silently drops the
``--allow-path`` roots and the auth token. This boots the real script and
checks that ``--allow-path`` reaches the handler.
"""

from __future__ import annotations

import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.parse
import urllib.request
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = APP_ROOT / "web" / "web_ui.py"


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


class TestWebScriptMode(unittest.TestCase):
    def test_allow_path_reaches_the_handler(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = root / "results" / "result_report.html"
            report.parent.mkdir(parents=True, exist_ok=True)
            report.write_text("<html>script-mode-ok</html>", encoding="utf-8")

            port = _free_port()
            proc = subprocess.Popen(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--no-browser",
                    "--port",
                    str(port),
                    "--allow-path",
                    str(root),
                ],
                cwd=str(APP_ROOT),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            try:
                deadline = time.time() + 40
                last_error: Exception | None = None
                started = False
                while time.time() < deadline:
                    if proc.poll() is not None:
                        output = proc.stdout.read() if proc.stdout else ""
                        self.fail(f"web_ui.py exited early: {output}")
                    try:
                        with urllib.request.urlopen(
                            f"http://127.0.0.1:{port}/tasks/count",
                            timeout=2,
                        ) as response:
                            if response.status == 200:
                                started = True
                                break
                    except Exception as exc:  # server not up yet
                        last_error = exc
                        time.sleep(0.3)
                self.assertTrue(started, f"server did not start: {last_error}")

                url = (
                    f"http://127.0.0.1:{port}/report?output="
                    + urllib.parse.quote(str(root))
                )
                with urllib.request.urlopen(url, timeout=5) as response:
                    self.assertEqual(response.status, 200)
                    self.assertIn(b"script-mode-ok", response.read())
            finally:
                proc.terminate()
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    proc.kill()


if __name__ == "__main__":
    unittest.main()
