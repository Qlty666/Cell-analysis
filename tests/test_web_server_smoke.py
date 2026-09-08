#!/usr/bin/env python3
"""End-to-end smoke tests for the real HTTP handler.

The other web tests call the request functions directly, so they never exercise
``Handler``'s routing. This module boots a real ``ThreadingHTTPServer`` on an
ephemeral loopback port and checks the routes and the security guards.
"""

from __future__ import annotations

import sys
import threading
import unittest
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

import web_handler  # noqa: E402
import web_ui  # noqa: E402


class TestWebServerSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._orig_token = web_ui.AUTH_TOKEN
        cls._orig_port = web_ui.SERVER_PORT
        web_ui.AUTH_TOKEN = None
        web_ui.SERVER_PORT = None
        ThreadingHTTPServer.allow_reuse_address = True
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), web_handler.Handler)
        cls.port = cls.server.server_address[1]
        cls.thread = threading.Thread(
            target=cls.server.serve_forever,
            daemon=True,
            name="web-ui-smoke",
        )
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=5)
        web_ui.AUTH_TOKEN = cls._orig_token
        web_ui.SERVER_PORT = cls._orig_port

    def _request(self, path: str, method: str = "GET", data: bytes | None = None):
        request = urllib.request.Request(
            f"http://127.0.0.1:{self.port}{path}",
            data=data,
            method=method,
        )
        try:
            with urllib.request.urlopen(request, timeout=15) as response:
                return response.status, response.read(), dict(response.headers)
        except urllib.error.HTTPError as exc:
            return exc.code, exc.read(), dict(exc.headers)

    def test_pages_and_status_endpoints_render(self):
        for path in (
            "/",
            "/guide",
            "/environment",
            "/datasets",
            "/dock",
            "/md-simulation",
            "/knockout",
            "/network",
            "/faers",
            "/validation",
            "/molecular-docking",
            "/full",
            "/results",
            "/tasks",
            "/result-guide",
            "/result-details",
            "/history",
            "/dock/history",
            "/molecular-docking/history",
            "/tasks/count",
            "/tasks/data",
            "/tasks/history",
        ):
            status, body, headers = self._request(path)
            self.assertEqual(status, 200, f"{path} -> {status}")
            self.assertTrue(body, path)
            if path.startswith("/tasks/"):
                self.assertIn("application/json", headers.get("Content-Type", ""))

    def test_unknown_path_returns_404(self):
        status, body, _ = self._request("/does-not-exist")
        self.assertEqual(status, 404)
        self.assertIn(b"not found", body)

    def test_nosniff_header_is_set(self):
        _, _, headers = self._request("/tasks/count")
        self.assertEqual(headers.get("X-Content-Type-Options"), "nosniff")

    def test_notifications_clear_accepts_empty_post(self):
        status, body, _ = self._request(
            "/tasks/notifications/clear", method="POST", data=b""
        )
        self.assertEqual(status, 200)
        self.assertIn(b"cleared", body)

    def test_full_file_rejects_disallowed_workdir(self):
        status, body, _ = self._request(
            "/full/file?workdir=C%3A%5CWindows&name=win.ini"
        )
        self.assertEqual(status, 403)
        self.assertIn(b"workdir not allowed", body)

    def test_report_rejects_disallowed_output(self):
        status, body, _ = self._request("/report?output=C%3A%5CWindows")
        self.assertEqual(status, 403)
        self.assertIn(b"output directory not allowed", body)

    def test_figure_without_job_returns_404(self):
        status, _, _ = self._request("/figure?job=missing&name=*.png")
        self.assertEqual(status, 404)

    def test_token_is_enforced_when_configured(self):
        """A non-loopback bind mints a token; the handler must require it."""
        web_ui.AUTH_TOKEN = "smoke-secret"
        try:
            status, body, _ = self._request("/tasks/count")
            self.assertEqual(status, 403)
            self.assertIn(b"unauthorized", body)

            status, _, headers = self._request("/tasks/count?token=smoke-secret")
            self.assertEqual(status, 200)
            self.assertIn("liverbio_token=smoke-secret", headers.get("Set-Cookie", ""))

            status, _, _ = self._request(
                "/tasks/count", method="GET", data=None
            )
            self.assertEqual(status, 403)
        finally:
            web_ui.AUTH_TOKEN = None


if __name__ == "__main__":
    unittest.main()
