#!/usr/bin/env python3
"""Regression tests for the web console security and resource boundaries."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

import web_ui  # noqa: E402


class _FakeProc:
    def __init__(self, returncode: int | None) -> None:
        self.returncode = returncode

    def poll(self):
        return self.returncode


class _FakeHandler:
    def __init__(self, path: str = "/", headers: dict | None = None) -> None:
        self.path = path
        self.headers = headers or {}


class TestWorkdirAllowlist(unittest.TestCase):
    def test_rejects_paths_outside_allowed_roots(self):
        with tempfile.TemporaryDirectory() as tmp:
            outside = Path(tmp) / "secret"
            outside.mkdir()
            self.assertIsNone(web_ui._workdir_allowed(str(outside)))

    def test_accepts_configured_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            web_ui.EXTRA_WORKDIR_ROOTS.append(root)
            try:
                nested = root / "run"
                nested.mkdir()
                self.assertEqual(web_ui._workdir_allowed(str(nested)), nested)
            finally:
                web_ui.EXTRA_WORKDIR_ROOTS.remove(root)

    def test_accepts_registered_job_workdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "job"
            workdir.mkdir()
            job_id = "test-allowlist"
            web_ui.FULL_JOBS[job_id] = {"workdir": workdir}
            try:
                self.assertEqual(
                    web_ui._workdir_allowed(str(workdir)),
                    workdir.resolve(),
                )
            finally:
                web_ui.FULL_JOBS.pop(job_id, None)


class TestModelPathAllowlist(unittest.TestCase):
    def test_rejects_model_outside_search_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            search = base / "search"
            search.mkdir()
            outside = base / "model.joblib"
            outside.write_bytes(b"x")
            original = web_ui.DATASET_SEARCH_DIR
            web_ui.DATASET_SEARCH_DIR = search
            try:
                self.assertIsNone(web_ui._model_path_allowed(outside))
                inside = search / "model.joblib"
                inside.write_bytes(b"x")
                self.assertEqual(
                    web_ui._model_path_allowed(inside),
                    inside.resolve(),
                )
            finally:
                web_ui.DATASET_SEARCH_DIR = original


class TestAuthToken(unittest.TestCase):
    def test_loopback_has_no_token(self):
        original_token = web_ui.AUTH_TOKEN
        original_port = web_ui.SERVER_PORT
        try:
            self.assertIsNone(web_ui._configure_runtime("127.0.0.1", 8000))
            self.assertTrue(web_ui._request_authorized(_FakeHandler()))
        finally:
            web_ui.AUTH_TOKEN = original_token
            web_ui.SERVER_PORT = original_port

    def test_remote_bind_requires_token(self):
        original_token = web_ui.AUTH_TOKEN
        original_port = web_ui.SERVER_PORT
        try:
            token = web_ui._configure_runtime("0.0.0.0", 8123)
            self.assertTrue(token)
            self.assertFalse(web_ui._request_authorized(_FakeHandler("/")))
            self.assertTrue(
                web_ui._request_authorized(
                    _FakeHandler("/", {"X-Auth-Token": token})
                )
            )
            self.assertTrue(
                web_ui._request_authorized(_FakeHandler(f"/?token={token}"))
            )
            self.assertTrue(
                web_ui._request_authorized(
                    _FakeHandler("/", {"Cookie": f"liverbio_token={token}"})
                )
            )
            self.assertFalse(
                web_ui._request_authorized(
                    _FakeHandler("/", {"X-Auth-Token": "wrong"})
                )
            )
        finally:
            web_ui.AUTH_TOKEN = original_token
            web_ui.SERVER_PORT = original_port


class TestOriginPortCheck(unittest.TestCase):
    def test_other_port_is_rejected_when_server_port_known(self):
        original_port = web_ui.SERVER_PORT
        try:
            web_ui.SERVER_PORT = 8000
            self.assertTrue(web_ui._origin_allowed("http://127.0.0.1:8000"))
            self.assertFalse(web_ui._origin_allowed("http://127.0.0.1:9999"))
            self.assertTrue(web_ui._origin_allowed(""))
        finally:
            web_ui.SERVER_PORT = original_port


class TestJobStorePruning(unittest.TestCase):
    def test_prune_keeps_running_jobs(self):
        store = web_ui.JOBS
        created: list[str] = []
        try:
            for index in range(web_ui.JOB_STORE_MAX_RECORDS + 5):
                job_id = f"prune-{index}"
                store[job_id] = {"proc": _FakeProc(0), "queued": False}
                created.append(job_id)
            running_id = "prune-running"
            store[running_id] = {"proc": _FakeProc(None), "queued": False}
            created.append(running_id)
            web_ui._prune_job_stores()
            self.assertLessEqual(len(store), web_ui.JOB_STORE_MAX_RECORDS + 1)
            self.assertIn(running_id, store)
        finally:
            for job_id in created:
                store.pop(job_id, None)


class TestLogTail(unittest.TestCase):
    def test_reads_only_the_tail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.log"
            path.write_text("A" * 5000 + "TAIL", encoding="utf-8")
            text = web_ui._log_tail(path, limit=8)
            self.assertEqual(text, "AAAA" + "TAIL")


class TestContentTypeTable(unittest.TestCase):
    def test_known_and_unknown_suffixes(self):
        self.assertEqual(web_ui._content_type(".PNG"), "image/png")
        self.assertEqual(web_ui._content_type(".json"), "application/json")
        self.assertEqual(
            web_ui._content_type(".unknown"), "application/octet-stream"
        )


class TestAtomicJsonWrite(unittest.TestCase):
    def test_write_replaces_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "data.json"
            path.write_text('{"old": true}', encoding="utf-8")
            web_ui._write_json_atomic(path, {"new": 1})
            self.assertEqual(
                path.read_text(encoding="utf-8").strip(),
                '{\n  "new": 1\n}',
            )
            self.assertFalse((Path(tmp) / "data.json.tmp").exists())


if __name__ == "__main__":
    unittest.main()
