#!/usr/bin/env python3
"""Tests for web UI status handling."""

from __future__ import annotations

import sys
import tempfile
import time
import unittest
from unittest import mock
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

from web_ui import FULL_JOBS, _full_status, start_full_job  # noqa: E402


class _FakeProc:
    def __init__(self, returncode: int) -> None:
        self.returncode = returncode

    def poll(self) -> int:
        return self.returncode


def _make_info(tmp: Path, returncode: int) -> dict:
    workdir = tmp / "work"
    log_path = workdir / "logs" / "web_full_test.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text(
        "[INFO] === stage 01 single_cell ===\n"
        "[INFO] full pipeline complete\n",
        encoding="utf-8",
    )
    marker_dir = workdir / "outputs" / "integration" / ".stages"
    marker_dir.mkdir(parents=True)
    for code, name in [
        ("01", "single_cell"),
        ("02", "key_targets"),
        ("03", "evidence"),
        ("04", "knockout_inputs"),
        ("05", "knockout"),
        ("06", "docking"),
        ("07", "report"),
    ]:
        (marker_dir / f"{code}_{name}.done").write_text(
            "done",
            encoding="utf-8",
        )
    return {
        "job_id": "test-job",
        "proc": _FakeProc(returncode),
        "workdir": workdir,
        "log": log_path,
        "accession": "GSE999999",
        "notified": True,
        "paused": False,
        "started": time.time(),
    }


class TestFullStatus(unittest.TestCase):
    def test_start_full_job_requires_output(self):
        with self.assertRaises(ValueError):
            start_full_job({"workdir": "somewhere", "output": ""})

    def test_start_full_job_requires_workdir(self):
        with self.assertRaises(ValueError):
            start_full_job({"workdir": "", "output": "somewhere"})

    def test_start_full_job_passes_workdir_to_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output = base / "out"
            workdir = base / "work"
            with mock.patch("web_ui._drain_full_queue"):
                result = start_full_job(
                    {
                        "output": [str(output)],
                        "workdir": [str(workdir)],
                    }
                )
            job_id = result["job"]
            try:
                cmd = FULL_JOBS[job_id]["cmd"]
                self.assertIn("--workdir", cmd)
                self.assertEqual(
                    cmd[cmd.index("--workdir") + 1],
                    str(workdir.resolve()),
                )
            finally:
                FULL_JOBS.pop(job_id, None)

    def test_completed_full_pipeline_uses_full_flow_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            status = _full_status(_make_info(Path(tmp), 0))
            self.assertTrue(status["ok"])
            self.assertFalse(status["paused"])
            self.assertEqual(status["stage"], "\u5168\u6d41\u7a0b\u5b8c\u6210")

    def test_pause_exit_code_is_reported_as_paused(self):
        with tempfile.TemporaryDirectory() as tmp:
            status = _full_status(_make_info(Path(tmp), 98))
            self.assertFalse(status["ok"])
            self.assertTrue(status["paused"])


if __name__ == "__main__":
    unittest.main()
