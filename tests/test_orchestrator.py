#!/usr/bin/env python3
"""Tests for the single-cell pipeline orchestrator helpers."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from pipeline import orchestrator  # noqa: E402


class TestOrchestratorHelpers(unittest.TestCase):
    def test_snapshot_r_script_matches_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            snapshot = orchestrator._snapshot_r_script(Path(tmp))
            source = (
                orchestrator.ROOT / "src" / "analysis" / "analysis_pipeline.R"
            )
            self.assertTrue(snapshot.exists())
            self.assertEqual(snapshot.read_bytes(), source.read_bytes())

    def test_cpu_seconds_is_non_negative_for_current_process(self):
        self.assertGreaterEqual(orchestrator._cpu_seconds(os.getpid()), 0)

    def test_cpu_seconds_returns_zero_for_missing_process(self):
        self.assertEqual(orchestrator._cpu_seconds(999_999_999), 0.0)

    def test_terminate_process_tree_kills_children(self):
        if orchestrator.psutil is None:
            self.skipTest("psutil not available")
        parent_code = (
            "import subprocess, sys, time\n"
            "child = subprocess.Popen("
            "[sys.executable, '-c', 'import time; time.sleep(60)'])\n"
            "print(child.pid, flush=True)\n"
            "time.sleep(60)\n"
        )
        proc = subprocess.Popen(
            [sys.executable, "-c", parent_code],
            stdout=subprocess.PIPE,
            text=True,
        )
        child_pid = int(proc.stdout.readline().strip())
        proc.stdout.close()
        self.assertTrue(orchestrator.psutil.pid_exists(child_pid))
        orchestrator._terminate_process_tree(proc, timeout=10)
        self.assertIsNotNone(proc.poll())
        self.assertFalse(orchestrator.psutil.pid_exists(child_pid))


if __name__ == "__main__":
    unittest.main()
