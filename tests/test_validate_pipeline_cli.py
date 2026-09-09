#!/usr/bin/env python3
"""CLI regression tests for scripts/validate_pipeline.py."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = APP_ROOT / "scripts" / "validate_pipeline.py"


class TestValidatePipelineCli(unittest.TestCase):
    def test_help_exits_without_running_the_pipeline(self):
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            cwd=str(APP_ROOT),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        self.assertEqual(
            proc.returncode,
            0,
            msg=f"--help failed:\n{proc.stderr or proc.stdout}",
        )
        output = proc.stdout + proc.stderr
        self.assertIn("usage", output.lower())
        self.assertIn("--keep-output", output)


if __name__ == "__main__":
    unittest.main()
