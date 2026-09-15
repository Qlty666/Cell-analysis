#!/usr/bin/env python3
"""CLI regression tests for scripts/validate_pipeline.py."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

APP_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = APP_ROOT / "scripts" / "validate_pipeline.py"
SCRIPTS = APP_ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_pipeline  # noqa: E402


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

    def test_runtime_failure_is_reported_as_missing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            fig_dir = out / "results" / "figures"
            fig_dir.mkdir(parents=True)
            for name in validate_pipeline.REQUIRED_FIGURES:
                (fig_dir / name).write_bytes(b"x")
            with patch.object(
                validate_pipeline,
                "valid_output",
                return_value=True,
            ), patch.object(
                validate_pipeline,
                "gsea_kegg_problem",
                return_value="GSEA KEGG unavailable",
            ):
                self.assertFalse(
                    validate_pipeline.verify_outputs(out, "GSE90001")
                )


if __name__ == "__main__":
    unittest.main()
