#!/usr/bin/env python3
"""Tests for the unified liverbio suite dispatcher."""

from __future__ import annotations

import io
import sys
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from liverbio_suite import cli  # noqa: E402


class TestLiverbioCli(unittest.TestCase):
    def test_project_version_matches_docking_version(self):
        docking_init = (
            cli.ROOT / "src" / "docking" / "__init__.py"
        ).read_text(encoding="utf-8")
        self.assertIn(cli.project_version(), docking_init)

    def test_all_entrypoints_exist(self):
        for name, path in cli.ENTRYPOINTS.items():
            with self.subTest(name=name):
                self.assertTrue(path.is_file(), path)

    def test_help_returns_zero(self):
        out = io.StringIO()
        with redirect_stdout(out):
            code = cli.main(["help"])
        self.assertEqual(code, 0)
        self.assertIn("liverbio <command>", out.getvalue())

    def test_version_returns_zero(self):
        out = io.StringIO()
        with redirect_stdout(out):
            code = cli.main(["version"])
        self.assertEqual(code, 0)
        self.assertIn("Liver Cancer Bioinformatics Suite", out.getvalue())

    def test_unknown_command_returns_two(self):
        err = io.StringIO()
        with redirect_stderr(err):
            code = cli.main(["not-a-command"])
        self.assertEqual(code, 2)
        self.assertIn("unknown command", err.getvalue())

    def test_forwards_arguments_to_expression_script(self):
        with mock.patch.object(cli.subprocess, "call", return_value=0) as call:
            code = cli.main(
                [
                    "expression",
                    "GSE125449",
                    "--output",
                    "../liver_cancer",
                ]
            )
        self.assertEqual(code, 0)
        cmd = call.call_args.args[0]
        self.assertEqual(cmd[0], sys.executable)
        self.assertEqual(
            Path(cmd[1]),
            cli.ENTRYPOINTS["expression"],
        )
        self.assertIn("GSE125449", cmd)

    def test_doctor_uses_environment_script(self):
        with mock.patch.object(cli, "run_script", return_value=0) as run:
            code = cli.main(["doctor", "pipeline"])
        self.assertEqual(code, 0)
        run.assert_called_once_with(cli.ENV_CHECKS["pipeline"], [])


if __name__ == "__main__":
    unittest.main()
