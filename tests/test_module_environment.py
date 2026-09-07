#!/usr/bin/env python3
"""Tests for the per-module environment installer."""

from __future__ import annotations

import io
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "launchers") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "launchers"))

from install_environment import (  # noqa: E402
    MODULES,
    check_module,
    install_module,
    main,
    resolve_module,
)


class TestModuleEnvironment(unittest.TestCase):
    def test_all_module_entrypoints_and_requirements_exist(self):
        for name, spec in MODULES.items():
            with self.subTest(module=name):
                for relative in spec["entrypoints"]:
                    if " " in relative:
                        continue
                    self.assertTrue(
                        (APP_ROOT / relative).is_file(),
                        f"{relative} should exist",
                    )
                for relative in spec["pip_requirements"]:
                    self.assertTrue(
                        (APP_ROOT / relative).is_file(),
                        f"{relative} should exist",
                    )

    def test_module_aliases_resolve(self):
        self.assertEqual(resolve_module("single"), "expression")
        self.assertEqual(resolve_module("pipeline"), "expression")
        self.assertEqual(resolve_module("dock"), "docking")
        self.assertEqual(resolve_module("all"), "full")

    def test_unknown_module_raises(self):
        with self.assertRaises(ValueError):
            resolve_module("not-a-module")

    def test_list_prints_all_modules(self):
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["list"])
        self.assertEqual(code, 0)
        for name in MODULES:
            self.assertIn(name, out.getvalue())

    def test_web_check_dry_run(self):
        self.assertEqual(check_module("web", dry_run=True), 0)

    def test_expression_install_dry_run(self):
        self.assertEqual(
            install_module(
                "expression",
                auto_install_r=False,
                dry_run=True,
            ),
            0,
        )

    def test_full_install_dry_run(self):
        self.assertEqual(
            install_module(
                "full",
                auto_install_r=False,
                dry_run=True,
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
