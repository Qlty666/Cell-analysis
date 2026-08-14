#!/usr/bin/env python3
"""Tests for the QC p-value comparison figure wiring."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

import web_ui  # noqa: E402


class TestQcPvalueFigure(unittest.TestCase):
    def test_figure_wired_into_web_ui(self):
        self.assertIn("fig_48_qc_pvalue_comparison.png", web_ui.FIGURE_NAMES)
        page = web_ui.render_page()
        self.assertIn("fig_48_qc_pvalue_comparison.png", page)
        self.assertIn("QC 质控差异度 P 值", page)


if __name__ == "__main__":
    unittest.main()
