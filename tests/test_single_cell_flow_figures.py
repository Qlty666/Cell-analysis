#!/usr/bin/env python3
"""Wiring tests for the WeChat single-cell flow figure additions."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

import web_ui  # noqa: E402
from report import generate_report  # noqa: E402


NEW_FIGURES = {
    "fig_49_qc_umi_feature_correlation.png",
    "fig_50_marker_ridgeplot.png",
    "fig_51_marker_stacked_violin.png",
}


class TestSingleCellFlowFigures(unittest.TestCase):
    def test_r_script_has_relation_qc_and_marker_views(self):
        script = (
            APP_ROOT / "src" / "analysis" / "analysis_pipeline.R"
        ).read_text(encoding="utf-8")
        for name in NEW_FIGURES:
            self.assertIn(name, script)
        self.assertIn(
            "fig_49_qc_umi_feature_correlation_stats.csv",
            script,
        )

    def test_web_and_report_wiring(self):
        for name in NEW_FIGURES:
            self.assertIn(name, web_ui.FIGURE_NAMES)
            self.assertIn(name, generate_report.FIGURE_GUIDE)
        template = (
            APP_ROOT
            / "web"
            / "templates"
            / "results_manifest_optimized.html"
        ).read_text(encoding="utf-8")
        for name in NEW_FIGURES:
            self.assertIn(name, template)
        details = json.loads(
            (APP_ROOT / "web" / "static" / "result_details.json").read_text(
                encoding="utf-8"
            )
        )
        names = {entry["name"] for entry in details["entries"]}
        for name in NEW_FIGURES:
            self.assertIn(name, names)
        self.assertIn(
            "fig_49_qc_umi_feature_correlation_stats.csv",
            names,
        )


if __name__ == "__main__":
    unittest.main()
