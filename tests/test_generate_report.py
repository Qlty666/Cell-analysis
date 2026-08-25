#!/usr/bin/env python3
"""Tests for the joint figure + data analysis report."""

import base64
import json
import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from report import generate_report  # noqa: E402

TINY_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk"
    "+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


class TestJointAnalysisReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.results = Path(self.tmp.name) / "results"
        (self.results / "figures").mkdir(parents=True)
        (self.results / "data").mkdir(parents=True)
        (self.results / "figures" / "fig_08_volcano.png").write_bytes(TINY_PNG)
        (self.results / "data" / "fig_08_deg_all.csv").write_text(
            "gene,avg_log2FC,p_val_adj\n"
            "APOC3,6.85,0\n"
            "RPL36A,-1.74,0.0001\n",
            encoding="utf-8",
        )
        (self.results / "data" / "fig_09_deg_significant.csv").write_text(
            "gene,avg_log2FC,p_val_adj\n"
            "APOC3,6.85,0\n",
            encoding="utf-8",
        )
        (self.results / "summary.json").write_text(
            json.dumps(
                {
                    "dataset": "GSE125449",
                    "title": "HCC vs iCCA",
                    "n_cells_raw": 100,
                    "n_cells_after_qc": 90,
                    "n_cells_after_doublet_removal": 85,
                    "deg_total": 5,
                    "deg_up": 3,
                    "deg_down": 2,
                    "top_degs": [
                        {"gene": "APOC3", "avg_log2FC": 6.85, "p_val_adj": 0}
                    ],
                }
            ),
            encoding="utf-8",
        )
        self.orig_res = generate_report.RES
        generate_report.RES = self.results
        generate_report.FIG = self.results / "figures"
        generate_report.DATA = self.results / "data"

    def tearDown(self):
        generate_report.RES = self.orig_res
        generate_report.FIG = self.orig_res / "figures"
        generate_report.DATA = self.orig_res / "data"
        self.tmp.cleanup()

    def test_main_writes_joint_reports(self):
        self.assertEqual(generate_report.main(), 0)
        html = (self.results / "result_report.html").read_text(encoding="utf-8")
        self.assertIn("结果图与结果数据联合分析", html)
        self.assertIn("fig_08_deg_all.csv", html)
        self.assertIn("APOC3", html)
        md = (self.results / "result_analysis_report.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("fig_08_volcano.png", md)
        self.assertIn("APOC3", md)
        payload = json.loads(
            (self.results / "result_analysis.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(payload["figures"]), 1)
        self.assertEqual(payload["figures"][0]["file"], "figures/fig_08_volcano.png")
        self.assertEqual(payload["figures"][0]["companion_data"][0]["file"],
                         "data/fig_08_deg_all.csv")

    def test_companion_matching_uses_actual_data_files(self):
        analyses = [
            {"kind": "figure", "rel": "figures/fig_08_volcano.png"},
            {"kind": "data", "rel": "data/fig_08_deg_all.csv"},
            {"kind": "data", "rel": "data/other.csv"},
        ]
        matched = generate_report.companion_data_analyses(
            "fig_08_volcano.png", analyses
        )
        self.assertEqual([item["rel"] for item in matched],
                         ["data/fig_08_deg_all.csv"])


if __name__ == "__main__":
    unittest.main()
