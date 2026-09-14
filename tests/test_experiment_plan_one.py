#!/usr/bin/env python3
"""Focused tests for the experiment-plan-one analysis package."""

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from experiment_plan_one.common import bh_fdr, split_gene_symbol
from experiment_plan_one.classify import classify_experiment_plan_results
from experiment_plan_one.pipeline import default_config, merge_config
from experiment_plan_one.targets import _parse_swiss_target_table, make_venn_figure


class TestExperimentPlanOne(unittest.TestCase):
    def test_split_gene_symbol(self):
        self.assertEqual(split_gene_symbol("EGFR"), "EGFR")
        self.assertEqual(split_gene_symbol("EGFR /// ERBB1"), "EGFR")
        self.assertEqual(split_gene_symbol("N/A"), "")

    def test_bh_fdr(self):
        adjusted = bh_fdr([0.01, 0.02, 0.5])
        self.assertEqual(len(adjusted), 3)
        self.assertTrue((adjusted >= 0).all())
        self.assertTrue((adjusted <= 1).all())
        self.assertLessEqual(adjusted[0], adjusted[-1])

    def test_swiss_target_parser_uses_genecards_anchor(self):
        page = """
        <table id="resultTable"><tr>
        <th>Target</th><th>Gene</th><th>UniProt</th><th>ChEMBL</th>
        <th>Class</th><th>Probability*</th></tr>
        <tr>
        <td>Epidermal growth factor receptor</td>
        <td><a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=EGFR">EGFR</a></td>
        <td><a href="https://www.uniprot.org/uniprot/P00533">P00533</a></td>
        <td>CHEMBL203</td><td>Kinase</td>
        <td><span>0.2234</span></td>
        </tr></table>
        """
        frame = _parse_swiss_target_table(page)
        self.assertEqual(frame.iloc[0]["gene"], "EGFR")
        self.assertAlmostEqual(frame.iloc[0]["probability"], 0.2234)

    def test_venn_figure_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "venn.png"
            make_venn_figure(
                {
                    "A": {"A1", "A2", "A3", "SHARED"},
                    "B": {"SHARED", "B1", "B2", "B3", "B4"},
                },
                output,
                title="test",
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 100)
            svg_path = output.with_suffix(".svg")
            self.assertTrue(svg_path.exists())
            svg = svg_path.read_text(encoding="utf-8")
            number_positions = {
                value: x
                for x, value in re.findall(
                    r'<text[^>]* x="([^"]+)"[^>]*>(\d+)</text>',
                    svg,
                )
            }
            self.assertEqual(set(number_positions), {"1", "3", "4"})
            self.assertEqual(len(set(number_positions.values())), 3)

    def test_config_merge(self):
        config = merge_config(default_config(), {"ml": {"seed": 7}})
        self.assertEqual(config["ml"]["seed"], 7)
        self.assertEqual(config["ml"]["cv_folds"], 5)

    def test_classify_results_creates_plan_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "01_compound_characterization"
            source.mkdir(parents=True)
            (source / "fig1a_workflow.png").write_bytes(b"png")
            result = classify_experiment_plan_results(root)
            classified = Path(result["classified_root"])
            self.assertTrue(
                (
                    classified
                    / "Figure1_化合物表征_靶点预测与通路富集"
                    / "Fig1a_研究全局流程图.png"
                ).exists()
            )
            self.assertTrue((classified / "分类清单.csv").exists())


if __name__ == "__main__":
    unittest.main()
