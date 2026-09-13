#!/usr/bin/env python3
"""Focused tests for the experiment-plan-one analysis package."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from experiment_plan_one.common import bh_fdr, split_gene_symbol
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
                {"A": {"EGFR", "STAT3"}, "B": {"STAT3", "TNF"}},
                output,
                title="test",
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 100)

    def test_config_merge(self):
        config = merge_config(default_config(), {"ml": {"seed": 7}})
        self.assertEqual(config["ml"]["seed"], 7)
        self.assertEqual(config["ml"]["cv_folds"], 5)


if __name__ == "__main__":
    unittest.main()
