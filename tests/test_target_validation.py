"""Tests for five-axis target validation scoring."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.target_validation import build_target_validation


class TestTargetValidation(unittest.TestCase):
    def test_builds_go_and_safety_adjusted_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            integration = root / "outputs" / "integration"
            integration.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "priority_rank": [1, 2],
                    "genetic_association": [0.9, 0.1],
                    "disease_association": [0.9, 0.2],
                    "direct_experimental": [0.9, 0.1],
                    "curated_target": [0.8, 0.1],
                    "clinical_precedent": [0.8, 0.0],
                    "structure": [0.9, 0.0],
                    "safety_risk": [0.0, 0.5],
                }
            ).to_csv(
                integration / "integrated_target_priority.csv",
                index=False,
            )
            result, summary = build_target_validation(root, integration)
            ranked = result.set_index("gene")
            self.assertEqual(ranked.loc["GENE1", "decision"], "GO")
            self.assertGreater(
                ranked.loc["GENE1", "adjusted_score"],
                ranked.loc["GENE2", "adjusted_score"],
            )
            self.assertGreater(
                ranked.loc["GENE2", "safety_penalty"],
                0,
            )
            self.assertEqual(summary["targets"], 2)


if __name__ == "__main__":
    unittest.main()
