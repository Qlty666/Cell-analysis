"""Tests for evidence-aware integrated target prioritisation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.target_priority import (
    build_target_priority,
    write_target_priority_summary,
)


class TestTargetPriority(unittest.TestCase):
    def test_combines_available_components_without_zeroing_missing_evidence(
        self,
    ):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidate_universe.csv"
            evidence = root / "evidence_priority.csv"
            knockout = root / "knockout.csv"
            output = root / "target_priority.csv"
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2", "GENE3"],
                    "deg_rank": [1, 2, 3],
                    "avg_log2fc": [3.0, 2.0, 1.0],
                    "p_val_adj": [0.001, 0.01, 0.05],
                    "advanced_priority_score": [0.8, 0.5, 0.2],
                }
            ).to_csv(candidates, index=False)
            pd.DataFrame(
                {
                    "target_symbol": ["GENE2", "GENE1"],
                    "priority_score": [0.90, 0.20],
                    "evidence_score": [0.95, 0.10],
                    "coverage_ratio": [0.60, 0.10],
                    "category_count": [3, 1],
                    "rank_stability": [0.95, 0.60],
                }
            ).to_csv(evidence, index=False)
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "target_score": [0.95],
                }
            ).to_csv(knockout, index=False)

            result = build_target_priority(
                candidates,
                evidence,
                output,
                knockout_csv=knockout,
            )
            ranked = result.set_index("gene")
            self.assertEqual(result.iloc[0]["gene"], "GENE2")
            self.assertEqual(ranked.loc["GENE1", "evidence_status"], "covered")
            self.assertEqual(ranked.loc["GENE3", "evidence_status"], "not_found")
            self.assertNotEqual(ranked.loc["GENE3", "integrated_score"], 0.0)
            summary = write_target_priority_summary(
                result,
                root / "target_priority_summary.json",
                evidence_hub_used=True,
                evidence_source_count=2,
            )
            self.assertEqual(
                summary["evidence_not_found"],
                1,
            )

    def test_empty_evidence_file_is_treated_as_unqueried(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidate_universe.csv"
            evidence = root / "evidence.csv"
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "deg_rank": [1],
                    "avg_log2fc": [2.0],
                    "p_val_adj": [0.01],
                }
            ).to_csv(candidates, index=False)
            evidence.write_text("", encoding="utf-8")
            result = build_target_priority(
                candidates,
                evidence,
                root / "priority.csv",
            )
            self.assertEqual(result.iloc[0]["evidence_status"], "not_queried")


if __name__ == "__main__":
    unittest.main()
