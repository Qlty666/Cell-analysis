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
            self.assertEqual(
                ranked.loc["GENE1", "disease_association"],
                18.0,
            )
            self.assertEqual(ranked.loc["GENE1", "druggability"], 16.0)
            self.assertEqual(ranked.loc["GENE1", "chemical_matter"], 18.0)
            self.assertGreater(
                ranked.loc["GENE2", "safety_penalty"],
                0,
            )
            self.assertEqual(summary["targets"], 2)
            self.assertIn("axis_availability", summary)

    def test_does_not_double_count_the_same_chemical_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            integration = root / "outputs" / "integration"
            integration.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "priority_rank": [1],
                    "genetic_association": [0.9],
                    "disease_association": [0.9],
                    "direct_experimental": [0.1],
                    "curated_target": [0.0],
                    "chembl_bioactivities": [50],
                    "clinical_precedent": [0.0],
                    "structure": [0.0],
                    "evidence_status": ["covered"],
                    "direct_experimental_sources": ["ChEMBL"],
                }
            ).to_csv(
                integration / "integrated_target_priority.csv",
                index=False,
            )
            result, _ = build_target_validation(root, integration)
            row = result.iloc[0]
            self.assertEqual(row["disease_association"], 18.0)
            self.assertEqual(row["druggability"], 0.0)
            self.assertEqual(row["chemical_matter"], 20.0)
            self.assertEqual(
                row["chemical_matter_evidence_status"],
                "available",
            )
            self.assertIn(
                "ChEMBL",
                row["chemical_matter_sources"],
            )

    def test_clinicaltrials_text_match_is_capped_and_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            integration = root / "outputs" / "integration"
            integration.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "priority_rank": [1],
                    "clinical_precedent": [1.0],
                    "clinical_precedent_sources": ["ClinicalTrialsGov"],
                }
            ).to_csv(
                integration / "integrated_target_priority.csv",
                index=False,
            )
            result, _ = build_target_validation(root, integration)
            row = result.iloc[0]
            self.assertEqual(row["clinical_precedent"], 10.0)
            self.assertEqual(
                row["clinical_evidence_quality"],
                "text_search_only",
            )
            self.assertTrue(row["clinical_evidence_cap_applied"])


if __name__ == "__main__":
    unittest.main()
