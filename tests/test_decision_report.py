"""Tests for conflict-aware final target decisions."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.decision_report import build_decision_report


class TestDecisionReport(unittest.TestCase):
    def test_downgrades_high_safety_risk_and_reports_action(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "evidence_hub").mkdir()
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "decision": ["GO", "CONDITIONAL_GO"],
                    "adjusted_score": [88.0, 65.0],
                    "safety_risk": [0.0, 0.8],
                }
            ).to_csv(out / "target_validation_scores.csv", index=False)
            (out / "evidence_hub" / "evidence_hub_summary.json").write_text(
                json.dumps({"status": "completed"}),
                encoding="utf-8",
            )
            (out / "external_validation_summary.json").write_text(
                json.dumps({"status": "completed", "auroc": 0.82}),
                encoding="utf-8",
            )
            (out / "omics_qc_summary.json").write_text(
                json.dumps({"status": "completed", "gate_passed": True}),
                encoding="utf-8",
            )
            (out / "structural_quality_summary.json").write_text(
                json.dumps(
                    {
                        "gates": {
                            "positive_control": True,
                            "replicate_consensus": True,
                            "md_rmsd_stability": True,
                        }
                    }
                ),
                encoding="utf-8",
            )
            frame, summary = build_decision_report(out)
            indexed = frame.set_index("gene")
            self.assertEqual(
                indexed.loc["GENE1", "action"],
                "PROCEED_VALIDATION",
            )
            self.assertEqual(
                indexed.loc["GENE2", "final_decision"],
                "REVIEW",
            )
            self.assertIn(
                "high_safety_risk",
                indexed.loc["GENE2", "conflict_flags"],
            )
            self.assertEqual(summary["targets"], 2)

    def test_platform_conflict_does_not_rewrite_target_action(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "decision": ["GO"],
                    "adjusted_score": [88.0],
                    "safety_risk": [0.0],
                }
            ).to_csv(out / "target_validation_scores.csv", index=False)
            frame, summary = build_decision_report(out)
            row = frame.iloc[0]
            self.assertEqual(row["target_action"], "PROCEED_VALIDATION")
            self.assertEqual(row["platform_action"], "PLATFORM_REVIEW")
            self.assertEqual(row["composite_action"], "REVIEW_CONFLICT")
            self.assertEqual(
                row["platform_conflicts"].split(";")[0],
                "evidence_hub_incomplete",
            )
            self.assertEqual(
                summary["platform_action"],
                "PLATFORM_REVIEW",
            )


if __name__ == "__main__":
    unittest.main()
