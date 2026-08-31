#!/usr/bin/env python3
"""Tests for ML/DL dataset search relevance ranking."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))
SCRIPTS_DIR = APP_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import dataset_search_ml  # noqa: E402


def _samples(n: int = 20) -> pd.DataFrame:
    rows = []
    for i in range(n):
        relevant = i % 2 == 0
        rows.append(
            {
                "disease": "liver cancer",
                "research_direction": "single cell RNA-seq",
                "accession": f"GSE{i + 1000}",
                "title": (
                    "HCC single cell RNA-seq tumor"
                    if relevant
                    else "Mouse kidney bulk tissue"
                ),
                "summary": (
                    "hepatocellular carcinoma scRNA-seq"
                    if relevant
                    else "mouse tissue"
                ),
                "organism": "Homo sapiens" if relevant else "Mus musculus",
                "label": 1 if relevant else 0,
            }
        )
    return pd.DataFrame(rows)


class TestDatasetSearchML(unittest.TestCase):
    def test_lexical_scores_and_rerank(self):
        rows = [
            {
                "accession": "GSE1",
                "title": "Mouse kidney",
                "summary": "mouse tissue",
            },
            {
                "accession": "GSE2",
                "title": "HCC single cell RNA-seq",
                "summary": "hepatocellular carcinoma",
            },
        ]
        scores = dataset_search_ml.lexical_scores(
            rows,
            "liver cancer",
            "single cell RNA-seq",
        )
        self.assertEqual(len(scores), 2)
        ranked = dataset_search_ml.rerank(
            rows,
            "liver cancer",
            "single cell RNA-seq",
        )
        self.assertEqual(ranked[0]["accession"], "GSE2")

    def test_train_and_evaluate(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            samples_csv = tmp_path / "samples.csv"
            _samples().to_csv(samples_csv, index=False)
            model_path = tmp_path / "model.joblib"
            dataset_search_ml.train(
                samples_csv,
                model_path,
                model_type="lr",
                seed=42,
            )
            self.assertTrue(model_path.exists())
            payload = dataset_search_ml.load_model(model_path)
            self.assertIn("vectorizer", payload)
            self.assertIn("model", payload)
            result = dataset_search_ml.evaluate(
                samples_csv,
                model_type="lr",
                seed=42,
                cv=2,
                allow_heuristic=True,
            )
            self.assertGreaterEqual(result["accuracy"], 0.0)
            self.assertIn("roc_auc", result)
            self.assertEqual(result["label_mode"], "heuristic")

    def test_train_requires_both_classes(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            samples_csv = tmp_path / "samples.csv"
            single_class = _samples().assign(label=1)
            single_class.to_csv(samples_csv, index=False)
            with self.assertRaises(ValueError):
                dataset_search_ml.train(
                    samples_csv,
                    tmp_path / "model.joblib",
                    model_type="lr",
                    seed=42,
                )

    def test_evaluate_requires_manual_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            samples_csv = Path(tmp) / "samples.csv"
            _samples().to_csv(samples_csv, index=False)
            with self.assertRaises(ValueError):
                dataset_search_ml.evaluate(
                    samples_csv,
                    model_type="lr",
                    seed=42,
                    cv=2,
                )


if __name__ == "__main__":
    unittest.main()
