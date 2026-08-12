#!/usr/bin/env python3
"""Tests for the random GEO dataset search validation script."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))
SCRIPTS_DIR = APP_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import validate_dataset_search  # noqa: E402


def _row(accession: str, title: str, summary: str) -> dict:
    return {
        "accession": accession,
        "disease": "",
        "research_direction": "",
        "title": title,
        "summary": summary,
        "organism": "Homo sapiens",
        "samples": "4",
        "platform": "GPL1",
        "date": "",
        "type": "",
        "url": "",
    }


class TestValidateDatasetSearch(unittest.TestCase):
    def test_is_relevant(self):
        row = _row(
            "GSE1",
            "HCC single cell RNA-seq",
            "hepatocellular carcinoma",
        )
        self.assertTrue(
            validate_dataset_search.is_relevant(
                row,
                "liver cancer",
                "single cell RNA-seq",
            )
        )
        row = _row("GSE2", "Mouse kidney", "mouse tissue")
        self.assertFalse(
            validate_dataset_search.is_relevant(
                row,
                "liver cancer",
                "single cell RNA-seq",
            )
        )

    def test_run_round_found_without_expansion(self):
        rows = [
            _row(
                "GSE1",
                "HCC single cell RNA-seq",
                "hepatocellular carcinoma",
            )
        ]
        with patch.object(
            validate_dataset_search.search_datasets,
            "search_datasets",
            return_value=rows,
        ):
            record = validate_dataset_search.run_round(
                "liver cancer",
                "single cell RNA-seq",
                max_results=5,
                expand=True,
            )
        self.assertTrue(record["found"])
        self.assertFalse(record["expanded"])
        self.assertEqual(record["first_accession"], "GSE1")

    def test_run_round_expands_when_combined_misses(self):
        irrelevant = [_row("GSE1", "Mouse kidney", "mouse tissue")]
        relevant = [
            _row(
                "GSE2",
                "Liver cancer scRNA-seq",
                "hepatocellular carcinoma",
            )
        ]
        calls = {"count": 0}

        def fake_search(query, max_results, disease, research_direction):
            calls["count"] += 1
            if calls["count"] == 1:
                return irrelevant
            return relevant

        with patch.object(
            validate_dataset_search.search_datasets,
            "search_datasets",
            side_effect=fake_search,
        ):
            record = validate_dataset_search.run_round(
                "liver cancer",
                "single cell RNA-seq",
                max_results=5,
                expand=True,
            )
        self.assertTrue(record["found"])
        self.assertTrue(record["expanded"])
        self.assertEqual(record["first_accession"], "GSE2")


if __name__ == "__main__":
    unittest.main()
