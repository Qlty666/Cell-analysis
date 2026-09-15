#!/usr/bin/env python3
"""Unit tests for the target-scoring validation script."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_new_features  # noqa: E402


class TestRequiredDatasetCount(unittest.TestCase):
    def test_full_run_keeps_twenty_dataset_requirement(self):
        self.assertEqual(
            validate_new_features.required_dataset_count(20, False),
            20,
        )
        self.assertEqual(
            validate_new_features.required_dataset_count(20, True),
            20,
        )

    def test_smoke_run_scales_to_requested_datasets(self):
        self.assertEqual(
            validate_new_features.required_dataset_count(1, True),
            1,
        )
        self.assertEqual(
            validate_new_features.required_dataset_count(1, False),
            2,
        )


class TestSkipBuild(unittest.TestCase):
    def test_limits_cached_datasets_and_avoids_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("alpha", "beta"):
                dataset_dir = root / name
                dataset_dir.mkdir()
                (dataset_dir / "dataset_summary.json").write_text(
                    json.dumps({"dataset": name, "samples": 2}),
                    encoding="utf-8",
                )
            result = {
                "dataset": "alpha",
                "samples": 2,
                "genes_scored": 1,
                "target_classes": {},
                "multidimensional_scoring": False,
                "manifest": True,
                "validation_files": [],
                "ko_seconds": 0.1,
                "validation_seconds": 0.1,
            }
            with patch.object(
                validate_new_features,
                "OUT_ROOT",
                root,
            ), patch.object(
                validate_new_features,
                "fetch_gene_map",
                side_effect=AssertionError("network call"),
            ), patch.object(
                validate_new_features,
                "fetch_gene_evidence",
                side_effect=AssertionError("network call"),
            ), patch.object(
                validate_new_features,
                "run_dataset",
                return_value=result,
            ) as run_dataset:
                code = validate_new_features.main(
                    ["--skip-build", "--skip-gse", "--max-studies", "1"]
                )
            self.assertEqual(code, 0)
            run_dataset.assert_called_once()


if __name__ == "__main__":
    unittest.main()
