#!/usr/bin/env python3
"""Tests for the random real-GSE full pipeline validation script."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_random_real_full_pipeline as runner  # noqa: E402


class TestOnlyAccessions(unittest.TestCase):
    def test_normalized_only_accessions_accepts_gse(self):
        self.assertEqual(
            runner.normalized_only_accessions(["gse125449", "GSE1"]),
            ["GSE125449", "GSE1"],
        )

    def test_normalized_only_accessions_rejects_path_traversal(self):
        with self.assertRaises(ValueError):
            runner.normalized_only_accessions(["../../outside", "GSE1"])


if __name__ == "__main__":
    unittest.main()
