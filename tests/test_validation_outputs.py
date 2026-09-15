#!/usr/bin/env python3
"""Tests for shared validation output checks."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from common.validation_outputs import gsea_kegg_problem  # noqa: E402


class TestGseaKeggProblem(unittest.TestCase):
    def _write_status(self, root: Path, payload) -> None:
        path = root / "results/data/06_enrichment/fig_21_gsea_kegg_status.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(payload), encoding="utf-8")

    def test_missing_status_is_a_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            problem = gsea_kegg_problem(Path(tmp))
        self.assertIn("fig_21_gsea_kegg_status.json", problem)

    def test_unavailable_status_is_a_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._write_status(
                root,
                {"status": "unavailable", "reason": "network failure"},
            )
            problem = gsea_kegg_problem(root)
        self.assertEqual(problem, "GSEA KEGG network failure")

    def test_ok_and_skipped_statuses_are_accepted(self):
        for status in ("ok", "skipped"):
            with self.subTest(status=status):
                with tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    self._write_status(root, {"status": status})
                    self.assertIsNone(gsea_kegg_problem(root))

    def test_invalid_json_is_a_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = (
                root
                / "results/data/06_enrichment/fig_21_gsea_kegg_status.json"
            )
            path.parent.mkdir(parents=True)
            path.write_text("{", encoding="utf-8")
            problem = gsea_kegg_problem(root)
        self.assertIn("invalid JSON", problem)


if __name__ == "__main__":
    unittest.main()
