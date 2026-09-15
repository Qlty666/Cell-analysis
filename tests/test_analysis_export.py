#!/usr/bin/env python3
"""Tests for exporting runs into a Codex analysis workspace."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from liverbio_suite import analysis_export  # noqa: E402


def _write_fake_run(source: Path) -> None:
    integration = source / "outputs" / "integration"
    (integration / ".stages").mkdir(parents=True)
    (integration / ".stages" / "08_report.done").write_text("{}", encoding="utf-8")
    (integration / "integration_summary.json").write_text(
        json.dumps({"single_cell": {"dataset": "GSE123456"}}),
        encoding="utf-8",
    )
    (source / "results").mkdir(parents=True)
    (source / "results" / "note.txt").write_text("ok", encoding="utf-8")


class TestAnalysisExport(unittest.TestCase):
    def test_safe_export_name_cannot_escape_directory(self):
        self.assertEqual(
            analysis_export._safe_export_name(r"..\..\outside"),
            "outside",
        )
        self.assertNotIn("/", analysis_export._safe_export_name("../evil"))
        self.assertNotIn(":", analysis_export._safe_export_name("C:\\evil"))

    def test_detect_accession_and_run_kind(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "run"
            _write_fake_run(source)
            self.assertEqual(analysis_export.detect_accession(source), "GSE123456")
            self.assertEqual(
                analysis_export.detect_run_kind(source),
                "full_pipeline",
            )

    def test_exports_without_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "run"
            _write_fake_run(source)
            analysis = base / "analysis"
            (analysis / "data").mkdir(parents=True)

            code = analysis_export.main(
                [
                    str(source),
                    "--analysis-root",
                    str(analysis),
                    "--no-inventory",
                ]
            )
            self.assertEqual(code, 0)
            destination = analysis / "data" / "imported_results" / "GSE123456"
            self.assertTrue((destination / "results" / "note.txt").exists())
            metadata = json.loads(
                (destination / "_source.json").read_text(encoding="utf-8")
            )
            self.assertEqual(metadata["dataset"], "GSE123456")
            self.assertEqual(metadata["run_kind"], "full_pipeline")

    def test_locates_run_by_accession_from_registry(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output_root = base / "roots" / "y2"
            source = output_root / "GSE123456"
            _write_fake_run(source)
            analysis = base / "analysis"
            registry_dir = analysis / "config"
            registry_dir.mkdir(parents=True)
            (registry_dir / "local_projects.json").write_text(
                json.dumps(
                    {
                        "projects": [
                            {
                                "output_roots": [
                                    {
                                        "path": str(output_root),
                                        "accessions": ["GSE123456"],
                                        "status": "completed",
                                    }
                                ]
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            (analysis / "data").mkdir(parents=True)

            code = analysis_export.main(
                ["GSE123456", "--analysis-root", str(analysis), "--dry-run"]
            )
            self.assertEqual(code, 0)
            destination = analysis / "data" / "imported_results" / "GSE123456"
            self.assertIn("GSE123456", destination.name)

    def test_dry_run_does_not_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "run"
            _write_fake_run(source)
            analysis = base / "analysis"
            (analysis / "data").mkdir(parents=True)

            code = analysis_export.main(
                [
                    str(source),
                    "--analysis-root",
                    str(analysis),
                    "--dry-run",
                ]
            )
            self.assertEqual(code, 0)
            destination = analysis / "data" / "imported_results" / "GSE123456"
            self.assertFalse(destination.exists())

    def test_requires_analysis_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "run"
            _write_fake_run(source)
            with mock.patch.object(
                analysis_export,
                "LOCAL_ANALYSIS_CONFIG",
                Path(tmp) / "analysis_workspace.json",
            ):
                code = analysis_export.main([str(source)])
            self.assertEqual(code, 2)

    def test_requires_run_or_accession(self):
        with mock.patch.object(
            analysis_export,
            "LOCAL_ANALYSIS_CONFIG",
            Path("missing_analysis_workspace.json"),
        ):
            code = analysis_export.main(["--dry-run"])
        self.assertEqual(code, 2)

    def test_remembers_analysis_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source = base / "run"
            _write_fake_run(source)
            analysis = base / "analysis"
            (analysis / "data").mkdir(parents=True)
            config_path = base / "config" / "analysis_workspace.json"

            with mock.patch.object(
                analysis_export,
                "LOCAL_ANALYSIS_CONFIG",
                config_path,
            ):
                code = analysis_export.main(
                    [
                        str(source),
                        "--analysis-root",
                        str(analysis),
                        "--remember-analysis-root",
                        "--no-inventory",
                    ]
                )
            self.assertEqual(code, 0)
            saved = json.loads(config_path.read_text(encoding="utf-8"))
            self.assertEqual(saved["analysis_root"], str(analysis.resolve()))


if __name__ == "__main__":
    unittest.main()
