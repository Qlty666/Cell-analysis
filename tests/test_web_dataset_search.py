#!/usr/bin/env python3
"""Tests for dataset search web integration."""

from __future__ import annotations

import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

import web_ui  # noqa: E402


def _fake_search_module():
    module = types.ModuleType("search_datasets")

    def build_query(disease, direction, query=None):
        if query:
            return query
        return " ".join(part for part in [disease, direction] if part)

    def search_fn(
        query,
        max_results=20,
        organism=None,
        keyword=None,
        disease="",
        research_direction="",
    ):
        return [
            {
                "accession": "GSE1",
                "disease": disease,
                "research_direction": research_direction,
                "title": "title",
                "summary": "summary",
                "organism": organism or "Homo sapiens",
                "samples": "3",
                "platform": "GPL1",
                "date": "2026",
                "type": "GSE",
                "url": "https://example.com/GSE1",
            }
        ]

    def write_outputs(rows, out_dir):
        out_dir.mkdir(parents=True, exist_ok=True)
        csv_path = out_dir / "dataset_search_results.csv"
        json_path = out_dir / "dataset_search_results.json"
        csv_path.write_text("accession\nGSE1\n", encoding="utf-8")
        json_path.write_text("[]", encoding="utf-8")
        return csv_path, json_path

    module.build_query = build_query
    module.search_datasets = search_fn
    module.write_outputs = write_outputs
    return module


class TestDatasetSearchRequest(unittest.TestCase):
    def test_requires_disease_or_query(self):
        with self.assertRaises(ValueError):
            web_ui.dataset_search_request({})

    def test_search_builds_query_and_returns_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            fake = _fake_search_module()
            with mock.patch.dict(sys.modules, {"search_datasets": fake}), mock.patch.object(
                web_ui, "DATASET_SEARCH_DIR", base / "search"
            ):
                result = web_ui.dataset_search_request(
                    {
                        "disease": ["liver cancer"],
                        "research_direction": ["single cell RNA-seq"],
                        "max_results": ["10"],
                    }
                )
            self.assertEqual(result["query"], "liver cancer single cell RNA-seq")
            self.assertEqual(result["count"], 1)
            self.assertFalse(result["model_applied"])
            self.assertTrue(result["csv_url"].startswith("/datasets/file?name="))
            self.assertTrue(
                (base / "search" / "dataset_search_results.csv").exists()
            )

    def test_model_rerank_applied_when_model_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            model_file = base / "model.joblib"
            model_file.write_bytes(b"model")
            fake_search = _fake_search_module()
            fake_ml = types.ModuleType("dataset_search_ml")

            def load_model(path):
                return {"path": str(path)}

            def rerank(rows, disease, direction, model=None):
                return [dict(row, relevance_score=0.9) for row in rows]

            fake_ml.load_model = load_model
            fake_ml.rerank = rerank
            with mock.patch.dict(
                sys.modules,
                {
                    "search_datasets": fake_search,
                    "dataset_search_ml": fake_ml,
                },
            ), mock.patch.object(web_ui, "DATASET_SEARCH_DIR", base / "search"):
                result = web_ui.dataset_search_request(
                    {
                        "disease": ["liver cancer"],
                        "model": [str(model_file)],
                    }
                )
            self.assertTrue(result["model_applied"])
            self.assertEqual(result["results"][0]["relevance_score"], 0.9)
            self.assertEqual(result["model_path"], str(model_file))

    def test_model_missing_raises(self):
        fake_search = _fake_search_module()
        with mock.patch.dict(sys.modules, {"search_datasets": fake_search}), mock.patch.object(
            web_ui, "DATASET_SEARCH_DIR", Path(".")
        ):
            with self.assertRaises(ValueError):
                web_ui.dataset_search_request(
                    {
                        "disease": ["liver cancer"],
                        "model": ["missing.joblib"],
                    }
                )


class _FakeThread:
    def __init__(self, target, args, daemon=False):
        self.target = target
        self.args = args
        self.daemon = daemon
        self.started = False

    def start(self):
        self.started = True


class TestDatasetDownload(unittest.TestCase):
    def test_requires_accession(self):
        with self.assertRaises(ValueError):
            web_ui.start_dataset_download({})

    def test_starts_background_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            captured = []

            def fake_thread(*args, **kwargs):
                thread = _FakeThread(*args, **kwargs)
                captured.append(thread)
                return thread

            with mock.patch.object(web_ui, "DATASET_SEARCH_DIR", base / "search"), mock.patch.object(
                web_ui.threading, "Thread", side_effect=fake_thread
            ):
                result = web_ui.start_dataset_download(
                    {
                        "accessions": ["gse1,GSE2"],
                        "download_root": [str(base / "out")],
                    }
                )
            self.assertIn("job", result)
            job = result["job"]
            thread = captured[0]
            self.assertTrue(thread.started)
            self.assertEqual([acc.upper() for acc in thread.args[1]], ["GSE1", "GSE2"])
            self.assertEqual(thread.args[2], (base / "out").resolve())
            status = web_ui.dataset_download_status(job)
            self.assertTrue(status["running"])
            with web_ui.DATASET_DOWNLOAD_LOCK:
                web_ui.DATASET_DOWNLOAD_JOBS.pop(job, None)

    def test_status_unknown_job_raises(self):
        with self.assertRaises(ValueError):
            web_ui.dataset_download_status("missing")

    def test_status_reads_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            log = base / "download.log"
            log.write_text("ok", encoding="utf-8")
            info = {
                "log": log,
                "running": False,
                "error": "",
                "results": {"GSE1": "ok"},
            }
            with web_ui.DATASET_DOWNLOAD_LOCK:
                web_ui.DATASET_DOWNLOAD_JOBS["j1"] = info
            try:
                status = web_ui.dataset_download_status("j1")
                self.assertFalse(status["running"])
                self.assertTrue(status["ok"])
                self.assertEqual(status["log"], "ok")
            finally:
                with web_ui.DATASET_DOWNLOAD_LOCK:
                    web_ui.DATASET_DOWNLOAD_JOBS.pop("j1", None)


class TestDatasetFileAndReport(unittest.TestCase):
    def test_dataset_file_path_prevents_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            search_dir = base / "search"
            search_dir.mkdir()
            (search_dir / "dataset_search_results.csv").write_text(
                "a",
                encoding="utf-8",
            )
            with mock.patch.object(web_ui, "DATASET_SEARCH_DIR", search_dir):
                self.assertEqual(
                    web_ui.dataset_file_path("dataset_search_results.csv"),
                    (search_dir / "dataset_search_results.csv").resolve(),
                )
                self.assertIsNone(web_ui.dataset_file_path("../secret.txt"))

    def test_single_report_path_only_serves_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            results = out / "results"
            results.mkdir(parents=True)
            (results / "result_report.html").write_text(
                "<html></html>",
                encoding="utf-8",
            )
            info = {"out": out}
            self.assertEqual(
                web_ui._single_report_path(info),
                (results / "result_report.html").resolve(),
            )
            self.assertEqual(
                web_ui._single_report_path(None, str(out)),
                (results / "result_report.html").resolve(),
            )
            info2 = {"out": Path(tmp) / "other"}
            self.assertIsNone(web_ui._single_report_path(info2))
            self.assertIsNone(web_ui._single_report_path(None, str(Path(tmp) / "other")))


if __name__ == "__main__":
    unittest.main()
