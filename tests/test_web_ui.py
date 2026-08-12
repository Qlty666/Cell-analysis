#!/usr/bin/env python3
"""Tests for web UI status handling."""

from __future__ import annotations

import sys
import tempfile
import time
import unittest
from unittest import mock
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

import web_ui as web_ui_module  # noqa: E402

from web_ui import (  # noqa: E402
    FULL_JOBS,
    HEARTBEAT_CLIENTS,
    HEARTBEAT_LOCK,
    _dock_file_path,
    _full_file_path,
    _full_status,
    _heartbeat_client_ids,
    _inject_heartbeat_script,
    _purge_stale_heartbeats,
    dock_results,
    full_results,
    register_heartbeat,
    start_full_job,
    unregister_heartbeat,
)


class _FakeProc:
    def __init__(self, returncode: int) -> None:
        self.returncode = returncode

    def poll(self) -> int:
        return self.returncode


class TestHeartbeatAutoShutdown(unittest.TestCase):
    def setUp(self) -> None:
        with HEARTBEAT_LOCK:
            HEARTBEAT_CLIENTS.clear()
            web_ui_module.HEARTBEAT_LAST_SEEN_AT = None

    def tearDown(self) -> None:
        with HEARTBEAT_LOCK:
            HEARTBEAT_CLIENTS.clear()
            web_ui_module.HEARTBEAT_LAST_SEEN_AT = None

    def test_heartbeat_script_injected_into_html(self):
        html = b"<html><body>test</body></html>"
        out = _inject_heartbeat_script(html)
        self.assertIn(b"navigator.sendBeacon", out)
        self.assertEqual(out.count(b"</body>"), 1)
        self.assertTrue(out.startswith(b"<html><body>test"))

    def test_heartbeat_script_not_injected_into_json(self):
        body = b'{"ok": true}'
        self.assertEqual(_inject_heartbeat_script(body), body)

    def test_register_and_unregister_clients(self):
        register_heartbeat("tab-a", now=100.0)
        register_heartbeat("tab-b", now=100.0)
        self.assertEqual(_heartbeat_client_ids(), {"tab-a", "tab-b"})
        unregister_heartbeat("tab-a")
        self.assertEqual(_heartbeat_client_ids(), {"tab-b"})

    def test_heartbeat_records_last_seen(self):
        self.assertIsNone(web_ui_module.HEARTBEAT_LAST_SEEN_AT)
        register_heartbeat("tab", now=100.0)
        self.assertEqual(web_ui_module.HEARTBEAT_LAST_SEEN_AT, 100.0)
        unregister_heartbeat("tab")
        self.assertEqual(web_ui_module.HEARTBEAT_LAST_SEEN_AT, 100.0)

    def test_stale_heartbeats_are_purged(self):
        register_heartbeat("old", now=10.0)
        register_heartbeat("new", now=100.0)
        _purge_stale_heartbeats(now=100.0, timeout=15.0)
        self.assertEqual(_heartbeat_client_ids(), {"new"})


def _make_info(tmp: Path, returncode: int) -> dict:
    workdir = tmp / "work"
    log_path = workdir / "logs" / "web_full_test.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text(
        "[INFO] === stage 01 single_cell ===\n"
        "[INFO] full pipeline complete\n",
        encoding="utf-8",
    )
    marker_dir = workdir / "outputs" / "integration" / ".stages"
    marker_dir.mkdir(parents=True)
    for code, name in [
        ("01", "single_cell"),
        ("02", "key_targets"),
        ("03", "evidence"),
        ("04", "knockout_inputs"),
        ("05", "knockout"),
        ("06", "docking"),
        ("07", "report"),
    ]:
        (marker_dir / f"{code}_{name}.done").write_text(
            "done",
            encoding="utf-8",
        )
    return {
        "job_id": "test-job",
        "proc": _FakeProc(returncode),
        "workdir": workdir,
        "log": log_path,
        "accession": "GSE999999",
        "notified": True,
        "paused": False,
        "started": time.time(),
    }


class TestFullStatus(unittest.TestCase):
    def test_start_full_job_requires_output(self):
        with self.assertRaises(ValueError):
            start_full_job({"workdir": "somewhere", "output": ""})

    def test_start_full_job_requires_workdir(self):
        with self.assertRaises(ValueError):
            start_full_job({"workdir": "", "output": "somewhere"})

    def test_start_full_job_passes_workdir_to_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output = base / "out"
            workdir = base / "work"
            with mock.patch("web_ui._drain_full_queue"):
                result = start_full_job(
                    {
                        "output": [str(output)],
                        "workdir": [str(workdir)],
                    }
                )
            job_id = result["job"]
            try:
                cmd = FULL_JOBS[job_id]["cmd"]
                self.assertIn("--workdir", cmd)
                self.assertEqual(
                    cmd[cmd.index("--workdir") + 1],
                    str(workdir.resolve()),
                )
            finally:
                FULL_JOBS.pop(job_id, None)

    def test_completed_full_pipeline_uses_full_flow_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            status = _full_status(_make_info(Path(tmp), 0))
            self.assertTrue(status["ok"])
            self.assertFalse(status["paused"])
            self.assertEqual(status["stage"], "\u5168\u6d41\u7a0b\u5b8c\u6210")

    def test_pause_exit_code_is_reported_as_paused(self):
        with tempfile.TemporaryDirectory() as tmp:
            status = _full_status(_make_info(Path(tmp), 98))
            self.assertFalse(status["ok"])
            self.assertTrue(status["paused"])


class TestResultDirectoryQueries(unittest.TestCase):
    def _make_full_workdir(self, tmp: Path) -> Path:
        workdir = tmp / "work"
        integration = workdir / "outputs" / "integration"
        integration.mkdir(parents=True)
        (integration / "key_genes.csv").write_text(
            "rank,gene\n1,EGFR\n", encoding="utf-8"
        )
        (integration / "gene_evidence.csv").write_text(
            "gene,uniprot\nEGFR,P00533\n", encoding="utf-8"
        )
        (integration / "docking_targets.csv").write_text(
            "gene,status\nEGFR,ok\n", encoding="utf-8"
        )
        (integration / "integration_summary.json").write_text(
            "{}", encoding="utf-8"
        )
        (integration / "integration_report.html").write_text(
            "<html></html>", encoding="utf-8"
        )

        ko_data = (
            workdir
            / "outputs"
            / "run_001"
            / "results"
            / "04_knockout"
            / "data"
        )
        ko_data.mkdir(parents=True)
        (ko_data / "fig_52_53_ranked_knockout.csv").write_text(
            "rank,gene\n1,EGFR\n", encoding="utf-8"
        )
        (ko_data / "fig_52_target_candidates.csv").write_text(
            "rank,gene\n1,EGFR\n", encoding="utf-8"
        )
        (ko_data / "knockout_report.md").write_text("# report", encoding="utf-8")

        val_data = (
            workdir
            / "outputs"
            / "run_001"
            / "results"
            / "05_validation"
            / "data"
        )
        val_data.mkdir(parents=True)
        (val_data / "validation_candidates.csv").write_text(
            "rank,gene\n1,EGFR\n", encoding="utf-8"
        )
        (val_data / "validation_plan.md").write_text(
            "# plan", encoding="utf-8"
        )

        gene_run = workdir / "work" / "EGFR" / "outputs" / "run_001"
        docked = gene_run / "docked"
        docked.mkdir(parents=True)
        (docked / "results.csv").write_text(
            "id,affinity\nL1,-8.0\n", encoding="utf-8"
        )
        figures = gene_run / "results" / "01_analysis" / "figures"
        figures.mkdir(parents=True)
        (figures / "fig_46_affinity_distribution.png").write_bytes(b"png")
        analysis_data = gene_run / "results" / "01_analysis" / "data"
        analysis_data.mkdir(parents=True)
        (analysis_data / "fig_46_47_ranked_results.csv").write_text(
            "rank,id,affinity\n1,L1,-8.0\n", encoding="utf-8"
        )
        (gene_run / "results" / "01_analysis" / "summary.json").write_text(
            "{}", encoding="utf-8"
        )
        return workdir

    def test_full_results_only_returns_result_figures_and_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = self._make_full_workdir(Path(tmp))
            data = full_results(workdir)
            files = data["files"]
            self.assertIn("key_genes.csv", files)
            self.assertIn(
                "outputs/run_001/results/04_knockout/data/"
                "fig_52_53_ranked_knockout.csv",
                files,
            )
            self.assertIn(
                "outputs/run_001/results/05_validation/data/"
                "validation_candidates.csv",
                files,
            )
            self.assertIn(
                "work/EGFR/outputs/run_001/docked/results.csv",
                files,
            )
            self.assertIn(
                "work/EGFR/outputs/run_001/results/01_analysis/figures/"
                "fig_46_affinity_distribution.png",
                files,
            )
            for excluded in [
                "integration_summary.json",
                "integration_report.html",
                "validation_plan.md",
                "summary.json",
                "knockout_report.md",
            ]:
                self.assertFalse(
                    any(name.endswith(excluded) for name in files),
                    excluded,
                )

    def test_full_file_path_allows_result_files_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = self._make_full_workdir(Path(tmp))
            result = _full_file_path(
                workdir,
                "outputs/run_001/results/04_knockout/data/"
                "fig_52_53_ranked_knockout.csv",
            )
            self.assertIsNotNone(result)
            self.assertTrue(result.is_file())
            self.assertIsNotNone(
                _full_file_path(
                    workdir,
                    "work/EGFR/outputs/run_001/docked/results.csv",
                )
            )
            self.assertIsNone(
                _full_file_path(
                    workdir,
                    "outputs/run_001/results/05_validation/data/"
                    "validation_plan.md",
                )
            )
            self.assertIsNone(
                _full_file_path(
                    workdir,
                    "outputs/integration/integration_summary.json",
                )
            )

    def test_full_results_includes_single_cell_result_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            workdir = self._make_full_workdir(base)
            single_cell = base / "single_cell"
            fig_dir = single_cell / "results" / "figures" / "01_qc"
            data_dir = single_cell / "results" / "data" / "01_qc"
            fig_dir.mkdir(parents=True)
            data_dir.mkdir(parents=True)
            (fig_dir / "fig_01_qc_raw_violin.png").write_bytes(b"png")
            (data_dir / "fig_01_qc_metrics.csv").write_text(
                "cell,nFeature_RNA\nC1,1000\n", encoding="utf-8"
            )
            (single_cell / "results" / "summary.json").write_text(
                "{}", encoding="utf-8"
            )
            context_dir = workdir / "outputs" / "integration" / ".stages"
            context_dir.mkdir(parents=True, exist_ok=True)
            (context_dir / "run_context.json").write_text(
                '{"single_cell_root": "%s"}'
                % str(single_cell).replace("\\", "\\\\"),
                encoding="utf-8",
            )

            data = full_results(workdir)
            files = data["files"]
            self.assertIn(
                "single_cell/results/figures/01_qc/fig_01_qc_raw_violin.png",
                files,
            )
            self.assertIn(
                "single_cell/results/data/01_qc/fig_01_qc_metrics.csv",
                files,
            )
            self.assertFalse(
                any(name.endswith("summary.json") for name in files)
            )
            self.assertIsNotNone(
                _full_file_path(
                    workdir,
                    "single_cell/results/figures/01_qc/fig_01_qc_raw_violin.png",
                )
            )
            self.assertIsNone(
                _full_file_path(
                    workdir,
                    "single_cell/results/summary.json",
                )
            )

    def test_dock_results_only_returns_result_figures_and_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            reports = out / "results"
            reports.mkdir(parents=True)
            analysis = reports / "01_analysis"
            (analysis / "data").mkdir(parents=True)
            (analysis / "figures").mkdir(parents=True)
            (analysis / "data" / "fig_46_47_ranked_results.csv").write_text(
                "rank,id,affinity\n1,L1,-8.0\n", encoding="utf-8"
            )
            (analysis / "figures" / "fig_46_affinity_distribution.png").write_bytes(
                b"png"
            )
            (analysis / "summary.json").write_text("{}", encoding="utf-8")
            (reports / "docking_report.html").write_text(
                "<html></html>", encoding="utf-8"
            )
            ml = reports / "03_ml" / "data"
            ml.mkdir(parents=True)
            (ml / "ml_model_info.json").write_text("{}", encoding="utf-8")
            (ml / "ml_model.joblib").write_bytes(b"model")
            docked = out / "docked"
            docked.mkdir(parents=True)
            (docked / "results.csv").write_text(
                "id,affinity\nL1,-8.0\n", encoding="utf-8"
            )

            info = {"output_dir": out}
            data = dock_results(info)
            files = data["files"]
            self.assertIn(
                "01_analysis/data/fig_46_47_ranked_results.csv",
                files,
            )
            self.assertIn(
                "01_analysis/figures/fig_46_affinity_distribution.png",
                files,
            )
            self.assertIn("docked/results.csv", files)
            for excluded in [
                "summary.json",
                "docking_report.html",
                "ml_model_info.json",
                "ml_model.joblib",
            ]:
                self.assertFalse(
                    any(name.endswith(excluded) for name in files),
                    excluded,
                )
            self.assertIsNotNone(
                _dock_file_path(
                    info,
                    "01_analysis/data/fig_46_47_ranked_results.csv",
                )
            )
            self.assertIsNotNone(
                _dock_file_path(info, "docked/results.csv")
            )
            self.assertIsNone(
                _dock_file_path(info, "01_analysis/summary.json")
            )
            self.assertIsNone(
                _dock_file_path(info, "03_ml/data/ml_model.joblib")
            )


if __name__ == "__main__":
    unittest.main()
