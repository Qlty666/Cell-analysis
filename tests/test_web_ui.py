#!/usr/bin/env python3
"""Tests for web UI status handling."""

from __future__ import annotations

import socket
import json
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
    DOCK_JOBS,
    FINISHED_NOTIFICATIONS,
    FULL_JOBS,
    HEARTBEAT_CLIENTS,
    HEARTBEAT_LOCK,
    NOTIFY_LOCK,
    _analysis_file_path,
    _dock_file_path,
    _full_file_path,
    _full_status,
    _single_status,
    _cleanup_stale_web_ui,
    _heartbeat_client_ids,
    _inject_heartbeat_script,
    _port_is_listening,
    _purge_stale_heartbeats,
    dock_results,
    full_results,
    register_heartbeat,
    run_faers_request,
    run_network_request,
    start_dock_job,
    start_full_job,
    unregister_heartbeat,
)


class _FakeProc:
    def __init__(self, returncode: int) -> None:
        self.returncode = returncode

    def poll(self) -> int:
        return self.returncode


class TestPortCleanup(unittest.TestCase):
    def test_port_listening_detection(self):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind(("127.0.0.1", 0))
            sock.listen()
            port = sock.getsockname()[1]
            self.assertTrue(_port_is_listening("127.0.0.1", port))
        self.assertFalse(_port_is_listening("127.0.0.1", 0))

    def test_cleanup_skips_when_port_free(self):
        with (
            mock.patch("web_ui._port_is_listening", return_value=False),
            mock.patch("web_ui._stop_stale_web_ui") as stop,
        ):
            self.assertTrue(_cleanup_stale_web_ui("127.0.0.1", 8000))
            stop.assert_not_called()

    def test_cleanup_stops_when_port_busy(self):
        with (
            mock.patch("web_ui._port_is_listening", return_value=True),
            mock.patch("web_ui._stop_stale_web_ui", return_value=True) as stop,
        ):
            self.assertTrue(_cleanup_stale_web_ui("127.0.0.1", 8000))
            stop.assert_called_once_with("127.0.0.1", 8000)


class TestSingleCellFigureList(unittest.TestCase):
    def test_top5_enrichment_figures_wired_into_web_ui(self):
        for name in ("fig_46_go_top5.png", "fig_47_kegg_top5.png"):
            self.assertIn(name, web_ui_module.FIGURE_NAMES)
        page = web_ui_module.render_page()
        self.assertIn("fig_46_go_top5.png", page)
        self.assertIn("GO BP 筛选后 Top5", page)
        self.assertIn("fig_47_kegg_top5.png", page)
        self.assertIn("KEGG 筛选后 Top5", page)

    def test_top5_enrichment_figures_expose_expected_styles(self):
        by_name = {item["file"]: item for item in web_ui_module.FIGURES}
        for name in ("fig_46_go_top5.png", "fig_47_kegg_top5.png"):
            self.assertIn("dotplot", by_name[name]["styles"])
            self.assertIn("barplot", by_name[name]["styles"])


class TestJobRecordPersistence(unittest.TestCase):
    def test_web_pages_do_not_clear_job_records_on_reload(self):
        web_ui_source = (APP_ROOT / "web" / "web_ui.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("clearJobRecordsOnReload", web_ui_source)
        self.assertNotIn("nav.type === 'reload'", web_ui_source)

        for name in (
            "full_page_template.html",
            "web_page_template.html",
            "dock_page_template.html",
            "tasks_template.html",
            "results_manifest_optimized.html",
        ):
            template = (APP_ROOT / "web" / "templates" / name).read_text(
                encoding="utf-8"
            )
            self.assertNotIn("clearJobRecordsOnReload", template)
            self.assertNotIn("nav.type === 'reload'", template)


class TestNetworkAndFaersWeb(unittest.TestCase):
    def test_network_request_runs_and_lists_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = workdir / "data" / "network"
            data_dir.mkdir(parents=True)
            (data_dir / "targets.csv").write_text(
                "gene,source\nALB,CTD\nGPC3,ChEMBL\nMMP9,CTD\n",
                encoding="utf-8",
            )
            (data_dir / "disease.csv").write_text(
                "gene\nALB\nGPC3\nEGFR\n",
                encoding="utf-8",
            )
            (data_dir / "ppi.tsv").write_text(
                "protein1\tprotein2\nALB\tGPC3\n",
                encoding="utf-8",
            )
            result = run_network_request(
                {
                    "net_workdir": [str(workdir)],
                    "net_compound_name": ["Test"],
                    "net_disease_name": ["Liver"],
                    "net_compound_targets": ["data/network/targets.csv"],
                    "net_disease_genes": ["data/network/disease.csv"],
                    "net_ppi": ["data/network/ppi.tsv"],
                }
            )
            self.assertEqual(result["summary"]["overlap_genes"], 2)
            self.assertTrue(result["summary"]["ppi_hub_scored"])
            self.assertTrue(
                any("compound_disease_overlap.csv" in f for f in result["files"])
            )
            out = Path(result["output_dir"])
            self.assertTrue((out / "figures" / "compound_disease_venn.png").exists())

    def test_faers_request_runs_and_lists_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = workdir / "data" / "faers"
            data_dir.mkdir(parents=True)
            (data_dir / "events.csv").write_text(
                "drug,event,count\nA,X,10\nA,Y,2\nB,X,2\nB,Y,20\n",
                encoding="utf-8",
            )
            result = run_faers_request(
                {
                    "faers_workdir": [str(workdir)],
                    "faers_input": ["data/faers/events.csv"],
                    "faers_min_count": ["3"],
                }
            )
            self.assertGreaterEqual(result["summary"]["signals"], 1)
            self.assertTrue(
                any("faers_signals.csv" in f for f in result["files"])
            )

    def test_analysis_file_path_rejects_escape(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            out_dir = workdir / "outputs" / "run_001" / "network_toxicology" / "data"
            out_dir.mkdir(parents=True)
            target = out_dir / "compound_disease_overlap.csv"
            target.write_text("gene,n_sources,sources\nALB,1,CTD\n", encoding="utf-8")
            self.assertEqual(
                _analysis_file_path(
                    str(workdir),
                    "data/compound_disease_overlap.csv",
                    "network",
                ),
                target.resolve(),
            )
            self.assertIsNone(
                _analysis_file_path(str(workdir), "../escape.csv", "network")
            )

    def test_web_templates_expose_new_sections(self):
        dock = (APP_ROOT / "web" / "templates" / "dock_page_template.html").read_text(
            encoding="utf-8"
        )
        results = (
            APP_ROOT / "web" / "templates" / "results_manifest_optimized.html"
        ).read_text(encoding="utf-8")
        guide = (APP_ROOT / "docs" / "result_figure_guide.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("网络毒理学分析", dock)
        self.assertIn("FAERS 不相称性信号检测", dock)
        self.assertIn("网络毒理学与 FAERS 信号", results)
        self.assertIn("### 6.3 网络毒理学与 FAERS 信号", guide)


class TestRecentWebIntegration(unittest.TestCase):
    def test_dock_job_passes_ml_model_options(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            with mock.patch("web_ui._drain_dock_queue"):
                result = start_dock_job(
                    {
                        "workdir": [str(workdir)],
                        "stage": ["ml-train"],
                        "model": ["lasso_svm"],
                        "training_csv": ["data/ml/training.csv"],
                        "label_column": ["active"],
                    }
                )
            job_id = result["job"]
            try:
                info = DOCK_JOBS[job_id]
                cfg_path = workdir / "config" / f"docking_web_{job_id}.json"
                cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
                self.assertEqual(cfg["ml"]["model"], "lasso_svm")
                self.assertEqual(
                    cfg["ml"]["training_csv"],
                    "data/ml/training.csv",
                )
                self.assertEqual(cfg["ml"]["label_column"], "active")
                self.assertEqual(info["stage"], "ml-train")
            finally:
                DOCK_JOBS.pop(job_id, None)

    def test_templates_expose_recent_controls(self):
        single = (APP_ROOT / "web" / "templates" / "web_page_template.html").read_text(
            encoding="utf-8"
        )
        full = (APP_ROOT / "web" / "templates" / "full_page_template.html").read_text(
            encoding="utf-8"
        )
        dock = (APP_ROOT / "web" / "templates" / "dock_page_template.html").read_text(
            encoding="utf-8"
        )
        self.assertIn('name="LIVER_ML_MODEL"', single)
        self.assertIn('name="ml_model"', full)
        self.assertIn('name="ppi_network_csv"', full)
        self.assertIn('name="depmap_csv"', full)
        self.assertIn('name="model"', dock)
        self.assertIn('name="training_csv"', dock)
        self.assertIn('name="ko_ppi"', dock)

    def test_figures_and_manifest_include_calibration(self):
        self.assertIn(
            "fig_45_ml_calibration_curve.png",
            web_ui_module.FIGURE_NAMES,
        )
        page = web_ui_module.render_page()
        self.assertIn("fig_45_ml_calibration_curve.png", page)
        results = (
            APP_ROOT / "web" / "templates" / "results_manifest_optimized.html"
        ).read_text(encoding="utf-8")
        self.assertIn("fig_45_ml_calibration_curve.png", results)
        self.assertIn("fig_24_ml_selected_features.csv", results)


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
        ("07", "cell_feedback"),
        ("08", "report"),
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
        "output": str(tmp / "single_cell"),
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

    def test_start_full_job_passes_new_optimization_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            output = base / "out"
            workdir = base / "work"
            with mock.patch("web_ui._drain_full_queue"):
                result = start_full_job(
                    {
                        "output": [str(output)],
                        "workdir": [str(workdir)],
                        "skip_qc_gate": ["1"],
                        "skip_differential_abundance": ["1"],
                        "dry_run": ["1"],
                        "ml_model": ["lasso_svm"],
                        "ppi_network_csv": ["data/network/string_edges.tsv"],
                        "depmap_csv": ["data/knockout/depmap.csv"],
                    }
                )
            job_id = result["job"]
            try:
                cmd = FULL_JOBS[job_id]["cmd"]
                self.assertIn("--skip-qc-gate", cmd)
                self.assertIn("--skip-differential-abundance", cmd)
                self.assertIn("--dry-run", cmd)
                self.assertIn("--ml-model", cmd)
                self.assertEqual(
                    cmd[cmd.index("--ml-model") + 1],
                    "lasso_svm",
                )
                self.assertIn("--ppi-network-csv", cmd)
                self.assertIn("--depmap-csv", cmd)
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

    def test_paused_notification_uses_readable_chinese_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            info = _make_info(Path(tmp), 98)
            info["notified"] = False
            with NOTIFY_LOCK:
                FINISHED_NOTIFICATIONS.clear()
            with mock.patch("web_ui._append_task_history"):
                status = _full_status(info)
            self.assertTrue(status["paused"])
            with NOTIFY_LOCK:
                items = list(FINISHED_NOTIFICATIONS)
                FINISHED_NOTIFICATIONS.clear()
            self.assertTrue(items)
            item = items[-1]
            self.assertEqual(item["page_label"], "全自动流水线")
            self.assertNotIn("鍏", item["title"])

    def test_paused_full_job_reports_single_cell_stage_from_r_log(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            info = _make_info(base, 98)
            sc = base / "single_cell"
            log_dir = sc / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            (log_dir / "pipeline_r.log").write_text(
                "[2026-08-17 16:52:47] start stage: 07_enrichment\n",
                encoding="utf-8",
            )
            marker_dir = (
                info["workdir"] / "outputs" / "integration" / ".stages"
            )
            for marker in marker_dir.glob("*.done"):
                marker.unlink()
            with (
                mock.patch("web_ui._append_task_history"),
                mock.patch("web_ui._drain_full_queue"),
            ):
                status = _full_status(info)
            self.assertTrue(status["paused"])
            self.assertEqual(status["stage"], "单细胞分析")

    def test_failed_full_job_reports_in_progress_stage_not_last_done(self):
        with tempfile.TemporaryDirectory() as tmp:
            info = _make_info(Path(tmp), 1)
            marker_dir = (
                info["workdir"] / "outputs" / "integration" / ".stages"
            )
            for code in ("07", "08"):
                for marker in marker_dir.glob(f"{code}_*.done"):
                    marker.unlink()
            info["log"].write_text(
                "[INFO] === stage 06 docking ===\n"
                "[INFO] === stage 07 cell_feedback ===\n"
                "RuntimeError: cell feedback R analysis failed\n"
                "Error in order(...) : argument 1 is not a vector\n",
                encoding="utf-8",
            )
            with mock.patch("web_ui._drain_full_queue"):
                status = _full_status(info)
            self.assertFalse(status["ok"])
            self.assertEqual(status["stage"], "细胞反馈")


class TestSingleCellStatus(unittest.TestCase):
    def test_paused_single_job_is_not_recorded_as_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            log_path = out / "logs" / "web_test.log"
            log_path.parent.mkdir(parents=True)
            log_path.write_text("", encoding="utf-8")
            info = {
                "job_id": "single-test",
                "proc": _FakeProc(98),
                "out": out,
                "log": log_path,
                "accession": "GSE123456",
                "species": "hs",
                "started": time.time(),
                "recorded": False,
                "notified": True,
            }
            with mock.patch("web_ui.record_job") as record:
                status = _single_status(info)
            self.assertTrue(status["paused"])
            self.assertFalse(status["ok"])
            record.assert_not_called()


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
        (integration / "qc_metrics.json").write_text(
            '{"qc_gate": {"status": "fail", "checks": []}}',
            encoding="utf-8",
        )
        (integration / "differential_abundance.csv").write_text(
            "celltype,n_cells,chi2,p_value,p_adjust,significant,direction\n"
            "T_cell,50,9.0,0.002,0.002,true,enriched_in_Tumor\n",
            encoding="utf-8",
        )
        feedback_data = integration / "cell_feedback" / "data"
        feedback_data.mkdir(parents=True)
        (feedback_data / "feedback_targets.csv").write_text(
            "gene,source,feedback_score,cell_support_score,top_celltype\n"
            "EGFR,knockout|docking,0.9,0.92,Hepatocyte\n",
            encoding="utf-8",
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
            self.assertIn("cell_feedback/data/feedback_targets.csv", files)
            self.assertEqual(data["cell_feedback"][0]["gene"], "EGFR")
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
            self.assertIsNotNone(
                _full_file_path(
                    workdir,
                    "cell_feedback/data/feedback_targets.csv",
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

    def test_full_results_includes_qc_and_differential_abundance(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = self._make_full_workdir(Path(tmp))
            data = full_results(workdir)
            self.assertEqual(
                data["qc_metrics"]["qc_gate"]["status"],
                "fail",
            )
            self.assertEqual(
                data["differential_abundance"][0]["celltype"],
                "T_cell",
            )
            self.assertIn(
                "differential_abundance.csv",
                data["files"],
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


class TestTemplatePolish(unittest.TestCase):
    def _read(self, name: str) -> str:
        path = APP_ROOT / "web" / "templates" / name
        return path.read_text(encoding="utf-8", errors="replace")

    def test_shared_css_exists(self):
        css = (APP_ROOT / "web" / "static" / "app.css").read_text(
            encoding="utf-8"
        )
        self.assertIn(".form-section", css)
        self.assertIn(".stat-card", css)

    def test_full_page_has_layout_and_form_helpers(self):
        html = self._read("full_page_template.html")
        for token in [
            "page-head",
            "form-section",
            "saveFormState",
            "restoreFormState",
            "resultStats",
            "skip_differential_abundance",
            "prefillDataType",
            "bulk",
        ]:
            self.assertIn(token, html)

    def test_single_cell_page_has_collapsible_sections(self):
        html = self._read("web_page_template.html")
        for token in ["form-section", "SINGLE_FORM_KEY", "saveFormState"]:
            self.assertIn(token, html)

    def test_dataset_page_groups_search_options(self):
        html = self._read("datasets_template.html")
        self.assertIn("form-section", html)
        self.assertIn("过滤与下载选项", html)
        self.assertIn("filterDataType", html)
        self.assertIn("applyFilters", html)
        self.assertIn("filterMinSamples", html)
        self.assertIn("filterMaxSamples", html)
        self.assertIn("filterStartDate", html)
        self.assertIn("filterEndDate", html)
        self.assertIn("filterPlatform", html)
        self.assertIn("filterType", html)
        self.assertIn("仅支持 single-cell", html)

    def test_tasks_page_has_stat_cards(self):
        html = self._read("tasks_template.html")
        for token in ["statTotal", "statRunning", "statQueued", "statPaused"]:
            self.assertIn(token, html)

    def test_results_page_has_filter_toolbar(self):
        html = self._read("results_manifest_optimized.html")
        self.assertIn("resultFilter", html)
        self.assertIn("filterResults", html)

    def test_results_manifest_includes_latest_figures(self):
        html = self._read("results_manifest_optimized.html")
        for token in [
            "figures/01_qc/fig_48_qc_pvalue_comparison.png",
            "data/01_qc/fig_48_qc_pvalue_comparison.csv",
            "figures/06_enrichment/fig_46_go_top5.png",
            "figures/06_enrichment/fig_47_kegg_top5.png",
        ]:
            self.assertIn(token, html)


if __name__ == "__main__":
    unittest.main()
