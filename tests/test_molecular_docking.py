#!/usr/bin/env python3
"""Tests for the standalone molecular docking board."""

from __future__ import annotations

import json
import logging
import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

from molecular_docking.cli import main as molecular_cli_main  # noqa: E402
from molecular_docking.config import (  # noqa: E402
    DEFAULT_CONFIG,
    load_config,
    save_config,
)
from molecular_docking.report import generate_report  # noqa: E402
from molecular_docking.pipeline import run_pipeline  # noqa: E402
from web_ui import (  # noqa: E402
    MOLECULAR_DOCK_JOBS,
    MOLECULAR_DOCK_QUEUE,
    NAV_HTML,
    molecular_docking_results,
    render_molecular_docking_page,
    start_molecular_docking_job,
)

LOG = logging.getLogger("test_molecular_docking")

_FAKE_VINA = r"""
import argparse

p = argparse.ArgumentParser()
p.add_argument("--receptor")
p.add_argument("--ligand")
p.add_argument("--out")
p.add_argument("--center_x", type=float)
p.add_argument("--center_y", type=float)
p.add_argument("--center_z", type=float)
p.add_argument("--size_x", type=float)
p.add_argument("--size_y", type=float)
p.add_argument("--size_z", type=float)
p.add_argument("--exhaustiveness", type=int)
p.add_argument("--num_modes", type=int)
p.add_argument("--energy_range", type=float)
p.add_argument("--cpu", type=int)
p.add_argument("--seed", type=int)
p.add_argument("--scoring")
p.add_argument("--flex", action="append")
args = p.parse_args()

with open(args.out, "w", encoding="utf-8") as fh:
    fh.write("MODEL 1\nATOM      1  C   LIG L   1       1.000   1.000   1.000\nENDMDL\n")

print("mode |   affinity | dist from best mode")
print("     | (kcal/mol) | rmsd l.b.| rmsd u.b.")
print("-----+------------+----------+----------")
print("   1        -8.3      0.000      0.000")
print("   2        -7.1      1.822      3.310")
"""


class TestMolecularDockingConfig(unittest.TestCase):
    def test_default_config_resolves_independent_workdir(self):
        cfg = load_config(DEFAULT_CONFIG)
        self.assertEqual(cfg.workdir, (APP_ROOT / "molecular_docking").resolve())
        self.assertEqual(cfg.get("docking", "engine"), "vina")
        self.assertEqual(cfg.get("docking", "exhaustiveness"), 8)

    def test_save_and_reload_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            workdir = tmp_path / "work"
            cfg = load_config(
                DEFAULT_CONFIG,
                {"workdir": str(workdir), "center": [1.0, 2.0, 3.0]},
            )
            cfg_path = workdir / "config" / "molecular_docking_config.json"
            save_config(cfg, cfg_path)
            reloaded = load_config(cfg_path)
            self.assertEqual(reloaded.workdir, workdir.resolve())
            self.assertEqual(reloaded.receptor_center(), [1.0, 2.0, 3.0])


class TestMolecularDockingCli(unittest.TestCase):
    def test_init_creates_workdir_skeleton(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cfg_path = tmp_path / "config.json"
            cfg_path.write_text(
                json.dumps(
                    {
                        "name": "molecular_docking",
                        "workdir": str(tmp_path / "work"),
                        "output_dir": "outputs/run_001",
                    }
                ),
                encoding="utf-8",
            )
            code = molecular_cli_main(
                ["init", "--config", str(cfg_path)]
            )
            self.assertEqual(code, 0)
            workdir = tmp_path / "work"
            self.assertTrue((workdir / "data" / "receptors").is_dir())
            self.assertTrue((workdir / "data" / "ligands").is_dir())
            self.assertTrue(
                (workdir / "config" / "molecular_docking_config.json").is_file()
            )


class TestMolecularDockingReport(unittest.TestCase):
    def test_report_lists_results_and_poses(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            cfg = load_config(
                DEFAULT_CONFIG,
                {"workdir": str(workdir), "cutoff": -7.0},
            )
            data_dir = cfg.analysis_dir() / "data"
            figures_dir = cfg.analysis_dir() / "figures"
            data_dir.mkdir(parents=True)
            figures_dir.mkdir(parents=True)
            (data_dir / "fig_46_47_ranked_results.csv").write_text(
                "id,smiles,affinity,mode,status\n"
                "L1,CCO,-8.3,1,ok\n",
                encoding="utf-8",
            )
            (cfg.analysis_dir() / "summary.json").write_text(
                json.dumps(
                    {
                        "total_docked": 1,
                        "hits": 1,
                        "best_affinity": -8.3,
                    }
                ),
                encoding="utf-8",
            )
            (figures_dir / "fig_46_affinity_distribution.png").write_bytes(b"png")
            docked = cfg.output_dir / "docked"
            docked.mkdir(parents=True)
            (docked / "L1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")

            report_path = generate_report(cfg, LOG)
            self.assertTrue(report_path.exists())
            text = report_path.read_text(encoding="utf-8")
            self.assertIn("分子对接报告", text)
            self.assertIn("L1", text)
            self.assertIn("L1.pdbqt", text)


class TestMolecularDockingPipeline(unittest.TestCase):
    def test_pipeline_resumes_after_preparation_stages(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            workdir = tmp_path / "work"
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "receptor": "data/receptors/receptor.pdbqt",
                    "executable": str(tmp_path / "fake_vina.py"),
                    "max_workers": 1,
                    "scoring": "",
                },
            )
            receptor = workdir / "data" / "receptors" / "receptor.pdbqt"
            receptor.parent.mkdir(parents=True)
            receptor.write_text(
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000\nTER\n",
                encoding="utf-8",
            )
            lig_dir = workdir / "data" / "ligands" / "prepared"
            lig_dir.mkdir(parents=True)
            lig_pdbqt = lig_dir / "L1.pdbqt"
            lig_pdbqt.write_text(
                "ATOM      1  C   LIG L   1       0.000   0.000   0.000\n",
                encoding="utf-8",
            )
            (lig_dir / "manifest.csv").write_text(
                "id,smiles,heavy_atoms,rotatable_bonds,pdbqt,status,error\n"
                f"L1,CCO,3,0,{lig_pdbqt},ok,\n",
                encoding="utf-8",
            )
            fake_vina = tmp_path / "fake_vina.py"
            fake_vina.write_text(_FAKE_VINA, encoding="utf-8")

            stage_dir = cfg.stage_dir()
            stage_dir.mkdir(parents=True)
            for code, name in [
                ("01", "prepare-receptor"),
                ("02", "prepare-ligands"),
            ]:
                (stage_dir / f"{code}_{name}.done").write_text(
                    "done",
                    encoding="utf-8",
                )

            run_pipeline(cfg, force=False)
            report = cfg.reports_dir() / "molecular_docking_report.html"
            self.assertTrue(report.exists())
            self.assertTrue((cfg.analysis_dir() / "summary.json").exists())
            for code, name in [
                ("03", "dock"),
                ("04", "analyze"),
                ("05", "redock"),
                ("06", "report"),
            ]:
                self.assertTrue((stage_dir / f"{code}_{name}.done").exists())


class TestMolecularDockingWeb(unittest.TestCase):
    def test_page_and_global_nav_expose_independent_entry(self):
        self.assertIn("/molecular-docking", NAV_HTML)
        page = render_molecular_docking_page()
        self.assertIn("独立处理受体", page)
        self.assertIn("/molecular-docking/start", page)
        self.assertIn("detectBox()", page)

    def test_results_include_report_and_pose_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            output_dir = workdir / "outputs" / "run_001"
            reports = output_dir / "results"
            data_dir = reports / "01_analysis" / "data"
            data_dir.mkdir(parents=True)
            (data_dir / "fig_46_47_ranked_results.csv").write_text(
                "id,smiles,affinity,status\nL1,CCO,-8.3,ok\n",
                encoding="utf-8",
            )
            (reports / "molecular_docking_report.html").write_text(
                "<h1>分子对接报告</h1>",
                encoding="utf-8",
            )
            docked = output_dir / "docked"
            docked.mkdir(parents=True)
            (docked / "L1.pdbqt").write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
            info = {
                "output_dir": output_dir,
                "stage": "pipeline",
                "started": 0,
            }
            result = molecular_docking_results(info)
            self.assertEqual(result["html_report"], "molecular_docking_report.html")
            self.assertIn("docked/L1.pdbqt", result["files"])
            self.assertEqual(result["rows"][0]["id"], "L1")

    def test_start_job_registers_independent_web_task(self):
        from unittest import mock

        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            MOLECULAR_DOCK_JOBS.clear()
            MOLECULAR_DOCK_QUEUE.clear()
            with mock.patch(
                "web_ui._drain_molecular_docking_queue"
            ):
                result = start_molecular_docking_job(
                    {
                        "workdir": [str(workdir)],
                        "receptor": ["data/receptors/receptor.pdb"],
                        "stage": ["pipeline"],
                    }
                )
            job_id = result["job"]
            self.assertIn(job_id, MOLECULAR_DOCK_JOBS)
            cfg_path = workdir / "config" / f"molecular_docking_web_{job_id}.json"
            self.assertTrue(cfg_path.exists())
            MOLECULAR_DOCK_JOBS.clear()
            MOLECULAR_DOCK_QUEUE.clear()


if __name__ == "__main__":
    unittest.main()
