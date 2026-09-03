#!/usr/bin/env python3
"""Tests for the GROMACS MD simulation module."""

from __future__ import annotations

import csv
import logging
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

from docking.config import load_config, save_config  # noqa: E402
from docking.md_simulation import (  # noqa: E402
    _gmx_env,
    _normalize_ligand_itp,
    _write_first_pdbqt_model,
    parse_xvg,
    run_md_simulation,
)

LOG = logging.getLogger("test_md_simulation")


def _make_workdir(tmp: Path) -> Path:
    workdir = tmp / "work"
    receptor = workdir / "data" / "receptors" / "receptor.pdb"
    receptor.parent.mkdir(parents=True)
    receptor.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000\n"
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000\n"
        "TER\n"
        "END\n",
        encoding="utf-8",
    )
    out = workdir / "outputs" / "run_001"
    data_dir = out / "results" / "01_analysis" / "data"
    data_dir.mkdir(parents=True)
    docked = out / "docked"
    docked.mkdir(parents=True)
    pose = docked / "L1.pdbqt"
    pose.write_text(
        "MODEL 1\n"
        "HETATM    1  C   LIG L   1       3.000   3.000   3.000\n"
        "ENDMDL\n"
        "MODEL 2\n"
        "HETATM    2  C   LIG L   1       9.000   9.000   9.000\n"
        "ENDMDL\n",
        encoding="utf-8",
    )
    (data_dir / "fig_48_diverse_hits.csv").write_text(
        f"rank,id,affinity,smiles,pose_file\n1,L1,-8.5,CCO,{pose}\n",
        encoding="utf-8",
    )
    return workdir


def _cfg(tmp: Path):
    cfg_path = tmp / "config.json"
    save_config(
        load_config(
            APP_ROOT / "config" / "docking_config.json",
            {"workdir": str(tmp)},
        ),
        cfg_path,
    )
    return load_config(cfg_path)


class TestMdSimulation(unittest.TestCase):
    def test_prepare_mode_writes_gromacs_inputs(self):
        with tempfile.TemporaryDirectory() as raw:
            workdir = _make_workdir(Path(raw))
            cfg = _cfg(workdir)
            summary = run_md_simulation(cfg, LOG, mode="prepare")
            self.assertEqual(summary["prepared"], 1)
            md_dir = cfg.md_dir()
            self.assertTrue((md_dir / "md_simulation_summary.json").exists())
            self.assertTrue((md_dir / "L1" / "complex.pdb").exists())
            self.assertTrue((md_dir / "L1" / "em.mdp").exists())
            rows = list(
                csv.DictReader(
                    (md_dir / "md_simulation_results.csv").open(
                        "r", newline="", encoding="utf-8"
                    )
                )
            )
            self.assertEqual(rows[0]["status"], "prepared")

    def test_first_model_pdbqt_is_extracted(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            src = root / "multi.pdbqt"
            out = root / "first.pdbqt"
            src.write_text(
                "MODEL 1\n"
                "ATOM      1  C   LIG L   1       1.000   1.000   1.000\n"
                "ENDMDL\n"
                "MODEL 2\n"
                "ATOM      2  C   LIG L   1       2.000   2.000   2.000\n"
                "ENDMDL\n",
                encoding="utf-8",
            )
            _write_first_pdbqt_model(src, out)
            text = out.read_text(encoding="utf-8")
            self.assertIn("1.000   1.000   1.000", text)
            self.assertNotIn("2.000   2.000   2.000", text)

    def test_parse_xvg_skips_comments(self):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "data.xvg"
            path.write_text(
                "# comment\n@ title\n"
                "0.000 0.100\n"
                "1.000 0.200\n",
                encoding="utf-8",
            )
            data = parse_xvg(path)
            self.assertEqual(len(data), 2)
            self.assertEqual(float(data[1, 1]), 0.2)

    def test_normalize_ligand_itp_keeps_comment_line(self):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "ligand.itp"
            path.write_text(
                "[ moleculetype ]\n"
                "; name nrexcl\n"
                "MOL 1\n"
                "[ atoms ]\n"
                "1 C 1 LIG C 1 0.0 12.011\n",
                encoding="utf-8",
            )
            _normalize_ligand_itp(path)
            text = path.read_text(encoding="utf-8")
            self.assertIn("; name nrexcl", text)
            self.assertIn("\nLIG 1", text)

    def test_gmx_env_sets_data_dir(self):
        with tempfile.TemporaryDirectory() as raw:
            cfg_path = Path(raw) / "config.json"
            cfg = load_config(
                APP_ROOT / "config" / "docking_config.json",
                {"md_gmx_data": raw},
            )
            save_config(cfg, cfg_path)
            cfg = load_config(cfg_path)
            env = _gmx_env("gmx", cfg)
            self.assertEqual(env["GMXDATA"], str(Path(raw).resolve()))

    def test_auto_mode_fails_clearly_without_tools(self):
        with tempfile.TemporaryDirectory() as raw:
            workdir = _make_workdir(Path(raw))
            cfg = _cfg(workdir)
            with mock.patch("docking.md_simulation.find_tool", return_value=None):
                with self.assertRaisesRegex(Exception, "GROMACS gmx"):
                    run_md_simulation(cfg, LOG, mode="auto")
            summary_path = cfg.md_dir() / "md_simulation_summary.json"
            self.assertTrue(summary_path.exists())


class TestMdWeb(unittest.TestCase):
    def test_dock_results_includes_md_section(self):
        import web_ui

        with tempfile.TemporaryDirectory() as raw:
            workdir = _make_workdir(Path(raw))
            cfg = _cfg(workdir)
            run_md_simulation(cfg, LOG, mode="prepare")
            data = web_ui.dock_results(
                {"output_dir": cfg.output_dir, "stage": "md-simulation"}
            )
            self.assertTrue(data["md"])
            self.assertEqual(data["md"]["summary"]["prepared"], 1)
            self.assertTrue(any("complex.pdb" not in name for name in data["files"]))

    def test_dock_page_exposes_md_controls(self):
        template = (
            APP_ROOT
            / "web"
            / "templates"
            / "dock_page_template.html"
        ).read_text(encoding="utf-8")
        for token in (
            "分子动力学模拟",
            'name="md_mode"',
            'name="md_prod_steps"',
            "startMdSimulation",
        ):
            self.assertIn(token, template)


if __name__ == "__main__":
    unittest.main()
