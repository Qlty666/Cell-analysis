#!/usr/bin/env python3
"""Tests for the GROMACS MD simulation module."""

from __future__ import annotations

import csv
import logging
import os
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from subprocess import CompletedProcess

import numpy as np

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

from docking.config import load_config, save_config  # noqa: E402
from docking.md_simulation import (  # noqa: E402
    _binding_contact_residues,
    _combine_hbond_series,
    _empty_md_metrics,
    _mean_rmsf_for_contact_residues,
    _stability_label,
    _store_mean_tail,
    _gmx_env,
    _normalize_ligand_itp,
    _split_external_command,
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
            self.assertEqual(rows[0]["rg_protein_mean_nm"], "")

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

    def test_parse_xvg_missing_file_returns_none(self):
        with tempfile.TemporaryDirectory() as raw:
            self.assertIsNone(parse_xvg(Path(raw) / "missing.xvg"))

    def test_analyze_collects_extended_gromacs_metrics(self):
        with tempfile.TemporaryDirectory() as raw:
            workdir = _make_workdir(Path(raw))
            cfg = _cfg(workdir)
            cfg.data["md_simulation"]["figures"] = False
            run_dir = cfg.md_dir() / "L1"
            run_dir.mkdir(parents=True)
            (run_dir / "md.tpr").write_bytes(b"tpr")
            (run_dir / "md_nojump.xtc").write_bytes(b"xtc")
            (run_dir / "index.ndx").write_text(
                "[ Protein ]\n1\n[ LIG ]\n2\n",
                encoding="utf-8",
            )
            (run_dir / "protein.gro").write_text(
                "protein\n1\n"
                "    1PRO      N    1   0.000   0.000   0.000\n"
                "  1.0   1.0   1.0\n",
                encoding="utf-8",
            )
            (run_dir / "ligand.gro").write_text(
                "ligand\n1\n"
                "    1LIG      C    1   0.100   0.100   0.100\n"
                "  1.0   1.0   1.0\n",
                encoding="utf-8",
            )

            def fake_run_command(cmd, **_kwargs):
                out = None
                for index, part in enumerate(cmd):
                    if part in ("-o", "-num") and index + 1 < len(cmd):
                        out = Path(cmd[index + 1])
                if out is None:
                    raise AssertionError("missing output flag in command")
                if "-res" in cmd:
                    out.write_text("1 0.1\n2 0.1\n", encoding="utf-8")
                else:
                    out.write_text(
                        "0.0 1.0\n10.0 1.0\n",
                        encoding="utf-8",
                    )
                return CompletedProcess(cmd, 0)

            with mock.patch(
                "docking.md_simulation.find_tool",
                return_value="gmx",
            ), mock.patch(
                "docking.md_simulation.run_command",
                side_effect=fake_run_command,
            ):
                from docking.md_simulation import _analyze_gromacs_output

                metrics = _analyze_gromacs_output(cfg, run_dir)
            self.assertAlmostEqual(float(metrics["rmsd_protein_mean_nm"]), 1.0)
            self.assertAlmostEqual(float(metrics["rg_protein_mean_nm"]), 1.0)
            self.assertAlmostEqual(float(metrics["sasa_protein_mean_nm2"]), 1.0)
            self.assertAlmostEqual(
                float(metrics["hbonds_protein_ligand_mean"]), 2.0
            )
            self.assertAlmostEqual(
                float(metrics["rmsf_contact_residue_mean_nm"]), 0.1
            )
            self.assertEqual(metrics["stability_label"], "stable")

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

    def test_store_mean_tail_records_last_half_stats(self):
        with tempfile.TemporaryDirectory() as raw:
            path = Path(raw) / "series.xvg"
            path.write_text(
                "0.0 0.1\n"
                "1.0 0.2\n"
                "2.0 0.3\n"
                "3.0 0.4\n"
                "4.0 0.5\n",
                encoding="utf-8",
            )
            metrics = _empty_md_metrics()
            _store_mean_tail(metrics, path, "demo", "_nm", fraction=0.4)
            self.assertAlmostEqual(float(metrics["demo_mean_nm"]), 0.3)
            self.assertAlmostEqual(float(metrics["demo_tail_mean_nm"]), 0.45)
            self.assertGreater(float(metrics["demo_tail_std_nm"]), 0.0)

    def test_binding_contact_residues_uses_gro_coordinates(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            protein = root / "protein.gro"
            protein.write_text(
                "protein\n2\n"
                "    1PRO      N    1   0.000   0.000   0.000\n"
                "    2PRO      N    2  10.000  10.000  10.000\n"
                "  2.0   2.0   2.0\n",
                encoding="utf-8",
            )
            ligand = root / "ligand.gro"
            ligand.write_text(
                "ligand\n1\n"
                "    1LIG      C    1   0.100   0.100   0.100\n"
                "  1.0   1.0   1.0\n",
                encoding="utf-8",
            )
            self.assertEqual(
                _binding_contact_residues(protein, ligand, cutoff_nm=0.6),
                [1],
            )

    def test_mean_rmsf_for_contact_residues_filters_rows(self):
        data = np.asarray(
            [[1.0, 0.10], [2.0, 0.50], [3.0, 0.90]],
            dtype=float,
        )
        self.assertAlmostEqual(
            float(_mean_rmsf_for_contact_residues(data, [1, 3])),
            0.5,
        )
        self.assertIsNone(_mean_rmsf_for_contact_residues(data, [9]))

    def test_stability_label_uses_last_half_noise(self):
        stable = _empty_md_metrics()
        stable.update(
            {
                "time_ns": 20.0,
                "rmsd_protein_tail_std_nm": 0.05,
                "rmsd_ligand_tail_std_nm": 0.08,
                "rg_protein_tail_std_nm": 0.02,
                "sasa_protein_tail_std_nm2": 0.4,
            }
        )
        self.assertEqual(_stability_label(stable), "stable")
        stable["rg_protein_tail_std_nm"] = 0.2
        self.assertEqual(_stability_label(stable), "review")
        self.assertEqual(_stability_label(_empty_md_metrics()), "insufficient")

    def test_combine_hbond_series_sums_both_donor_directions(self):
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            (root / "a.xvg").write_text(
                "0.0 2.0\n1.0 3.0\n",
                encoding="utf-8",
            )
            (root / "b.xvg").write_text(
                "0.0 1.0\n1.0 1.0\n",
                encoding="utf-8",
            )
            out = _combine_hbond_series(root, ["a.xvg", "b.xvg"], "combined.xvg")
            self.assertIsNotNone(out)
            data = parse_xvg(out)
            self.assertEqual(float(data[0, 1]), 3.0)
            self.assertEqual(float(data[1, 1]), 4.0)

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


class TestExternalCommandParsing(unittest.TestCase):
    @unittest.skipUnless(os.name == "nt", "Windows path parsing")
    def test_windows_path_is_not_mangled(self):
        command = _split_external_command(
            "\"C:\\Program Files\\gmx_MMPBSA.exe\" {run_dir}"
        )
        self.assertEqual(
            command,
            [
                r"C:\Program Files\gmx_MMPBSA.exe",
                "{run_dir}",
            ],
        )

    def test_argument_list_is_preserved(self):
        command = _split_external_command(["gmx_MMPBSA", "--version"])
        self.assertEqual(command, ["gmx_MMPBSA", "--version"])


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

    def test_md_page_exposes_md_controls(self):
        template = (
            APP_ROOT
            / "web"
            / "templates"
            / "md_simulation_page_template.html"
        ).read_text(encoding="utf-8")
        for token in (
            "分子动力学模拟",
            'name="md_mode"',
            'name="md_prod_steps"',
            "startMdSimulation",
        ):
            self.assertIn(token, template)
        dock = (
            APP_ROOT / "web" / "templates" / "dock_page_template.html"
        ).read_text(encoding="utf-8")
        self.assertNotIn('id="mdForm"', dock)


if __name__ == "__main__":
    unittest.main()
