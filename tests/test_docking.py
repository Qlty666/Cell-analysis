#!/usr/bin/env python3
"""Unit tests for the docking pipeline (no RDKit/Meeko required)."""

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

from docking.analysis import analyze_results  # noqa: E402
from docking.config import load_config, save_config  # noqa: E402
from docking.docking import build_vina_command, run_docking  # noqa: E402
from docking import ml as docking_ml  # noqa: E402
from docking.receptor import _sanitize_receptor_models  # noqa: E402
from docking.utils import parse_vina_affinities, safe_name  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_docking")


class TestConfig(unittest.TestCase):
    def test_default_workdir(self):
        cfg = load_config(DEFAULT_CONFIG)
        self.assertEqual(cfg.workdir, (APP_ROOT / "dock").resolve())
        self.assertEqual(len(cfg.receptor_center()), 3)

    def test_overrides(self):
        cfg = load_config(
            DEFAULT_CONFIG,
            {"center": [1.0, 2.0, 3.0], "exhaustiveness": 16},
        )
        self.assertEqual(cfg.receptor_center(), [1.0, 2.0, 3.0])
        self.assertEqual(cfg.get("docking", "exhaustiveness"), 16)


class TestUtils(unittest.TestCase):
    def test_parse_vina_affinities(self):
        text = (
            "mode |   affinity | dist from best mode\n"
            "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
            "-----+------------+----------+----------\n"
            "   1        -8.3      0.000      0.000\n"
            "   2        -7.1      1.822      3.310\n"
        )
        modes = parse_vina_affinities(text)
        self.assertEqual(modes[0]["affinity"], -8.3)
        self.assertEqual(modes[1]["mode"], 2)

    def test_safe_name(self):
        self.assertEqual(safe_name("ABC_123"), "ABC_123")
        self.assertEqual(safe_name("x/y z"), "x_y_z")
        self.assertEqual(safe_name(""), "ligand")


class TestCommandBuild(unittest.TestCase):
    def test_build_vina_command(self):
        cfg = load_config(
            DEFAULT_CONFIG,
            {
                "center": [10.0, 20.0, 30.0],
                "size": [25.0, 25.0, 25.0],
                "executable": "fake_vina.py",
            },
        )
        cmd = build_vina_command(
            cfg,
            "fake_vina.py",
            "receptor.pdbqt",
            "ligand.pdbqt",
            "out.pdbqt",
        )
        self.assertIn("--center_x", cmd)
        self.assertIn("10.0", cmd)
        self.assertTrue(Path(cmd[0]).name.lower().startswith("python"))


class TestReceptorSanitize(unittest.TestCase):
    def test_keeps_only_first_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "receptor.pdbqt"
            path.write_text(
                "MODEL 1\n"
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
                "ENDMDL\n"
                "MODEL 2\n"
                "ATOM      2  N   ALA A   2       1.000   1.000   1.000\n"
                "ENDMDL\n",
                encoding="utf-8",
            )
            _sanitize_receptor_models(path)
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("MODEL", text)
            self.assertIn("ALA A   1", text)
            self.assertNotIn("ALA A   2", text)

    def test_leaves_single_model_pdbqt_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "receptor.pdbqt"
            body = "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
            path.write_text(body, encoding="utf-8")
            _sanitize_receptor_models(path)
            self.assertEqual(path.read_text(encoding="utf-8"), body)


class TestEndToEndDockAndAnalyze(unittest.TestCase):
    def test_dock_and_analyze(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            workdir = tmp_path / "work"
            (workdir / "data" / "receptors").mkdir(parents=True)
            (workdir / "data" / "ligands" / "prepared").mkdir(parents=True)
            receptor = workdir / "data" / "receptors" / "receptor.pdbqt"
            receptor.write_text(
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
                "TER\n",
                encoding="utf-8",
            )
            lig_pdbqt = workdir / "data" / "ligands" / "prepared" / "L1.pdbqt"
            lig_pdbqt.write_text("ATOM      1  C   LIG L   1       0.000   0.000   0.000\n", encoding="utf-8")
            manifest = workdir / "data" / "ligands" / "prepared" / "manifest.csv"
            manifest.write_text(
                "id,smiles,heavy_atoms,rotatable_bonds,pdbqt,status,error\n"
                f"L1,CCO,3,0,{lig_pdbqt},ok,\n",
                encoding="utf-8",
            )

            fake_vina = tmp_path / "fake_vina.py"
            fake_vina.write_text(
                _FAKE_VINA,
                encoding="utf-8",
            )

            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "outdir": "outputs/run_001",
                    "receptor": "data/receptors/receptor.pdbqt",
                    "ligand": "data/ligands/library.sdf",
                    "executable": str(fake_vina),
                    "max_workers": 1,
                    "scoring": "",
                },
            )
            cfg_path = workdir / "config" / "test.json"
            save_config(cfg, cfg_path)
            cfg = load_config(cfg_path)

            run_docking(cfg, LOG)
            results = cfg.results_path()
            self.assertTrue(results.exists())
            rows = results.read_text(encoding="utf-8").splitlines()
            self.assertTrue(any("ok" in row and "-8.3" in row for row in rows))

            analyze_results(cfg, LOG)
            ranked = cfg.analysis_dir() / "data" / "fig_46_47_ranked_results.csv"
            self.assertTrue(ranked.exists())
            hits = cfg.analysis_dir() / "data" / "fig_47_top_hits.csv"
            self.assertTrue(hits.exists())
            self.assertTrue((cfg.analysis_dir() / "summary.json").exists())


class TestDockResumeFresh(unittest.TestCase):
    def test_resume_false_does_not_accumulate_old_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            workdir = tmp_path / "work"
            (workdir / "data" / "receptors").mkdir(parents=True)
            (workdir / "data" / "ligands" / "prepared").mkdir(parents=True)
            receptor = workdir / "data" / "receptors" / "receptor.pdbqt"
            receptor.write_text(
                "ATOM      1  N   ALA A   1       0.000   0.000   0.000\n"
                "TER\n",
                encoding="utf-8",
            )
            lig_pdbqt = workdir / "data" / "ligands" / "prepared" / "L1.pdbqt"
            lig_pdbqt.write_text(
                "ATOM      1  C   LIG L   1       0.000   0.000   0.000\n",
                encoding="utf-8",
            )
            manifest = workdir / "data" / "ligands" / "prepared" / "manifest.csv"
            manifest.write_text(
                "id,smiles,heavy_atoms,rotatable_bonds,pdbqt,status,error\n"
                f"L1,CCO,3,0,{lig_pdbqt},ok,\n",
                encoding="utf-8",
            )
            fake_vina = tmp_path / "fake_vina.py"
            fake_vina.write_text(_FAKE_VINA, encoding="utf-8")
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "outdir": "outputs/run_001",
                    "receptor": "data/receptors/receptor.pdbqt",
                    "ligand": "data/ligands/library.sdf",
                    "executable": str(fake_vina),
                    "max_workers": 1,
                    "scoring": "",
                    "resume": False,
                },
            )
            run_docking(cfg, LOG)
            run_docking(cfg, LOG)
            rows = list(csv.DictReader(cfg.results_path().open("r", newline="")))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["id"], "L1")


class TestMLTorchFallback(unittest.TestCase):
    def test_torch_fallback_returns_sklearn_mlp(self):
        from sklearn.neural_network import MLPClassifier

        class _FailingTorch:
            def __init__(self, *args, **kwargs):
                raise RuntimeError("torch unavailable")

        with mock.patch.object(docking_ml, "_TorchMLP", _FailingTorch):
            model = docking_ml._build_model(
                "torch",
                "classification",
                hidden=16,
                epochs=5,
                random_state=42,
            )
        self.assertIsInstance(model, MLPClassifier)


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


if __name__ == "__main__":
    unittest.main()
