#!/usr/bin/env python3
"""Validate the docking pipeline with a fake Vina executable."""

import logging
import sys
import tempfile
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking.analysis import analyze_results  # noqa: E402
from docking.config import load_config, save_config  # noqa: E402
from docking.docking import run_docking  # noqa: E402
from docking.utils import setup_logging  # noqa: E402


def main() -> int:
    log = setup_logging(verbose=True)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        workdir = tmp_path / "work"
        (workdir / "data" / "receptors").mkdir(parents=True)
        (workdir / "data" / "ligands" / "prepared").mkdir(parents=True)
        receptor = workdir / "data" / "receptors" / "receptor.pdbqt"
        receptor.write_text(
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000\nTER\n",
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
            Path(__file__).resolve().parent / "config" / "docking_config.json",
            {
                "workdir": str(workdir),
                "receptor": "data/receptors/receptor.pdbqt",
                "executable": str(fake_vina),
                "max_workers": 1,
                "scoring": "",
            },
        )
        save_config(cfg, workdir / "config" / "validate.json")
        cfg = load_config(workdir / "config" / "validate.json")

        run_docking(cfg, log)
        summary = analyze_results(cfg, log)
        log.info("validation passed: %s", summary)
        print(f"Validation passed. Output: {cfg.output_dir}")
        return 0


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
    sys.exit(main())
