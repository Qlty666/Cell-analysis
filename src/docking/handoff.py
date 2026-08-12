#!/usr/bin/env python3
"""Export docking hits to MD and external docking workflow templates."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from .config import ResolvedConfig
from .utils import DockingError


def export_md(cfg: ResolvedConfig, log) -> Path:
    source = _pick_results(cfg)
    top_n = int(cfg.get("md", "top_n", 10))
    df = pd.read_csv(source, dtype={"id": str}).head(top_n)
    out = cfg.md_export_dir()
    poses = out / "poses"
    amber = out / "amber"
    gromacs = out / "gromacs"
    for folder in (poses, amber, gromacs):
        folder.mkdir(parents=True, exist_ok=True)

    copied = 0
    for _, row in df.iterrows():
        pose = row.get("pose_file", "")
        if pose and Path(pose).exists():
            shutil.copyfile(pose, poses / f"{row['id']}.pdbqt")
            copied += 1

    _write_amber_templates(amber)
    _write_gromacs_templates(gromacs)
    _write_md_readme(out, df, copied)
    log.info(
        "MD handoff exported: %s hits, %s poses -> %s",
        len(df),
        copied,
        out,
    )
    return out


def export_external(cfg: ResolvedConfig, log) -> Path:
    out = cfg.output_dir / "external"
    out.mkdir(parents=True, exist_ok=True)
    receptor = cfg.receptor_output()
    manifest = cfg.manifest_path()
    _write_external_readme(
        out / "unidock_pro" / "README.md",
        "UniDock-Pro",
        [
            "GPU classical docking / similarity search / hybrid docking",
            f"receptor: {receptor}",
            f"ligand manifest: {manifest}",
            "run_unidock --receptor receptor.pdbqt --ligand ligand_index.txt "
            "--search_box center size --search_mode classical",
        ],
    )
    _write_external_readme(
        out / "hdock" / "README.md",
        "HDOCK",
        [
            "protein-protein or protein-nucleic acid local docking",
            f"receptor: {cfg.receptor_input()}",
            "HDOCKlite receptor.pdb ligand.pdb",
            "outputs: hdock.out, models.pdb",
        ],
    )
    _write_external_readme(
        out / "haddock" / "README.md",
        "HADDOCK",
        [
            "information-driven docking with active/passive residues",
            "requires receptor.pdb, ligand.pdb and AIR restraints",
            "configure active residues from RCSB/UniProt evidence",
            "run haddock2.5 -> cluster ranking -> top models",
        ],
    )
    log.info("external docking templates exported -> %s", out)
    return out


def _pick_results(cfg: ResolvedConfig) -> Path:
    reports = cfg.reports_dir()
    for name in (
        "fig_46_47_ranked_results.csv",
        "fig_48_diverse_hits.csv",
        "fig_47_top_hits.csv",
    ):
        path = reports / "01_analysis" / "data" / name
        if path.exists():
            return path
    if cfg.results_path().exists():
        return cfg.results_path()
    raise DockingError("no docking results found; run dock/analyze first")


def _write_amber_templates(amber: Path) -> None:
    (amber / "tleap.in").write_text(
        "source leaprc.protein.ff19SB\n"
        "source leaprc.water.tip3p\n"
        "source leaprc.gaff2\n"
        "loadamberparams ligand.frcmod\n"
        "complex = loadpdb complex.pdb\n"
        "solvateoct complex TIP3PBOX 12.0\n"
        "addions complex Na+ 0\n"
        "addions complex Cl- 0\n"
        "saveamberparm complex complex.prmtop complex.rst7\n"
        "quit\n",
        encoding="utf-8",
    )
    for name, content in {
        "run_min.sh": (
            "#!/bin/bash\n"
            "# Minimization\n"
            "pmemd.cuda -O -i min.in -o min.out -p complex.prmtop "
            "-c complex.rst7 -r min.rst7 -ref complex.rst7\n"
        ),
        "run_equil.sh": (
            "#!/bin/bash\n"
            "# NVT then NPT equilibration\n"
            "pmemd.cuda -O -i nvt.in -o nvt.out -p complex.prmtop "
            "-c min.rst7 -r nvt.rst7 -ref min.rst7\n"
            "pmemd.cuda -O -i npt.in -o npt.out -p complex.prmtop "
            "-c nvt.rst7 -r npt.rst7 -ref nvt.rst7\n"
        ),
        "run_prod.sh": (
            "#!/bin/bash\n"
            "# Production MD\n"
            "pmemd.cuda -O -i prod.in -o prod.out -p complex.prmtop "
            "-c npt.rst7 -r prod.rst7 -x prod.nc\n"
        ),
        "analyze_cpptraj.in": (
            "parm complex.prmtop\n"
            "trajin prod.nc\n"
            "rmsd protein out rmsd.dat :1-100&!@H=\n"
            "rmsd ligand out ligand_rmsd.dat :LIG&!@H=\n"
            "run\n"
        ),
    }.items():
        (amber / name).write_text(content, encoding="utf-8")


def _write_gromacs_templates(gromacs: Path) -> None:
    mdp = {
        "em.mdp": (
            "integrator = steep\n"
            "nsteps = 5000\n"
            "emtol = 1000\n"
            "cutoff-scheme = Verlet\n"
        ),
        "nvt.mdp": (
            "integrator = md\n"
            "nsteps = 50000\n"
            "tcoupl = v-rescale\n"
            "tc-grps = Protein_LIG Water_and_ions\n"
        ),
        "npt.mdp": (
            "integrator = md\n"
            "nsteps = 50000\n"
            "pcoupl = Parrinello-Rahman\n"
            "refcoord-scaling = com\n"
        ),
        "md.mdp": (
            "integrator = md\n"
            "nsteps = 500000\n"
            "tcoupl = v-rescale\n"
            "pcoupl = Parrinello-Rahman\n"
        ),
    }
    for name, content in mdp.items():
        (gromacs / name).write_text(content, encoding="utf-8")
    (gromacs / "run_gmx.bat").write_text(
        "@echo off\r\n"
        "gmx pdb2gmx -f complex.pdb -o complex.gro -ff amber99sb-ildn "
        "-water tip3p\r\n"
        "gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt cubic\r\n"
        "gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top\r\n"
        "gmx grompp -f em.mdp -c solv.gro -p topol.top -o em.tpr\r\n"
        "gmx mdrun -deffnm em\r\n",
        encoding="utf-8",
    )


def _write_md_readme(out: Path, df: pd.DataFrame, copied: int) -> None:
    lines = [
        "# MD Handoff",
        "",
        f"Exported {len(df)} top hits, {copied} pose files.",
        "",
        "| rank | id | affinity |",
        "|---|---|---|",
    ]
    for _, row in df.iterrows():
        lines.append(f"| {row.get('rank', '')} | {row['id']} | {row.get('affinity', '')} |")
    lines += [
        "",
        "## Amber",
        "1. Prepare ligand with antechamber and write ligand.frcmod.",
        "2. Place complex.pdb in amber/ and run tleap.",
        "3. Run run_min.sh -> run_equil.sh -> run_prod.sh.",
        "4. Analyze with cpptraj using analyze_cpptraj.in.",
        "",
        "## GROMACS",
        "1. Convert receptor PDB and ligand to a single complex.pdb.",
        "2. Generate ligand topology (CGenFF/GAFF) and merge into topol.top.",
        "3. Run run_gmx.bat or the equivalent Linux commands.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines), encoding="utf-8")


def _write_external_readme(path: Path, title: str, bullets: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {title}", ""]
    lines += [f"- {item}" for item in bullets]
    path.write_text("\n".join(lines), encoding="utf-8")
