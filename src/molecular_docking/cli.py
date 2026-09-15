"""Command line interface for the standalone molecular docking board."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from docking.box import detect_and_update_config
from docking.cli import _add_common
from docking.environment import check_environment
from docking.utils import DockingError, setup_logging, write_json

from .config import DEFAULT_CONFIG, load_config, save_config
from .pipeline import run_pipeline


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_molecular_docking",
        description="Standalone molecular docking board (AutoDock Vina)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name, help_text in [
        ("init", "create the molecular docking workdir skeleton"),
        ("pipeline", "run all docking stages with resume support"),
        ("prepare-receptor", "prepare receptor PDBQT"),
        ("prepare-ligands", "prepare ligand library PDBQT files"),
        ("dock", "run AutoDock Vina over prepared ligands"),
        ("analyze", "rank docking results and build result files"),
        ("redock", "re-dock top hits with higher exhaustiveness"),
        ("report", "generate the molecular docking HTML report"),
        ("detect-box", "detect docking box from cocrystal ligand and save config"),
        ("check-env", "check Python packages and external docking tools"),
    ]:
        sub = subparsers.add_parser(name, help=help_text)
        _add_common(
            sub,
            extended=False,
            config_default=DEFAULT_CONFIG,
            config_help="path to molecular docking config JSON",
        )

    args = parser.parse_args(argv)
    log = setup_logging(verbose=args.verbose)
    overrides = {
        key: value
        for key, value in {
            "workdir": args.workdir,
            "outdir": args.outdir,
            "receptor": args.receptor,
            "ligand": args.ligand,
            "center": args.center,
            "size": args.size,
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "energy_range": args.energy_range,
            "cpu": args.cpu,
            "max_workers": args.max_workers,
            "seed": args.seed,
            "max_ligands": args.max_ligands,
            "cutoff": args.cutoff,
            "moderate_cutoff": args.moderate_cutoff,
            "top_n": args.top_n,
            "executable": args.executable,
            "scoring": args.scoring,
        }.items()
        if value is not None
    }
    cfg = load_config(args.config, overrides)

    if args.command == "check-env":
        failed = 0
        print()
        print("Molecular docking environment check")
        print("=" * 72)
        for check in check_environment():
            ok = check["ok"]
            failed += 0 if ok else 1
            print(f"[{'OK  ' if ok else 'FAIL'}] {check['name']:<28} {check['detail']}")
        print("=" * 72)
        return 1 if failed else 0

    if args.command == "init":
        for path in [
            cfg.workdir / "data" / "receptors",
            cfg.workdir / "data" / "ligands",
            cfg.output_dir,
            cfg.logs_dir(),
        ]:
            path.mkdir(parents=True, exist_ok=True)
        save_config(cfg, cfg.workdir / "config" / "molecular_docking_config.json")
        log.info("initialized molecular docking workdir: %s", cfg.workdir)
        return 0

    if args.command == "pipeline":
        run_pipeline(cfg, force=args.force, start_stage=args.start_stage)
        return 0

    if args.command == "prepare-receptor":
        from docking.receptor import prepare_receptor

        prepare_receptor(cfg, log)
    elif args.command == "prepare-ligands":
        from docking.ligands import prepare_ligands

        prepare_ligands(cfg, log)
    elif args.command == "dock":
        from docking.docking import run_docking

        run_docking(cfg, log)
    elif args.command == "analyze":
        from docking.analysis import analyze_results

        analyze_results(cfg, log)
    elif args.command == "redock":
        from docking.redock import run_redock

        run_redock(cfg, log)
    elif args.command == "report":
        from .report import generate_report

        generate_report(cfg, log)
    elif args.command == "detect-box":
        _detect_box(cfg, log)
    else:
        parser.error(f"unknown command: {args.command}")

    log.info("%s complete", args.command)
    return 0


def _detect_box(cfg, log):
    center, size, mode = detect_and_update_config(
        cfg,
        log,
        save=save_config,
    )
    write_json(
        cfg.output_dir / "detect_box_result.json",
        {"center": center, "size": size, "mode": mode},
    )
    return center, size, mode


if __name__ == "__main__":
    sys.exit(main())
