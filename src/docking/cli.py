#!/usr/bin/env python3
"""Command line interface for the docking pipeline."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import (
    analysis,
    box,
    docking,
    evidence,
    handoff,
    knockout,
    ligands,
    ml,
    network_toxicology,
    pipeline,
    receptor,
    redock,
    report,
    signal_detection,
    validation,
)
from .config import load_config, save_config
from .environment import check_environment
from .utils import setup_logging

APP_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_docking",
        description="Liver cancer virtual screening pipeline (AutoDock Vina)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    for name, help_text in [
        ("pipeline", "run all docking stages with resume support"),
        ("prepare-receptor", "prepare receptor PDBQT"),
        ("prepare-ligands", "prepare ligand library PDBQT files"),
        ("dock", "run AutoDock Vina over prepared ligands"),
        ("analyze", "rank docking results and build reports"),
        ("evidence", "collect target/ligand evidence from public databases"),
        ("ml-train", "train ML/DL rescoring model"),
        ("ml-predict", "apply trained ML/DL model to docking results"),
        ("export-md", "export top hits to Amber/GROMACS templates"),
        ("export-external", "export UniDock/HDOCK/HADDOCK templates"),
        ("redock", "re-dock top hits with higher exhaustiveness"),
        ("report", "generate HTML summary report"),
        ("virtual-knockout", "score virtual gene knockouts from expression/DepMap data"),
        ("network", "compound-disease overlap, PPI hub and C-T-P-D network"),
        ("faers", "FAERS-style disproportionality signal detection"),
        ("export-validation", "export wet-lab validation plan for top targets"),
        ("cell-feedback", "re-run single-cell analysis from knockout/docking results"),
        ("detect-box", "detect docking box from cocrystal ligand and save config"),
        ("check-env", "check Python packages and external tools"),
        ("check-cadd", "check CADD workflow skills and ML libraries"),
        ("init", "create the workdir skeleton"),
    ]:
        sub = subparsers.add_parser(name, help=help_text)
        _add_common(sub)

    args = parser.parse_args(argv)
    log = setup_logging(verbose=args.verbose)
    overrides = {
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
        "top_n": args.top_n,
        "executable": args.executable,
        "scoring": args.scoring,
        "model": args.model,
        "label_column": args.label_column,
        "training_csv": args.training_csv,
        "epochs": args.epochs,
        "hidden_size": args.hidden_size,
        "uniprot": args.uniprot,
        "pdb": args.pdb,
        "chembl_target": args.chembl_target,
        "target_name": args.target_name,
        "ligand_name": args.ligand_name,
        "ligand_smiles": args.ligand_smiles,
        "max_items": args.max_items,
        "expression_csv": args.expression_csv,
        "metadata_csv": args.metadata_csv,
        "depmap_csv": args.depmap_csv,
        "prognosis_csv": args.prognosis_csv,
        "druggability_csv": args.druggability_csv,
        "off_target_csv": args.off_target_csv,
        "ppi_network_csv": args.ppi_network_csv,
        "cell_type_column": args.cell_type_column,
        "group_column": args.group_column,
        "case_label": args.case_label,
        "normal_label": args.normal_label,
        "ko_top_n": args.ko_top_n,
        "validation_top_n": args.validation_top_n,
        "compound_name": args.compound_name,
        "disease_name": args.disease_name,
        "compound_targets_csv": args.compound_targets_csv,
        "disease_genes_csv": args.disease_genes_csv,
        "disease_gene_column": args.disease_gene_column,
        "network_output_dir": args.network_output_dir,
        "faers_input": args.faers_input,
        "faers_drug_column": args.faers_drug_column,
        "faers_event_column": args.faers_event_column,
        "faers_count_column": args.faers_count_column,
        "faers_min_count": args.faers_min_count,
    }
    overrides = {key: value for key, value in overrides.items() if value is not None}
    cfg = load_config(args.config, overrides)

    if args.command == "check-env":
        return 0 if print_environment(check_environment()) else 1
    if args.command == "check-cadd":
        _print_cadd()
        return 0
    if args.command == "init":
        _init_workdir(cfg, log)
        return 0
    if args.command == "pipeline":
        pipeline.run_pipeline(cfg, force=args.force, start_stage=args.start_stage)
        return 0

    if args.command == "prepare-receptor":
        receptor.prepare_receptor(cfg, log)
    elif args.command == "prepare-ligands":
        ligands.prepare_ligands(cfg, log)
    elif args.command == "dock":
        docking.run_docking(cfg, log)
    elif args.command == "analyze":
        analysis.analyze_results(cfg, log)
    elif args.command == "evidence":
        evidence.gather_evidence(cfg, log)
    elif args.command == "ml-train":
        ml.train_ml(
            cfg,
            log,
            model_type=args.model,
            label_column=args.label_column,
            training_csv=args.training_csv,
        )
    elif args.command == "ml-predict":
        ml.predict_ml(cfg, log)
    elif args.command == "export-md":
        handoff.export_md(cfg, log)
    elif args.command == "export-external":
        handoff.export_external(cfg, log)
    elif args.command == "redock":
        redock.run_redock(cfg, log)
    elif args.command == "report":
        report.generate_report(cfg, log)
    elif args.command == "virtual-knockout":
        knockout.run_knockout(cfg, log)
    elif args.command == "network":
        network_toxicology.run_network_toxicology(cfg, log)
    elif args.command == "faers":
        signal_detection.run_faers(cfg, log)
    elif args.command == "export-validation":
        validation.export_validation(cfg, log)
    elif args.command == "cell-feedback":
        from pipeline.cell_feedback import run_cell_feedback

        single_cell_root = getattr(args, "single_cell_root", "")
        if not single_cell_root:
            parser.error("cell-feedback requires --single-cell-root")
        summary = run_cell_feedback(
            cfg.workdir,
            Path(single_cell_root).resolve(),
            top_n=args.feedback_top_n or 12,
            max_features=args.feedback_max_features or 8,
            timeout_seconds=args.feedback_timeout or 3600,
            species=getattr(args, "feedback_species", "hs") or "hs",
        )
        log.info(
            "cell feedback %s: %s matched genes",
            summary.get("status", "unknown"),
            summary.get("genes_matched", 0),
        )
    elif args.command == "detect-box":
        box.detect_and_update_config(cfg, log)
    else:
        parser.error(f"unknown command: {args.command}")
    log.info("%s complete", args.command)
    return 0


def _add_common(sub: argparse.ArgumentParser) -> None:
    sub.add_argument(
        "--config",
        default=str(DEFAULT_CONFIG),
        help="path to docking config JSON/YAML",
    )
    sub.add_argument("--workdir", help="override working directory")
    sub.add_argument("--outdir", help="override output directory")
    sub.add_argument("--receptor", help="override receptor input file")
    sub.add_argument("--ligand", help="override ligand library file")
    sub.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    sub.add_argument("--size", nargs=3, type=float, metavar=("X", "Y", "Z"))
    sub.add_argument("--exhaustiveness", type=int)
    sub.add_argument("--num-modes", type=int)
    sub.add_argument("--energy-range", type=float)
    sub.add_argument("--cpu", type=int)
    sub.add_argument("--max-workers", type=int)
    sub.add_argument("--seed", type=int)
    sub.add_argument("--max-ligands", type=int)
    sub.add_argument("--cutoff", type=float)
    sub.add_argument("--top-n", type=int)
    sub.add_argument("--executable", help="Vina executable or script path")
    sub.add_argument("--scoring", help="Vina scoring function, e.g. vina/vinardo")
    sub.add_argument("--model", help="ML model: rf, gbm, mlp, lasso_svm or torch")
    sub.add_argument("--label-column", help="training label column")
    sub.add_argument("--training-csv", help="training CSV path")
    sub.add_argument("--epochs", type=int, help="MLP/torch training epochs")
    sub.add_argument("--hidden-size", type=int, help="MLP hidden layer size")
    sub.add_argument("--uniprot", help="UniProt accession for evidence")
    sub.add_argument("--pdb", help="PDB ID for evidence")
    sub.add_argument("--chembl-target", help="ChEMBL target ID")
    sub.add_argument("--target-name", help="target name for evidence report")
    sub.add_argument("--ligand-name", help="ligand name for ChEBI/PubChem")
    sub.add_argument("--ligand-smiles", help="ligand SMILES for PubChem")
    sub.add_argument("--max-items", type=int, help="max records per database")
    sub.add_argument(
        "--expression-csv",
        help="expression matrix or long-format CSV for virtual knockout",
    )
    sub.add_argument("--metadata-csv", help="sample metadata CSV (sample + group)")
    sub.add_argument("--depmap-csv", help="DepMap CRISPR gene effect CSV")
    sub.add_argument(
        "--prognosis-csv",
        help="prognosis CSV (gene + hazard ratio) for target scoring",
    )
    sub.add_argument(
        "--druggability-csv",
        help="druggability CSV (gene + known ligands/structures/assays)",
    )
    sub.add_argument(
        "--off-target-csv",
        help="off-target CSV (gene + paralogs/safety concern)",
    )
    sub.add_argument(
        "--cell-type-column",
        help="cell-type column in metadata for specificity scoring",
    )
    sub.add_argument("--group-column", help="group/condition column in metadata")
    sub.add_argument("--case-label", help="case/tumor group label")
    sub.add_argument("--normal-label", help="normal/control group label")
    sub.add_argument("--ko-top-n", type=int, help="top N genes for knockout report")
    sub.add_argument(
        "--ppi-network-csv",
        help="STRING-style PPI edge table used by virtual knockout and network",
    )
    sub.add_argument(
        "--validation-top-n",
        type=int,
        help="top N genes for the wet-lab validation plan",
    )
    sub.add_argument("--compound-name", help="compound name for network toxicology")
    sub.add_argument("--disease-name", help="disease name for network toxicology")
    sub.add_argument(
        "--compound-targets-csv",
        help="compound-target gene table for network toxicology",
    )
    sub.add_argument(
        "--disease-genes-csv",
        help="disease gene list for network toxicology",
    )
    sub.add_argument(
        "--disease-gene-column",
        help="gene column in the disease gene list",
    )
    sub.add_argument(
        "--network-output-dir",
        help="output directory for network toxicology results",
    )
    sub.add_argument(
        "--faers-input",
        help="FAERS-style event CSV (drug, event, optional count)",
    )
    sub.add_argument("--faers-drug-column", help="drug column in FAERS CSV")
    sub.add_argument("--faers-event-column", help="event column in FAERS CSV")
    sub.add_argument("--faers-count-column", help="count column in FAERS CSV")
    sub.add_argument("--faers-min-count", type=int, help="minimum count for signals")
    sub.add_argument(
        "--single-cell-root",
        help="single-cell output root containing results/checkpoints",
    )
    sub.add_argument("--feedback-top-n", type=int, help="top N feedback genes")
    sub.add_argument(
        "--feedback-max-features",
        type=int,
        help="max genes shown in UMAP/DotPlot feedback figures",
    )
    sub.add_argument(
        "--feedback-timeout",
        type=int,
        help="cell feedback R analysis timeout in seconds",
    )
    sub.add_argument(
        "--feedback-species",
        choices=["hs", "mm"],
        default="hs",
        help="species used by cell feedback enrichment (default: hs)",
    )
    sub.add_argument("--force", action="store_true", help="rerun stages from scratch")
    sub.add_argument("--start-stage", default=None, help="stage code to start from")
    sub.add_argument("--verbose", action="store_true")


def _init_workdir(cfg, log) -> None:
    for path in [
        cfg.workdir / "data" / "receptors",
        cfg.workdir / "data" / "ligands",
        cfg.output_dir,
        cfg.logs_dir(),
    ]:
        path.mkdir(parents=True, exist_ok=True)
    save_config(cfg, cfg.workdir / "config" / "docking_config.json")
    log.info("initialized docking workdir: %s", cfg.workdir)


def _print_environment(checks: list[dict]) -> None:
    print()
    print("Docking environment check")
    print("=" * 72)
    failed = 0
    for check in checks:
        status = "OK  " if check["ok"] else "FAIL"
        if not check["ok"]:
            failed += 1
        print(f"[{status}] {check['name']:<28} {check['detail']}")
    print("=" * 72)
    if failed:
        print("Missing items. Run launchers\\install_dock_dependencies.bat first.")
        print("AutoDock Vina can be installed via conda-forge: autodock-vina")
    else:
        print("Environment ready")
    print()


def print_environment(checks: list[dict]) -> bool:
    """Print environment checks and return whether every check passed."""
    _print_environment(checks)
    return all(check["ok"] for check in checks)


def _print_cadd() -> None:
    import importlib.util

    print()
    print("CADD workflow check")
    print("=" * 72)
    for name, script in evidence.SKILL_SCRIPTS.items():
        ok = script.exists()
        print(f"[{'OK  ' if ok else 'FAIL'}] {name:<24} {script}")
    for package in ["sklearn", "torch", "joblib", "rdkit"]:
        ok = importlib.util.find_spec(package) is not None
        print(f"[{'OK  ' if ok else 'FAIL'}] {package:<24} installed")
    print("=" * 72)
    print()


if __name__ == "__main__":
    sys.exit(main())
