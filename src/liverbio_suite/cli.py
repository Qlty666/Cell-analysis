#!/usr/bin/env python3
"""One command for the expression, docking, full-pipeline and web tools.

This is a thin dispatcher: each suite feature keeps its own tested CLI, and
``liverbio`` only forwards arguments to that feature.  Existing launch scripts
and web entry points remain unchanged.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

ENTRYPOINTS = {
    "expression": ROOT / "scripts" / "run_pipeline.py",
    "docking": ROOT / "scripts" / "run_docking.py",
    "full": ROOT / "scripts" / "run_full_pipeline.py",
    "datasets": ROOT / "scripts" / "search_datasets.py",
    "web": ROOT / "web" / "web_ui.py",
    "install-skills": ROOT / "scripts" / "install_codex_skills.py",
}

ENV_CHECKS = {
    "pipeline": ROOT / "launchers" / "check_pipeline_environment.py",
    "docking": ROOT / "launchers" / "check_dock_environment.py",
    "dock": ROOT / "launchers" / "check_dock_environment.py",
}

FEATURE_SUMMARY = (
    "expression: single-cell / bulk RNA-seq / microarray expression pipeline\n"
    "docking:    virtual screening, knockout, network and FAERS commands\n"
    "full:       integrated expression-to-docking pipeline\n"
    "datasets:   GEO / BioStudies / Expression Atlas search and download\n"
    "web:        local web console for the whole suite\n"
    "doctor:     environment check (pipeline, docking, or all)\n"
    "install-skills: copy project Codex skills into the user skill root\n"
)

USAGE = """Liver Cancer Bioinformatics Suite

Usage:
  liverbio <command> [feature arguments ...]
  liverbio help
  liverbio version

Commands:
%s
Examples:
  liverbio expression GSE125449 --output ../liver_cancer --species auto
  liverbio full --accession GSE125449 --output ../liver_cancer
  liverbio docking pipeline --config config/docking_config.json
  liverbio datasets --disease "liver cancer" --max-results 20
  liverbio web --page full
  liverbio doctor
  liverbio install-skills
""" % FEATURE_SUMMARY


def project_root() -> Path:
    """Absolute project root that contains scripts, web and src."""
    return ROOT


def project_version() -> str:
    """Read the version from src/docking/__init__.py without importing more."""
    init_path = ROOT / "src" / "docking" / "__init__.py"
    try:
        text = init_path.read_text(encoding="utf-8")
    except OSError:
        return "unknown"
    match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', text)
    return match.group(1) if match else "unknown"


def print_usage() -> None:
    print(USAGE)


def run_script(script: Path, args: list[str]) -> int:
    if not script.is_file():
        print(f"ERROR: entry point not found: {script}", file=sys.stderr)
        return 1
    cmd = [sys.executable, str(script), *args]
    return subprocess.call(cmd, cwd=ROOT)


def run_doctor(kind: str) -> int:
    if kind == "all":
        checks = list(ENV_CHECKS.values())
    else:
        script = ENV_CHECKS.get(kind)
        if script is None:
            print("doctor accepts: pipeline, docking, or all", file=sys.stderr)
            return 2
        checks = [script]
    results = [run_script(script, []) for script in checks]
    return 0 if all(code == 0 for code in results) else 1


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("help", "--help", "-h"):
        print_usage()
        return 0
    if args[0] in ("version", "--version", "-V"):
        print(f"Liver Cancer Bioinformatics Suite {project_version()}")
        return 0
    if args[0] == "doctor":
        kind = args[1] if len(args) > 1 else "all"
        return run_doctor(kind)
    script = ENTRYPOINTS.get(args[0])
    if script is None:
        print(f"ERROR: unknown command: {args[0]}", file=sys.stderr)
        print_usage()
        return 2
    return run_script(script, args[1:])


if __name__ == "__main__":
    sys.exit(main())
