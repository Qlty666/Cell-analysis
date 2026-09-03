#!/usr/bin/env python3
"""Install project Codex skills into the current user's skill root."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKILL_NAMES = [
    "liver-expression-analysis",
    "liver-virtual-screening",
    "liver-full-pipeline",
    "liver-dataset-search",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Install liver cancer suite skills into Codex"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing skill of the same name",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="print the skills that would be installed and exit",
    )
    args = parser.parse_args(argv)

    codex_home = Path(os.environ.get("CODEX_HOME", "")).resolve() if os.environ.get("CODEX_HOME") else Path.home() / ".codex"
    skills_root = codex_home / "skills"
    if args.list:
        for name in SKILL_NAMES:
            print(f"{name} -> {skills_root / name}")
        return 0

    installed: list[str] = []
    skipped: list[str] = []
    failed: list[str] = []
    skills_root.mkdir(parents=True, exist_ok=True)
    for name in SKILL_NAMES:
        source = ROOT / "skills" / name
        target = skills_root / name
        if not source.is_dir():
            failed.append(name)
            continue
        if target.exists() and not args.force:
            skipped.append(name)
            continue
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
        installed.append(name)

    for name in installed:
        print(f"installed: {name} -> {skills_root / name}")
    for name in skipped:
        print(f"skipped (already installed, use --force to update): {name}")
    for name in failed:
        print(f"failed: missing source skill {ROOT / 'skills' / name}", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
