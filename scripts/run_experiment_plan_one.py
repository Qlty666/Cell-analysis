#!/usr/bin/env python3
"""Run the 6PPD-Q / NAFLD experiment-plan-one analysis."""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from experiment_plan_one.pipeline import (  # noqa: E402
    STAGES,
    ExperimentPlanOne,
    load_config,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the local 6PPD-Q / NAFLD pipeline defined in 实验方案一: "
            "compound targets, disease targets, PPI, ML, single-cell, docking "
            "and GROMACS preparation."
        )
    )
    parser.add_argument(
        "--output-root",
        default=None,
        help="result root; defaults to output_root in the JSON config",
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "experiment_plan_one.json"),
        help="pipeline JSON configuration",
    )
    parser.add_argument(
        "--stage",
        action="append",
        choices=STAGES,
        help="run only this stage; repeat for multiple stages",
    )
    parser.add_argument(
        "--start-stage",
        choices=STAGES,
        help="start at this stage and run all following stages",
    )
    parser.add_argument(
        "--skip-stage",
        action="append",
        choices=STAGES,
        default=[],
        help="skip a stage; repeat for multiple stages",
    )
    parser.add_argument("--force", action="store_true", help="rerun completed stages")
    parser.add_argument("--dry-run", action="store_true", help="print the selected stages only")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def resolve_stages(args: argparse.Namespace) -> list[str]:
    if args.stage:
        selected = list(dict.fromkeys(args.stage))
    elif args.start_stage:
        selected = list(STAGES[STAGES.index(args.start_stage) :])
    else:
        selected = list(STAGES)
    skipped = set(args.skip_stage)
    return [stage for stage in selected if stage not in skipped]


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        config = load_config(Path(args.config))
        output_value = args.output_root or config.get("output_root") or "../y1"
        output_root = Path(str(output_value)).expanduser()
        if not output_root.is_absolute():
            output_root = (ROOT / output_root).resolve()
        stages = resolve_stages(args)
        print(f"experiment plan one root: {output_root}")
        print("stages: " + ", ".join(stages))
        if args.dry_run:
            return 0
        runner = ExperimentPlanOne(
            output_root,
            config,
            force=args.force,
            verbose=args.verbose,
        )
        runner.run(stages)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
