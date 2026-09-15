#!/usr/bin/env python3
"""Print environment readiness for the docking pipeline."""

import argparse
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking.cli import print_environment  # noqa: E402
from docking.environment import check_environment  # noqa: E402


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
    )
    parser.add_argument(
        "--skip-ml",
        action="store_true",
        help="ignore optional ML packages (scikit-learn, joblib, torch)",
    )
    args = parser.parse_args()
    checks = check_environment()
    if args.skip_ml:
        checks = [
            item
            for item in checks
            if item["name"] not in ("scikit-learn", "joblib", "torch")
        ]
    sys.exit(0 if print_environment(checks) else 1)
