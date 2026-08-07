#!/usr/bin/env python3
"""Print environment readiness for the docking pipeline."""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking.cli import _print_environment  # noqa: E402
from docking.environment import check_environment  # noqa: E402


if __name__ == "__main__":
    _print_environment(check_environment())
