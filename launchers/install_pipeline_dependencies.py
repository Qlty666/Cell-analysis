#!/usr/bin/env python3
"""Install all R dependencies required by the pipeline."""

import subprocess
import sys
import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from common.env import find_rscript  # noqa: E402

# Never let a child installer hang forever.
R_DEPS_TIMEOUT_SECONDS = 3600
PIP_TIMEOUT_SECONDS = 3600
CHECK_TIMEOUT_SECONDS = 900


def _cmd_text(cmd: list[str]) -> str:
    return " ".join(str(part) for part in cmd)


def _run_logged(cmd: list[str], timeout: int, **kwargs) -> subprocess.CompletedProcess:
    print(f"> {_cmd_text(cmd)} (timeout {timeout}s)", flush=True)
    return subprocess.run(cmd, timeout=timeout, **kwargs)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="")
    args = parser.parse_args()

    rscript = find_rscript()
    if not rscript:
        print("Rscript not found. Install R >= 4.5 first: https://www.r-project.org/")
        return 1

    print(f"Using Rscript: {rscript}")
    env = os.environ.copy()
    if args.target:
        target = Path(args.target).expanduser().resolve()
        target.mkdir(parents=True, exist_ok=True)
        env["R_LIBS_USER"] = str(target)
        env["PIP_TARGET"] = str(target)
        env["PYTHONPATH"] = str(target) + os.pathsep + env.get("PYTHONPATH", "")
        print(f"Install target: {target}")
    try:
        result = _run_logged(
            [rscript, str(ROOT / "src" / "analysis" / "install_deps.R")],
            R_DEPS_TIMEOUT_SECONDS,
            cwd=ROOT,
            env=env,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except subprocess.TimeoutExpired:
        print("R dependency installation timed out.")
        return 124
    if result.returncode != 0:
        print("R dependency installation failed.")
        return 1

    if args.target:
        try:
            pip = _run_logged(
                [
                    sys.executable,
                    "-m",
                    "pip",
                    "install",
                    "--target",
                    str(Path(args.target).expanduser().resolve()),
                    "-r",
                    str(ROOT / "requirements.txt"),
                ],
                PIP_TIMEOUT_SECONDS,
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
        except subprocess.TimeoutExpired:
            print("Python dependency installation timed out.")
            return 124
        if pip.returncode != 0:
            print("Python dependency installation failed.")
            return 1

    print("Verifying environment...")
    try:
        check = _run_logged(
            [
                sys.executable,
                str(ROOT / "launchers" / "check_pipeline_environment.py"),
            ],
            CHECK_TIMEOUT_SECONDS,
            cwd=ROOT,
            env=env,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except subprocess.TimeoutExpired:
        print("environment check timed out.")
        return 124
    return check.returncode


if __name__ == "__main__":
    sys.exit(main())
