#!/usr/bin/env python3
"""Install all R dependencies required by the pipeline."""

import shutil
import subprocess
import sys
import argparse
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def find_rscript() -> str | None:
    if shutil.which("Rscript"):
        return shutil.which("Rscript")
    if shutil.which("Rscript.exe"):
        return shutil.which("Rscript.exe")
    for base in [
        Path(r"C:\Program Files\R"),
        Path(r"C:\Program Files\Microsoft\R Open"),
    ]:
        if base.exists():
            candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
            if candidates:
                return str(candidates[0])
    return None


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
    result = subprocess.run(
        [rscript, str(ROOT / "src" / "analysis" / "install_deps.R")],
        cwd=ROOT,
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        print("R dependency installation failed.")
        return 1

    if args.target:
        pip = subprocess.run(
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
            cwd=ROOT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if pip.returncode != 0:
            print("Python dependency installation failed.")
            return 1

    print("Verifying environment...")
    check = subprocess.run(
        [sys.executable, str(ROOT / "launchers" / "check_pipeline_environment.py")],
        cwd=ROOT,
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return check.returncode


if __name__ == "__main__":
    sys.exit(main())
