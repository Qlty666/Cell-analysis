#!/usr/bin/env python3
"""Install the Python dependencies used by the docking pipeline."""

import subprocess
import sys
import zipfile
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIREMENTS = ROOT / "requirements_dock.txt"
TARGET = ""


def main() -> int:
    global TARGET
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="")
    args = parser.parse_args()
    TARGET = args.target
    if TARGET:
        Path(TARGET).expanduser().resolve().mkdir(parents=True, exist_ok=True)
        print(f"Install target: {TARGET}")
    print("Installing docking Python dependencies (RDKit, Meeko, ...)")
    if _run_pip(["install", "-r", str(REQUIREMENTS)]) != 0:
        return 1
    if not _ensure_autodocktools():
        return 1
    print()
    print("Python dependencies ready.")
    print("AutoDock Vina should be placed in dock/tools/vina.exe or added to PATH.")
    print("Run: python scripts\\run_docking.py check-env")
    return 0


def _run_pip(args: list[str]) -> int:
    cmd = [sys.executable, "-m", "pip", *args]
    if TARGET:
        cmd += ["--target", str(Path(TARGET).expanduser().resolve())]
    return subprocess.run(cmd, cwd=ROOT).returncode


def _ensure_autodocktools() -> bool:
    tools = ROOT / "dock" / "tools"
    src = tools / "AutoDockTools_py3"
    zip_path = tools / "autodocktools_py3.zip"
    if not src.exists():
        if not zip_path.exists():
            print(
                "AutoDockTools_py3 source and zip both missing; "
                "prepare_receptor4 unavailable"
            )
            return False
        print("Extracting AutoDockTools_py3...")
        with zipfile.ZipFile(zip_path) as archive:
            tools_resolved = tools.resolve()
            for member in archive.infolist():
                member_name = member.filename.replace("\\", "/")
                target = (tools_resolved / member_name).resolve()
                if not target.is_relative_to(tools_resolved):
                    raise RuntimeError(f"unsafe zip member: {member.filename}")
            archive.extractall(tools)
        extracted = tools / "AutoDockTools_py3-master"
        if extracted.exists() and not src.exists():
            extracted.rename(src)
    if not src.exists():
        print("AutoDockTools_py3 source still missing after extraction")
        return False
    print("Installing AutoDockTools_py3 (prepare_receptor4)...")
    return _run_pip(["install", "versioneer"]) == 0 and _run_pip(
        ["install", str(src)]
    ) == 0


if __name__ == "__main__":
    sys.exit(main())
