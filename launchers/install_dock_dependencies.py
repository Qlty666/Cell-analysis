#!/usr/bin/env python3
"""Install the Python dependencies used by the docking pipeline."""

import subprocess
import sys
import zipfile
import argparse
import json
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIREMENTS = ROOT / "requirements_dock.txt"
TARGET = ""
ADT_ZIP_URL = (
    "https://codeload.github.com/Valdes-Tresanco-MS/AutoDockTools_py3/"
    "zip/refs/heads/master"
)
VINA_FALLBACK_URL = (
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/"
    "v1.2.7/vina_1.2.7_win.exe"
)


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
    if not _ensure_vina():
        return 1
    print()
    print("Python dependencies ready.")
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
            print("AutoDockTools_py3 source and zip missing; downloading...")
            if not _download(ADT_ZIP_URL, zip_path):
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


def _download(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 liver-cancer-pipeline-env"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as resp, dest.open(
            "wb"
        ) as handle:
            while True:
                chunk = resp.read(1024 * 256)
                if not chunk:
                    break
                handle.write(chunk)
    except Exception as exc:
        print(f"download failed: {exc}", file=sys.stderr)
        if dest.exists():
            dest.unlink()
        return False
    return dest.is_file()


def _latest_vina_win_url() -> str:
    api = "https://api.github.com/repos/ccsb-scripps/AutoDock-Vina/releases/latest"
    request = urllib.request.Request(
        api,
        headers={"User-Agent": "Mozilla/5.0 liver-cancer-pipeline-env"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as resp:
            payload = json.loads(resp.read().decode("utf-8", "replace"))
        for asset in payload.get("assets", []):
            name = asset.get("name", "")
            if (
                name.endswith("_win.exe")
                and "_split_" not in name
                and name.startswith("vina_")
            ):
                return asset.get("browser_download_url") or VINA_FALLBACK_URL
    except Exception as exc:
        print(f"could not query Vina release metadata: {exc}", file=sys.stderr)
    return VINA_FALLBACK_URL


def _ensure_vina() -> bool:
    tools = ROOT / "dock" / "tools"
    target = tools / "vina.exe"
    if target.is_file():
        print(f"AutoDock Vina: {target}")
        return True
    if shutil.which("vina") or shutil.which("vina.exe"):
        print("AutoDock Vina found in PATH.")
        return True
    if not sys.platform.startswith("win"):
        print(
            "AutoDock Vina is not present. Install autodock-vina from conda-forge.",
            file=sys.stderr,
        )
        return False
    tmp = target.with_suffix(".exe.download")
    print("AutoDock Vina missing; downloading Windows binary...")
    if not _download(_latest_vina_win_url(), tmp):
        return False
    tmp.replace(target)
    print(f"AutoDock Vina: {target}")
    return True


if __name__ == "__main__":
    sys.exit(main())
