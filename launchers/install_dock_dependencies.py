#!/usr/bin/env python3
"""Install the Python dependencies used by the docking pipeline."""

import subprocess
import sys
import zipfile
import argparse
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIREMENTS = ROOT / "requirements_dock.txt"
TARGET = ""

# Pinned AutoDockTools_py3 release (GitHub tag, resolved 2026-06).
ADT_REPO = "Valdes-Tresanco-MS/AutoDockTools_py3"
ADT_TAG = "1.5.7.post1"
ADT_ZIP_URL = (
    f"https://codeload.github.com/{ADT_REPO}/zip/refs/tags/{ADT_TAG}"
)
# Pinned AutoDock Vina release. GitHub publishes no checksum for these assets,
# so the exact release URL is the integrity control here.
VINA_VERSION = "1.2.7"
VINA_FALLBACK_URL = (
    "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/"
    f"v{VINA_VERSION}/vina_{VINA_VERSION}_win.exe"
)
PIP_TIMEOUT_SECONDS = 3600
DOWNLOAD_TIMEOUT_SECONDS = 300


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


def _cmd_text(cmd: list[str]) -> str:
    return " ".join(str(part) for part in cmd)


def _run_pip(args: list[str], timeout: int = PIP_TIMEOUT_SECONDS) -> int:
    cmd = [sys.executable, "-m", "pip", *args]
    if TARGET:
        cmd += ["--target", str(Path(TARGET).expanduser().resolve())]
    print(f"> {_cmd_text(cmd)} (timeout {timeout}s)")
    try:
        return subprocess.run(cmd, cwd=ROOT, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        print(
            f"pip timed out after {timeout}s: {_cmd_text(cmd)}",
            file=sys.stderr,
        )
        return 124


def _extracted_adt_dir(tools: Path) -> Path | None:
    preferred = tools / f"AutoDockTools_py3-{ADT_TAG}"
    if preferred.is_dir():
        return preferred
    for candidate in sorted(tools.glob("AutoDockTools_py3-*")):
        if candidate.is_dir():
            return candidate
    return None


def _ensure_autodocktools() -> bool:
    tools = ROOT / "dock" / "tools"
    src = tools / "AutoDockTools_py3"
    zip_path = tools / "autodocktools_py3.zip"
    ref_marker = tools / "autodocktools_py3.ref"
    if not src.exists():
        cached_ref = (
            ref_marker.read_text(encoding="utf-8").strip()
            if ref_marker.is_file()
            else ""
        )
        if zip_path.exists() and cached_ref != ADT_TAG:
            print(
                "discarding cached AutoDockTools_py3 archive "
                f"(ref {cached_ref or 'unknown'} != {ADT_TAG})"
            )
            zip_path.unlink()
        if not zip_path.exists():
            print(
                "AutoDockTools_py3 source and zip missing; "
                f"downloading tag {ADT_TAG}..."
            )
            if not _download(ADT_ZIP_URL, zip_path):
                return False
            ref_marker.write_text(ADT_TAG, encoding="utf-8")
        print("Extracting AutoDockTools_py3...")
        with zipfile.ZipFile(zip_path) as archive:
            tools_resolved = tools.resolve()
            for member in archive.infolist():
                member_name = member.filename.replace("\\", "/")
                target = (tools_resolved / member_name).resolve()
                if not target.is_relative_to(tools_resolved):
                    raise RuntimeError(f"unsafe zip member: {member.filename}")
            archive.extractall(tools)
        extracted = _extracted_adt_dir(tools)
        if extracted is not None and not src.exists():
            extracted.rename(src)
    if not src.exists():
        print("AutoDockTools_py3 source still missing after extraction")
        return False
    print("Installing AutoDockTools_py3 (prepare_receptor4)...")
    return _run_pip(["install", "versioneer"]) == 0 and _run_pip(
        ["install", str(src)]
    ) == 0


def _download(url: str, dest: Path) -> bool:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 liver-cancer-pipeline-env"},
    )
    try:
        with urllib.request.urlopen(
            request, timeout=DOWNLOAD_TIMEOUT_SECONDS
        ) as resp, dest.open("wb") as handle:
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
    if not dest.is_file() or dest.stat().st_size == 0:
        print(f"download produced no data: {dest}", file=sys.stderr)
        if dest.exists():
            dest.unlink()
        return False
    return True


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
    print(
        f"AutoDock Vina missing; downloading pinned Windows binary "
        f"v{VINA_VERSION}..."
    )
    if not _download(VINA_FALLBACK_URL, tmp):
        print(
            "could not download AutoDock Vina. Download it manually from "
            "https://github.com/ccsb-scripps/AutoDock-Vina/releases and place "
            f"it at {target}.",
            file=sys.stderr,
        )
        return False
    tmp.replace(target)
    print(f"AutoDock Vina: {target}")
    return True


if __name__ == "__main__":
    sys.exit(main())
