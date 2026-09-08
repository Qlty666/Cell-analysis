#!/usr/bin/env python3
"""Create a clean source package for installing the suite on a new computer.

The package contains only files tracked by git plus the generated
``README_FIRST.txt`` guide, so local outputs, caches, downloaded tools and
runtime logs are not carried to the new machine.
"""

from __future__ import annotations

import argparse
import datetime
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
VERSION_PATTERN = re.compile(r'__version__\s*=\s*["\']([^"\']+)["\']')
PROJECT_NAME = "Cell-analysis-portable"
EXCLUDE_DIRS = {
    ".git",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    "__pycache__",
    "data_cache",
    "logs",
    "results",
    "dock/tools",
    "dock/outputs",
    "dock/work",
    "molecular_docking",
    "portable",
    "venv",
    ".venv",
}
# Never ship credentials or private keys, even in the non-git fallback path.
SECRET_SUFFIXES = {".pem", ".p12", ".key", ".pfx"}
SECRET_NAMES = {".env"}


def project_version() -> str:
    init_path = ROOT / "src" / "docking" / "__init__.py"
    try:
        text = init_path.read_text(encoding="utf-8")
    except OSError:
        return "unknown"
    match = VERSION_PATTERN.search(text)
    return match.group(1) if match else "unknown"


def source_files() -> list[Path]:
    """Return clean source files using git when the repository is present."""
    if (ROOT / ".git").exists():
        proc = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=60,
        )
        if proc.returncode == 0:
            paths = [
                ROOT / relative
                for relative in proc.stdout.split("\0")
                if relative and (ROOT / relative).is_file()
            ]
            return sorted(paths)

    known_root_files = [
        ".gitattributes",
        ".gitignore",
        "AGENTS.md",
        "LICENSE",
        "NEW_COMPUTER_SETUP.md",
        "README.md",
        "VIRTUAL_SCREENING_REQUIREMENTS.md",
        "check_new_computer.bat",
        "check_new_computer.sh",
        "environment_dock.yml",
        "liverbio.bat",
        "package_for_new_computer.bat",
        "package_for_new_computer.sh",
        "requirements.txt",
        "requirements_dock.txt",
        "setup_new_computer.bat",
        "setup_new_computer.sh",
    ]
    paths = [ROOT / name for name in known_root_files if (ROOT / name).is_file()]
    for directory in (
        "config",
        "docs",
        "launchers",
        "scripts",
        "skills",
        "src",
        "tests",
        "web",
    ):
        for path in sorted((ROOT / directory).rglob("*")):
            if path.is_file() and not _excluded(path):
                paths.append(path)
    return sorted(paths)


def _excluded(path: Path) -> bool:
    relative = path.relative_to(ROOT).as_posix()
    parts = relative.split("/")
    if set(parts) & {
        "__pycache__",
        ".pytest_cache",
        ".mypy_cache",
        ".ruff_cache",
    }:
        return True
    if any(part.endswith(".egg-info") for part in parts):
        return True
    if any(
        relative == item or relative.startswith(item + "/")
        for item in EXCLUDE_DIRS
    ):
        return True
    if path.suffix.lower() in {".pyc", ".pyo", ".log"}:
        return True
    if path.suffix.lower() in SECRET_SUFFIXES:
        return True
    name = path.name.lower()
    if name in SECRET_NAMES or any(
        name.startswith(secret + ".") for secret in SECRET_NAMES
    ):
        return True
    return False


def first_readme_text() -> str:
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"""This is a clean source package generated on {now}.

Quick start on Windows:
  1. Extract the archive anywhere.
  2. Run setup_new_computer.bat (installs Python packages, R when missing,
     docking tools, project skills and optional ML packages).
  3. Run check_new_computer.bat to verify.
  4. Use liverbio.bat or python scripts\\liverbio.py.

Linux/macOS quick start:
  python3 launchers/install_environment.py install full --with-ml
  python3 launchers/install_environment.py check full

More details: NEW_COMPUTER_SETUP.md and README.md
"""


def build_package(
    output: Path,
    dry_run: bool = False,
    force: bool = False,
) -> Path:
    files = source_files()
    if not files:
        raise SystemExit("no source files found to package")
    output.parent.mkdir(parents=True, exist_ok=True)
    archive_path = output.resolve()
    if archive_path.exists() and not force and not dry_run:
        raise SystemExit(
            f"output archive already exists: {archive_path}\n"
            "Pass --force to overwrite it."
        )
    print(f"Packaging {len(files)} tracked/source files -> {archive_path}")
    if dry_run:
        return archive_path

    with zipfile.ZipFile(
        archive_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
    ) as archive:
        for path in files:
            archive.write(path, path.relative_to(ROOT).as_posix())
        archive.writestr("README_FIRST.txt", first_readme_text())

    size_mb = archive_path.stat().st_size / 1024 / 1024
    print(f"Package ready: {archive_path} ({size_mb:.2f} MB)")
    return archive_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default=str(ROOT / "portable" / f"{PROJECT_NAME}_v{project_version()}.zip"),
        help="output zip path",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the planned package without writing it",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite the output zip when it already exists",
    )
    args = parser.parse_args(argv)
    build_package(
        Path(args.output),
        dry_run=args.dry_run,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
