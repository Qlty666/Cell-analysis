#!/usr/bin/env python3
"""Install or verify dependencies for one functional module.

Examples:
  python launchers\\install_environment.py install expression
  python launchers\\install_environment.py install docking --with-ml
  python launchers\\install_environment.py check full
  python launchers\\install_environment.py list
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from common.env import find_rscript  # noqa: E402


MODULES = {
    "expression": {
        "label": "Expression analysis",
        "aliases": ("single", "single-cell", "pipeline"),
        "entrypoints": (
            "scripts/run_pipeline.py",
            "launchers/run_GSE125449.bat",
            "launchers/run_pipeline_prompt.bat",
        ),
        "pip_requirements": ("requirements.txt",),
        "pip_packages": (),
        "r_deps": True,
        "dock_installer": False,
        "skills": False,
        "check": "pipeline",
    },
    "datasets": {
        "label": "Dataset search",
        "aliases": ("dataset", "search"),
        "entrypoints": (
            "scripts/search_datasets.py",
            "scripts/dataset_search_ml.py",
        ),
        "pip_requirements": ("requirements.txt",),
        "pip_packages": ("joblib",),
        "r_deps": False,
        "dock_installer": False,
        "skills": False,
        "check": "datasets",
    },
    "docking": {
        "label": "Virtual screening / docking",
        "aliases": ("dock", "virtual-screening", "cadd"),
        "entrypoints": ("scripts/run_docking.py",),
        "pip_requirements": ("requirements.txt",),
        "pip_packages": ("joblib",),
        "r_deps": False,
        "dock_installer": True,
        "skills": False,
        "check": "dock",
    },
    "molecular-docking": {
        "label": "Standalone molecular docking",
        "aliases": ("molecular_docking", "standalone-docking"),
        "entrypoints": (
            "scripts/run_molecular_docking.py",
            "launchers/run_molecular_docking.bat",
        ),
        "pip_requirements": (),
        "pip_packages": (),
        "r_deps": False,
        "dock_installer": True,
        "skills": False,
        "check": "dock",
    },
    "md": {
        "label": "Molecular dynamics (GROMACS)",
        "aliases": ("molecular-dynamics", "gromacs"),
        "entrypoints": ("scripts/run_docking.py md-simulation",),
        "pip_requirements": (),
        "pip_packages": (),
        "r_deps": False,
        "dock_installer": True,
        "skills": False,
        "check": "md",
    },
    "full": {
        "label": "Full integrated pipeline",
        "aliases": ("full-pipeline", "all", "everything"),
        "entrypoints": (
            "scripts/run_full_pipeline.py",
            "launchers/run_full_pipeline.bat",
        ),
        "pip_requirements": ("requirements.txt",),
        "pip_packages": ("joblib",),
        "r_deps": True,
        "dock_installer": True,
        "skills": True,
        "check": "both",
    },
    "web": {
        "label": "Web console",
        "aliases": ("web-ui", "ui"),
        "entrypoints": (
            "web/web_ui.py",
            "launchers/run_web_ui.bat",
        ),
        "pip_requirements": (),
        "pip_packages": (),
        "r_deps": False,
        "dock_installer": False,
        "skills": False,
        "check": "web",
    },
    "skills": {
        "label": "Project Codex skills",
        "aliases": ("codex-skills",),
        "entrypoints": ("scripts/install_codex_skills.py",),
        "pip_requirements": (),
        "pip_packages": (),
        "r_deps": False,
        "dock_installer": False,
        "skills": True,
        "check": "skills",
    },
}

ALIASES: dict[str, str] = {}
for _name, _spec in MODULES.items():
    ALIASES[_name] = _name
    for _alias in _spec["aliases"]:
        ALIASES[_alias] = _name


def resolve_module(name: str) -> str:
    key = ALIASES.get(str(name).strip().lower())
    if key is None:
        choices = ", ".join(sorted(MODULES))
        raise ValueError(f"unknown module {name!r}; available: {choices}")
    return key


def print_list() -> None:
    print("Functional modules and their install/check names:")
    for name in sorted(MODULES):
        spec = MODULES[name]
        aliases = ", ".join(spec["aliases"])
        print(f"  {name:<20} {spec['label']} (aliases: {aliases})")


def _cmd_text(cmd: list[str]) -> str:
    return " ".join(str(part) for part in cmd)


def _run(cmd: list[str], dry_run: bool, env: dict[str, str] | None = None) -> int:
    print(f"> {_cmd_text(cmd)}")
    if dry_run:
        return 0
    return subprocess.call(cmd, cwd=ROOT, env=env)


def _target_env(target: str) -> dict[str, str]:
    env = os.environ.copy()
    if target:
        target_path = str(Path(target).expanduser().resolve())
        env["PYTHONPATH"] = target_path + os.pathsep + env.get("PYTHONPATH", "")
        env["R_LIBS_USER"] = target_path
    return env


def _pip_install(args: list[str], target: str, dry_run: bool) -> int:
    if target:
        target_path = Path(target).expanduser().resolve()
        target_path.mkdir(parents=True, exist_ok=True)
        args = list(args) + ["--target", str(target_path)]
    return _run(
        [sys.executable, "-m", "pip", "install", *args],
        dry_run,
    )


def _run_project_script(
    relative: str,
    args: list[str],
    target: str,
    dry_run: bool,
) -> int:
    script = ROOT / relative
    if not script.is_file():
        print(f"missing project script: {script}", file=sys.stderr)
        return 1
    return _run(
        [sys.executable, str(script), *args],
        dry_run,
        env=_target_env(target),
    )


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
        return False
    return dest.is_file()


def _bootstrap_r(auto_install: bool, dry_run: bool) -> str | None:
    found = find_rscript()
    if found:
        return found
    if dry_run:
        print("Rscript was not found; dry-run assumes R is available.")
        return "Rscript.exe"
    if not auto_install:
        print(
            "Rscript not found. Install R 4.5+ from https://cran.r-project.org/ "
            "or rerun with R in PATH.",
            file=sys.stderr,
        )
        return None
    if os.name != "nt":
        print(
            "Automatic R installation is only implemented on Windows. "
            "Install R 4.5+ manually.",
            file=sys.stderr,
        )
        return None

    base_url = "https://cran.r-project.org/bin/windows/base/"
    try:
        request = urllib.request.Request(
            base_url,
            headers={"User-Agent": "Mozilla/5.0 liver-cancer-pipeline-env"},
        )
        with urllib.request.urlopen(request, timeout=60) as resp:
            html = resp.read().decode("utf-8", "replace")
    except Exception as exc:
        print(f"could not query CRAN R releases: {exc}", file=sys.stderr)
        return None
    match = re.search(r'href="(R-([0-9.]+)-win\.exe)"', html)
    if not match:
        print("could not find the latest R Windows installer on CRAN", file=sys.stderr)
        return None
    filename = match.group(1)
    version = match.group(2)
    installer_url = base_url + filename
    installers_dir = ROOT / "data_cache" / "installers"
    installer_path = installers_dir / filename
    if not installer_path.is_file():
        if not _download(installer_url, installer_path):
            return None

    local = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    install_dir = local / "Programs" / "R" / f"R-{version}"
    print(f"Installing R {version} to {install_dir}")
    code = _run(
        [
            str(installer_path),
            "/VERYSILENT",
            "/SUPPRESSMSGBOXES",
            "/NORESTART",
            f"/DIR={install_dir}",
        ],
        dry_run=False,
    )
    if code != 0:
        print(f"R installer exited with code {code}", file=sys.stderr)
        return None
    rscript = install_dir / "bin" / "Rscript.exe"
    return str(rscript) if rscript.is_file() else find_rscript()


def _install_r_packages(target: str, auto_install_r: bool, dry_run: bool) -> bool:
    rscript = _bootstrap_r(auto_install_r, dry_run)
    if not rscript:
        return False
    args = ["--target", target] if target else []
    return _run_project_script(
        "launchers/install_pipeline_dependencies.py",
        args,
        target,
        dry_run,
    ) == 0


def _install_dock_tools(target: str, dry_run: bool) -> bool:
    args = ["--target", target] if target else []
    return _run_project_script(
        "launchers/install_dock_dependencies.py",
        args,
        target,
        dry_run,
    ) == 0


def _install_skills(dry_run: bool) -> bool:
    return _run_project_script(
        "scripts/install_codex_skills.py",
        [],
        "",
        dry_run,
    ) == 0


def install_module(
    name: str,
    target: str = "",
    auto_install_r: bool = True,
    with_ml: bool = False,
    dry_run: bool = False,
) -> int:
    key = resolve_module(name)
    spec = MODULES[key]
    print(f"Installing environment for: {spec['label']} ({key})")

    for relative in spec["pip_requirements"]:
        if _pip_install(["-r", str(ROOT / relative)], target, dry_run) != 0:
            return 1
    if spec["pip_packages"] and _pip_install(
        list(spec["pip_packages"]),
        target,
        dry_run,
    ) != 0:
        return 1
    if with_ml and _pip_install(["joblib", "torch"], target, dry_run) != 0:
        return 1
    if spec["r_deps"] and not _install_r_packages(
        target,
        auto_install_r,
        dry_run,
    ):
        return 1
    if spec["dock_installer"] and not _install_dock_tools(target, dry_run):
        return 1
    if spec["skills"] and not _install_skills(dry_run):
        return 1
    if key == "md":
        print("NOTE: GROMACS gmx is not downloaded automatically.")
        print("Install gmx separately, then run the MD module.")
    return check_module(
        key,
        target=target,
        dry_run=dry_run,
        skip_ml=not with_ml,
    )


def _check_python_packages(
    packages: tuple[str, ...],
    target: str,
    dry_run: bool,
) -> int:
    code = (
        "import importlib.util, sys\n"
        "ok = True\n"
        "for name in sys.argv[1:]:\n"
        "    found = importlib.util.find_spec(name) is not None\n"
        "    print(('[OK ] ' if found else '[FAIL] ') + name)\n"
        "    ok = ok and found\n"
        "raise SystemExit(0 if ok else 1)\n"
    )
    return _run(
        [sys.executable, "-c", code, *packages],
        dry_run,
        env=_target_env(target),
    )


def check_module(
    name: str,
    target: str = "",
    dry_run: bool = False,
    skip_ml: bool = False,
) -> int:
    key = resolve_module(name)
    spec = MODULES[key]
    mode = spec["check"]
    if dry_run:
        print(f"Would check environment for: {spec['label']} ({key})")
    if mode == "pipeline":
        return _run_project_script(
            "launchers/check_pipeline_environment.py",
            [],
            target,
            dry_run,
        )
    if mode == "dock":
        args = ["--skip-ml"] if skip_ml else []
        return _run_project_script(
            "launchers/check_dock_environment.py",
            args,
            target,
            dry_run,
        )
    if mode == "md":
        args = ["--skip-ml"] if skip_ml else []
        result = _run_project_script(
            "launchers/check_dock_environment.py",
            args,
            target,
            dry_run,
        )
        print("MD additionally requires gmx and a ligand topology tool.")
        return result
    if mode == "datasets":
        return _check_python_packages(
            ("numpy", "pandas", "sklearn", "joblib"),
            target,
            dry_run,
        )
    if mode == "web":
        if not (ROOT / "web" / "web_ui.py").is_file():
            print("web/web_ui.py is missing", file=sys.stderr)
            return 1
        print(f"Python {sys.version.split()[0]}; web entry point found.")
        return 0
    if mode == "skills":
        missing = [
            relative
            for relative in (
                "skills/liver-expression-analysis",
                "skills/liver-virtual-screening",
                "skills/liver-full-pipeline",
                "skills/liver-dataset-search",
            )
            if not (ROOT / relative).is_dir()
        ]
        for item in missing:
            print(f"[FAIL] {item}", file=sys.stderr)
        print("OK: all project skill sources found" if not missing else "skills missing")
        return 0 if not missing else 1
    if mode == "both":
        first = check_module(
            "expression",
            target=target,
            dry_run=dry_run,
            skip_ml=skip_ml,
        )
        if first != 0:
            return first
        return check_module(
            "docking",
            target=target,
            dry_run=dry_run,
            skip_ml=skip_ml,
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=("install", "check", "list"),
        help="action to run",
    )
    parser.add_argument(
        "module",
        nargs="?",
        help="module name, e.g. expression, docking, md, full",
    )
    parser.add_argument("--target", default="", help="pip/R target directory")
    parser.add_argument(
        "--auto-install-r",
        action="store_true",
        default=True,
        help="download and install R on Windows when Rscript is missing",
    )
    parser.add_argument(
        "--no-auto-install-r",
        action="store_false",
        dest="auto_install_r",
        help="do not download R automatically",
    )
    parser.add_argument(
        "--with-ml",
        action="store_true",
        help="also install ML/DL packages such as torch",
    )
    parser.add_argument(
        "--skip-ml",
        action="store_true",
        help="ignore optional ML packages when checking",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print commands without executing them",
    )
    args = parser.parse_args(argv)

    if args.command == "list":
        print_list()
        return 0
    if not args.module:
        parser.error("module is required for install/check")
    try:
        module = resolve_module(args.module)
    except ValueError as exc:
        parser.error(str(exc))
        return 2
    if args.command == "check":
        return check_module(
            module,
            target=args.target,
            dry_run=args.dry_run,
            skip_ml=args.skip_ml,
        )
    return install_module(
        module,
        target=args.target,
        auto_install_r=args.auto_install_r,
        with_ml=args.with_ml,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    sys.exit(main())
