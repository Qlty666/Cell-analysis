"""Export completed LiverBio runs into a Codex analysis workspace.

The analysis workspace keeps third-party pipeline outputs under
``data/imported_results`` and registers them with
``scripts/identify_imported_results.ps1``.  This module implements that
handoff for runs produced by this project.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCAL_ANALYSIS_CONFIG = ROOT / "config" / "analysis_workspace.json"

SKIP_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".git",
}


def project_version() -> str:
    """Read the suite version from the docking package without importing."""
    init_path = ROOT / "src" / "docking" / "__init__.py"
    try:
        text = init_path.read_text(encoding="utf-8")
    except OSError:
        return "unknown"
    match = re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', text)
    return match.group(1) if match else "unknown"


def git_revision() -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        value = result.stdout.strip()
        return value or "unknown"
    except Exception:
        return "unknown"


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _string_value(data: dict, *keys: str) -> str:
    for key in keys:
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()
    return ""


def detect_accession(source: Path) -> str:
    """Return the dataset accession recorded in a LiverBio run directory."""
    candidates: list[tuple[Path, tuple[str, ...]]] = []

    integration = source / "outputs" / "integration"
    candidates.append(
        (integration / "integration_summary.json", ("dataset", "accession"))
    )
    candidates.append(
        (integration / "run_manifest.json", ("accession", "dataset"))
    )
    candidates.append(
        (source / "outputs" / "results" / "summary.json", ("dataset", "accession"))
    )
    candidates.append(
        (source / "results" / "summary.json", ("dataset", "accession"))
    )

    for path, keys in candidates:
        accession = _string_value(_read_json(path), *keys)
        if accession:
            return accession.upper()

    nested = _read_json(integration / "run_manifest.json")
    single_cell = nested.get("parameters", {}).get("single_cell", {})
    if isinstance(single_cell, dict):
        accession = _string_value(single_cell, "dataset", "accession")
        if accession:
            return accession.upper()

    summary = _read_json(source / "outputs" / "integration" / "integration_summary.json")
    single_cell = summary.get("single_cell", {})
    if isinstance(single_cell, dict):
        accession = _string_value(single_cell, "dataset", "accession")
        if accession:
            return accession.upper()

    for folder in ("data", "outputs" / "data"):
        manifest_dir = source / folder
        if manifest_dir.is_dir():
            for manifest in sorted(manifest_dir.glob("*_manifest.json")):
                accession = _string_value(_read_json(manifest), "accession", "dataset")
                if accession:
                    return accession.upper()
    return ""


def detect_run_kind(source: Path) -> str:
    """Distinguish full-pipeline, expression-only and generic result roots."""
    integration = source / "outputs" / "integration"
    if integration.is_dir() and any(integration.rglob("*.done")):
        return "full_pipeline"
    if integration.is_dir():
        return "integration"
    if (source / "outputs" / "results" / "pipeline_complete.json").exists():
        return "expression"
    if (source / "results" / "pipeline_complete.json").exists():
        return "expression"
    return "other"


def _as_list(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, dict):
        return [value]
    return []


def _output_root_infos(analysis_root: Path) -> list[dict]:
    """Collect output roots from analysis-workspace registry files."""
    infos: list[dict] = []
    registry = analysis_root / "config" / "local_projects.json"
    if registry.exists():
        data = _read_json(registry)
        for project in _as_list(data.get("projects") or data.get("project")):
            infos.extend(_as_list(project.get("output_roots")))
        infos.extend(_as_list(data.get("output_roots")))

    config = resolve_analysis_config()
    infos.extend(_as_list(config.get("output_roots")))

    env_root = os.environ.get("LIVER_OUTPUT_ROOT")
    if env_root:
        infos.append({"path": env_root, "status": "environment"})
    return infos


def _source_score(path: Path) -> int:
    score = 0
    if (path / "outputs" / "integration").is_dir():
        score += 10
    if (path / "outputs" / "integration" / "integration_summary.json").exists():
        score += 5
    if (path / "outputs" / "results" / "pipeline_complete.json").exists():
        score += 3
    if (path / "results" / "pipeline_complete.json").exists():
        score += 2
    return score


def locate_run_source(accession: str, analysis_root: Path) -> Path | None:
    """Find a run directory for an accession under known output roots."""
    accession = accession.upper()
    candidates: list[tuple[int, float, Path]] = []
    for info in _output_root_infos(analysis_root):
        root_value = info.get("path") if isinstance(info, dict) else None
        if not root_value:
            continue
        root = Path(str(root_value)).expanduser()
        if not root.is_dir():
            continue
        accessions = info.get("accessions") if isinstance(info, dict) else []
        matches = (
            accession in {str(item).upper() for item in _as_list(accessions)}
            if isinstance(accessions, (list, tuple, set))
            else False
        )
        subdir = root / accession
        possible = [subdir] if subdir.is_dir() else []
        if not possible and matches:
            possible.append(root)
        for path in possible:
            if path.is_dir() and _source_score(path) > 0:
                age = path.stat().st_mtime
                candidates.append((_source_score(path), age, path))

    if not candidates:
        return None
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return candidates[0][2]


def resolve_analysis_config() -> dict:
    data: dict = {}
    if LOCAL_ANALYSIS_CONFIG.exists():
        data = _read_json(LOCAL_ANALYSIS_CONFIG)
    try:
        env = json.loads(os.environ.get("LIVER_ANALYSIS_CONFIG", "{}"))
        if isinstance(env, dict):
            data.update(env)
    except Exception:
        pass
    return data


def _sync_file(source: Path, target: Path, stats: dict) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    source_stat = source.stat()
    if target.exists():
        target_stat = target.stat()
        same_size = target_stat.st_size == source_stat.st_size
        same_mtime = abs(target_stat.st_mtime - source_stat.st_mtime) < 2.0
        if same_size and same_mtime:
            stats["skipped"] += 1
            return
    shutil.copy2(source, target)
    stats["copied"] += 1


def _sync_tree(source: Path, target: Path, stats: dict) -> None:
    if not target.exists():
        target.mkdir(parents=True, exist_ok=True)
    for item in sorted(source.iterdir()):
        if item.name in SKIP_NAMES:
            continue
        dest = target / item.name
        if item.is_dir():
            _sync_tree(item, dest, stats)
        elif item.is_file():
            _sync_file(item, dest, stats)


def write_source_metadata(
    source: Path,
    destination: Path,
    analysis_root: Path,
    accession: str,
    kind: str,
    copied: int,
    skipped: int,
) -> Path:
    metadata = {
        "exporter": "liverbio-analysis-export",
        "source_root": str(source.resolve()),
        "analysis_root": str(analysis_root.resolve()),
        "destination": str(destination.resolve()),
        "dataset": accession,
        "run_kind": kind,
        "files_copied": copied,
        "files_skipped": skipped,
        "project_version": project_version(),
        "git_revision": git_revision(),
        "exported_at": datetime.now().isoformat(timespec="seconds"),
    }
    out_path = destination / "_source.json"
    out_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return out_path


def refresh_inventory(
    analysis_root: Path,
    inventory_script: str | None,
) -> tuple[int, str]:
    script = Path(inventory_script or "").resolve()
    if not inventory_script:
        default = analysis_root / "scripts" / "identify_imported_results.ps1"
        script = default if default.exists() else Path()
    if not script.exists():
        return 2, f"inventory script not found: {script}"
    shell = shutil.which("pwsh.exe") or shutil.which("powershell.exe")
    if not shell:
        return 2, "PowerShell executable not found on PATH"
    command = [
        shell,
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-File",
        str(script),
        "-Folder",
        str(analysis_root / "data" / "imported_results"),
    ]
    result = subprocess.run(
        command,
        cwd=str(analysis_root),
        capture_output=True,
        text=True,
        check=False,
    )
    output = (result.stdout or "") + (result.stderr or "")
    return result.returncode, output.strip()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="liverbio analysis-export",
        description=(
            "Export a finished run into a Codex analysis workspace at "
            "data/imported_results and refresh its inventory"
        ),
    )
    parser.add_argument(
        "run",
        nargs="?",
        help="finished run root or dataset accession (for example GSE235863)",
    )
    parser.add_argument("--source", help="finished run root")
    parser.add_argument("--analysis-root", help="analysis workspace root")
    parser.add_argument(
        "--remember-analysis-root",
        action="store_true",
        help="save --analysis-root into config/analysis_workspace.json",
    )
    parser.add_argument("--name", help="destination folder name")
    parser.add_argument(
        "--inventory-script",
        help="identify_imported_results.ps1 path (default under analysis root)",
    )
    parser.add_argument(
        "--no-inventory",
        action="store_true",
        help="copy files without regenerating _inventory.json",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show what would be copied without writing files",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    run_value = (args.source or args.run or "").strip()
    if not run_value:
        print(
            "a run directory or dataset accession is required "
            "(for example: analysis-export GSE235863)",
            file=sys.stderr,
        )
        return 2

    config = resolve_analysis_config()
    analysis_root_value = (
        args.analysis_root
        or os.environ.get("LIVER_ANALYSIS_ROOT")
        or config.get("analysis_root")
    )
    if not analysis_root_value:
        print(
            "analysis root is required; pass --analysis-root, set "
            "LIVER_ANALYSIS_ROOT, or write config/analysis_workspace.json",
            file=sys.stderr,
        )
        return 2
    analysis_root = Path(analysis_root_value).expanduser().resolve()
    if not analysis_root.is_dir():
        print(f"analysis workspace not found: {analysis_root}", file=sys.stderr)
        return 1

    if args.remember_analysis_root:
        saved = {
            "analysis_root": str(analysis_root),
        }
        if args.inventory_script:
            saved["inventory_script"] = args.inventory_script
        LOCAL_ANALYSIS_CONFIG.parent.mkdir(parents=True, exist_ok=True)
        LOCAL_ANALYSIS_CONFIG.write_text(
            json.dumps(saved, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        print("remembered analysis root: %s" % LOCAL_ANALYSIS_CONFIG)

    source_candidate = Path(run_value)
    if source_candidate.is_dir():
        source = source_candidate.expanduser().resolve()
    else:
        source = locate_run_source(run_value, analysis_root)
        if source is None:
            print(
                "could not locate a run for %s; pass an existing run "
                "directory or register its output root" % run_value,
                file=sys.stderr,
            )
            return 2
    inventory_script = args.inventory_script or config.get("inventory_script")

    accession = detect_accession(source)
    fallback = run_value.upper() if not Path(run_value).is_dir() else accession
    name = (args.name or accession or fallback or source.name).strip().replace("/", "_")
    destination = analysis_root / "data" / "imported_results" / name
    kind = detect_run_kind(source)

    stats = {"copied": 0, "skipped": 0}
    if args.dry_run:
        for item in sorted(source.rglob("*")):
            if any(part in SKIP_NAMES for part in item.parts):
                continue
            if item.is_file():
                rel = item.relative_to(source)
                target = destination / rel
                if target.exists() and target.stat().st_size == item.stat().st_size:
                    stats["skipped"] += 1
                else:
                    stats["copied"] += 1
        print(
            "dry run: source=%s destination=%s kind=%s accession=%s "
            "files_to_copy=%s files_unchanged=%s"
            % (
                source,
                destination,
                kind,
                accession or name,
                stats["copied"],
                stats["skipped"],
            )
        )
        return 0

    _sync_tree(source, destination, stats)
    metadata_path = write_source_metadata(
        source,
        destination,
        analysis_root,
        accession or name,
        kind,
        stats["copied"],
        stats["skipped"],
    )
    print(
        "exported %s -> %s (%s files copied, %s unchanged)"
        % (source, destination, stats["copied"], stats["skipped"])
    )
    print("source metadata: %s" % metadata_path)

    if args.no_inventory:
        return 0
    code, output = refresh_inventory(analysis_root, inventory_script)
    if output:
        print(output)
    if code != 0:
        print("inventory refresh failed with exit code %s" % code, file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
