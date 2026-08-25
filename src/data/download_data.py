"""Data acquisition for the GSE125449 liver cancer single-cell dataset."""

import re
from pathlib import Path
import shutil
import subprocess

try:
    from .geo_downloader import canonical_accession, ensure_geo_dataset
except ImportError:
    from geo_downloader import canonical_accession, ensure_geo_dataset


def find_curl() -> str:
    found = shutil.which("curl.exe") or shutil.which("curl")
    if found:
        return found
    raise RuntimeError("curl not found")


def ensure_data(cfg: dict, root: Path, log) -> None:
    raw_dir = root / cfg["raw_dir"]
    raw_dir.mkdir(parents=True, exist_ok=True)

    for item in cfg["downloads"]:
        out = raw_dir / item["file"]
        min_bytes = int(item.get("min_bytes", 0))
        if out.exists() and out.stat().st_size >= min_bytes:
            continue
        if out.exists():
            out.unlink()
        log(f"downloading {item['file']}")
        cmd = [
            find_curl(),
            "-L",
            "--ssl-no-revoke",
            "--retry", "5",
            "--retry-delay", "3",
            "--fail",
            "--silent",
            "--show-error",
            "--max-time", "600",
            "-o", str(out),
            item["url"],
        ]
        subprocess.run(cmd, check=True)
        if not out.exists() or out.stat().st_size < min_bytes:
            raise RuntimeError(f"Downloaded file too small or missing: {item['file']}")

    log("data files ready")


def ensure_data_for_accession(accession: str, cfg: dict, root: Path, log) -> None:
    acc = canonical_accession(accession)
    if acc == cfg.get("dataset", "").upper():
        ensure_data(cfg, root, log)
    elif re.fullmatch(r"GSE\d+", acc):
        ensure_geo_dataset(acc, root, log)
    elif re.fullmatch(r"(?:E-[A-Z0-9]+-\d+|S-BSST\d+)", acc):
        try:
            from .biostudies_downloader import ensure_biostudies_dataset
        except ImportError:
            from biostudies_downloader import ensure_biostudies_dataset
        ensure_biostudies_dataset(acc, root, log)
    else:
        raise RuntimeError(
            "dataset accession must look like GSE125449, E-MTAB-1234, "
            "or S-BSST123"
        )
