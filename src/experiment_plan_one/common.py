"""Shared I/O, metadata and plotting helpers for experiment plan one."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import logging
import math
import os
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

LOG = logging.getLogger("experiment_plan_one")

USER_AGENT = "liver-cancer-bioinformatics/1.6 experiment-plan-one"
GEO_RAW_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"

# Okabe-Ito subset with at least 3:1 contrast on white. The palette is used
# only as a starting point; every figure still needs direct labels or another
# redundant cue because color alone is not an accessible encoding.
FIGURE_PALETTE = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#000000",
)
STATUS_COLORS = {
    "completed": "#009E73",
    "cached": "#0072B2",
    "empty": "#D55E00",
    "no_results": "#D55E00",
    "unavailable": "#7A7A7A",
    "disabled": "#7A7A7A",
    "failed": "#B22222",
    "timeout": "#B22222",
    "not_configured": "#7A7A7A",
}


def configure_logging(log_path: Path | None = None, verbose: bool = False) -> None:
    """Configure console and optional file logging once."""
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(
            level=level,
            format="%(asctime)s [%(levelname)s] %(message)s",
        )
    root.setLevel(level)
    logging.getLogger("fontTools.subset").setLevel(logging.WARNING)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        targets = {Path(handler.baseFilename) for handler in root.handlers if hasattr(handler, "baseFilename")}
        if log_path.resolve() not in targets:
            handler = logging.FileHandler(log_path, encoding="utf-8")
            handler.setFormatter(
                logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
            )
            root.addHandler(handler)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _content_length(response: Any) -> int | None:
    value = response.headers.get("Content-Length")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def download_file(
    url: str,
    destination: Path,
    *,
    expected_size: int | None = None,
    timeout: int = 120,
    retries: int = 4,
) -> Path:
    """Download a URL with progress logging and optional byte-size validation."""
    destination = destination.resolve()
    ensure_dir(destination.parent)
    if destination.exists() and destination.stat().st_size > 0:
        if expected_size is None or destination.stat().st_size == expected_size:
            LOG.info("using cached %s", destination)
            return destination
        LOG.warning(
            "cached file has unexpected size: %s (%s != %s)",
            destination,
            destination.stat().st_size,
            expected_size,
        )

    partial = destination.with_suffix(destination.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        start = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": USER_AGENT}
        if start:
            headers["Range"] = f"bytes={start}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", 200)
                if start and status != 206:
                    LOG.warning("server did not honor resume for %s; restarting", url)
                    partial.unlink(missing_ok=True)
                    start = 0
                total = _content_length(response)
                if total is not None and status == 206:
                    total += start
                mode = "ab" if start else "wb"
                next_log = 0
                with partial.open(mode) as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                        start += len(chunk)
                        if start >= next_log:
                            if total:
                                LOG.info(
                                    "download %s: %.1f%% (%s/%s bytes)",
                                    destination.name,
                                    100.0 * start / total,
                                    start,
                                    total,
                                )
                            else:
                                LOG.info("download %s: %s bytes", destination.name, start)
                            next_log = start + (
                                max(total // 10, 25 * 1024 * 1024)
                                if total
                                else 25 * 1024 * 1024
                            )
            if expected_size is not None and partial.stat().st_size != expected_size:
                raise ValueError(
                    f"incomplete download {destination.name}: "
                    f"{partial.stat().st_size} != {expected_size}"
                )
            os.replace(partial, destination)
            LOG.info("downloaded %s", destination)
            return destination
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if (
                isinstance(exc, urllib.error.HTTPError)
                and exc.code == 416
                and start
                and attempt < retries
            ):
                LOG.warning("resume offset is beyond server range; restarting download")
                partial.unlink(missing_ok=True)
                continue
            LOG.warning(
                "download attempt %s/%s failed for %s: %s",
                attempt,
                retries,
                url,
                exc,
            )
            if attempt < retries:
                time.sleep(min(30, 2**attempt))
    raise RuntimeError(f"failed to download {url}: {last_error}")


def get_json(url: str, timeout: int = 60, retries: int = 3) -> Any:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"request failed: {url}: {last_error}")


def post_form(url: str, data: dict[str, Any], timeout: int = 120) -> tuple[str, str]:
    encoded = urllib.parse.urlencode(data).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=encoded,
        headers={
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.geturl(), response.read().decode("utf-8", "replace")


def _safe_tar_members(archive: tarfile.TarFile, destination: Path) -> Iterable[tarfile.TarInfo]:
    root = destination.resolve()
    for member in archive.getmembers():
        target = (destination / member.name).resolve()
        if root != target and root not in target.parents:
            raise RuntimeError(f"unsafe archive member: {member.name}")
        if member.issym() or member.islnk():
            link_target = (target.parent / member.linkname).resolve()
            if root != link_target and root not in link_target.parents:
                raise RuntimeError(f"unsafe archive link: {member.name} -> {member.linkname}")
        yield member


def extract_tar(archive_path: Path, destination: Path) -> Path:
    ensure_dir(destination)
    marker = destination / ".extract_complete"
    if marker.exists():
        LOG.info("using extracted archive %s", destination)
        return destination
    LOG.info("extracting %s -> %s", archive_path.name, destination)
    with tarfile.open(archive_path, "r:*") as archive:
        archive.extractall(destination, members=_safe_tar_members(archive, destination))
    marker.write_text("ok\n", encoding="utf-8")
    return destination


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Parse a GEO series matrix into genes x samples and sample metadata."""
    sample_fields: dict[str, list[list[str]]] = {}
    series_fields: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
                continue
            if not line.startswith("!"):
                continue
            cells = next(csv.reader([line], delimiter="\t", quotechar='"'))
            key = cells[0].lstrip("!")
            values = [cell.strip() for cell in cells[1:]]
            payload = [value for value in values if value != ""]
            if key.startswith("Sample_"):
                sample_fields.setdefault(key, []).append(payload)
            else:
                series_fields.setdefault(key, []).extend(payload)

    if not table_lines:
        raise ValueError(f"series matrix table not found: {path}")
    matrix = pd.read_csv(
        io.StringIO("\n".join(table_lines)),
        sep="\t",
        index_col=0,
        low_memory=False,
    )
    sample_ids = (sample_fields.get("Sample_geo_accession") or [[]])[0]
    if not sample_ids:
        sample_ids = [str(column) for column in matrix.columns]
    if len(sample_ids) != matrix.shape[1]:
        raise ValueError(
            f"sample count mismatch in {path}: {len(sample_ids)} metadata vs "
            f"{matrix.shape[1]} matrix columns"
        )

    metadata_rows: list[dict[str, Any]] = []
    for index, sample_id in enumerate(sample_ids):
        row: dict[str, Any] = {"sample_id": sample_id}
        for key, occurrences in sample_fields.items():
            if key == "Sample_geo_accession":
                continue
            values = [
                occurrence[index]
                for occurrence in occurrences
                if index < len(occurrence) and occurrence[index] != ""
            ]
            row[key] = "; ".join(values)
            if key == "Sample_characteristics_ch1":
                for value in values:
                    if ":" in value:
                        feature, feature_value = value.split(":", 1)
                        normalized = _slug(feature)
                        row[normalized] = feature_value.strip()
        metadata_rows.append(row)
    metadata = pd.DataFrame(metadata_rows).set_index("sample_id")
    matrix.columns = sample_ids
    return matrix, metadata, {"series": series_fields}


def _slug(value: str) -> str:
    text = value.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "field"


def parse_platform_soft(path: Path) -> dict[str, str]:
    """Return probe ID -> HGNC symbol from a GEO platform SOFT archive."""
    mapping: dict[str, str] = {}
    in_table = False
    header: list[str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
                header = None
                continue
            if line == "!platform_table_end":
                in_table = False
                continue
            if not in_table:
                continue
            cells = next(csv.reader([line], delimiter="\t", quotechar='"'))
            if header is None:
                header = cells
                continue
            if len(cells) < len(header):
                cells.extend([""] * (len(header) - len(cells)))
            row = dict(zip(header, cells))
            probe = (
                row.get("ID")
                or row.get("Probe Set ID")
                or row.get("Probe_Set_ID")
                or ""
            ).strip()
            symbol_value = ""
            for column in header:
                normalized = column.strip().lower().replace("_", " ")
                if normalized in {
                    "gene symbol",
                    "gene symbols",
                    "symbol",
                    "hgnc symbol",
                }:
                    symbol_value = row.get(column, "")
                    break
            if not symbol_value:
                for column in header:
                    if "gene symbol" in column.lower():
                        symbol_value = row.get(column, "")
                        break
            symbol = split_gene_symbol(symbol_value)
            if probe and symbol:
                mapping[probe] = symbol
    if not mapping:
        raise ValueError(f"no probe-to-symbol mapping found in {path}")
    return mapping


def split_gene_symbol(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() in {"nan", "none", "---"}:
        return ""
    text = re.split(r"///|//|;|,", text, maxsplit=1)[0].strip()
    text = re.sub(r"\s+", "", text)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", text):
        return ""
    return text.upper()


def collapse_probes(
    matrix: pd.DataFrame,
    probe_to_symbol: dict[str, str],
    *,
    log: logging.Logger | None = None,
) -> pd.DataFrame:
    """Collapse probes to symbols, retaining the probe with highest variance."""
    logger = log or LOG
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    numeric.index = numeric.index.astype(str)
    rows: list[tuple[str, float, pd.Series]] = []
    for probe, row in numeric.iterrows():
        symbol = probe_to_symbol.get(probe) or split_gene_symbol(probe)
        if not symbol:
            continue
        values = row.dropna()
        if values.empty:
            continue
        variance = float(values.var(ddof=0)) if len(values) > 1 else 0.0
        rows.append((symbol, variance, row))
    if not rows:
        raise ValueError("no probes could be mapped to gene symbols")
    selected: dict[str, tuple[float, pd.Series]] = {}
    for symbol, variance, row in rows:
        existing = selected.get(symbol)
        if existing is None or variance > existing[0]:
            selected[symbol] = (variance, row)
    collapsed = pd.DataFrame(
        {symbol: values for symbol, (_, values) in selected.items()}
    ).T
    collapsed.index.name = "gene"
    logger.info(
        "mapped %s probes to %s unique gene symbols",
        len(rows),
        collapsed.shape[0],
    )
    return collapsed


def standardize_expression(matrix: pd.DataFrame, *, already_log: bool = False) -> pd.DataFrame:
    numeric = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if not already_log:
        numeric = np.log2(numeric.clip(lower=0) + 1.0)
    return numeric


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    means = numeric.mean(axis=1)
    std = numeric.std(axis=1, ddof=0).replace(0, np.nan)
    return numeric.sub(means, axis=0).div(std, axis=0).fillna(0.0)


def save_expression(matrix: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    matrix.to_csv(path, index=True, compression="gzip" if path.suffix == ".gz" else None)


def read_expression(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    return frame


def figure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "figure.facecolor": "white",
            "savefig.dpi": 600,
            "savefig.transparent": False,
            "savefig.facecolor": "white",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 8.5,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.1,
            "lines.markersize": 4,
            "patch.linewidth": 0.7,
            "savefig.pad_inches": 0.035,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.axisbelow": True,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def save_figure(fig: plt.Figure, path: Path, *, tight: bool = True) -> Path:
    ensure_dir(path.parent)
    if tight:
        fig.tight_layout()
    save_kwargs = {
        "dpi": 600,
        "bbox_inches": "tight",
        "facecolor": "white",
        "metadata": {
            "Software": "experiment_plan_one",
            "Title": path.stem,
        },
    }
    fig.savefig(path, **save_kwargs)
    if path.suffix.lower() == ".png":
        with Image.open(path) as image:
            if image.mode == "RGBA":
                alpha = image.getchannel("A")
                flattened = Image.new("RGB", image.size, "white")
                flattened.paste(image.convert("RGB"), mask=alpha)
                flattened.save(
                    path,
                    format="PNG",
                    dpi=(600, 600),
                    optimize=True,
                )
        for suffix in (".pdf", ".svg"):
            fig.savefig(
                path.with_suffix(suffix),
                bbox_inches="tight",
                facecolor="white",
            )
    plt.close(fig)
    return path


figure_style()


def bh_fdr(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    result = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return result
    order = np.argsort(p[valid])
    sorted_p = p[valid][order]
    n = len(sorted_p)
    adjusted = sorted_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    valid_indices = np.flatnonzero(valid)
    result[valid_indices[order]] = adjusted
    return result


def benjamini_hochberg(frame: pd.DataFrame, p_column: str = "p_value") -> pd.DataFrame:
    output = frame.copy()
    output["fdr"] = bh_fdr(output[p_column])
    output["neg_log10_fdr"] = -np.log10(output["fdr"].clip(lower=1e-300))
    return output


def slug(value: str) -> str:
    return _slug(value)


def release_memory(*objects: Any) -> None:
    """Best-effort garbage collection used by large single-cell stages."""
    import gc

    del objects
    gc.collect()
