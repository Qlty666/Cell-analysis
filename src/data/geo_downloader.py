"""Generic GEO supplementary-file downloader."""

from pathlib import Path
import concurrent.futures
import json
import gzip
import logging
import re
import shutil
import subprocess
import time
import tarfile

from common.http import DEFAULT_USER_AGENT, HttpError, http_download, http_get

try:
    from .series_matrix import series_matrix_to_tsv
except ImportError:  # pragma: no cover - direct script execution fallback
    from series_matrix import series_matrix_to_tsv

logger = logging.getLogger(__name__)

USER_AGENT = DEFAULT_USER_AGENT

CACHE_ROOT = Path(__file__).resolve().parents[2] / "data_cache"

BULK_COUNT_TABLE_RE = re.compile(
    r"(?:^|/)GSM\d+_[^/]+\.(?:txt|tsv|csv)(?:\.gz)?$",
    re.IGNORECASE,
)

BULK_MATRIX_NAME_RE = re.compile(
    r"bulk[-_ ]?rna|bulk[-_ ]?seq|_bulk_|bulk_count|"
    r"raw_counts|count_matrix|expr_matrix|expression_matrix",
    re.IGNORECASE,
)

ARCHIVE_SUFFIXES = (
    ".tar",
    ".tar.gz",
    ".tgz",
    ".tar.bz2",
    ".tbz2",
    ".tar.xz",
    ".txz",
)

SINGLE_CELL_MATRIX_RE = re.compile(
    r"(?:^|[/_. -])filtered[-_. ]?feature[-_. ]?bc[-_. ]?matrix"
    r"|(?:^|[/_. -])raw[-_. ]?feature[-_. ]?bc[-_. ]?matrix"
    r"|(?:^|[/_. -])single[-_. ]?cell|singlecell"
    r"|scrna|scrnaseq|cellranger|10xgenomics|seurat|singlecellexperiment",
    re.IGNORECASE,
)


def _files_look_bulk(files: dict) -> bool:
    matrices = files.get("matrix") or []
    barcodes = files.get("barcodes") or []
    genes = files.get("genes") or []
    if any(BULK_COUNT_TABLE_RE.search(name) for name in matrices):
        return True
    if any(BULK_MATRIX_NAME_RE.search(name) for name in matrices):
        return True
    if any(
        name.lower().endswith((".rds", ".h5", ".h5ad", ".loom"))
        for name in matrices
    ):
        return False
    return bool(matrices) and not barcodes and not genes and not any(
        SINGLE_CELL_MATRIX_RE.search(name) for name in matrices
    )


def _files_look_single_cell(files: dict) -> bool:
    matrices = files.get("matrix") or []
    if any(
        name.lower().endswith((".rds", ".h5", ".h5ad", ".loom"))
        for name in matrices
    ):
        return True
    if files.get("barcodes") or files.get("genes"):
        return True
    return any(
        SINGLE_CELL_MATRIX_RE.search(name) for name in matrices
    )


_CELL_BARCODE_RE = re.compile(r"^[ACGTN]{8,}(?:-[0-9]+)?$", re.IGNORECASE)


def _matrix_header_looks_single_cell(path: Path) -> bool:
    try:
        opener = gzip.open if str(path).lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            first_line = fh.readline()
        if not first_line:
            return False
        fields = [
            field.strip('"').strip()
            for field in re.split(r"[,\t]", first_line.rstrip("\r\n"))
            if field.strip('"').strip()
        ]
        if not fields:
            return False
        return bool(
            _CELL_BARCODE_RE.match(fields[0])
            or (len(fields) >= 2 and _CELL_BARCODE_RE.match(fields[1]))
        )
    except Exception as exc:  # noqa: BLE001 - unreadable file means "not single cell"
        logger.warning("could not inspect matrix header %s: %s", path, exc)
        return False


def _series_matrix_says_single_cell(path: Path) -> bool:
    try:
        opener = gzip.open if str(path).lower().endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            text = fh.read(50000)
        low = text.lower()
        return any(
            token in low
            for token in (
                "single-cell",
                "single cell",
                "singlecell",
                "scrna",
                "10x genomics",
                "cell ranger",
                "smart-seq",
            )
        )
    except Exception as exc:  # noqa: BLE001 - unreadable file means "not single cell"
        logger.warning("could not read series matrix %s: %s", path, exc)
        return False


def _matrix_files_look_single_cell(files: dict, base_dir: Path) -> bool:
    for rel in files.get("matrix") or []:
        if _matrix_header_looks_single_cell(base_dir / rel):
            return True
    return False


# Ensembl-style organism codes. Aliases use word boundaries so that words such
# as "generated" or "collaborative" are not misread as "rat".
ORGANISM_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"\bhomo\s+sapiens\b", "hs"),
    (r"\bhuman\b", "hs"),
    (r"\bmus\s+musculus\b", "mm"),
    (r"\bmouse\b", "mm"),
    (r"\brattus\s+norvegicus\b", "rn"),
    (r"\brat\b", "rn"),
    (r"\bdanio\s+rerio\b", "dr"),
    (r"\bzebrafish\b", "dr"),
    (r"\bdrosophila\s+melanogaster\b", "dm"),
    (r"\bsaccharomyces\s+cerevisiae\b", "sc"),
    (r"\bcaenorhabditis\s+elegans\b", "ce"),
    (r"\bmacaca\s+mulatta\b", "mmul"),
    (r"\bsus\s+scrofa\b", "ss"),
    (r"\bcanis\s+lupus\s+familiaris\b", "cf"),
)


def detect_organism_code(*texts: str) -> str | None:
    """Return the Ensembl-style code for the first recognized organism.

    ``None`` means no supported organism was mentioned; callers must record
    that explicitly instead of silently assuming human.
    """
    haystack = " ".join(text for text in texts if text).lower()
    if not haystack:
        return None
    for pattern, code in ORGANISM_PATTERNS:
        if re.search(pattern, haystack):
            return code
    return None


def _is_archive(name: str) -> bool:
    return name.lower().endswith(ARCHIVE_SUFFIXES)


def _curl() -> str:
    found = shutil.which("curl.exe") or shutil.which("curl")
    if found:
        return found
    raise RuntimeError("curl not found")


def _fetch(url: str) -> str:
    """Fetch a small text/HTML resource, preferring the shared HTTP helper."""
    try:
        return http_get(
            url,
            timeout=120,
            retries=3,
            backoff=2.0,
            user_agent=USER_AGENT,
        ).decode("utf-8", "replace")
    except HttpError as exc:
        logger.warning(
            "urllib fetch failed for %s (%s); retrying with curl", url, exc
        )
    result = subprocess.run(
        [
            _curl(),
            "-L",
            "--ssl-no-revoke",
            "-A", USER_AGENT,
            "--silent",
            "--show-error",
            "--max-time", "120",
            url,
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True,
    )
    return result.stdout


def _download(url: str, out: Path, log, force: bool = False) -> None:
    if not force and out.exists() and out.stat().st_size > 0:
        log(f"{out.name} already exists; skipping download")
        return
    for attempt in range(1, 4):
        log(f"downloading {out.name} (attempt {attempt}/3)")
        try:
            subprocess.run(
                [
                    _curl(),
                    "-L",
                    "--ssl-no-revoke",
                    "-A", USER_AGENT,
                    "-C", "-",
                    "--retry", "5",
                    "--retry-delay", "3",
                    "--fail",
                    "--silent",
                    "--show-error",
                    "--max-time", "5400",
                    "-o", str(out),
                    url,
                ],
                check=True,
            )
            if not out.exists() or out.stat().st_size == 0:
                raise RuntimeError("downloaded file is empty")
            return
        except (subprocess.CalledProcessError, OSError, RuntimeError) as exc:
            if attempt == 3:
                # curl keeps resume support for multi-GB GEO archives; the
                # shared helper is the fallback when curl is missing/blocked.
                logger.warning(
                    "curl download failed for %s (%s); falling back to urllib",
                    url,
                    exc,
                )
                http_download(
                    url,
                    out,
                    timeout=5400,
                    retries=3,
                    backoff=3.0,
                    user_agent=USER_AGENT,
                    log=log,
                )
                return
            delay = 3 * attempt
            log(f"download attempt {attempt} failed ({exc}); retrying in {delay}s")
            time.sleep(delay)


def series_prefix(accession: str) -> str:
    acc = normalize_accession(accession)
    digits = acc[3:]
    if len(digits) <= 3:
        return "GSE" + digits
    return "GSE" + digits[: len(digits) - 3] + "nnn"


def canonical_accession(accession: str) -> str:
    acc = str(accession or "").strip().upper()
    match = re.fullmatch(r"E-GEOD-(\d+)", acc)
    if match:
        return "GSE" + match.group(1)
    return acc


def normalize_accession(accession: str) -> str:
    acc = accession.strip().upper()
    if not re.fullmatch(r"GSE\d+", acc):
        raise ValueError("GSE accession must look like GSE125449")
    return acc


def _select_files(names: list[str]) -> dict:
    barcodes = []
    genes = []
    matrices = []
    metadata = []
    bulk_matrices = []
    single_cell_matrices = []

    for name in names:
        low = name.lower()
        if _is_archive(name):
            continue
        if any(
            token in low
            for token in (
                "normalized",
                "tpm",
                "fpkm",
                "rpkm",
                "cpm",
                "natural_log",
                "series_matrix",
            )
        ):
            continue
        if re.search(
            r"samples|metadata|cellinfo|cell_info|phenotype|celltype|annotation",
            low,
        ):
            metadata.append(name)
        elif "barcode" in low:
            barcodes.append(name)
        elif "genes" in low or "features" in low:
            genes.append(name)
        elif re.search(
            r"matrix|\.mtx|counts?|read_counts|umi_counts|"
            r"rna[-_ ]?seq|rnaseq|expression|\.rds$|"
            r"\.h5ad$|\.h5$|\.loom$",
            low,
        ) or BULK_COUNT_TABLE_RE.search(name):
            if BULK_COUNT_TABLE_RE.search(name) or BULK_MATRIX_NAME_RE.search(name):
                bulk_matrices.append(name)
            else:
                single_cell_matrices.append(name)

    # Prefer bulk count tables when a series contains both sample-level and
    # single-cell matrices (e.g. a bulk cohort plus scRNA-seq subsets). Mixing
    # them produces matrices with incompatible column semantics.
    if bulk_matrices:
        matrices = bulk_matrices
    else:
        matrices = single_cell_matrices

    return {
        "matrix": matrices,
        "barcodes": barcodes,
        "genes": genes,
        "metadata": metadata,
        "bulk": _files_look_bulk(
            {"matrix": matrices, "barcodes": barcodes, "genes": genes}
        ),
    }


def _expand_archive_files(files: dict, base_dir: Path, log) -> dict:
    """Replace archive entries in a manifest with files extracted from them."""
    expanded = {
        "matrix": list(files.get("matrix", [])),
        "barcodes": list(files.get("barcodes", [])),
        "genes": list(files.get("genes", [])),
        "metadata": list(files.get("metadata", [])),
        "series_matrices": list(files.get("series_matrices", [])),
    }
    archives = []
    for group in ("matrix", "barcodes", "genes", "metadata", "series_matrices"):
        for name in files.get(group, []):
            if _is_archive(name) and name not in archives:
                archives.append(name)
    if not archives:
        return files

    extract_dir = base_dir / "_extracted"
    for name in archives:
        archive_path = base_dir / name
        if not archive_path.exists() or archive_path.stat().st_size == 0:
            continue
        log(f"expanding archive {name}")
        extract_dir.mkdir(parents=True, exist_ok=True)
        _extract_archive(archive_path, extract_dir)
        inner = _select_files(_walk_relative(extract_dir))
        for group in ("matrix", "barcodes", "genes", "metadata"):
            expanded[group] = [
                item for item in expanded[group] if item != name
            ] + [f"_extracted/{item}" for item in inner[group]]

    for group, values in expanded.items():
        seen = set()
        unique = []
        for value in values:
            if value not in seen:
                seen.add(value)
                unique.append(value)
        expanded[group] = unique
    return expanded


def _refresh_manifest_mode(manifest: dict, base_dir: Path | None = None) -> dict:
    files = manifest.get("files")
    if isinstance(files, dict):
        if manifest.get("single_cell_hint"):
            manifest["mode"] = "single_cell"
        elif _files_look_bulk(files):
            looks_single_cell = (
                base_dir is not None
                and _matrix_files_look_single_cell(
                    files,
                    base_dir,
                )
            )
            if looks_single_cell:
                manifest["single_cell_hint"] = True
                manifest["mode"] = "single_cell"
            else:
                manifest["mode"] = "bulk"
        elif _files_look_single_cell(files):
            manifest["mode"] = "single_cell"
        else:
            manifest["mode"] = "generic"
    return manifest


def _single_cell_matrix_relatives(base_dir: Path) -> list[str]:
    """Find raw single-cell matrix files in a download directory."""
    if base_dir is None or not base_dir.exists():
        return []
    found: list[str] = []
    suffixes = (".h5ad", ".h5ad.gz", ".h5", ".loom", ".rds")
    for path in sorted(base_dir.rglob("*")):
        if not path.is_file():
            continue
        name = path.name.lower()
        if not name.endswith(suffixes):
            continue
        if re.search(r"bulk[-_ ]?rna|bulk[-_ ]?seq|_bulk_", name):
            continue
        rel = path.relative_to(base_dir).as_posix()
        if rel not in found:
            found.append(rel)
    return found


def _prefer_single_cell_matrices(manifest: dict, base_dir: Path | None) -> dict:
    """Keep .h5ad/.h5/.loom matrices when a series is single-cell hinted.

    A GEO series may bundle bulk count tables and single-cell matrices. The
    old downloader deliberately picked bulk tables first, which caused
    datasets such as GSE235863 to be analyzed as a handful of pseudo-cells.
    When the series metadata clearly says single-cell, retain the actual
    single-cell matrices so downstream code can use their cell metadata.
    """
    files = manifest.get("files")
    if not isinstance(files, dict) or not manifest.get("single_cell_hint"):
        return manifest
    suffixes = (".h5ad", ".h5ad.gz", ".h5", ".loom", ".rds")
    listed = [
        name
        for name in files.get("matrix") or []
        if name.lower().endswith(suffixes)
        and not re.search(r"bulk[-_ ]?rna|bulk[-_ ]?seq|_bulk_", name, re.IGNORECASE)
    ]
    sc_from_disk = _single_cell_matrix_relatives(base_dir)
    sc_names = list(dict.fromkeys(listed + sc_from_disk))
    if not sc_names:
        return manifest
    files["matrix"] = sc_names
    files["barcodes"] = []
    files["genes"] = []
    manifest["mode"] = "single_cell"
    return manifest


def _download_files(urls: dict, raw_dir: Path, log) -> dict:
    downloaded = {}
    for category in ["matrix", "barcodes", "genes", "metadata"]:
        downloaded[category] = []
        files = urls[category]
        for name in files:
            out = raw_dir / Path(name).name
            _download(urls["_base"] + name, out, log)
            downloaded[category].append(name)
    return downloaded


def _convert_downloaded(downloaded: dict, raw_dir: Path) -> dict:
    try:
        from .h5_converter import convert_h5ad, convert_loom
    except ImportError:
        convert_h5ad = None
        convert_loom = None

    def convert_one(name: str) -> dict:
        path = raw_dir / name
        prefix = Path(name).parent.as_posix()
        if prefix == ".":
            prefix = ""
        low = name.lower()
        if low.endswith(".h5ad.gz"):
            if convert_h5ad is None:
                raise RuntimeError(
                    "h5py/scipy are required to convert .h5ad inputs; "
                    "run: python -m pip install h5py scipy"
                )
            from .h5_converter import convert_h5ad_gz

            result = convert_h5ad_gz(path)
            return {
                "matrix": f"{prefix}/{result['matrix']}" if prefix else result["matrix"],
                "barcodes": (
                    f"{prefix}/{result['barcodes']}" if prefix else result["barcodes"]
                ),
                "genes": f"{prefix}/{result['genes']}" if prefix else result["genes"],
            }
        if low.endswith(".h5ad"):
            if convert_h5ad is None:
                raise RuntimeError(
                    "h5py/scipy are required to convert .h5ad inputs; "
                    "run: python -m pip install h5py scipy"
                )
            result = convert_h5ad(path)
            return {
                "matrix": f"{prefix}/{result['matrix']}" if prefix else result["matrix"],
                "barcodes": (
                    f"{prefix}/{result['barcodes']}" if prefix else result["barcodes"]
                ),
                "genes": f"{prefix}/{result['genes']}" if prefix else result["genes"],
            }
        elif name.lower().endswith(".loom"):
            if convert_loom is None:
                raise RuntimeError(
                    "h5py/scipy are required to convert .loom inputs; "
                    "run: python -m pip install h5py scipy"
                )
            result = convert_loom(path)
            return {
                "matrix": f"{prefix}/{result['matrix']}" if prefix else result["matrix"],
                "barcodes": (
                    f"{prefix}/{result['barcodes']}" if prefix else result["barcodes"]
                ),
                "genes": f"{prefix}/{result['genes']}" if prefix else result["genes"],
            }
        return {"matrix": name, "barcodes": [], "genes": []}

    names = downloaded["matrix"]
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=min(4, max(1, len(names)))
    ) as pool:
        results = list(pool.map(convert_one, names))

    matrix = []
    barcodes = []
    genes = []
    for result in results:
        matrix.append(result["matrix"])
        barcodes.extend(
            result["barcodes"]
            if isinstance(result["barcodes"], list)
            else [result["barcodes"]]
        )
        genes.extend(
            result["genes"]
            if isinstance(result["genes"], list)
            else [result["genes"]]
        )
    downloaded["matrix"] = matrix
    downloaded["barcodes"] = barcodes or downloaded["barcodes"]
    downloaded["genes"] = genes or downloaded["genes"]
    return downloaded


def _extract_archive(archive: Path, dest: Path) -> None:
    dest = dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive, "r:*") as tf:
            for member in tf.getmembers():
                if member.issym() or member.islnk():
                    raise RuntimeError(f"unsafe tar member: {member.name}")
                member_path = member.name.replace("\\", "/")
                target = (dest / member_path).resolve()
                if not target.is_relative_to(dest):
                    raise RuntimeError(f"unsafe tar member: {member.name}")
            tf.extractall(dest)
    except (tarfile.TarError, EOFError) as exc:
        raise RuntimeError(
            f"archive is corrupted or truncated: {archive.name} ({exc})"
        ) from exc


def _extract_geo_archive(
    archive_path: Path,
    extract_dir: Path,
    url: str,
    log,
) -> None:
    """Extract a GEO archive, re-downloading once when it is corrupt."""
    try:
        with tarfile.open(archive_path, "r:*") as tf:
            tf.getmembers()
        valid = True
    except (tarfile.TarError, EOFError, OSError):
        valid = False
    if not valid:
        log(f"{archive_path.name} is corrupted or truncated; re-downloading")
        shutil.rmtree(extract_dir, ignore_errors=True)
        archive_path.unlink(missing_ok=True)
        _download(url, archive_path, log, force=True)
    elif extract_dir.exists() and any(extract_dir.iterdir()):
        log(f"{archive_path.name} already extracted; skipping extraction")
        return
    extract_dir.mkdir(parents=True, exist_ok=True)
    _extract_archive(archive_path, extract_dir)


def _walk_relative(directory: Path) -> list[str]:
    return [
        p.relative_to(directory).as_posix()
        for p in directory.rglob("*")
        if p.is_file()
    ]


def ensure_geo_dataset(accession: str, root: Path, log) -> dict:
    acc = normalize_accession(accession)
    prefix = series_prefix(acc)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/suppl/"
    raw_dir = root / "data" / "raw" / acc
    raw_dir.mkdir(parents=True, exist_ok=True)

    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = data_dir / f"{acc}_manifest.json"

    cache_dir = CACHE_ROOT / acc
    cache_manifest = cache_dir / "manifest.json"
    if (
        not manifest_path.exists()
        and cache_manifest.exists()
    ):
        cached = json.loads(cache_manifest.read_text(encoding="utf-8"))
        cached["files"] = _expand_archive_files(
            cached.get("files", {}),
            cache_dir,
            log,
        )
        cached = _refresh_manifest_mode(cached, cache_dir)
        cached = _prefer_single_cell_matrices(cached, cache_dir)
        all_files = []
        for group in ["matrix", "barcodes", "genes", "metadata", "series_matrices"]:
            all_files.extend(cached.get("files", {}).get(group, []))
        if all_files and all(
            (cache_dir / rel).exists()
            for rel in all_files
        ):
            for rel in all_files:
                dest = raw_dir / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(cache_dir / rel, dest)
            manifest_path.write_text(
                json.dumps(cached, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            cache_manifest.write_text(
                json.dumps(cached, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            log(f"using cached download for {acc}")
            return cached

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["files"] = _expand_archive_files(
            manifest.get("files", {}),
            raw_dir,
            log,
        )
        manifest = _refresh_manifest_mode(manifest, raw_dir)
        manifest = _prefer_single_cell_matrices(manifest, raw_dir)
        all_files = []
        for group in ["matrix", "barcodes", "genes", "metadata", "series_matrices"]:
            all_files.extend(manifest.get("files", {}).get(group, []))
        if all_files and all(
            (raw_dir / rel).exists() and (raw_dir / rel).stat().st_size > 0
            for rel in all_files
        ):
            log(f"{acc} already downloaded; skipping download")
            manifest_path.write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            return manifest

    html = _fetch(base)
    names = [
        re.sub(r"^.*/", "", m)
        for m in re.findall(r'href="([^"]+)"', html)
        if not m.endswith("/")
    ]
    selected = _select_files(names)

    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/matrix/"
    series_files = []
    try:
        matrix_html = _fetch(matrix_url)
        series_files = [
            re.sub(r"^.*/", "", m)
            for m in re.findall(r'href="([^"]+)"', matrix_html)
            if "series_matrix" in m and not m.endswith("/")
        ]
    except subprocess.CalledProcessError:
        log("series matrix directory not available; continuing without it")

    all_urls = {
        "_base": base,
        "matrix": selected["matrix"],
        "barcodes": selected["barcodes"],
        "genes": selected["genes"],
        "metadata": selected["metadata"],
    }
    downloaded = _download_files(all_urls, raw_dir, log)
    downloaded = _convert_downloaded(downloaded, raw_dir)
    downloaded = _expand_archive_files(downloaded, raw_dir, log)
    downloaded["bulk"] = bool(selected["bulk"])

    archive_names = [
        n for n in names if re.search(r"\.tar(\.gz)?$", n, re.IGNORECASE)
    ]
    if not downloaded["matrix"] and archive_names:
        archive_name = Path(archive_names[0]).name
        archive_path = raw_dir / archive_name
        extract_dir = raw_dir / "_extracted"
        _extract_geo_archive(
            archive_path,
            extract_dir,
            base + archive_names[0],
            log,
        )
        inner = _select_files(_walk_relative(extract_dir))
        downloaded["matrix"] = [
            "_extracted/" + name for name in inner["matrix"]
        ]
        downloaded["barcodes"] = [
            "_extracted/" + name for name in inner["barcodes"]
        ]
        downloaded["genes"] = [
            "_extracted/" + name for name in inner["genes"]
        ]
        downloaded["metadata"] = [
            "_extracted/" + name for name in inner["metadata"]
        ]
        downloaded["bulk"] = bool(inner["bulk"])
        downloaded = _convert_downloaded(downloaded, raw_dir)
        log(f"found {len(downloaded['matrix'])} matrix files inside archive")
    if not downloaded["matrix"] and not series_files:
        raise RuntimeError(
            f"No count matrix files found for {acc}; "
            "this dataset may require additional manual configuration."
        )

    series_paths = []
    for name in series_files:
        safe_name = Path(name).name
        out = raw_dir / safe_name
        _download(matrix_url + name, out, log)
        series_paths.append(safe_name)

    series_matrix_fallback = False
    if not downloaded["matrix"] and series_paths:
        for series_name in series_paths:
            source = raw_dir / series_name
            destination = raw_dir / (
                Path(series_name).name.replace(".txt.gz", "")
                + ".expression.tsv.gz"
            )
            try:
                series_matrix_to_tsv(source, destination)
            except (OSError, ValueError) as exc:
                log(f"could not parse series matrix {series_name}: {exc}")
                continue
            downloaded["matrix"].append(destination.name)
            series_matrix_fallback = True
        if series_matrix_fallback:
            downloaded["bulk"] = True
            log(
                f"using {len(downloaded['matrix'])} series-matrix expression "
                "table(s) because no supplementary count matrix was found"
            )

    organism_texts: list[str] = []
    for name in series_paths:
        try:
            with gzip.open(raw_dir / name, "rt", encoding="utf-8", errors="replace") as fh:
                organism_texts.append(fh.read(20000))
        except (OSError, EOFError) as exc:
            logger.warning(
                "could not read series matrix %s for organism detection: %s",
                raw_dir / name,
                exc,
            )
    organism = detect_organism_code(*organism_texts) or "unknown"

    single_cell_hint = any(
        _series_matrix_says_single_cell(raw_dir / name)
        for name in series_paths
    )
    if not single_cell_hint:
        single_cell_hint = _matrix_files_look_single_cell(downloaded, raw_dir)

    if not downloaded["matrix"]:
        raise RuntimeError(
            f"No count matrix files found for {acc}; "
            "this dataset may require additional manual configuration."
        )

    manifest = {
        "accession": acc,
        "mode": "generic",
        "single_cell_hint": bool(single_cell_hint),
        "organism": organism,
        "files": {
            "matrix": downloaded["matrix"],
            "barcodes": downloaded["barcodes"],
            "genes": downloaded["genes"],
            "metadata": downloaded["metadata"],
            "series_matrices": series_paths,
        },
    }
    manifest = _refresh_manifest_mode(manifest, raw_dir)
    manifest = _prefer_single_cell_matrices(manifest, raw_dir)
    if series_matrix_fallback and not any(
        str(name).lower().endswith((".h5ad", ".h5ad.gz", ".h5", ".loom", ".rds"))
        for name in manifest["files"]["matrix"]
    ):
        manifest["mode"] = "bulk"
        manifest["single_cell_hint"] = False
        manifest["data_type"] = "microarray_or_normalized"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    for rel in (
        manifest["files"]["matrix"]
        + manifest["files"]["barcodes"]
        + manifest["files"]["genes"]
        + manifest["files"]["metadata"]
        + manifest["files"]["series_matrices"]
    ):
        src = raw_dir / rel
        dest = cache_dir / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if src.exists():
            shutil.copy2(src, dest)
    cache_manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    log(f"manifest written: {manifest_path}")
    return manifest
