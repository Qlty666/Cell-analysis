"""Generic GEO supplementary-file downloader."""

from pathlib import Path
import json
import gzip
import re
import shutil
import subprocess
import tarfile
import urllib.request
import shutil

try:
    from .h5_converter import convert_h5ad, convert_loom
except ImportError:
    from h5_converter import convert_h5ad, convert_loom

CACHE_ROOT = Path(__file__).resolve().parents[2] / "data_cache"


def _curl() -> str:
    found = shutil.which("curl.exe") or shutil.which("curl")
    if found:
        return found
    raise RuntimeError("curl not found")


def _fetch(url: str) -> str:
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.read().decode("utf-8", "replace")
    except Exception:
        result = subprocess.run(
            [
                _curl(),
                "-L",
                "--ssl-no-revoke",
                "-A", "Mozilla/5.0",
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


def _download(url: str, out: Path, log) -> None:
    log(f"downloading {out.name}")
    subprocess.run(
        [
            _curl(),
            "-L",
            "--ssl-no-revoke",
            "-A", "Mozilla/5.0",
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


def series_prefix(accession: str) -> str:
    acc = normalize_accession(accession)
    digits = acc[3:]
    if len(digits) <= 3:
        return "GSE" + digits
    return "GSE" + digits[: len(digits) - 3] + "nnn"


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

    for name in names:
        low = name.lower()
        if "normalized" in low or "tpm" in low or "natural_log" in low:
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
        elif re.search(r"matrix|\.mtx|counts?|read_counts|umi_counts|\.rds$", low):
            matrices.append(name)

    return {
        "matrix": matrices,
        "barcodes": barcodes,
        "genes": genes,
        "metadata": metadata,
    }


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
    matrix = []
    barcodes = []
    genes = []
    for name in downloaded["matrix"]:
        path = raw_dir / name
        if name.lower().endswith(".h5ad"):
            result = convert_h5ad(path)
            matrix.append(result["matrix"])
            barcodes.append(result["barcodes"])
            genes.append(result["genes"])
        elif name.lower().endswith(".loom"):
            result = convert_loom(path)
            matrix.append(result["matrix"])
            barcodes.append(result["barcodes"])
            genes.append(result["genes"])
        else:
            matrix.append(name)
    downloaded["matrix"] = matrix
    downloaded["barcodes"] = barcodes or downloaded["barcodes"]
    downloaded["genes"] = genes or downloaded["genes"]
    return downloaded


def _extract_archive(archive: Path, dest: Path) -> None:
    dest = dest.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tf:
        for member in tf.getmembers():
            if member.issym() or member.islnk():
                raise RuntimeError(f"unsafe tar member: {member.name}")
            member_path = member.name.replace("\\", "/")
            target = (dest / member_path).resolve()
            if not target.is_relative_to(dest):
                raise RuntimeError(f"unsafe tar member: {member.name}")
        tf.extractall(dest)


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
            log(f"using cached download for {acc}")
            return cached

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        all_files = []
        for group in ["matrix", "barcodes", "genes", "metadata", "series_matrices"]:
            all_files.extend(manifest.get("files", {}).get(group, []))
        if all_files and all(
            (raw_dir / rel).exists() and (raw_dir / rel).stat().st_size > 0
            for rel in all_files
        ):
            log(f"{acc} already downloaded; skipping download")
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

    if not downloaded["matrix"]:
        archive_names = [
            n for n in names if re.search(r"\.tar(\.gz)?$", n, re.IGNORECASE)
        ]
        if not archive_names:
            raise RuntimeError(
                f"No count matrix files found for {acc}; "
                "this dataset may require additional manual configuration."
            )
        archive_name = Path(archive_names[0]).name
        archive_path = raw_dir / archive_name
        _download(base + archive_names[0], archive_path, log)
        extract_dir = raw_dir / "_extracted"
        extract_dir.mkdir(parents=True, exist_ok=True)
        _extract_archive(archive_path, extract_dir)
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
        downloaded = _convert_downloaded(downloaded, raw_dir)
        log(f"found {len(downloaded['matrix'])} matrix files inside archive")

    series_paths = []
    for name in series_files:
        safe_name = Path(name).name
        out = raw_dir / safe_name
        _download(matrix_url + name, out, log)
        series_paths.append(safe_name)

    organism = "hs"
    for name in series_paths:
        try:
            with gzip.open(raw_dir / name, "rt", encoding="utf-8", errors="replace") as fh:
                text = fh.read(20000)
            if "mus musculus" in text.lower():
                organism = "mm"
                break
        except Exception:
            continue

    if not downloaded["matrix"]:
        raise RuntimeError(
            f"No count matrix files found for {acc}; "
            "this dataset may require additional manual configuration."
        )

    manifest = {
        "accession": acc,
        "mode": "generic",
        "organism": organism,
        "files": {
            "matrix": downloaded["matrix"],
            "barcodes": downloaded["barcodes"],
            "genes": downloaded["genes"],
            "metadata": downloaded["metadata"],
            "series_matrices": series_paths,
        },
    }
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
