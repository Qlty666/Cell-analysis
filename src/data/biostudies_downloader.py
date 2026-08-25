"""Download processed expression files from EBI BioStudies/ArrayExpress."""

from __future__ import annotations

import json
import re
import shutil
import time
import urllib.request
from pathlib import Path

from . import geo_downloader as gd

BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"
CACHE_ROOT = Path(__file__).resolve().parents[2] / "data_cache"
USER_AGENT = "Mozilla/5.0 (liver-cancer-pipeline; dataset-download)"

ACCESSION_RE = re.compile(r"(?:E-[A-Z0-9]+-\d+|S-BSST\d+)", re.IGNORECASE)


def normalize_accession(accession: str) -> str:
    acc = accession.strip().upper()
    if not ACCESSION_RE.fullmatch(acc):
        raise ValueError(
            "BioStudies accession must look like E-MTAB-1234 or S-BSST123"
        )
    return acc


def _http_get_json(url: str, timeout: int = 120, retries: int = 3) -> dict:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": USER_AGENT},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt < retries - 1:
                time.sleep(1 + attempt * 2)
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"failed to fetch: {url}")


def _safe_rel(path: str) -> str:
    rel = str(path or "").replace("\\", "/").lstrip("/")
    if not rel or ".." in Path(rel).parts:
        raise RuntimeError(f"unsafe BioStudies file path: {rel}")
    return rel


def _attributes(section: dict) -> dict:
    result: dict[str, str] = {}
    for item in section.get("attributes") or []:
        if isinstance(item, dict) and item.get("name"):
            result[str(item["name"])] = str(item.get("value") or "")
    return result


def _study_metadata(accession: str) -> dict:
    data = _http_get_json(f"{BIOSTUDIES_API}/studies/{accession}")
    attributes = _attributes(data)
    section = data.get("section") or {}
    section_attributes = _attributes(section)
    return {
        "root_path": attributes.get("RootPath") or accession,
        "title": section_attributes.get("Title") or attributes.get("Title") or "",
        "organism": section_attributes.get("Organism") or "",
        "description": section_attributes.get("Description") or "",
        "study_type": section_attributes.get("Study type") or "",
    }


def _list_files(accession: str) -> list[dict]:
    try:
        data = _http_get_json(f"{BIOSTUDIES_API}/studies/{accession}/files")
        items = data.get("items") or []
        if items:
            return items
    except Exception:  # noqa: BLE001
        pass
    data = _http_get_json(f"{BIOSTUDIES_API}/studies/{accession}")
    items: list[dict] = []

    def walk(node) -> None:
        if isinstance(node, dict):
            files = node.get("files")
            if isinstance(files, list):
                section = str(node.get("accno") or node.get("type") or "")
                for item in files:
                    if isinstance(item, dict) and item.get("path"):
                        items.append(
                            {
                                "path": _safe_rel(str(item["path"])),
                                "Size": item.get("size", 0),
                                "Section": section,
                                "Description": str(item.get("Description") or ""),
                            }
                        )
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    walk(data)
    return items


def _select_processed_files(items: list[dict]) -> dict:
    matrix: list[str] = []
    barcodes: list[str] = []
    genes: list[str] = []
    metadata: list[str] = []
    for item in items:
        path = _safe_rel(str(item.get("path") or ""))
        low = path.lower()
        name = Path(path).name.lower()
        section = str(item.get("Section") or "").lower()
        if re.search(
            r"\.(fastq|fq|bam|cram|bcf|vcf|vcf\.gz|tif|tiff|png|jpe?g|pdf|"
            r"html?|xlsx?|docx?|pptx?)(\.gz)?$",
            low,
        ):
            continue
        if re.search(r"raw[-_ ]?data|mage-tab|idf|sdrf", section):
            continue
        if gd._is_archive(path):
            if re.search(r"processed|count|matrix|expression", low) or "processed" in section:
                matrix.append(path)
            continue
        if re.search(
            r"samples|sample_info|sample_metadata|metadata|cellinfo|"
            r"phenotype|celltype|annotation|sdrf|idf|design|factors",
            low,
        ):
            metadata.append(path)
        elif "barcode" in low:
            barcodes.append(path)
        elif re.search(r"(^|[/_. -])genes?($|[/_. -])|features", low):
            genes.append(path)
        elif re.search(
            r"\.(mtx(\.gz)?|h5ad|rds|loom|h5|txt(\.gz)?|tsv(\.gz)?|csv(\.gz)?)$",
            low,
        ) or re.search(r"count|matrix|expression|tpm|fpkm|rpkm|normaliz|processed", low):
            matrix.append(path)
    return {
        "matrix": matrix,
        "barcodes": barcodes,
        "genes": genes,
        "metadata": metadata,
    }


def _download_files(
    selected: dict,
    root_path: str,
    raw_dir: Path,
    log,
) -> dict:
    downloaded: dict = {
        "matrix": [],
        "barcodes": [],
        "genes": [],
        "metadata": [],
    }
    for category in ("matrix", "barcodes", "genes", "metadata"):
        for rel in selected.get(category, []):
            rel = _safe_rel(rel)
            out = raw_dir / rel
            out.parent.mkdir(parents=True, exist_ok=True)
            url = f"https://www.ebi.ac.uk/biostudies/files/{root_path}/{rel}"
            gd._download(url, out, log)
            downloaded[category].append(rel)
    return downloaded


def _manifest_files_exist(manifest: dict, raw_dir: Path) -> bool:
    files = manifest.get("files") or {}
    all_files: list[str] = []
    for group in ("matrix", "barcodes", "genes", "metadata", "series_matrices"):
        all_files.extend(files.get(group) or [])
    return bool(all_files) and all(
        (raw_dir / rel).exists() and (raw_dir / rel).stat().st_size > 0
        for rel in all_files
    )


def ensure_biostudies_dataset(
    accession: str,
    root: Path,
    log,
) -> dict:
    acc = normalize_accession(accession)
    raw_dir = root / "data" / "raw" / acc
    raw_dir.mkdir(parents=True, exist_ok=True)
    data_dir = root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = data_dir / f"{acc}_manifest.json"
    cache_dir = CACHE_ROOT / acc
    cache_manifest = cache_dir / "manifest.json"

    if not manifest_path.exists() and cache_manifest.exists():
        cached = json.loads(cache_manifest.read_text(encoding="utf-8"))
        cached["files"] = gd._expand_archive_files(
            cached.get("files", {}),
            cache_dir,
            log,
        )
        cached = gd._refresh_manifest_mode(cached, cache_dir)
        if _manifest_files_exist(cached, cache_dir):
            for group in ("matrix", "barcodes", "genes", "metadata", "series_matrices"):
                for rel in cached.get("files", {}).get(group, []):
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
        manifest["files"] = gd._expand_archive_files(
            manifest.get("files", {}),
            raw_dir,
            log,
        )
        manifest = gd._refresh_manifest_mode(manifest, raw_dir)
        if _manifest_files_exist(manifest, raw_dir):
            log(f"{acc} already downloaded; skipping download")
            manifest_path.write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            return manifest

    metadata = _study_metadata(acc)
    items = _list_files(acc)
    selected = _select_processed_files(items)
    if not selected["matrix"]:
        raise RuntimeError(
            f"No processed count matrix files found for {acc}; "
            "this ArrayExpress/BioStudies dataset may only contain raw or "
            "non-expression files. Prefer a GEO GSE dataset or download "
            "files manually."
        )

    downloaded = _download_files(selected, metadata["root_path"], raw_dir, log)
    downloaded = gd._convert_downloaded(downloaded, raw_dir)
    downloaded = gd._expand_archive_files(downloaded, raw_dir, log)
    if not downloaded["matrix"]:
        raise RuntimeError(
            f"No usable count matrix files found after download for {acc}."
        )

    single_cell_hint = gd._matrix_files_look_single_cell(downloaded, raw_dir)
    organism = "hs"
    if re.search(
        r"mus musculus",
        f"{metadata['organism']} {metadata['description']}",
        re.IGNORECASE,
    ):
        organism = "mm"

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
            "series_matrices": [],
        },
    }
    manifest = gd._refresh_manifest_mode(manifest, raw_dir)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    cache_dir.mkdir(parents=True, exist_ok=True)
    for group in ("matrix", "barcodes", "genes", "metadata", "series_matrices"):
        for rel in manifest["files"].get(group, []):
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
