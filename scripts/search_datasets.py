#!/usr/bin/env python3
"""Search public expression datasets and optionally download matching series.

Currently searches NCBI GEO, EBI ArrayExpress/BioStudies and EBI Expression
Atlas. Only datasets with runnable processed expression files are marked as
pipeline-ready.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1"
ATLAS_API = "https://www.ebi.ac.uk/gxa/json/experiments"
ATLAS_CACHE = APP_ROOT / "data_cache" / "dataset_search" / "expression_atlas_index.json"
ATLAS_CACHE_TTL_SECONDS = 7 * 24 * 3600
USER_AGENT = "Mozilla/5.0 (liver-cancer-pipeline; dataset-search)"

SUPPORTED_DATABASES = ("geo", "biostudies", "atlas")
DATABASE_LABELS = {
    "geo": "GEO",
    "biostudies": "ArrayExpress/BioStudies",
    "atlas": "Expression Atlas",
}

CSV_COLUMNS = [
    "accession",
    "database",
    "disease",
    "research_direction",
    "title",
    "summary",
    "data_type",
    "organism",
    "samples",
    "platform",
    "date",
    "type",
    "url",
    "quality_score",
    "run_supported",
    "relevance_score",
]

DISEASE_SYNONYMS = {
    "liver cancer": [
        "hepatocellular carcinoma",
        "hcc",
        "hepatic",
    ],
    "breast cancer": [
        "breast cancer",
        "breast carcinoma",
        "brca",
    ],
    "lung cancer": ["lung cancer", "nsclc", "lung carcinoma"],
    "colorectal cancer": ["colorectal cancer", "colon cancer", "crc"],
    "gastric cancer": ["gastric cancer", "stomach cancer"],
    "pancreatic cancer": [
        "pancreatic cancer",
        "pdac",
        "pancreatic ductal adenocarcinoma",
    ],
    "kidney cancer": [
        "renal cell carcinoma",
        "renal carcinoma",
        "rcc",
    ],
    "ovarian cancer": ["ovarian cancer", "ovarian carcinoma"],
    "melanoma": ["melanoma", "skin cancer"],
    "glioma": ["glioma", "glioblastoma", "gbm", "brain tumor"],
    "alzheimer's disease": ["alzheimer", "alzheimer's disease"],
    "diabetes": ["diabetes", "type 2 diabetes"],
    "covid-19": ["covid-19", "sars-cov-2", "covid 19"],
    "rheumatoid arthritis": ["rheumatoid arthritis"],
    "atherosclerosis": ["atherosclerosis", "atherosclerotic"],
}

DIRECTION_SYNONYMS = {
    "single cell RNA-seq": [
        "single cell",
        "single-cell",
        "scrna",
        "scrna-seq",
        "rna-seq",
    ],
    "bulk RNA-seq": ["bulk rna-seq", "rna-seq", "bulk rna"],
    "spatial transcriptomics": [
        "spatial transcriptomic",
        "spatial",
        "visium",
        "stereo-seq",
    ],
    "tumor microenvironment": [
        "tumor microenvironment",
        "tumour microenvironment",
        "microenvironment",
        "immune infiltration",
    ],
    "drug resistance": [
        "drug resistance",
        "chemoresistance",
        "resistance",
    ],
    "immunotherapy": [
        "immunotherapy",
        "immune checkpoint",
        "pd-1",
        "pd-l1",
        "car-t",
    ],
    "DNA methylation": ["dna methylation", "methylation", "epigenetic"],
    "copy number variation": ["copy number", "cnv"],
    "prognosis": ["prognosis", "prognostic", "survival"],
    "metastasis": [
        "metastasis",
        "metastatic",
        "invasion",
        "migration",
    ],
    "cell differentiation": ["differentiation", "cell fate", "lineage"],
    "cellular senescence": ["senescence", "aging", "ageing"],
}


def build_expanded_query(disease: str, direction: str) -> str:
    disease_terms = [disease] + _mapping_values(
        DISEASE_SYNONYMS,
        disease,
    )
    direction_terms = [direction] + _mapping_values(
        DIRECTION_SYNONYMS,
        direction,
    )
    parts: list[str] = []
    if any(disease_terms):
        parts.append(
            "(" + " OR ".join(f'"{term}"' for term in disease_terms) + ")"
        )
    if any(direction_terms):
        parts.append(
            "(" + " OR ".join(f'"{term}"' for term in direction_terms) + ")"
        )
    return " AND ".join(parts)


def _mapping_values(mapping: dict, key: str) -> list[str]:
    target = key.lower()
    for name, values in mapping.items():
        if name.lower() == target:
            return values
    return []


def _to_int(value) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def _parse_date(value) -> tuple[int, int, int] | None:
    text = str(value or "")
    match = re.search(r"(\d{4})[-/](\d{1,2})[-/](\d{1,2})", text)
    if match:
        return tuple(int(part) for part in match.groups())
    match = re.search(r"(\d{4})", text)
    if match:
        return (int(match.group(1)), 1, 1)
    return None


def _http_get(url: str, timeout: int = 90, retries: int = 3) -> str:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt < retries - 1:
                time.sleep(1 + attempt * 2)
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"failed to fetch: {url}")


def _http_get_json(url: str, timeout: int = 90, retries: int = 3) -> dict:
    return json.loads(_http_get(url, timeout=timeout, retries=retries))


SPECIES_TERMS = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Danio rerio",
    "Drosophila melanogaster",
    "Saccharomyces cerevisiae",
    "Caenorhabditis elegans",
    "Macaca mulatta",
    "Sus scrofa",
    "Canis lupus familiaris",
]


def _detect_organism(text: str) -> str:
    low = text.lower()
    for term in SPECIES_TERMS:
        if term.lower() in low:
            return term
    return ""


def _detect_study_type(text: str) -> str:
    patterns = [
        r"spatial transcriptomics[^.;\n]*",
        r"transcription profiling by [^.;\n]*",
        r"transcriptomics by [^.;\n]*",
        r"genome binding/occupancy profiling[^.;\n]*",
        r"methylation profiling[^.;\n]*",
        r"protein expression profiling[^.;\n]*",
        r"expression profiling by [^.;\n]*",
        r"proteomics[^.;\n]*",
    ]
    for pattern in patterns:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return match.group(0).strip()
    return ""


def infer_data_type(text: str) -> str:
    haystack = text.lower()
    single_cell_terms = (
        "single cell",
        "single-cell",
        "scrna",
        "scrna-seq",
        "10x genomics",
        "cell ranger",
        "smart-seq",
        "spatial transcriptomic",
        "spatial transcriptomics",
    )
    if any(token in haystack for token in single_cell_terms):
        return "single-cell"
    if any(
        token in haystack
        for token in ("bulk rna", "bulk-rna", "rna-seq", "rnaseq", "bulk")
    ):
        return "bulk"
    return "other"


def _samples_int(row: dict) -> int:
    try:
        return int(float(str(row.get("samples") or 0)))
    except (TypeError, ValueError):
        return 0


def quality_score(row: dict) -> float:
    score = 0.0
    database = str(row.get("database", ""))
    if database == "Expression Atlas":
        score += 0.25
    elif database == "ArrayExpress/BioStudies":
        score += 0.10
    organism = str(row.get("organism", "")).lower()
    if organism in ("homo sapiens", "human", "mus musculus", "mouse"):
        score += 0.10
    if row.get("data_type") == "single-cell":
        score += 0.25
    elif row.get("data_type") == "bulk":
        score += 0.15
    samples = _samples_int(row)
    if samples >= 10:
        score += 0.15
    elif samples >= 3:
        score += 0.05
    type_text = str(row.get("type", "")).lower()
    if re.search(
        r"expression profiling|transcription profiling|rna-seq|rnaseq|transcriptomics",
        type_text,
    ):
        score += 0.15
    if row.get("run_supported"):
        score += 0.10
    return round(min(score, 1.0), 4)


def canonical_row_accession(accession: str) -> str:
    acc = str(accession or "").strip().upper()
    match = re.fullmatch(r"E-GEOD-(\d+)", acc)
    if match:
        return "GSE" + match.group(1)
    return acc


def _normalize_row(row: dict) -> dict:
    normalized = dict(row)
    normalized.setdefault("accession", "")
    normalized.setdefault("database", "GEO")
    normalized.setdefault("disease", "")
    normalized.setdefault("research_direction", "")
    normalized.setdefault("title", "")
    normalized.setdefault("summary", "")
    if not normalized.get("data_type"):
        normalized["data_type"] = infer_data_type(
            str(normalized.get("title", ""))
            + "\n"
            + str(normalized.get("summary", ""))
        )
    normalized.setdefault("organism", "")
    normalized.setdefault("samples", "")
    normalized.setdefault("platform", "")
    normalized.setdefault("date", "")
    normalized.setdefault("type", "")
    normalized.setdefault("url", "")
    normalized.setdefault("quality_score", 0.0)
    normalized.setdefault("run_supported", False)
    return normalized


def esearch(query: str, max_results: int = 20) -> list[str]:
    params = urllib.parse.urlencode(
        {
            "db": "gds",
            "term": query,
            "retmax": max(max_results, 1),
            "retmode": "xml",
        }
    )
    root = ET.fromstring(_http_get(f"{EUTILS}/esearch.fcgi?{params}"))
    return [node.text for node in root.findall(".//IdList/Id") if node.text]


def _doc_summary(doc: ET.Element) -> dict:
    row: dict[str, str] = {}
    for item in doc.findall("Item"):
        name = item.get("Name")
        if not name:
            continue
        if item.get("Type") == "List":
            if name == "Samples":
                row["n_samples"] = str(
                    len(item.findall("Item[@Name='Sample']"))
                )
            continue
        text = (item.text or "").strip()
        if text:
            row[name] = text
    return row


def esummary(uid: str) -> dict:
    params = urllib.parse.urlencode(
        {
            "db": "gds",
            "id": uid,
            "retmode": "xml",
        }
    )
    root = ET.fromstring(_http_get(f"{EUTILS}/esummary.fcgi?{params}"))
    doc = root.find(".//DocSum")
    return _doc_summary(doc) if doc is not None else {}


def esummary_many(uids: list[str]) -> dict[str, dict]:
    if not uids:
        return {}
    params = urllib.parse.urlencode(
        {
            "db": "gds",
            "id": ",".join(uids),
            "retmode": "xml",
        }
    )
    root = ET.fromstring(_http_get(f"{EUTILS}/esummary.fcgi?{params}"))
    rows: dict[str, dict] = {}
    for doc in root.findall(".//DocSum"):
        uid = doc.findtext("Id") or ""
        if uid:
            rows[uid] = _doc_summary(doc)
    return rows


def to_row(
    raw: dict,
    disease: str = "",
    research_direction: str = "",
    database: str = "GEO",
    run_supported: bool | None = None,
) -> dict:
    accession = str(raw.get("Accession", ""))
    match = re.search(r"(GSE\d+)", accession)
    gse = match.group(1) if match else accession
    title = str(raw.get("Title") or raw.get("title") or "")
    summary = str(raw.get("Summary") or raw.get("summary") or "")
    data_type = infer_data_type(f"{title}\n{summary}")
    if run_supported is None:
        run_supported = bool(re.fullmatch(r"GSE\d+", canonical_row_accession(gse)))
    row = {
        "accession": gse,
        "database": database,
        "disease": disease,
        "research_direction": research_direction,
        "title": title,
        "summary": summary,
        "data_type": data_type,
        "organism": str(raw.get("taxon", "")),
        "samples": str(
            raw.get("n_samples")
            or raw.get("Samples")
            or ""
        ),
        "platform": str(raw.get("GPL", "")),
        "date": str(raw.get("PDAT", "")),
        "type": str(
            raw.get("Type")
            or raw.get("gdsType")
            or raw.get("entryType")
            or ""
        ),
        "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}",
        "quality_score": 0.0,
        "run_supported": bool(run_supported),
    }
    row["quality_score"] = quality_score(row)
    return row


def filter_rows(
    rows: list[dict],
    organism: str | None = None,
    keyword: str | None = None,
    data_type: str | None = None,
    min_samples: int | None = None,
    max_samples: int | None = None,
    start_date: str | None = None,
    end_date: str | None = None,
    platform: str | None = None,
    dataset_type: str | None = None,
) -> list[dict]:
    def matches(row: dict) -> bool:
        if organism and organism.lower() not in row["organism"].lower():
            return False
        if keyword:
            haystack = " ".join(
                [row["accession"], row["title"], row["summary"]]
            ).lower()
            if keyword.lower() not in haystack:
                return False
        if data_type and row.get("data_type", "") != data_type:
            return False
        sample_count = _to_int(row.get("samples"))
        if min_samples is not None and (
            sample_count is None or sample_count < min_samples
        ):
            return False
        if max_samples is not None and (
            sample_count is None or sample_count > max_samples
        ):
            return False
        row_date = _parse_date(row.get("date"))
        start = _parse_date(start_date) if start_date else None
        end = _parse_date(end_date) if end_date else None
        if row_date is not None:
            if start is not None and row_date < start:
                return False
            if end is not None and row_date > end:
                return False
        elif start_date or end_date:
            return False
        if (
            platform
            and platform.lower() not in str(row.get("platform", "")).lower()
        ):
            return False
        if (
            dataset_type
            and dataset_type.lower() not in str(row.get("type", "")).lower()
        ):
            return False
        return True

    return [row for row in rows if matches(row)]


def search_geo(
    query: str,
    max_results: int = 20,
    disease: str = "",
    research_direction: str = "",
) -> list[dict]:
    ids = esearch(query, max_results)
    rows: list[dict] = []
    summaries = (
        {ids[0]: esummary(ids[0])}
        if len(ids) == 1
        else esummary_many(ids)
    )
    for uid in ids:
        rows.append(
            to_row(
                summaries.get(uid, {}),
                disease=disease,
                research_direction=research_direction,
                database="GEO",
            )
        )
    time.sleep(0.34)
    return rows


def search_biostudies(
    query: str,
    max_results: int = 20,
    disease: str = "",
    research_direction: str = "",
) -> list[dict]:
    params = urllib.parse.urlencode(
        {
            "query": query,
            "collection": "ArrayExpress",
            "pageSize": max(max_results, 1),
        }
    )
    data = _http_get_json(f"{BIOSTUDIES_API}/search?{params}")
    rows: list[dict] = []
    for hit in data.get("hits") or []:
        accession = str(hit.get("accession") or "").strip().upper()
        if not accession:
            continue
        title = str(hit.get("title") or "")
        content = str(hit.get("content") or "")
        summary = content
        if title and title.lower() in summary.lower():
            idx = summary.lower().index(title.lower())
            summary = (summary[:idx] + summary[idx + len(title):]).strip()
        summary = re.sub(rf"^{re.escape(accession)}\s*", "", summary).strip()
        try:
            files = int(hit.get("files") or 0)
        except (TypeError, ValueError):
            files = 0
        row = _normalize_row(
            {
                "accession": accession,
                "database": "ArrayExpress/BioStudies",
                "disease": disease,
                "research_direction": research_direction,
                "title": title,
                "summary": summary,
                "data_type": infer_data_type(f"{title}\n{content}"),
                "organism": _detect_organism(content),
                "samples": "",
                "platform": "",
                "date": str(hit.get("release_date") or ""),
                "type": _detect_study_type(content),
                "url": f"https://www.ebi.ac.uk/biostudies/studies/{accession}",
                "run_supported": files > 0,
            }
        )
        row["quality_score"] = quality_score(row)
        rows.append(row)
    return rows


def _load_atlas_index(force_refresh: bool = False) -> list[dict]:
    if not force_refresh and ATLAS_CACHE.exists():
        try:
            cached = json.loads(ATLAS_CACHE.read_text(encoding="utf-8"))
            fetched_at = float(cached.get("fetched_at") or 0)
            if time.time() - fetched_at < ATLAS_CACHE_TTL_SECONDS:
                return cached.get("experiments") or []
        except (OSError, ValueError, TypeError):
            pass
    data = _http_get_json(ATLAS_API)
    experiments = data.get("experiments") or []
    ATLAS_CACHE.parent.mkdir(parents=True, exist_ok=True)
    ATLAS_CACHE.write_text(
        json.dumps(
            {"fetched_at": time.time(), "experiments": experiments},
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    return experiments


def search_atlas(
    query: str,
    max_results: int = 20,
    disease: str = "",
    research_direction: str = "",
) -> list[dict]:
    disease_terms = [term.lower() for term in ([disease] + _mapping_values(DISEASE_SYNONYMS, disease)) if term]
    direction_terms = [
        term.lower()
        for term in (
            [research_direction]
            + _mapping_values(DIRECTION_SYNONYMS, research_direction)
        )
        if term
    ]
    rows: list[dict] = []
    for experiment in _load_atlas_index():
        accession = str(experiment.get("experimentAccession") or "").strip().upper()
        description = str(experiment.get("experimentDescription") or "")
        if not accession:
            continue
        text = f"{accession} {description}".lower()
        disease_hit = not disease_terms or any(term in text for term in disease_terms)
        direction_hit = not direction_terms or any(term in text for term in direction_terms)
        if not (disease_hit and direction_hit):
            continue
        technology = experiment.get("technologyType") or []
        if isinstance(technology, str):
            technology = [technology]
        row = _normalize_row(
            {
                "accession": accession,
                "database": "Expression Atlas",
                "disease": disease,
                "research_direction": research_direction,
                "title": description,
                "summary": description,
                "data_type": infer_data_type(description),
                "organism": str(experiment.get("species") or ""),
                "samples": str(experiment.get("numberOfAssays") or ""),
                "platform": "",
                "date": str(experiment.get("lastUpdate") or experiment.get("loadDate") or ""),
                "type": " | ".join(
                    part
                    for part in [
                        str(experiment.get("experimentType") or ""),
                        " | ".join(str(item) for item in technology),
                    ]
                    if part
                ),
                "url": f"https://www.ebi.ac.uk/gxa/experiments/{accession}",
                "run_supported": True,
            }
        )
        row["quality_score"] = quality_score(row)
        rows.append(row)
    return rows


def _prefer_row(new: dict, existing: dict) -> bool:
    new_db = str(new.get("database", ""))
    existing_db = str(existing.get("database", ""))
    if new_db == "GEO" and existing_db != "GEO":
        return True
    if existing_db == "GEO" and new_db != "GEO":
        return False
    for field in ("organism", "samples", "summary"):
        if str(new.get(field, "")) and not str(existing.get(field, "")):
            return True
        if str(existing.get(field, "")) and not str(new.get(field, "")):
            return False
    return bool(new.get("run_supported")) and not bool(existing.get("run_supported"))


def merge_rows(rows: list[dict]) -> list[dict]:
    merged: dict[str, dict] = {}
    for row in rows:
        accession = canonical_row_accession(row.get("accession", ""))
        if not accession:
            continue
        key = accession.upper()
        existing = merged.get(key)
        if existing is None or _prefer_row(row, existing):
            merged[key] = dict(_normalize_row(row), accession=accession)
    result = list(merged.values())
    result.sort(
        key=lambda item: (
            -float(item.get("quality_score") or 0.0),
            -_samples_int(item),
            str(item.get("date", "")),
        )
    )
    return result


def search_datasets(
    query: str,
    max_results: int = 20,
    organism: str | None = None,
    keyword: str | None = None,
    data_type: str | None = None,
    min_samples: int | None = None,
    max_samples: int | None = None,
    start_date: str | None = None,
    end_date: str | None = None,
    platform: str | None = None,
    dataset_type: str | None = None,
    disease: str = "",
    research_direction: str = "",
    databases: list[str] | None = None,
) -> list[dict]:
    if databases is None:
        databases = list(SUPPORTED_DATABASES)
    selected = [str(item).strip().lower() for item in databases if str(item).strip()]
    unknown = [item for item in selected if item not in SUPPORTED_DATABASES]
    if unknown:
        raise ValueError(f"unknown databases: {', '.join(unknown)}")
    rows: list[dict] = []
    if "geo" in selected:
        rows.extend(search_geo(query, max_results, disease, research_direction))
    if "biostudies" in selected:
        rows.extend(search_biostudies(query, max_results, disease, research_direction))
    if "atlas" in selected:
        rows.extend(search_atlas(query, max_results, disease, research_direction))
    rows = filter_rows(
        rows,
        organism=organism,
        keyword=keyword,
        data_type=data_type,
        min_samples=min_samples,
        max_samples=max_samples,
        start_date=start_date,
        end_date=end_date,
        platform=platform,
        dataset_type=dataset_type,
    )
    rows = merge_rows(rows)
    return rows[:max_results]


def build_query(
    disease: str,
    research_direction: str,
    query: str | None = None,
) -> str:
    if query:
        return query
    parts = [disease, research_direction]
    return " ".join(part.strip() for part in parts if part.strip())


def write_outputs(rows: list[dict], out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "dataset_search_results.csv"
    json_path = out_dir / "dataset_search_results.json"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return csv_path, json_path


def download_accessions(
    accessions: list[str],
    root: Path,
    log=print,
) -> dict[str, str]:
    from data.geo_downloader import canonical_accession, ensure_geo_dataset
    from data.biostudies_downloader import ensure_biostudies_dataset

    results: dict[str, str] = {}
    for accession in accessions:
        try:
            acc = canonical_accession(accession)
            if re.fullmatch(r"GSE\d+", acc):
                ensure_geo_dataset(acc, root, log)
            elif re.fullmatch(r"(?:E-[A-Z0-9]+-\d+|S-BSST\d+)", acc):
                ensure_biostudies_dataset(acc, root, log)
            else:
                raise ValueError(
                    "dataset accession must look like GSE125449, "
                    "E-MTAB-1234, or S-BSST123"
                )
            results[accession] = "ok"
        except Exception as exc:  # noqa: BLE001
            results[accession] = f"error: {exc}"
    return results


def _print_table(rows: list[dict]) -> None:
    for index, row in enumerate(rows, start=1):
        print(f"{index:>3}. {row['accession']} | {row['title']}")
        print(f"     organism: {row['organism']} | samples: {row['samples']}")
        print(f"     {row['url']}")


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(
        description=(
            "Search expression datasets via NCBI GEO, EBI ArrayExpress/"
            "BioStudies and Expression Atlas, then optionally download "
            "matching runnable series."
        )
    )
    parser.add_argument(
        "--databases",
        default=",".join(SUPPORTED_DATABASES),
        help=(
            "comma-separated sources: geo, biostudies, atlas "
            f"(default: {', '.join(SUPPORTED_DATABASES)})"
        ),
    )
    parser.add_argument("--query", help="raw search term across selected databases")
    parser.add_argument(
        "--disease",
        help="disease name, e.g. 'liver cancer'",
    )
    parser.add_argument(
        "--research-direction",
        help="research direction, e.g. 'single cell RNA-seq'",
    )
    parser.add_argument("--max-results", type=int, default=20)
    parser.add_argument(
        "--organism",
        help="filter by organism substring, e.g. 'Homo sapiens'",
    )
    parser.add_argument(
        "--keyword",
        help="filter by extra keyword in accession/title/summary",
    )
    parser.add_argument(
        "--data-type",
        choices=["single-cell", "bulk", "other"],
        help="filter by inferred data type",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        help="minimum number of samples",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        help="maximum number of samples",
    )
    parser.add_argument(
        "--start-date",
        help="earliest publication/update date, e.g. 2024-01-01",
    )
    parser.add_argument(
        "--end-date",
        help="latest publication/update date, e.g. 2025-12-31",
    )
    parser.add_argument(
        "--platform",
        help="filter by platform GPL id or platform name substring",
    )
    parser.add_argument(
        "--dataset-type",
        help="filter by GEO dataset type substring, e.g. 'high throughput sequencing'",
    )
    parser.add_argument(
        "--output",
        default=str(APP_ROOT / "data_cache" / "dataset_search"),
        help="directory for search result CSV/JSON",
    )
    parser.add_argument(
        "--download",
        help="comma-separated GSE accessions to download",
    )
    parser.add_argument(
        "--download-top",
        type=int,
        default=0,
        help="download the first N search results",
    )
    parser.add_argument(
        "--download-root",
        default=str(APP_ROOT.parent / "liver_cancer"),
        help="root directory used by ensure_geo_dataset",
    )
    parser.add_argument(
        "--model",
        default="",
        help="ML/DL relevance model (joblib) for reranking results",
    )
    args = parser.parse_args()
    if not args.query and not args.disease:
        parser.error("provide --query or --disease")
    databases = [
        item.strip().lower()
        for item in args.databases.split(",")
        if item.strip()
    ]
    query = build_query(
        args.disease or "",
        args.research_direction or "",
        args.query,
    )

    rows = search_datasets(
        query,
        max_results=args.max_results,
        organism=args.organism,
        keyword=args.keyword,
        data_type=args.data_type,
        min_samples=args.min_samples,
        max_samples=args.max_samples,
        start_date=args.start_date,
        end_date=args.end_date,
        platform=args.platform,
        dataset_type=args.dataset_type,
        disease=args.disease or "",
        research_direction=args.research_direction or "",
        databases=databases,
    )
    if args.model:
        from dataset_search_ml import load_model, rerank

        model = load_model(Path(args.model))
        rows = rerank(
            rows,
            args.disease or "",
            args.research_direction or "",
            model=model,
        )
        print(f"ML reranking applied: {Path(args.model).name}")
    if not rows:
        print("No datasets matched the search criteria.")
        return 0

    csv_path, json_path = write_outputs(rows, Path(args.output))
    print(f"Search query: {query}")
    print(f"Databases: {', '.join(DATABASE_LABELS.get(item, item) for item in databases)}")
    print(f"Matched {len(rows)} datasets.")
    print(f"CSV: {csv_path}")
    print(f"JSON: {json_path}")
    _print_table(rows)

    accessions: list[str] = []
    if args.download:
        accessions.extend(
            acc.strip().upper()
            for acc in args.download.split(",")
            if acc.strip()
        )
    if args.download_top:
        accessions.extend(
            row["accession"] for row in rows[: args.download_top]
        )
    accessions = list(dict.fromkeys(accessions))
    if accessions:
        print(f"Downloading {len(accessions)} dataset(s)...")
        results = download_accessions(accessions, Path(args.download_root))
        out_dir = Path(args.output)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "download_results.json").write_text(
            json.dumps(results, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        failed = [
            acc for acc, status in results.items()
            if status != "ok"
        ]
        if failed:
            print("Failed downloads:")
            for acc in failed:
                print(f"  {acc}: {results[acc]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
