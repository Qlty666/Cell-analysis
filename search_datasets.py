#!/usr/bin/env python3
"""Search GEO datasets and optionally download matching series."""

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

APP_ROOT = Path(__file__).resolve().parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "Mozilla/5.0 (liver-cancer-pipeline; dataset-search)"

CSV_COLUMNS = [
    "accession",
    "disease",
    "research_direction",
    "title",
    "summary",
    "organism",
    "samples",
    "platform",
    "date",
    "type",
    "url",
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
) -> dict:
    accession = str(raw.get("Accession", ""))
    match = re.search(r"(GSE\d+)", accession)
    gse = match.group(1) if match else accession
    return {
        "accession": gse,
        "disease": disease,
        "research_direction": research_direction,
        "title": str(raw.get("Title") or raw.get("title") or ""),
        "summary": str(raw.get("Summary") or raw.get("summary") or ""),
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
    }


def filter_rows(
    rows: list[dict],
    organism: str | None = None,
    keyword: str | None = None,
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
        return True

    return [row for row in rows if matches(row)]


def search_datasets(
    query: str,
    max_results: int = 20,
    organism: str | None = None,
    keyword: str | None = None,
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
            )
        )
    time.sleep(0.34)
    return filter_rows(rows, organism, keyword)


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
    from data.geo_downloader import ensure_geo_dataset

    results: dict[str, str] = {}
    for accession in accessions:
        try:
            ensure_geo_dataset(accession, root, log)
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
            "Search GEO datasets via NCBI E-utilities and optionally "
            "download matching GSE series."
        )
    )
    parser.add_argument("--query", help="raw GEO search term")
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
        disease=args.disease or "",
        research_direction=args.research_direction or "",
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
