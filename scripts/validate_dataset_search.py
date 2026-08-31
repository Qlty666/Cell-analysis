#!/usr/bin/env python3
"""Randomly validate GEO dataset search relevance across diseases."""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
APP_ROOT = SCRIPTS_DIR.parent
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

import search_datasets  # noqa: E402

DISEASES = search_datasets.DISEASE_SYNONYMS
DIRECTIONS = search_datasets.DIRECTION_SYNONYMS


def _mapping_values(mapping: dict, key: str) -> list[str]:
    return search_datasets._mapping_values(mapping, key)


def is_relevant(row: dict, disease: str, direction: str) -> bool:
    text = " ".join(
        [str(row.get("accession", "")), str(row.get("title", "")), str(row.get("summary", ""))]
    ).lower()
    disease_terms = [disease.lower()] + [
        term.lower()
        for term in _mapping_values(DISEASES, disease)
    ]
    direction_terms = [direction.lower()] + [
        term.lower()
        for term in _mapping_values(DIRECTIONS, direction)
    ]
    disease_hit = any(term in text for term in disease_terms)
    direction_hit = any(term in text for term in direction_terms)
    return disease_hit and direction_hit


def relevant_rows(rows: list[dict], disease: str, direction: str) -> list[dict]:
    return [
        row for row in rows
        if is_relevant(row, disease, direction)
    ]


def _labeled_samples(
    rows: list[dict],
    disease: str,
    direction: str,
    manual_labels: dict | None = None,
) -> list[dict]:
    samples: list[dict] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        accession = str(row.get("accession", ""))
        key = (
            accession,
            disease.lower(),
            direction.lower(),
        )
        if key in seen:
            continue
        seen.add(key)
        accession_upper = accession.upper()
        manual_key = (accession_upper, disease.lower(), direction.lower())
        manual_label = manual_labels.get(manual_key) if manual_labels else None
        label_source = "manual" if manual_label is not None else "heuristic"
        label = (
            manual_label
            if manual_label is not None
            else (1 if is_relevant(row, disease, direction) else 0)
        )
        samples.append(
            {
                "disease": disease,
                "research_direction": direction,
                "accession": accession,
                "title": str(row.get("title", "")),
                "summary": str(row.get("summary", "")),
                "organism": str(row.get("organism", "")),
                "label": label,
                "label_source": label_source,
            }
        )
    return samples


def run_round(
    disease: str,
    direction: str,
    max_results: int,
    expand: bool,
    fallback_max_results: int = 15,
    model: dict | None = None,
    top_size: int | None = None,
    fetch_size: int | None = None,
    databases: list[str] | None = None,
    manual_labels: dict | None = None,
) -> dict:
    started = time.time()
    top_size = top_size or max_results
    fetch_size = fetch_size or max_results
    query = search_datasets.build_query(disease, direction)
    combined = search_datasets.search_datasets(
        query,
        max_results=fetch_size,
        disease=disease,
        research_direction=direction,
        databases=databases,
    )
    if model is not None:
        from dataset_search_ml import rerank

        combined = rerank(
            combined,
            disease,
            direction,
            model=model,
        )[:top_size]
    else:
        combined = combined[:top_size]
    hits = relevant_rows(combined, disease, direction)
    samples = _labeled_samples(
        combined,
        disease,
        direction,
        manual_labels,
    )
    expanded = False
    if not hits and expand:
        fallback_queries = [
            ("disease", search_datasets.build_expanded_query(disease, "")),
            ("direction", search_datasets.build_expanded_query("", direction)),
            (
                "combined_expanded",
                search_datasets.build_expanded_query(disease, direction),
            ),
        ]
        for _label, fallback_query in fallback_queries:
            if hits:
                break
            fallback_rows = search_datasets.search_datasets(
                fallback_query,
                max_results=max(fallback_max_results, fetch_size),
                disease=disease,
                research_direction=direction,
                databases=databases,
            )
            if model is not None:
                from dataset_search_ml import rerank

                fallback_rows = rerank(
                    fallback_rows,
                    disease,
                    direction,
                    model=model,
                )[:top_size]
            else:
                fallback_rows = fallback_rows[:top_size]
            samples.extend(
                _labeled_samples(
                    fallback_rows,
                    disease,
                    direction,
                    manual_labels,
                )
            )
            fallback_hits = relevant_rows(fallback_rows, disease, direction)
            if fallback_hits:
                expanded = True
                hits = fallback_hits

    first = hits[0] if hits else {}
    manual_round = any(
        sample.get("label_source") == "manual" for sample in samples
    )
    manual_found = any(
        sample.get("label_source") == "manual"
        and int(sample.get("label") or 0) == 1
        for sample in samples
    )
    return {
        "disease": disease,
        "research_direction": direction,
        "query": query,
        "combined_results": len(combined),
        "relevant_combined": len(
            relevant_rows(combined, disease, direction)
        ),
        "expanded": expanded,
        "found": manual_found if manual_round else bool(hits),
        "found_source": "manual" if manual_round else "heuristic",
        "first_accession": str(first.get("accession", "")),
        "first_title": str(first.get("title", "")),
        "elapsed_seconds": round(time.time() - started, 2),
        "samples": samples,
    }


def write_report(records: list[dict], summary: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"validation_{len(records)}_rounds.csv"
    json_path = out_dir / f"validation_{len(records)}_rounds.json"
    columns = [
        "disease",
        "research_direction",
        "query",
        "combined_results",
        "relevant_combined",
        "expanded",
        "found",
        "found_source",
        "first_accession",
        "first_title",
        "elapsed_seconds",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(records)
    json_path.write_text(
        json.dumps(
            {"summary": summary, "rounds": records},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Validation CSV: {csv_path}")
    print(f"Validation JSON: {json_path}")


def write_training_samples(
    samples: list[dict],
    out_dir: Path,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "training_samples.csv"
    columns = [
        "disease",
        "research_direction",
        "accession",
        "title",
        "summary",
        "organism",
        "label",
        "label_source",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(samples)
    return path


def main() -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(
        description="Run repeated random GEO dataset search validation."
    )
    parser.add_argument("--rounds", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--max-results", type=int, default=5)
    parser.add_argument("--fallback-max-results", type=int, default=15)
    parser.add_argument(
        "--model",
        default="",
        help="ML/DL relevance model (joblib) for reranking",
    )
    parser.add_argument(
        "--rerank-top",
        type=int,
        default=0,
        help="keep top N rows after ML reranking (default: --max-results)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(APP_ROOT / "data_cache" / "dataset_search"),
    )
    parser.add_argument(
        "--databases",
        default="geo",
        help=(
            "comma-separated dataset sources: geo, biostudies, atlas "
            "(default: geo for GEO-only validation)"
        ),
    )
    parser.add_argument(
        "--no-expand",
        action="store_true",
        help="disable disease-only/direction-only fallback search",
    )
    parser.add_argument(
        "--manual-labels",
        default="",
        help=(
            "CSV with accession,disease,research_direction,label columns "
            "from manual review; labels override heuristic labels"
        ),
    )
    args = parser.parse_args()
    model = None
    if args.model:
        from dataset_search_ml import load_model

        model = load_model(Path(args.model))
    databases = [
        item.strip().lower()
        for item in args.databases.split(",")
        if item.strip()
    ]

    manual_labels: dict[tuple[str, str, str], int] = {}
    if args.manual_labels:
        manual_path = Path(args.manual_labels)
        with manual_path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                accession = str(row.get("accession", "")).strip().upper()
                disease_key = str(row.get("disease", "")).strip().lower()
                direction_key = str(
                    row.get("research_direction", "")
                ).strip().lower()
                try:
                    label = int(float(str(row.get("label", ""))))
                except (TypeError, ValueError):
                    continue
                if not accession or not disease_key or not direction_key:
                    continue
                manual_labels[(accession, disease_key, direction_key)] = (
                    1 if label else 0
                )
        print(f"Loaded {len(manual_labels)} manual labels from {manual_path}")

    rng = random.Random(args.seed)
    pairs = [
        (disease, direction)
        for disease in DISEASES
        for direction in DIRECTIONS
    ]
    if len(pairs) >= args.rounds:
        selected = rng.sample(pairs, args.rounds)
    else:
        selected = [rng.choice(pairs) for _ in range(args.rounds)]

    records: list[dict] = []
    all_samples: list[dict] = []
    for index, (disease, direction) in enumerate(selected, start=1):
        record = run_round(
            disease,
            direction,
            max_results=args.max_results,
            expand=not args.no_expand,
            fallback_max_results=args.fallback_max_results,
            model=model,
            top_size=args.rerank_top or args.max_results,
            fetch_size=(
                args.max_results * 3
                if model is not None
                else args.max_results
            ),
            databases=databases,
            manual_labels=manual_labels,
        )
        all_samples.extend(record.pop("samples", []))
        records.append(record)
        if index == 1 or index % 10 == 0 or index == args.rounds:
            print(
                f"[{index}/{args.rounds}] {disease} + {direction} -> "
                f"found={record['found']} expanded={record['expanded']} "
                f"({record['elapsed_seconds']}s)"
            )

    found = sum(1 for record in records if record["found"])
    expanded = sum(1 for record in records if record["expanded"])
    total_results = sum(record["combined_results"] for record in records)
    manual_labeled = sum(
        1
        for sample in all_samples
        if sample.get("label_source") == "manual"
    )
    manual_rounds = sum(
        1 for record in records if record.get("found_source") == "manual"
    )
    manual_found = sum(
        1
        for record in records
        if record.get("found_source") == "manual" and record["found"]
    )
    summary = {
        "rounds": len(records),
        "seed": args.seed,
        "label_mode": "manual" if manual_labeled else "heuristic",
        "manual_labeled_samples": manual_labeled,
        "manual_rounds": manual_rounds,
        "manual_found_rate": (
            round(manual_found / manual_rounds, 4) if manual_rounds else None
        ),
        "found": found,
        "found_rate": round(found / len(records), 4) if records else 0.0,
        "expanded_rounds": expanded,
        "total_combined_results": total_results,
        "avg_combined_results": round(
            total_results / len(records), 2
        ) if records else 0.0,
        "finished_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    write_report(records, summary, Path(args.output_dir))
    training_path = write_training_samples(
        all_samples,
        Path(args.output_dir),
    )
    print(f"Training samples: {training_path} ({len(all_samples)} rows)")
    if manual_rounds:
        print(
            f"Manual-labeled rounds: {manual_found}/{manual_rounds} found "
            f"({summary['manual_found_rate']:.1%})"
        )
    if manual_labeled == 0:
        print(
            "WARNING: training labels are heuristic synonyms; ML metrics "
            "will not validate real relevance without --manual-labels"
        )
    print(
        f"Summary: {found}/{len(records)} rounds found relevant datasets "
        f"({summary['found_rate']:.1%}); expanded {expanded} rounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
