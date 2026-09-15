#!/usr/bin/env python3
"""Collect, score and export target evidence from versioned sources."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from evidence import EvidenceContext, EvidenceHub, SQLiteEvidenceStore  # noqa: E402


def _comma_values(value: str | None) -> list[str]:
    if not value:
        return []
    return [
        item.strip()
        for item in value.replace("\n", ",").split(",")
        if item.strip()
    ]


def _json_object(value: str | None, label: str) -> dict:
    if not value:
        return {}
    try:
        payload = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} must be valid JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must be a JSON object")
    return payload


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect and score target evidence from multiple databases."
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "evidence_sources.json"),
    )
    parser.add_argument("--database", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--pubchem-cid", default="")
    parser.add_argument("--chembl-id", default="")
    parser.add_argument("--smiles", default="")
    parser.add_argument("--isomeric-smiles", default="")
    parser.add_argument("--inchi-key", default="")
    parser.add_argument("--disease", default="")
    parser.add_argument("--disease-id", default="")
    parser.add_argument("--targets", default="")
    parser.add_argument("--ensembl-map", default="")
    parser.add_argument("--max-records", type=int, default=None)
    parser.add_argument("--timeout", type=int, default=None)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--benchmark-positive", default="")
    parser.add_argument("--benchmark-negative", default="")
    parser.add_argument("--benchmark-top-n", type=int, default=20)
    parser.add_argument("--benchmark-source", default="")
    parser.add_argument("--benchmark-version", default="")
    parser.add_argument("--benchmark-independent", action="store_true")
    parser.add_argument("--benchmark-exclude-sources", default="")
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    import logging

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    try:
        config_path = Path(args.config).expanduser().resolve()
        config = json.loads(config_path.read_text(encoding="utf-8"))
        output = Path(args.output).expanduser().resolve()
        database = (
            Path(args.database).expanduser().resolve()
            if args.database
            else output / "evidence.sqlite"
        )
        limits = config.get("limits") or {}
        context = EvidenceContext(
            compound={
                "pubchem_cid": args.pubchem_cid,
                "chembl_id": args.chembl_id,
                "canonical_smiles": args.smiles,
                "isomeric_smiles": args.isomeric_smiles,
                "inchi_key": args.inchi_key,
            },
            disease={
                "name": args.disease,
                "id": args.disease_id,
            },
            target_symbols=_comma_values(args.targets),
            ensembl_ids={
                str(key).upper(): str(value)
                for key, value in _json_object(
                    args.ensembl_map,
                    "--ensembl-map",
                ).items()
            },
            cache_dir=output,
            max_records_per_source=int(
                args.max_records
                if args.max_records is not None
                else limits.get("max_records_per_source", 1000)
            ),
            timeout_seconds=int(
                args.timeout
                if args.timeout is not None
                else limits.get("timeout_seconds", 120)
            ),
            allow_network=not args.offline,
            source_options={
                name: dict(options or {})
                for name, options in (config.get("sources") or {}).items()
            },
        )
        if args.strict:
            config["strict"] = True
        with SQLiteEvidenceStore(database) as store:
            hub = EvidenceHub.from_config(
                database,
                config,
                store=store,
            )
            status = (
                store.source_runs()
                if args.export_only
                else hub.collect(context)
            )
            scored = hub.score(
                target_symbols=context.target_symbols or None,
                benchmark_positives=_comma_values(args.benchmark_positive),
                benchmark_negatives=_comma_values(args.benchmark_negative),
                benchmark_top_n=args.benchmark_top_n,
                benchmark_source=args.benchmark_source,
                benchmark_version=args.benchmark_version,
                benchmark_independent=args.benchmark_independent,
                benchmark_exclude_sources=_comma_values(
                    args.benchmark_exclude_sources
                ),
            )
            paths = hub.export(
                output,
                context=context,
                source_status=status,
                scored=scored,
            )
        print(
            json.dumps(
                {
                    "status": "completed",
                    "database": str(database),
                    "outputs": paths,
                    "benchmark": scored.get("benchmark"),
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
