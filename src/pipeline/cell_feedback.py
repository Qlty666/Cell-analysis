#!/usr/bin/env python3
"""Re-analyze single-cell data using virtual knockout and docking results.

The full pipeline currently flows from single-cell analysis into target
prioritization and docking. This module closes the loop by merging the
knockout scores and docking hits into a compact feedback manifest, then
re-loads the Seurat object through an R helper to score target expression,
cell-type enrichment, module activity and refined target support.
"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
from pathlib import Path

import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent.parent
log = logging.getLogger("cell_feedback")

try:
    from common.env import require_rscript as _require_rscript
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from common.env import require_rscript as _require_rscript

FEEDBACK_COLUMNS = [
    "feedback_rank",
    "gene",
    "source",
    "target_class",
    "target_score",
    "knockout_score",
    "docking_status",
    "docking_hits",
    "best_affinity",
    "pdb_id",
    "feedback_score",
    "key_rank",
    "avg_log2fc",
]


def _find_rscript() -> str:
    return _require_rscript(
        "Rscript not found; install R 4.x before cell feedback"
    )


def _num(value, default: float = 0.0) -> float:
    try:
        if value is None or pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _int_value(value, default: int = 0) -> int:
    return int(round(_num(value, float(default))))


def _str_value(value, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    return str(value).strip()


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _read_ko_frame(workdir: Path) -> pd.DataFrame:
    candidates = [
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "data"
        / "fig_52_53_ranked_knockout.csv",
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "data"
        / "fig_52_target_candidates.csv",
    ]
    for path in candidates:
        if path.exists():
            try:
                frame = pd.read_csv(path)
                if not frame.empty and "gene" in frame.columns:
                    return frame
            except Exception as exc:  # noqa: BLE001
                log.warning("could not read knockout CSV %s: %s", path, exc)
    return pd.DataFrame()


def _read_docking_frame(workdir: Path) -> pd.DataFrame:
    path = workdir / "outputs" / "integration" / "docking_targets.csv"
    if not path.exists():
        return pd.DataFrame()
    try:
        frame = pd.read_csv(path)
        return frame if not frame.empty else pd.DataFrame()
    except Exception as exc:  # noqa: BLE001
        log.warning("could not read docking CSV %s: %s", path, exc)
        return pd.DataFrame()


def _read_key_genes(workdir: Path) -> pd.DataFrame:
    path = workdir / "outputs" / "integration" / "key_genes.csv"
    if not path.exists():
        return pd.DataFrame()
    try:
        frame = pd.read_csv(path)
        return frame if not frame.empty else pd.DataFrame()
    except Exception as exc:  # noqa: BLE001
        log.warning("could not read key genes CSV %s: %s", path, exc)
        return pd.DataFrame()


def build_feedback_manifest(
    workdir: Path,
    top_n: int = 12,
    ko_csv: Path | None = None,
    docking_csv: Path | None = None,
    key_genes_csv: Path | None = None,
) -> pd.DataFrame:
    """Merge knockout and docking results into one prioritised manifest."""
    workdir = Path(workdir).resolve()
    records: dict[str, dict] = {}

    ko_frame = (
        pd.read_csv(ko_csv)
        if ko_csv is not None and ko_csv.exists()
        else _read_ko_frame(workdir)
    )
    for _, row in ko_frame.iterrows():
        gene = _str_value(row.get("gene"))
        if not gene:
            continue
        rec = records.setdefault(gene, {"gene": gene, "source": ""})
        rec["source"] = _merge_source(rec.get("source", ""), "knockout")
        rec["target_score"] = max(
            _num(rec.get("target_score")),
            _num(row.get("target_score")),
        )
        rec["knockout_score"] = max(
            _num(rec.get("knockout_score")),
            _num(row.get("knockout_score")),
        )
        if row.get("rank") is not None and pd.notna(row.get("rank")):
            rec["knockout_rank"] = min(
                _int_value(rec.get("knockout_rank"), 999999),
                _int_value(row.get("rank"), 999999),
            )
        if row.get("target_class") is not None and pd.notna(row.get("target_class")):
            rec["target_class"] = _str_value(row.get("target_class"))

    docking_frame = (
        pd.read_csv(docking_csv)
        if docking_csv is not None and docking_csv.exists()
        else _read_docking_frame(workdir)
    )
    for _, row in docking_frame.iterrows():
        gene = _str_value(row.get("gene"))
        if not gene:
            continue
        rec = records.setdefault(gene, {"gene": gene, "source": ""})
        rec["source"] = _merge_source(rec.get("source", ""), "docking")
        rec["docking_status"] = _str_value(row.get("status"))
        rec["docking_hits"] = max(
            _int_value(rec.get("docking_hits")),
            _int_value(row.get("hits")),
        )
        rec["best_affinity"] = _str_value(row.get("best_affinity"))
        rec["pdb_id"] = _str_value(row.get("pdb_id"))

    key_frame = (
        pd.read_csv(key_genes_csv)
        if key_genes_csv is not None and key_genes_csv.exists()
        else _read_key_genes(workdir)
    )
    if not records and not key_frame.empty:
        for _, row in key_frame.iterrows():
            gene = _str_value(row.get("gene"))
            if not gene:
                continue
            rec = records.setdefault(gene, {"gene": gene, "source": ""})
            rec["source"] = _merge_source(rec.get("source", ""), "key_target")
            rec["key_rank"] = _int_value(row.get("rank"))
            rec["avg_log2fc"] = _num(row.get("avg_log2fc"))
            rec["target_class"] = "key_target"

    if not records:
        return pd.DataFrame(columns=FEEDBACK_COLUMNS)

    rows = []
    for gene, rec in records.items():
        target_score = _num(rec.get("target_score"), 0.5)
        knockout_score = _num(rec.get("knockout_score"), 0.5)
        base = target_score if rec.get("target_score") is not None else knockout_score
        bonus = 0.0
        hits = _int_value(rec.get("docking_hits"))
        if rec.get("docking_status") == "ok" and hits > 0:
            bonus = min(0.20, hits * 0.04)
        rows.append(
            {
                "gene": gene,
                "source": rec.get("source", ""),
                "target_class": rec.get("target_class", ""),
                "target_score": target_score,
                "knockout_score": knockout_score,
                "docking_status": rec.get("docking_status", ""),
                "docking_hits": hits,
                "best_affinity": rec.get("best_affinity", ""),
                "pdb_id": rec.get("pdb_id", ""),
                "feedback_score": max(0.0, min(1.0, base + bonus)),
                "key_rank": _int_value(rec.get("key_rank")),
                "avg_log2fc": _num(rec.get("avg_log2fc")),
            }
        )

    frame = pd.DataFrame(rows)
    frame = frame.sort_values(
        ["feedback_score", "docking_hits", "target_score", "knockout_score"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    frame.insert(0, "feedback_rank", pd.RangeIndex(1, len(frame) + 1))
    frame = frame.head(top_n).reset_index(drop=True)
    frame["feedback_rank"] = pd.RangeIndex(1, len(frame) + 1)
    return frame.reindex(columns=FEEDBACK_COLUMNS)


def _merge_source(current: str, new: str) -> str:
    parts = [p for p in current.split("|") if p]
    if new not in parts:
        parts.append(new)
    return "|".join(parts)


def _run_rscript(
    script: Path,
    args: list[str],
    timeout: int,
    cwd: Path,
) -> subprocess.CompletedProcess:
    cmd = [_find_rscript(), str(script), *map(str, args)]
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


def run_cell_feedback(
    workdir: Path,
    single_cell_root: Path,
    top_n: int = 12,
    max_features: int = 8,
    timeout_seconds: int = 3600,
    species: str = "hs",
) -> dict:
    """Write the feedback manifest and run the Seurat-level R analysis."""
    workdir = Path(workdir).resolve()
    single_cell_root = Path(single_cell_root).resolve()
    out_dir = workdir / "outputs" / "integration" / "cell_feedback"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = build_feedback_manifest(workdir, top_n=top_n)
    manifest_path = out_dir / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    if manifest.empty:
        summary = {
            "status": "skipped",
            "reason": "no knockout, docking or key-gene results found",
        }
        _write_json(out_dir / "cell_feedback_summary.json", summary)
        log.info("cell feedback skipped: %s", summary["reason"])
        return summary

    script = APP_ROOT / "src" / "pipeline" / "cell_feedback.R"
    if not script.exists():
        raise RuntimeError(f"cell feedback R script missing: {script}")

    log.info(
        "running cell feedback: %s genes, %s features, single-cell root %s",
        len(manifest),
        max_features,
        single_cell_root,
    )
    proc = _run_rscript(
        script,
        [single_cell_root, out_dir, top_n, max_features, species],
        timeout_seconds,
        APP_ROOT,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "cell feedback R analysis failed:\n"
            + (proc.stderr or proc.stdout)[-3000:]
        )

    summary_path = out_dir / "cell_feedback_summary.json"
    if not summary_path.exists():
        raise RuntimeError(
            "cell feedback R analysis did not write "
            f"{summary_path}\n{proc.stdout[-2000:]}"
        )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        summary = {"status": "failed", "reason": "invalid summary JSON"}
    summary["manifest"] = str(manifest_path)
    _write_json(summary_path, summary)
    log.info(
        "cell feedback complete: status=%s matched=%s",
        summary.get("status"),
        summary.get("genes_matched"),
    )
    return summary
