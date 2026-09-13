"""QC metric collection and gate evaluation for the integrated pipeline."""

from __future__ import annotations

import logging
import math
from pathlib import Path

import pandas as pd

from docking.utils import write_json

from .stage_paths import _integration_dir, _read_json

log = logging.getLogger("full_pipeline")


def _as_float(value) -> float | None:
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError) as exc:
        log.debug("non-numeric value %r: %s", value, exc)
        return None


def _min_affinity_text(values) -> str:
    """Return the numeric minimum of an affinity column as text."""
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return str(numeric.min()) if not numeric.empty else ""


def collect_qc_metrics(single_cell_root: Path, workdir: Path) -> dict:
    """Aggregate single-cell and downstream metrics for the QC gate."""
    data_dir = single_cell_root / "results" / "data"
    metrics: dict = {
        "single_cell": _read_json(single_cell_root / "results" / "summary.json"),
        "qc_thresholds": [],
        "doublet_rate": None,
        "doublet_rate_by_sample": [],
        "pseudobulk_used": None,
        "pseudobulk_warning": "",
        "knockout": {},
        "docking": {},
    }

    qc_threshold_csv = data_dir / "01_qc" / "fig_01_qc_thresholds.csv"
    if qc_threshold_csv.exists():
        try:
            threshold_frame = pd.read_csv(qc_threshold_csv)
            metrics["qc_thresholds"] = threshold_frame.to_dict(
                orient="records"
            )
        except (OSError, ValueError) as exc:
            log.warning(
                "could not read QC threshold CSV %s: %s",
                qc_threshold_csv,
                exc,
            )

    doublet_csv = (
        data_dir
        / "08_publication"
        / "fig_29_doublet_rate_by_sample.csv"
    )
    if not doublet_csv.exists():
        doublet_csv = (
            data_dir
            / "08_publication"
            / "fig_29_doublet_rate_sample.csv"
        )
    if doublet_csv.exists():
        try:
            frame = pd.read_csv(doublet_csv)
            rates = pd.to_numeric(
                frame.get("doublet_rate"),
                errors="coerce",
            )
            metrics["doublet_rate_by_sample"] = frame.to_dict(
                orient="records"
            )
            if rates.notna().any():
                metrics["doublet_rate"] = float(rates.mean())
        except (OSError, ValueError) as exc:
            log.warning(
                "could not read doublet-rate CSV %s: %s",
                doublet_csv,
                exc,
            )
    else:
        doublet_tbl = (
            data_dir / "02_doublets" / "fig_02_doublet_results.csv"
        )
        if doublet_tbl.exists():
            try:
                frame = pd.read_csv(doublet_tbl)
                if (
                    "doublet_call" in frame.columns
                    and "sample" in frame.columns
                ):
                    grouped = (
                        frame.groupby("sample")["doublet_call"]
                        .agg(
                            n_cells="count",
                            n_doublets=lambda series: int(
                                (
                                    series.astype(str).str.lower()
                                    == "doublet"
                                ).sum()
                            ),
                        )
                        .reset_index()
                    )
                    grouped["doublet_rate"] = (
                        grouped["n_doublets"] / grouped["n_cells"]
                    )
                    metrics["doublet_rate_by_sample"] = grouped.to_dict(
                        orient="records"
                    )
                    rates = pd.to_numeric(
                        grouped["doublet_rate"],
                        errors="coerce",
                    )
                    if rates.notna().any():
                        metrics["doublet_rate"] = float(rates.mean())
            except (OSError, ValueError) as exc:
                log.warning(
                    "could not read doublet results CSV %s: %s",
                    doublet_tbl,
                    exc,
                )

    warn_path = data_dir / "pseudobulk_warning.txt"
    if warn_path.exists():
        try:
            metrics["pseudobulk_warning"] = warn_path.read_text(
                encoding="utf-8",
                errors="replace",
            ).strip()
            metrics["pseudobulk_used"] = False
        except OSError as exc:
            log.warning(
                "could not read pseudobulk warning %s: %s",
                warn_path,
                exc,
            )

    integration = _integration_dir(workdir)
    ko_summary = _read_json(integration / "knockout_summary.json")
    metrics["knockout"] = ko_summary.get("knockout") or {}
    metrics["docking"] = _read_json(
        integration / "docking_summary.json"
    )
    return metrics


def evaluate_qc_gate(metrics: dict, config: dict) -> dict:
    """Turn collected metrics into an explicit pass/warn/fail gate."""
    if not config.get("enabled", True):
        return {
            "status": "skipped",
            "checks": [],
            "summary": "QC gate disabled by configuration",
        }

    single_cell = metrics.get("single_cell") or {}
    checks: list[dict] = []
    fail_on_missing = bool(
        config.get("fail_on_missing_metrics", False)
    )

    def check(
        name: str,
        value,
        threshold,
        message: str,
        missing_message: str,
        below_is_bad: bool = True,
    ) -> None:
        number = _as_float(threshold)
        if number is None or number == 0:
            return
        current = _as_float(value)
        if current is None:
            checks.append(
                {
                    "name": name,
                    "level": "fail" if fail_on_missing else "warn",
                    "ok": False,
                    "message": missing_message,
                }
            )
            return
        ok = current >= number if below_is_bad else current <= number
        checks.append(
            {
                "name": name,
                "level": "fail" if not ok else "pass",
                "ok": ok,
                "message": message.format(
                    value=current,
                    threshold=number,
                ),
            }
        )

    check(
        "min_cells_after_qc",
        single_cell.get("n_cells_after_qc"),
        config.get("min_cells_after_qc"),
        "cells after QC = {value:.0f} (min {threshold:.0f})",
        "n_cells_after_qc metric is missing",
    )
    check(
        "min_genes",
        single_cell.get("n_genes"),
        config.get("min_genes"),
        "genes = {value:.0f} (min {threshold:.0f})",
        "n_genes metric is missing",
    )
    check(
        "min_deg_genes",
        single_cell.get("deg_total"),
        config.get("min_deg_genes"),
        "DEGs = {value:.0f} (min {threshold:.0f})",
        "deg_total metric is missing",
    )
    check(
        "max_doublet_rate",
        metrics.get("doublet_rate"),
        config.get("max_doublet_rate"),
        "doublet rate = {value:.3f} (max {threshold:.3f})",
        "doublet rate metric is missing",
        below_is_bad=False,
    )

    if (
        config.get("require_pseudobulk")
        and metrics.get("pseudobulk_used") is False
    ):
        checks.append(
            {
                "name": "require_pseudobulk",
                "level": "fail",
                "ok": False,
                "message": metrics.get("pseudobulk_warning")
                or "DE fell back to cells-as-replicates",
            }
        )

    failed = [
        check_item
        for check_item in checks
        if check_item["level"] == "fail" and not check_item["ok"]
    ]
    warned = [
        check_item
        for check_item in checks
        if check_item["level"] == "warn" and not check_item["ok"]
    ]
    status = "fail" if failed else ("warn" if warned else "pass")
    return {
        "status": status,
        "checks": checks,
        "summary": (
            f"{len(failed)} failed, {len(warned)} warned"
            if checks
            else "no QC thresholds configured"
        ),
    }


def write_qc_metrics(
    workdir: Path,
    single_cell_root: Path,
    qc_config: dict,
) -> dict:
    """Collect and evaluate QC metrics, writing qc_metrics.json."""
    metrics = collect_qc_metrics(single_cell_root, workdir)
    metrics["qc_gate"] = evaluate_qc_gate(metrics, qc_config)
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(out_dir / "qc_metrics.json", metrics)
    return metrics
