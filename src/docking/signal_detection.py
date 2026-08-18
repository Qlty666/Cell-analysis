#!/usr/bin/env python3
"""Disproportionality signal detection for pharmacovigilance-style tables.

Implements the four FAERS-style signal measures referenced in the reviewed
articles: ROR, PRR, BCPNN IC and EBGM. The BCPNN and EBGM values use the
standard closed-form approximations; they are suitable for screening and
should be treated as approximations for exploratory use.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import DockingError, write_json


def _safe(x: float, epsilon: float = 0.5) -> float:
    return x + epsilon


def _ci_lower_upper(point: float, se: float) -> tuple[float, float]:
    if point <= 0 or not np.isfinite(point) or not np.isfinite(se):
        return np.nan, np.nan
    log_point = math.log(point)
    lower = math.exp(log_point - 1.96 * se)
    upper = math.exp(log_point + 1.96 * se)
    return lower, upper


def detect_signals(
    df: pd.DataFrame,
    drug_column: str = "drug",
    event_column: str = "event",
    count_column: str | None = None,
    min_count: int = 3,
) -> pd.DataFrame:
    """Build a 2x2 disproportionality table for every drug-event pair."""
    if df.empty:
        raise DockingError("FAERS event table is empty")
    if drug_column not in df.columns or event_column not in df.columns:
        raise DockingError(
            f"FAERS table must contain '{drug_column}' and '{event_column}' columns"
        )
    count_column = count_column or next(
        (
            c
            for c in df.columns
            if str(c).lower() in {"count", "n", "cases", "reports"}
        ),
        None,
    )
    count = (
        df[count_column]
        if count_column and count_column in df.columns
        else pd.Series(1, index=df.index)
    )
    count = pd.to_numeric(count, errors="coerce").fillna(1.0)
    work = pd.DataFrame(
        {
            "drug": df[drug_column].astype(str).str.strip(),
            "event": df[event_column].astype(str).str.strip(),
            "count": count,
        }
    )
    work = work[(work["drug"] != "") & (work["event"] != "")]
    grouped = work.groupby(["drug", "event"], as_index=False)["count"].sum()
    total = float(grouped["count"].sum())
    if total <= 0:
        raise DockingError("FAERS event table has no report counts")

    drug_totals = grouped.groupby("drug", as_index=False)["count"].transform("sum")
    event_totals = grouped.groupby("event", as_index=False)["count"].transform("sum")
    out = grouped.copy()
    out["a"] = out["count"].astype(float)
    out["b"] = (drug_totals - out["count"]).astype(float)
    out["c"] = (event_totals - out["count"]).astype(float)
    out["d"] = (total - drug_totals - event_totals + out["count"]).astype(float)

    a = out["a"]
    b = out["b"]
    c = out["c"]
    d = out["d"]
    n = total
    ror = (a * d) / (b * c)
    ror_se = np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    out["ror"] = ror
    out["ror_lower"], out["ror_upper"] = zip(
        *[_ci_lower_upper(p, se) for p, se in zip(ror, ror_se)]
    )

    prr = (a / (a + b)) / (c / (c + d))
    prr_se = np.sqrt(1.0 / a - 1.0 / (a + b) + 1.0 / c - 1.0 / (c + d))
    out["prr"] = prr
    out["prr_lower"], out["prr_upper"] = zip(
        *[_ci_lower_upper(p, se) for p, se in zip(prr, prr_se)]
    )

    chi_sq = (
        n
        * np.power(a * d - b * c, 2.0)
        / ((a + b) * (c + d) * (a + c) * (b + d))
    )
    out["chi_square"] = chi_sq

    expected = ((a + b) * (a + c)) / n
    ic = np.log2((a + 0.5) / (expected + 0.5))
    ic025 = ic - 3.3 / np.sqrt(a)
    out["ic"] = ic
    out["ic_lower"] = ic025

    rr = a / expected
    ebgm = rr
    ebgm05 = np.exp(np.log(rr + 1e-9) - 1.96 * np.sqrt(1.0 / a + 1.0 / expected))
    out["ebgm"] = ebgm
    out["ebgm05"] = ebgm05

    out["signal_ror"] = (a >= min_count) & (out["ror_lower"] > 1.0)
    out["signal_prr"] = (a >= min_count) & (out["prr"] >= 2.0) & (chi_sq >= 4.0)
    out["signal_bcpnn"] = (a >= min_count) & (ic025 > 0.0)
    out["signal_ebgm"] = (a >= min_count) & (ebgm05 > 1.0)
    out["signal"] = (
        out["signal_ror"]
        & out["signal_prr"]
        & (out["signal_bcpnn"] | out["signal_ebgm"])
    )
    out = out.sort_values(
        ["signal", "ror_lower", "a"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return out


def run_faers(cfg, log) -> dict:
    """Run FAERS-style disproportionality signal detection from config."""
    section = cfg.data.get("faers", {}) or {}
    input_csv = section.get("input_csv")
    if not input_csv:
        raise DockingError("faers.input_csv is required")
    path = cfg._resolve(input_csv, cfg.workdir)
    if not path.exists():
        raise DockingError(f"FAERS event table not found: {path}")
    df = pd.read_csv(path)
    signals = detect_signals(
        df,
        drug_column=section.get("drug_column") or "drug",
        event_column=section.get("event_column") or "event",
        count_column=section.get("count_column"),
        min_count=int(section.get("min_count", 3)),
    )
    out_dir = cfg._resolve(
        section.get("output_dir") or "outputs/run_001/faers",
        cfg.workdir,
    )
    (out_dir / "data").mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "data" / "faers_signals.csv"
    signals.to_csv(csv_path, index=False)

    html_path = out_dir / "data" / "faers_signals.html"
    top = signals.head(100)
    rows = "".join(
        f"<tr><td>{row['drug']}</td><td>{row['event']}</td>"
        f"<td>{int(row['a'])}</td><td>{row['ror']:.2f}</td>"
        f"<td>{row['prr']:.2f}</td><td>{row['ic']:.2f}</td>"
        f"<td>{row['ebgm']:.2f}</td><td>{row['signal']}</td></tr>"
        for row in top.to_dict("records")
    )
    html = f"""<!doctype html>
<html lang="zh"><head><meta charset="utf-8"><title>FAERS signals</title></head>
<body><h1>FAERS Disproportionality Signals</h1>
<p>{len(signals)} drug-event pairs, {int(signals['signal'].sum())} signals</p>
<table border="1" cellpadding="4" cellspacing="0">
<thead><tr><th>Drug</th><th>Event</th><th>Count</th><th>ROR</th>
<th>PRR</th><th>IC</th><th>EBGM</th><th>Signal</th></tr></thead>
<tbody>{rows}</tbody></table></body></html>
"""
    html_path.write_text(html, encoding="utf-8")
    summary = {
        "pairs": int(len(signals)),
        "signals": int(signals["signal"].sum()),
        "min_count": int(section.get("min_count", 3)),
        "top_drugs": (
            signals[signals["signal"]]
            .groupby("drug")["a"]
            .sum()
            .sort_values(ascending=False)
            .head(10)
            .to_dict()
        ),
        "output_csv": str(csv_path),
        "output_html": str(html_path),
    }
    write_json(out_dir / "faers_summary.json", summary)
    log.info(
        "FAERS signal detection complete: %s pairs, %s signals -> %s",
        summary["pairs"],
        summary["signals"],
        out_dir,
    )
    return summary
