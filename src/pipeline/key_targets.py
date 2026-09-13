"""Key-gene ranking from single-cell differential-expression outputs."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from docking.utils import write_json

from .differential import _bh_adjust
from .errors import IntegrationError

log = logging.getLogger("full_pipeline")

DEFAULT_GENE_BLACKLIST = [
    r"^RPL",
    r"^RPS",
    r"^MRPL",
    r"^MRPS",
    r"^MT-",
    r"^MTRNR",
    r"^SNORD",
    r"^SCGB",
    r"^IGH",
    r"^IGK",
    r"^IGL",
    r"^TRA",
    r"^TRB",
    r"^TRG",
    r"^HLA-D",
    r"^LINC",
    r"^RP[0-9]",
    r"^AC[0-9]",
    r"^AL[0-9]",
]


def extract_key_genes(
    single_cell_root: Path,
    out_dir: Path,
    top_n: int = 50,
    keep_all: bool = False,
    blacklist_patterns: list[str] | None = None,
    advanced_priority_csv: str | Path | None = None,
) -> pd.DataFrame:
    """Rank significant DEGs into a compact key-gene table."""
    data_dir = single_cell_root / "results" / "data"
    significant_path = data_dir / "05_deg" / "fig_09_deg_significant.csv"
    all_path = data_dir / "05_deg" / "fig_08_deg_all.csv"
    deg_path = significant_path if significant_path.exists() else all_path
    if not deg_path.exists():
        raise IntegrationError(f"DEG table not found under {data_dir}")

    frame = pd.read_csv(deg_path)
    if frame.empty and deg_path == significant_path and all_path.exists():
        log.warning(
            "significant DEG table is empty (%s); "
            "falling back to the full DEG table",
            significant_path,
        )
        deg_path = all_path
        frame = pd.read_csv(all_path)
    if frame.empty:
        raise IntegrationError(f"DEG table is empty: {deg_path}")

    rename = {
        "avg_log2FC": "avg_log2fc",
        "log2FoldChange": "avg_log2fc",
        "pvalue": "p_val",
        "padj": "p_val_adj",
    }
    frame = frame.rename(columns=rename)
    # Duplicate source columns would otherwise return DataFrames on lookup.
    frame = frame.loc[:, ~frame.columns.duplicated()]
    if "gene" not in frame.columns:
        raise IntegrationError(f"DEG table has no gene column: {deg_path}")

    if "p_val_adj" not in frame.columns:
        if "p_val" not in frame.columns:
            raise IntegrationError(
                f"DEG table has neither p_val_adj nor p_val: {deg_path}"
            )
        frame["p_val_adj"] = _bh_adjust(
            pd.to_numeric(frame["p_val"], errors="coerce").fillna(1.0)
        )
        log.warning(
            "DEG table %s has no p_val_adj column; computed "
            "BH-adjusted p-values from p_val",
            deg_path,
        )

    if "significant" in frame.columns:
        flag = (
            frame["significant"]
            .astype(str)
            .str.strip()
            .str.upper()
            .isin(["TRUE", "1", "YES"])
        )
        frame = frame[flag]
    if "direction" in frame.columns and not keep_all:
        frame = frame[
            frame["direction"].astype(str).str.strip().isin(["Up", "Down"])
        ]

    if "avg_log2fc" not in frame.columns:
        raise IntegrationError(f"DEG table has no log2FC column: {deg_path}")

    frame = frame.copy()
    frame["gene"] = frame["gene"].astype(str)
    frame["avg_log2fc"] = pd.to_numeric(
        frame["avg_log2fc"],
        errors="coerce",
    )
    frame["p_val_adj"] = pd.to_numeric(
        frame["p_val_adj"],
        errors="coerce",
    )
    frame["abs_log2fc"] = frame["avg_log2fc"].abs()
    advanced_priority = str(
        advanced_priority_csv
        or os.environ.get("LIVER_ADVANCED_PRIORITY_CSV", "")
    ).strip()
    if advanced_priority:
        priority_path = Path(advanced_priority).expanduser()
        if not priority_path.exists():
            log.warning(
                "advanced priority CSV not found; keeping default DEG "
                "ranking: %s",
                priority_path,
            )
        else:
            try:
                priority = pd.read_csv(priority_path)
                gene_col = next(
                    (
                        column
                        for column in ("gene", "symbol", "target")
                        if column in priority.columns
                    ),
                    None,
                )
                if gene_col and "priority_score" in priority.columns:
                    priority = priority[
                        [gene_col, "priority_score"]
                    ].rename(columns={gene_col: "gene"})
                    priority["gene"] = priority["gene"].astype(str)
                    priority["priority_score"] = pd.to_numeric(
                        priority["priority_score"],
                        errors="coerce",
                    )
                    frame = frame.merge(
                        priority.drop_duplicates("gene", keep="first"),
                        on="gene",
                        how="left",
                    )
                    frame["advanced_priority_score"] = frame[
                        "priority_score"
                    ]
                    frame = frame.drop(columns=["priority_score"])
                else:
                    log.warning(
                        "advanced priority CSV has no gene/priority_score "
                        "columns: %s",
                        priority_path,
                    )
            except Exception as exc:
                log.warning(
                    "could not read advanced priority CSV %s: %s",
                    priority_path,
                    exc,
                )
    if "advanced_priority_score" in frame.columns:
        frame = frame.sort_values(
            ["advanced_priority_score", "p_val_adj", "abs_log2fc"],
            ascending=[False, True, False],
            na_position="last",
        ).reset_index(drop=True)
    else:
        frame = frame.sort_values(
            ["p_val_adj", "abs_log2fc"],
            ascending=[True, False],
            na_position="last",
        ).reset_index(drop=True)
    frame = frame.drop_duplicates(
        subset=["gene"],
        keep="first",
    ).reset_index(drop=True)
    frame["deg_rank"] = np.arange(1, len(frame) + 1)

    total_before = len(frame)
    if not keep_all:
        patterns = blacklist_patterns or DEFAULT_GENE_BLACKLIST
        if patterns:
            expr = re.compile("|".join(patterns), re.IGNORECASE)
            frame = frame[~frame["gene"].str.match(expr)]
    frame = frame.reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)

    ml_path = data_dir / "07_ml" / "fig_24_ml_feature_importance.csv"
    if ml_path.exists():
        try:
            ml = pd.read_csv(ml_path, index_col=0)
            ml.index = ml.index.astype(str)
            ml_values = ml.iloc[:, 0].astype(float)
            frame["ml_importance"] = (
                frame["gene"].map(ml_values).fillna(0.0)
            )
        except Exception as exc:
            log.warning(
                "could not read ML feature importance %s: %s",
                ml_path,
                exc,
            )
            frame["ml_importance"] = 0.0
    else:
        frame["ml_importance"] = 0.0

    out_cols = [
        "rank",
        "gene",
        "direction",
        "avg_log2fc",
        "p_val_adj",
        "pct.1",
        "pct.2",
        "deg_rank",
        "ml_importance",
        "advanced_priority_score",
    ]
    for col in out_cols:
        if col not in frame.columns:
            frame[col] = ""
    frame = frame.head(top_n)[out_cols].reset_index(drop=True)
    frame["source"] = "DEG"

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "key_genes.csv"
    frame.to_csv(csv_path, index=False)
    summary = {
        "deg_table": str(deg_path),
        "deg_total": int(total_before),
        "after_blacklist": int(len(frame)),
        "top_n": int(top_n),
        "keep_all": bool(keep_all),
        "output_csv": str(csv_path),
    }
    write_json(out_dir / "key_genes_summary.json", summary)
    log.info(
        "key targets: %s genes kept from %s DEGs -> %s",
        len(frame),
        total_before,
        csv_path,
    )
    return frame
