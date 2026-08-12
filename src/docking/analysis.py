#!/usr/bin/env python3
"""Rank docking results and generate hit reports and figures."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import ResolvedConfig
from .utils import DockingError, write_json


def analyze_results(cfg: ResolvedConfig, log):
    results_path = cfg.results_path()
    if not results_path.exists():
        raise DockingError(f"docking results not found: {results_path}")

    df = pd.read_csv(results_path, dtype={"id": str})
    if df.empty:
        raise DockingError("docking results CSV is empty")

    df["affinity"] = pd.to_numeric(df.get("affinity"), errors="coerce")
    ok = df[
        (df.get("status", "") == "ok") & df["affinity"].notna()
    ].copy()
    ok = ok.sort_values("affinity", na_position="last")
    ok["rank"] = np.arange(1, len(ok) + 1)

    cutoff = float(cfg.get("analysis", "cutoff", -7.0))
    top_n = int(cfg.get("analysis", "top_n", 100))
    hits = ok[ok["affinity"] <= cutoff].copy()
    top = hits.head(top_n) if not hits.empty else ok.head(top_n)

    reports_dir = cfg.analysis_dir()
    figures_dir = reports_dir / "figures"
    data_dir = reports_dir / "data"
    reports_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    ok.to_csv(data_dir / "fig_46_47_ranked_results.csv", index=False)
    hits.to_csv(data_dir / "fig_47_top_hits.csv", index=False)

    diverse = top
    if cfg.get("analysis", "diversity", True) and len(top) > 1:
        diverse = select_diverse(top, cfg, log)
        diverse.to_csv(data_dir / "fig_48_diverse_hits.csv", index=False)

    if cfg.get("analysis", "figures", True):
        try:
            make_figures(ok, hits, top, diverse, figures_dir, cutoff)
        except Exception as exc:
            log.warning("figure generation failed: %s", exc)

    try:
        import openpyxl  # noqa: F401

        with pd.ExcelWriter(
            data_dir / "docking_results.xlsx", engine="openpyxl"
        ) as xw:
            ok.to_excel(xw, sheet_name="all", index=False)
            hits.to_excel(xw, sheet_name="hits", index=False)
            top.to_excel(xw, sheet_name="top", index=False)
            diverse.to_excel(xw, sheet_name="diverse", index=False)
    except Exception as exc:
        log.warning("Excel report skipped: %s", exc)

    summary = {
        "total_docked": int(len(ok)),
        "hits": int(len(hits)),
        "top_n": int(len(top)),
        "diverse": int(len(diverse)),
        "cutoff": cutoff,
        "best_affinity": float(ok["affinity"].min()) if len(ok) else None,
        "median_affinity": float(ok["affinity"].median()) if len(ok) else None,
        "reports_dir": str(reports_dir),
    }
    write_json(reports_dir / "summary.json", summary)
    log.info(
        "analysis complete: %s hits (cutoff %s), top %s",
        len(hits),
        cutoff,
        len(top),
    )
    return summary


def select_diverse(frame: pd.DataFrame, cfg: ResolvedConfig, log) -> pd.DataFrame:
    if "smiles" not in frame.columns or frame.empty:
        return frame
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.DataStructs import BulkTanimotoSimilarity
    except ImportError:
        log.warning("RDKit not installed; diversity selection skipped")
        return frame

    mols = [Chem.MolFromSmiles(str(s)) for s in frame["smiles"]]
    valid = [i for i, mol in enumerate(mols) if mol is not None]
    if len(valid) < 2:
        return frame
    fps = [
        AllChem.GetMorganFingerprintAsBitVect(mols[i], 2, 2048)
        for i in valid
    ]
    k = min(int(cfg.get("analysis", "top_n", 100)), len(valid))
    selected = [valid[0]]
    remaining = [i for i in valid if i != valid[0]]
    while len(selected) < k and remaining:
        scores = []
        for i in remaining:
            sims = BulkTanimotoSimilarity(fps[i], [fps[j] for j in selected])
            scores.append((min(sims), i))
        _, best_idx = max(scores)
        selected.append(best_idx)
        remaining.remove(best_idx)
    return frame.iloc[selected].copy()


def make_figures(
    ok: pd.DataFrame,
    hits: pd.DataFrame,
    top: pd.DataFrame,
    diverse: pd.DataFrame,
    figures_dir: Path,
    cutoff: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.hist(ok["affinity"], bins=40, color="#4c7bb8", edgecolor="white")
    ax.axvline(cutoff, color="#c0392b", linestyle="--", linewidth=1.5)
    ax.set_xlabel("Affinity (kcal/mol)")
    ax.set_ylabel("Ligand count")
    ax.set_title("Docking affinity distribution")
    fig.tight_layout()
    fig.savefig(figures_dir / "fig_46_affinity_distribution.png", dpi=150)
    plt.close(fig)

    top20 = top.head(20).iloc[::-1]
    if not top20.empty:
        fig, ax = plt.subplots(
            figsize=(8, max(3.0, len(top20) * 0.32))
        )
        ax.barh(top20["id"].astype(str), top20["affinity"], color="#2e7d32")
        ax.set_xlabel("Affinity (kcal/mol)")
        ax.set_title("Top ranked docking hits")
        fig.tight_layout()
        fig.savefig(figures_dir / "fig_47_top_hits.png", dpi=150)
        plt.close(fig)

    if not diverse.empty:
        fig, ax = plt.subplots(figsize=(8, max(3.0, len(diverse) * 0.32)))
        ax.barh(
            diverse["id"].astype(str),
            diverse["affinity"],
            color="#8e44ad",
        )
        ax.set_xlabel("Affinity (kcal/mol)")
        ax.set_title("Diverse hit selection")
        fig.tight_layout()
        fig.savefig(figures_dir / "fig_48_diverse_hits.png", dpi=150)
        plt.close(fig)
