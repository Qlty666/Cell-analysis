#!/usr/bin/env python3
"""Mendelian randomisation and colocalisation for local GWAS/eQTL files.

The module accepts pre-exported summary statistics so it can run offline and
reproducibly. R packages are used when available; otherwise a transparent
Python IVW/weighted-median/MR-Egger fallback is written.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

APP_ROOT = Path(__file__).resolve().parents[2]
LOG = logging.getLogger("mr_coloc")


def _resolve(base: Path, value: str | Path | None) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _read_summary(spec: dict, base: Path, role: str) -> pd.DataFrame:
    path = _resolve(base, spec.get("file"))
    if path is None or not path.exists():
        raise FileNotFoundError(f"{role} summary statistics not found: {path}")
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(path, sep=sep)
    lower = {str(column).strip().lower(): str(column) for column in frame.columns}

    def pick(value, aliases: list[str], default: str) -> str:
        if value:
            return str(value)
        for alias in aliases:
            if alias in lower:
                return lower[alias]
        return default

    required = {
        "snp": pick(spec.get("snp"), ["snp", "rsid", "id"], "SNP"),
        "beta": pick(spec.get("beta"), ["beta", "b", "effect"], "beta"),
        "se": pick(spec.get("se"), ["se", "stderr", "standard_error"], "se"),
        "effect_allele": pick(
            spec.get("effect_allele"),
            ["effect_allele", "ea", "a1", "allele1"],
            "effect_allele",
        ),
        "other_allele": pick(
            spec.get("other_allele"),
            ["other_allele", "oa", "a2", "allele2"],
            "other_allele",
        ),
        "pval": pick(spec.get("pval"), ["pval", "p", "p_value", "pvalue"], "pval"),
    }
    missing = [column for column in required.values() if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{role} summary is missing columns {missing}; "
            f"available columns: {list(frame.columns)}"
        )
    output = pd.DataFrame(
        {
            "snp": frame[required["snp"]].astype(str),
            "beta": pd.to_numeric(frame[required["beta"]], errors="coerce"),
            "se": pd.to_numeric(frame[required["se"]], errors="coerce"),
            "effect_allele": frame[required["effect_allele"]]
            .astype(str)
            .str.upper(),
            "other_allele": frame[required["other_allele"]]
            .astype(str)
            .str.upper(),
            "pval": pd.to_numeric(frame[required["pval"]], errors="coerce"),
        }
    )
    if "eaf" in spec:
        output["eaf"] = pd.to_numeric(frame[spec["eaf"]], errors="coerce")
    for optional in ("chromosome", "position", "rsid"):
        key = spec.get(optional)
        if key and key in frame.columns:
            output[optional] = frame[key]
    output = output.dropna(subset=["snp", "beta", "se", "pval"])
    output = output[
        np.isfinite(output["beta"])
        & np.isfinite(output["se"])
        & (output["se"] > 0)
        & np.isfinite(output["pval"])
    ]
    output = output.drop_duplicates("snp", keep="first")
    return output.reset_index(drop=True)


def harmonise(exposure: pd.DataFrame, outcome: pd.DataFrame) -> pd.DataFrame:
    merged = exposure.merge(
        outcome,
        on="snp",
        suffixes=("_exposure", "_outcome"),
        how="inner",
    )
    if merged.empty:
        raise ValueError("no overlapping SNPs between exposure and outcome")
    exp_effect = merged["effect_allele_exposure"].astype(str).str.upper()
    exp_other = merged["other_allele_exposure"].astype(str).str.upper()
    out_effect = merged["effect_allele_outcome"].astype(str).str.upper()
    out_other = merged["other_allele_outcome"].astype(str).str.upper()
    same = (exp_effect == out_effect) & (exp_other == out_other)
    flipped = (exp_effect == out_other) & (exp_other == out_effect)
    merged["allele_status"] = np.where(
        same,
        "aligned",
        np.where(flipped, "flipped", "incompatible"),
    )
    merged.loc[flipped, "beta_outcome"] = -pd.to_numeric(
        merged.loc[flipped, "beta_outcome"],
        errors="coerce",
    )
    merged["palindromic"] = (
        merged["effect_allele_exposure"].isin(["A", "T"])
        & merged["other_allele_exposure"].isin(["A", "T"])
    ) | (
        merged["effect_allele_exposure"].isin(["C", "G"])
        & merged["other_allele_exposure"].isin(["C", "G"])
    )
    has_eaf = {"eaf_exposure", "eaf_outcome"}.issubset(merged.columns)
    merged["strand_status"] = "not_palindromic"
    merged["harmonisation_keep"] = merged["allele_status"].isin(
        ["aligned", "flipped"]
    )
    missing_eaf = int(0)
    ambiguous_palindromic = int(0)
    if has_eaf:
        eaf_exposure = pd.to_numeric(
            merged["eaf_exposure"],
            errors="coerce",
        )
        eaf_outcome = pd.to_numeric(
            merged["eaf_outcome"],
            errors="coerce",
        )
        eaf_outcome_aligned = eaf_outcome.where(~flipped, 1.0 - eaf_outcome)
        ambiguous = (
            eaf_exposure.between(0.42, 0.58)
            | eaf_outcome_aligned.between(0.42, 0.58)
            | (
                (eaf_exposure > 0.5)
                != (eaf_outcome_aligned > 0.5)
            )
            | eaf_exposure.isna()
            | eaf_outcome_aligned.isna()
        )
        merged.loc[merged["palindromic"] & ambiguous, "strand_status"] = (
            "ambiguous_palindromic"
        )
        merged.loc[
            merged["palindromic"] & ambiguous,
            "harmonisation_keep",
        ] = False
        ambiguous_palindromic = int(
            (merged["palindromic"] & ambiguous).sum()
        )
    else:
        missing_eaf = int(
            (merged["palindromic"] & merged["harmonisation_keep"]).sum()
        )
        merged.loc[
            merged["palindromic"] & merged["harmonisation_keep"],
            "strand_status",
        ] = "ambiguous_no_eaf"
        merged.loc[
            merged["palindromic"] & merged["harmonisation_keep"],
            "harmonisation_keep",
        ] = False
    merged["f_statistic"] = (
        pd.to_numeric(merged["beta_exposure"], errors="coerce")
        / pd.to_numeric(merged["se_exposure"], errors="coerce")
    ) ** 2
    merged["ratio"] = (
        pd.to_numeric(merged["beta_outcome"], errors="coerce")
        / pd.to_numeric(merged["beta_exposure"], errors="coerce")
    )
    incompatible = int(
        (~merged["allele_status"].isin(["aligned", "flipped"])).sum()
    )
    result = merged[
        merged["harmonisation_keep"]
        & np.isfinite(merged["ratio"])
    ].reset_index(drop=True)
    result.attrs["incompatible_snps"] = incompatible
    result.attrs["ambiguous_palindromic_snps"] = ambiguous_palindromic
    result.attrs["missing_eaf_palindromic_snps"] = missing_eaf
    return result


def _distance_clump(
    frame: pd.DataFrame,
    distance_kb: float,
) -> tuple[pd.DataFrame, str, int]:
    """Greedy distance pruning when no LD reference panel is configured."""
    if "chromosome" not in frame.columns or "position" not in frame.columns:
        return frame, "not_available_no_position", 0
    ordered = frame.sort_values(
        ["pval", "snp"],
        ascending=[True, True],
    ).copy()
    ordered["position"] = pd.to_numeric(
        ordered["position"],
        errors="coerce",
    )
    ordered = ordered[
        ordered["chromosome"].notna()
        & np.isfinite(ordered["position"])
    ]
    selected_indices: list[int] = []
    distance_bp = max(1.0, float(distance_kb) * 1000.0)
    for index, row in ordered.iterrows():
        chromosome = str(row.get("chromosome", ""))
        position = float(row["position"])
        if any(
            str(ordered.loc[selected, "chromosome"]) == chromosome
            and abs(float(ordered.loc[selected, "position"]) - position)
            < distance_bp
            for selected in selected_indices
        ):
            continue
        selected_indices.append(index)
    return (
        ordered.loc[selected_indices].reset_index(drop=True),
        "distance_pruning_no_ld_reference",
        int(len(ordered) - len(selected_indices)),
    )


def _plink_clump(
    frame: pd.DataFrame,
    clump_config: dict,
) -> tuple[pd.DataFrame, dict]:
    """Use PLINK for LD clumping when a local reference panel is configured."""
    executable = clump_config.get("plink_executable")
    bfile = clump_config.get("bfile")
    if not executable or not bfile:
        return frame, {
            "enabled": bool(clump_config.get("enabled", False)),
            "status": "not_configured",
            "method": "distance_pruning_no_ld_reference",
        }
    resolved_executable = shutil.which(str(executable)) or str(executable)
    if not Path(resolved_executable).exists() and shutil.which(resolved_executable) is None:
        return frame, {
            "enabled": True,
            "status": "failed",
            "method": "plink",
            "reason": f"PLINK executable not found: {executable}",
        }
    with tempfile.TemporaryDirectory(prefix="mr_clump_") as tmp:
        tmp_dir = Path(tmp)
        snp_file = tmp_dir / "instruments.txt"
        out_prefix = tmp_dir / "clump"
        snp_file.write_text(
            "\n".join(frame["snp"].astype(str)) + "\n",
            encoding="utf-8",
        )
        command = [
            resolved_executable,
            "--bfile",
            str(bfile),
            "--clump",
            str(snp_file),
            "--clump-p1",
            str(clump_config.get("p1", 1.0)),
            "--clump-p2",
            str(clump_config.get("p2", 1.0)),
            "--clump-r2",
            str(clump_config.get("r2", 0.001)),
            "--clump-kb",
            str(clump_config.get("distance_kb", 10000)),
            "--out",
            str(out_prefix),
        ]
        proc = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=int(clump_config.get("timeout_seconds", 1800)),
        )
        clumped_path = out_prefix.with_suffix(".clumped")
        if proc.returncode != 0 or not clumped_path.exists():
            return frame, {
                "enabled": True,
                "status": "failed",
                "method": "plink",
                "reason": (proc.stderr or proc.stdout)[-2000:],
            }
        clumped = pd.read_csv(clumped_path, sep=r"\s+")
        keep_snps = set(clumped["SNP"].astype(str))
        result = frame[frame["snp"].astype(str).isin(keep_snps)].copy()
        return result.reset_index(drop=True), {
            "enabled": True,
            "status": "completed",
            "method": "plink_ld_clump",
            "r2": float(clump_config.get("r2", 0.001)),
            "distance_kb": float(clump_config.get("distance_kb", 10000)),
            "removed": int(len(frame) - len(result)),
        }


def _ivw(frame: pd.DataFrame) -> dict:
    beta_x = frame["beta_exposure"].to_numpy(dtype=float)
    beta_y = frame["beta_outcome"].to_numpy(dtype=float)
    se_y = frame["se_outcome"].to_numpy(dtype=float)
    weights = 1.0 / np.square(se_y)
    estimate = float(np.sum(beta_x * beta_y * weights) / np.sum(np.square(beta_x) * weights))
    se = float(1.0 / math.sqrt(np.sum(np.square(beta_x) * weights)))
    p_value = float(2.0 * stats.norm.sf(abs(estimate / se)))
    return {
        "method": "inverse_variance_weighted",
        "estimate": estimate,
        "se": se,
        "pvalue": p_value,
        "ci_lower": estimate - 1.96 * se,
        "ci_upper": estimate + 1.96 * se,
        "nsnp": int(len(frame)),
    }


def _weighted_median(frame: pd.DataFrame) -> dict:
    ratios = frame["ratio"].to_numpy(dtype=float)
    beta_x = frame["beta_exposure"].to_numpy(dtype=float)
    se_y = frame["se_outcome"].to_numpy(dtype=float)
    weights = np.square(beta_x) / np.square(se_y)
    order = np.argsort(ratios)
    ratios = ratios[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cutoff = cumulative[-1] / 2.0
    estimate = float(ratios[np.searchsorted(cumulative, cutoff)])
    # Standard error via bootstrap to avoid a misleading normal approximation.
    rng = np.random.default_rng(42)
    boot = []
    for _ in range(500):
        idx = rng.integers(0, len(ratios), len(ratios))
        values = ratios[idx]
        w = weights[idx]
        order_boot = np.argsort(values)
        values = values[order_boot]
        w = w[order_boot]
        boot.append(float(values[np.searchsorted(np.cumsum(w), np.sum(w) / 2.0)]))
    se = float(np.std(boot, ddof=1)) if len(boot) > 2 else float("nan")
    p_value = (
        float(2.0 * stats.norm.sf(abs(estimate / se)))
        if math.isfinite(se) and se > 0
        else float("nan")
    )
    return {
        "method": "weighted_median",
        "estimate": estimate,
        "se": se,
        "pvalue": p_value,
        "ci_lower": estimate - 1.96 * se if math.isfinite(se) else float("nan"),
        "ci_upper": estimate + 1.96 * se if math.isfinite(se) else float("nan"),
        "nsnp": int(len(frame)),
    }


def _egger(frame: pd.DataFrame) -> dict:
    if len(frame) < 3:
        return {
            "method": "mr_egger",
            "estimate": float("nan"),
            "se": float("nan"),
            "pvalue": float("nan"),
            "intercept": float("nan"),
            "intercept_pvalue": float("nan"),
            "nsnp": int(len(frame)),
            "reason": "at least three instruments required",
        }
    beta_x = frame["beta_exposure"].to_numpy(dtype=float)
    beta_y = frame["beta_outcome"].to_numpy(dtype=float)
    se_y = frame["se_outcome"].to_numpy(dtype=float)
    weights = 1.0 / np.square(se_y)
    x = np.column_stack([np.ones(len(beta_x)), beta_x])
    weighted_x = x * np.sqrt(weights)[:, None]
    weighted_y = beta_y * np.sqrt(weights)
    coefficients, _, _, _ = np.linalg.lstsq(weighted_x, weighted_y, rcond=None)
    intercept, slope = coefficients
    residual = weighted_y - weighted_x @ coefficients
    dof = max(1, len(beta_x) - 2)
    sigma2 = float(np.sum(np.square(residual)) / dof)
    covariance = sigma2 * np.linalg.inv(weighted_x.T @ weighted_x)
    se_slope = float(math.sqrt(max(covariance[1, 1], 0.0)))
    se_intercept = float(math.sqrt(max(covariance[0, 0], 0.0)))
    return {
        "method": "mr_egger",
        "estimate": float(slope),
        "se": se_slope,
        "pvalue": float(2.0 * stats.t.sf(abs(slope / se_slope), dof)),
        "ci_lower": float(slope - stats.t.ppf(0.975, dof) * se_slope),
        "ci_upper": float(slope + stats.t.ppf(0.975, dof) * se_slope),
        "intercept": float(intercept),
        "intercept_pvalue": float(2.0 * stats.t.sf(abs(intercept / se_intercept), dof)),
        "nsnp": int(len(frame)),
    }


def _heterogeneity(frame: pd.DataFrame, ivw_estimate: float) -> dict:
    ratios = frame["ratio"].to_numpy(dtype=float)
    beta_x = frame["beta_exposure"].to_numpy(dtype=float)
    se_y = frame["se_outcome"].to_numpy(dtype=float)
    weights = np.square(beta_x) / np.square(se_y)
    q_value = float(np.sum(weights * np.square(ratios - ivw_estimate)))
    dof = max(1, len(frame) - 1)
    p_value = float(stats.chi2.sf(q_value, dof))
    return {
        "cochran_q": q_value,
        "q_df": int(dof),
        "q_pvalue": p_value,
        "i_squared": float(max(0.0, (q_value - dof) / q_value) * 100.0)
        if q_value > 0
        else 0.0,
    }


def _leave_one_out(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for snp in frame["snp"].astype(str):
        subset = frame[frame["snp"].astype(str) != snp]
        if len(subset) < 2:
            continue
        result = _ivw(subset)
        rows.append(
            {
                "snp_removed": snp,
                "estimate": result["estimate"],
                "se": result["se"],
                "pvalue": result["pvalue"],
            }
        )
    return pd.DataFrame(rows)


def run_r_backend(
    harmonised_path: Path,
    config_path: Path,
    out_dir: Path,
) -> dict:
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    script = APP_ROOT / "src" / "analysis" / "R" / "mr_coloc.R"
    if rscript is None or not script.exists():
        return {"status": "skipped", "reason": "Rscript or mr_coloc.R unavailable"}
    proc = subprocess.run(
        [rscript, str(script), str(harmonised_path), str(config_path), str(out_dir)],
        cwd=APP_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=int(_read_json(config_path).get("timeout_seconds", 3600)),
    )
    (out_dir / "mr_r.log").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or ""),
        encoding="utf-8",
    )
    if proc.returncode != 0:
        return {
            "status": "failed",
            "returncode": proc.returncode,
            "reason": proc.stderr[-2000:],
        }
    path = out_dir / "mr_r_summary.json"
    return _read_json(path) if path.exists() else {"status": "completed"}


def run(config: dict, out_dir: Path, skip_r: bool = False) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    base = Path(str(config.get("_config_dir") or Path.cwd())).resolve()
    exposure_spec = config.get("exposure", {})
    outcome_spec = config.get("outcome", {})
    exposure = _read_summary(exposure_spec, base, "exposure")
    outcome = _read_summary(outcome_spec, base, "outcome")
    p_threshold = float(config.get("p_threshold", 5e-6))
    exposure = exposure[exposure["pval"] <= p_threshold].copy()
    if exposure.empty:
        raise ValueError(
            f"no exposure instruments passed p <= {p_threshold}; "
            "relax p_threshold or provide pre-selected SNPs"
        )
    clump_config = config.get("clump", {}) or {}
    clump_summary: dict
    if clump_config.get("enabled", True):
        if clump_config.get("bfile"):
            exposure, clump_summary = _plink_clump(exposure, clump_config)
        else:
            exposure, clump_method, removed = _distance_clump(
                exposure,
                float(clump_config.get("distance_kb", 10000)),
            )
            clump_summary = {
                "enabled": True,
                "status": (
                    "completed"
                    if clump_method != "not_available_no_position"
                    else "not_available"
                ),
                "method": clump_method,
                "removed": int(removed),
                "distance_kb": float(
                    clump_config.get("distance_kb", 10000)
                ),
            }
    else:
        clump_summary = {"enabled": False, "status": "disabled", "removed": 0}
    harmonised = harmonise(exposure, outcome)
    if len(harmonised) < 2:
        raise ValueError("fewer than two valid instruments after harmonisation")
    harmonised.to_csv(out_dir / "harmonised_instruments.csv", index=False)

    ivw = _ivw(harmonised)
    median = _weighted_median(harmonised)
    egger = _egger(harmonised)
    methods = pd.DataFrame([ivw, median, egger])
    methods.to_csv(out_dir / "mr_methods.csv", index=False)
    heterogeneity = _heterogeneity(harmonised, ivw["estimate"])
    heterogeneity["mean_f_statistic"] = float(
        harmonised["f_statistic"].mean()
    )
    heterogeneity["weak_instrument_count"] = int(
        (harmonised["f_statistic"] < 10).sum()
    )
    _write_json(out_dir / "mr_sensitivity.json", heterogeneity)
    _leave_one_out(harmonised).to_csv(
        out_dir / "mr_leave_one_out.csv",
        index=False,
    )

    config_path = out_dir / "resolved_mr_config.json"
    _write_json(config_path, config)
    r_summary = (
        {"status": "skipped", "reason": "disabled by --skip-r"}
        if skip_r
        else run_r_backend(
            out_dir / "harmonised_instruments.csv",
            config_path,
            out_dir,
        )
    )
    summary = {
        "status": "completed",
        "n_exposure_snps": int(len(exposure)),
        "n_harmonised_snps": int(len(harmonised)),
        "clumping": clump_summary,
        "incompatible_snps": int(
            harmonised.attrs.get("incompatible_snps", 0)
        ),
        "ambiguous_palindromic_snps": int(
            harmonised.attrs.get("ambiguous_palindromic_snps", 0)
        ),
        "missing_eaf_palindromic_snps": int(
            harmonised.attrs.get("missing_eaf_palindromic_snps", 0)
        ),
        "p_threshold": p_threshold,
        "primary": ivw,
        "methods": methods.to_dict(orient="records"),
        "sensitivity": heterogeneity,
        "r_backend": r_summary,
        "outputs": {
            "harmonised": str(out_dir / "harmonised_instruments.csv"),
            "methods": str(out_dir / "mr_methods.csv"),
            "sensitivity": str(out_dir / "mr_sensitivity.json"),
            "leave_one_out": str(out_dir / "mr_leave_one_out.csv"),
        },
    }
    _write_json(out_dir / "mr_summary.json", summary)
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Local Mendelian randomisation and colocalisation analysis."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--skip-r", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    config_path = Path(args.config).expanduser().resolve()
    out_dir = Path(args.output).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        config = _read_json(config_path)
        config["_config_dir"] = str(config_path.parent)
        summary = run(config, out_dir, skip_r=args.skip_r)
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return 0
    except Exception as exc:
        LOG.exception("MR/colocalisation failed: %s", exc)
        _write_json(
            out_dir / "mr_summary.json",
            {"status": "failed", "reason": str(exc)},
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
