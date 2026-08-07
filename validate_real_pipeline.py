#!/usr/bin/env python3
"""Run 10 real liver-disease GEO datasets through the full pipeline."""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

DATASETS = [
    ("GSE125449", "HCC vs iCCA", "hs"),
    ("GSE149614", "HCC tumor vs normal", "hs"),
    ("GSE140228", "HCC/CC tissues", "hs"),
    ("GSE295950", "MASH vs normal", "mm"),
    ("GSE264335", "HBV vs PHH", "hs"),
    ("GSE239964", "Liver fibrosis vs control", "mm"),
    ("GSE320340", "Alcohol-associated liver disease", "mm"),
    ("GSE291325", "FALD vs normal liver", "hs"),
    ("GSE330471", "HCC KLF2 vs vector", "mm"),
    ("GSE146409", "Liver tumor microenvironment", "hs"),
]

REQUIRED_FIGURES = [
    "fig_01_qc_raw_violin.png",
    "fig_01_qc_filtered_violin.png",
    "fig_02_doublet_scores.png",
    "fig_03_umap_clusters.png",
    "fig_04_umap_condition.png",
    "fig_05_umap_annotation.png",
    "fig_06_dotplot_markers.png",
    "fig_07_annotation_confusion_heatmap.png",
    "fig_08_volcano.png",
    "fig_09_deg_heatmap.png",
    "fig_10_go_up.png",
    "fig_11_go_down.png",
    "fig_12_kegg_up.png",
    "fig_13_kegg_down.png",
    "fig_14_pca.png",
    "fig_15_elbow.png",
    "fig_16_featureplot_markers.png",
    "fig_17_marker_violin.png",
    "fig_18_celltype_proportion.png",
    "fig_19_condition_proportion.png",
    "fig_20_gsea_go.png",
    "fig_21_gsea_kegg.png",
    "fig_22_go_network.png",
    "fig_23_kegg_network.png",
    "fig_24_ml_feature_importance.png",
    "fig_25_ml_shap.png",
]


def verify_outputs(out: Path, accession: str) -> bool:
    problems = []
    for rel in [
        "results/pipeline_complete.json",
        "results/summary.json",
        "results/result_report.html",
    ]:
        if not (out / rel).exists():
            problems.append(rel)
    fig_dir = out / "results" / "figures"
    for name in REQUIRED_FIGURES:
        if not (fig_dir / name).exists():
            problems.append(name)
    if problems:
        print(f"[{accession}] missing: {', '.join(problems[:20])}")
        return False
    return True


def run_one(accession: str, out: Path, species: str, timeout: int) -> bool:
    print(f"Running pipeline for {accession}")
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "run_pipeline.py"),
            accession,
            "--output",
            str(out),
            "--species",
            species,
        ],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )
    if result.returncode != 0:
        print(f"[{accession}] pipeline failed")
        return False
    if not verify_outputs(out, accession):
        return False
    print(f"[{accession}] PASS")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--only", default=None)
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--timeout", type=int, default=5400)
    args = parser.parse_args()

    validation_root = ROOT / "validation_real"
    validation_root.mkdir(parents=True, exist_ok=True)
    all_ok = True

    for accession, label, species in DATASETS:
        if args.only and accession != args.only.upper():
            continue
        out = validation_root / accession
        print(f"[{accession}] {label}")
        if (out / "results" / "pipeline_complete.json").exists() and verify_outputs(
            out, accession
        ):
            print(f"[{accession}] already passed, skipping")
            continue
        ok = run_one(accession, out, species, args.timeout)
        if not ok:
            all_ok = False
            if args.only:
                return 1

    if all_ok:
        print("All requested real datasets passed.")
        if not args.keep and validation_root.exists():
            shutil.rmtree(validation_root)
            print("Validation output cleaned.")
    else:
        print("Some datasets failed; validation output kept for inspection.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
