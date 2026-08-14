#!/usr/bin/env python3
"""Run 10 synthetic single-cell datasets through the full pipeline."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from data.validation_generator import generate_dataset  # noqa: E402

DATASETS = [
    ("GSE90001", "mtx_single"),
    ("GSE90002", "barcode_csv_single"),
    ("GSE90003", "gene_csv_single"),
    ("GSE90004", "mtx_multi"),
    ("GSE90005", "barcode_csv_multi"),
    ("GSE90006", "gene_csv_multi"),
    ("GSE90007", "mtx_no_sample"),
    ("GSE90008", "mixed_gene_mtx"),
    ("GSE90009", "gene_csv_single"),
    ("GSE90010", "barcode_csv_single"),
]

REQUIRED_FIGURES = [
    "fig_01_qc_raw_violin.png",
    "fig_01_qc_filtered_violin.png",
    "fig_48_qc_pvalue_comparison.png",
    "fig_02_doublet_scores.png",
    "fig_03_umap_clusters.png",
    "fig_04_umap_condition.png",
    "fig_05_umap_annotation.png",
    "fig_06_dotplot_markers.png",
    "fig_07_annotation_confusion_heatmap.png",
    "fig_08_volcano.png",
    "fig_09_deg_heatmap.png",
    "fig_09_deg_horizontal_violin.png",
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
    "fig_46_go_top5.png",
    "fig_47_kegg_top5.png",
    "fig_24_ml_feature_importance.png",
    "fig_25_ml_shap.png",
]


def verify_outputs(out: Path, accession: str) -> bool:
    def find_figure(fig_dir: Path, name: str) -> bool:
        return any(p.is_file() for p in fig_dir.rglob(name))

    problems = []
    if not (out / "results" / "pipeline_complete.json").exists():
        problems.append("pipeline_complete.json")
    if not (out / "results" / "summary.json").exists():
        problems.append("summary.json")
    if not (out / "results" / "result_report.html").exists():
        problems.append("result_report.html")
    fig_dir = out / "results" / "figures"
    ml_summary = out / "results" / "data" / "07_ml" / "ml_model_summary.json"
    ml_skipped = False
    if ml_summary.exists():
        try:
            ml_skipped = (
                json.loads(ml_summary.read_text(encoding="utf-8")).get("status")
                == "skipped"
            )
        except Exception:
            ml_skipped = False
    for name in REQUIRED_FIGURES:
        if ml_skipped and name in {
            "fig_24_ml_feature_importance.png",
            "fig_25_ml_shap.png",
        }:
            continue
        if not find_figure(fig_dir, name):
            problems.append(name)
    if problems:
        print(f"[{accession}] missing outputs: {', '.join(problems)}")
        return False
    return True


def main() -> int:
    validation_root = ROOT / "validation_output"
    validation_root.mkdir(parents=True, exist_ok=True)
    success = False

    try:
        for accession, kind in DATASETS:
            out = validation_root / accession
            if (
                (out / "results" / "pipeline_complete.json").exists()
                and verify_outputs(out, accession)
            ):
                print(f"[{accession}] already passed, skipping")
                continue
            if out.exists():
                shutil.rmtree(out)
            print(f"Generating {accession} ({kind})")
            generate_dataset(
                accession,
                validation_root,
                kind,
                n_cells=1200,
                n_genes=3000,
                n_samples=4,
                seed=sum(ord(c) for c in accession),
            )

        for accession, kind in DATASETS:
            out = validation_root / accession
            if (
                (out / "results" / "pipeline_complete.json").exists()
                and verify_outputs(out, accession)
            ):
                continue
            print(f"Running pipeline for {accession} ({kind})")
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "run_pipeline.py"),
                    accession,
                    "--output",
                    str(out),
                    "--species",
                    "hs",
                    "--skip-download",
                    "--skip-deps",
                ],
                cwd=ROOT,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=1800,
            )
            if result.returncode != 0:
                print(f"[{accession}] pipeline failed")
                return 1
            if not verify_outputs(out, accession):
                return 1
            print(f"[{accession}] PASS")

        print("All 10 validation datasets passed.")
        success = True
        return 0
    finally:
        if success and validation_root.exists():
            shutil.rmtree(validation_root)
            print("Validation output cleaned.")


if __name__ == "__main__":
    sys.exit(main())
