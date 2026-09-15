#!/usr/bin/env python3
"""Wet-lab validation handoff for prioritized targets.

Reads the ranked virtual knockout results and generates a staged validation
plan (cell lines -> organoids -> drug dose-response -> animal -> PDX) plus a
candidate CSV that can be handed directly to a wet-lab team.
"""

from __future__ import annotations

import pandas as pd

from .config import ResolvedConfig
from .provenance import write_run_manifest
from .utils import DockingError, write_json

PHASES = [
    (
        "phase_1_cell_line",
        "CRISPR 敲除 / 过表达（细胞系或原代细胞）",
        "增殖、凋亡、分化表型；至少两个独立 sgRNA",
    ),
    (
        "phase_2_organoid",
        "类器官 / 共培养体系",
        "微环境重塑、旁分泌信号、免疫共培养效应",
    ),
    (
        "phase_3_drug_response",
        "候选药物剂量反应",
        "IC50/EC50、靶点结合实验、关键通路验证",
    ),
    (
        "phase_4_animal",
        "动物模型",
        "药效学、安全性、药代动力学评价",
    ),
    (
        "phase_5_pdx",
        "PDX / 临床样本",
        "个体化疗效差异与耐药监测验证",
    ),
]

LIVER_CELL_LINES = [
    "HepG2",
    "Huh7",
    "Hep3B",
    "PLC/PRF/5",
    "MHCC97H",
    "SNU-449",
]


def export_validation(cfg: ResolvedConfig, log) -> dict:
    """Generate wet-lab validation plan files from ranked knockout results."""
    ko_dir = cfg.knockout_dir()
    ranked = ko_dir / "data" / "fig_52_53_ranked_knockout.csv"
    if not ranked.exists():
        raise DockingError(
            f"ranked knockout CSV not found: {ranked}; "
            "run virtual-knockout first"
        )
    frame = pd.read_csv(ranked)
    if frame.empty:
        raise DockingError(f"ranked knockout CSV is empty: {ranked}")

    top_n = min(int(cfg.get("validation", "top_n", 10)), len(frame))
    top = frame.head(top_n)
    out_dir = cfg.validation_dir()
    out_data = out_dir / "data"
    out_data.mkdir(parents=True, exist_ok=True)

    cols = [
        c
        for c in [
            "rank",
            "gene",
            "target_class",
            "target_score",
            "knockout_score",
            "reversal_score",
            "pathway_score",
            "specificity_score",
            "prognosis_score",
            "druggability_score",
            "off_target_paralogs",
            "safety_concern",
        ]
        if c in top.columns
    ]
    candidates = top[cols].copy()
    for phase, title, readout in PHASES:
        candidates[f"{phase}_assay"] = title
        candidates[f"{phase}_readout"] = readout
    candidates_csv = out_data / "validation_candidates.csv"
    candidates.to_csv(candidates_csv, index=False)

    plan_md = out_data / "validation_plan.md"
    plan_md.write_text(_render_plan(candidates), encoding="utf-8")

    summary = {
        "candidates": int(len(candidates)),
        "top_n": top_n,
        "classes": (
            candidates["target_class"].value_counts().to_dict()
            if "target_class" in candidates.columns
            else {}
        ),
        "validation_candidates_csv": str(candidates_csv),
        "validation_plan_md": str(plan_md),
    }
    write_json(out_data / "summary.json", summary)
    manifest = write_run_manifest(
        out_data,
        cfg,
        "export-validation",
        {"ranked_knockout_csv": ranked},
        {"top_n": top_n},
    )
    summary["manifest"] = str(manifest)
    log.info(
        "validation plan exported: %s candidates -> %s",
        len(candidates),
        out_dir,
    )
    return summary


def _render_plan(candidates: pd.DataFrame) -> str:
    lines = [
        "# 湿实验验证方案",
        "",
        "> 本方案由 AI 靶点评分结果自动生成，用于指导从计算靶点到实验验证"
        "的闭环；每一级实验通过后才进入下一级。",
        "",
        "## 验证阶段",
        "",
        "| 阶段 | 实验 | 核心指标 |",
        "| --- | --- | --- |",
    ]
    for phase, title, readout in PHASES:
        lines.append(f"| {phase} | {title} | {readout} |")
    lines += [
        "",
        "## 建议细胞模型",
        "",
        ", ".join(LIVER_CELL_LINES),
        "",
        "## 候选靶点",
        "",
        "| 排名 | 基因 | 类别 | 靶点得分 | 敲除得分 | 安全性标记 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for _, row in candidates.iterrows():
        lines.append(
            "| {rank} | {gene} | {cls} | {target:.3f} | {ko:.3f} | {safe:.0f} |".format(
                rank=row.get("rank", ""),
                gene=row.get("gene", ""),
                cls=row.get("target_class", "unknown"),
                target=float(row.get("target_score", 0) or 0),
                ko=float(row.get("knockout_score", 0) or 0),
                safe=float(row.get("safety_concern", 0) or 0),
            )
        )
    lines += [
        "",
        "## 质量要求",
        "",
        "- 细胞系阶段：至少两个独立 sgRNA，设置对照，重复三次。",
        "- 类器官阶段：同时记录肿瘤细胞、免疫细胞与基质细胞的组成变化。",
        "- 药物阶段：报告剂量反应曲线、靶点结合证据与通路标志物变化。",
        "- 动物阶段：报告药效、体重/脏器毒性及血浆药物浓度。",
        "- PDX 阶段：关联入组患者的分子分型与疗效，检验个体化差异。",
        "",
    ]
    return "\n".join(lines)
