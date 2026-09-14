"""Classify generated results into the five figures of experiment plan one."""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from PIL import Image, UnidentifiedImageError

from .common import ensure_dir, write_json

CLASSIFIED_DIR = "按方案分类"


@dataclass(frozen=True)
class PanelAlias:
    figure: str
    panel: str
    description: str
    source: str
    destination: str
    status: str = "available"


PLAN_GROUPS: dict[str, tuple[str, ...]] = {
    "Figure1_化合物表征_靶点预测与通路富集": (
        "01_compound_characterization",
        "02_disease_targets",
        "03_intersection_ppi/enrichment",
    ),
    "Figure2_PPI网络与枢纽基因初步筛选": (
        "02b_evidence",
        "03_intersection_ppi",
        "04_bulk_training/fig2f_candidate_gene_heatmap.png",
        "04_bulk_training/fig2f_candidate_gene_heatmap.csv",
        "04_bulk_training/fig2g_gse49541_candidate_boxplots.png",
        "04_bulk_training/fig2g_gse49541_candidate_boxplot_stats.csv",
        "04_bulk_training/fig2h_gse164441_candidate_boxplots.png",
        "04_bulk_training/fig2h_gse164441_candidate_boxplot_stats.csv",
    ),
    "Figure3_机器学习模型构建与SHAP核心特征": (
        "04_bulk_training",
        "05_machine_learning",
    ),
    "Figure4_单细胞图谱_细胞通讯与虚拟扰动": (
        "06_single_cell_mouse",
        "07_single_cell_human",
    ),
    "Figure5_分子对接与分子动力学模拟": (
        "08_docking",
        "09_md_mmpbsa",
    ),
}

PANEL_ALIASES: tuple[PanelAlias, ...] = (
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "a",
        "研究全局流程图",
        "01_compound_characterization/fig1a_workflow.png",
        "Fig1a_研究全局流程图.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "b",
        "6PPD-Q 2D 化学结构",
        "01_compound_characterization/fig1b_compound_2d.png",
        "Fig1b_6PPD-Q_2D化学结构.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "c",
        "6PPD-Q 3D 结构",
        "01_compound_characterization/fig1c_compound_3d.png",
        "Fig1c_6PPD-Q_3D结构.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "c",
        "理化性质表",
        "01_compound_characterization/fig1c_physicochemical_properties.png",
        "Fig1c_理化性质表.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "d",
        "化合物靶点来源覆盖",
        "01_compound_characterization/fig1d_compound_target_source_counts.png",
        "Fig1d_化合物靶点来源覆盖.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "e",
        "疾病靶点来源状态",
        "02_disease_targets/fig1e_disease_target_source_status.png",
        "Fig1e_疾病靶点来源状态.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "f",
        "化合物-疾病靶点交集",
        "03_intersection_ppi/fig1f_compound_disease_venn.png",
        "Fig1f_化合物与疾病靶点交集.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "g",
        "KEGG 富集 Top 10",
        "03_intersection_ppi/enrichment/fig1g_kegg_enrichment_bar.png",
        "Fig1g_KEGG富集_Top10.png",
    ),
    PanelAlias(
        "Figure1_化合物表征_靶点预测与通路富集",
        "h",
        "GO 富集 BP/CC/MF",
        "03_intersection_ppi/enrichment/fig1h_go_enrichment_bubble.png",
        "Fig1h_GO富集_BP_CC_MF.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "a",
        "STRING PPI 网络",
        "03_intersection_ppi/fig2a_string_network.png",
        "Fig2a_STRING_PPI网络.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "b",
        "Cytoscape 风格模块网络",
        "03_intersection_ppi/fig2b_cytoscape_module_network.png",
        "Fig2b_PPI模块网络.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "c",
        "Degree 排名 Top 20",
        "03_intersection_ppi/fig2c_degree_top20.png",
        "Fig2c_Degree排名_Top20.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "d",
        "Betweenness 排名",
        "03_intersection_ppi/fig2d_betweenness_top20.png",
        "Fig2d_Betweenness排名.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "e",
        "MCC 与 Degree 交集",
        "03_intersection_ppi/fig2e_mcc_degree_venn.png",
        "Fig2e_MCC与Degree交集.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "f",
        "GSE89632 候选基因热图",
        "04_bulk_training/fig2f_candidate_gene_heatmap.png",
        "Fig2f_GSE89632_候选基因热图.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "g",
        "GSE49541 纤维化验证",
        "04_bulk_training/fig2g_gse49541_candidate_boxplots.png",
        "Fig2g_GSE49541_纤维化验证.png",
    ),
    PanelAlias(
        "Figure2_PPI网络与枢纽基因初步筛选",
        "h",
        "GSE164441 验证",
        "04_bulk_training/fig2h_gse164441_candidate_boxplots.png",
        "Fig2h_GSE164441_验证.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "a",
        "多模型 AUC 热图",
        "05_machine_learning/fig3a_multimodel_auc_heatmap.png",
        "Fig3a_多模型AUC热图.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "b",
        "最优模型训练 ROC",
        "05_machine_learning/fig3b_best_model_roc_training.png",
        "Fig3b_最优模型训练ROC.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "c",
        "GSE49541 外部验证 ROC",
        "05_machine_learning/fig3c_gse49541_external_roc.png",
        "Fig3c_GSE49541_外部验证ROC.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "d",
        "GSE164441 外部验证 ROC",
        "05_machine_learning/fig3d_gse164441_external_roc.png",
        "Fig3d_GSE164441_外部验证ROC.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "d",
        "GSE135251 补充 NAFLD 验证 ROC",
        "05_machine_learning/fig3d_supplementary_gse135251_roc.png",
        "Fig3d_补充_GSE135251_NAFLD验证ROC.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "e",
        "校准曲线",
        "05_machine_learning/fig3e_calibration_curve.png",
        "Fig3e_校准曲线.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "f",
        "SHAP 特征重要性",
        "05_machine_learning/fig3f_shap_importance.png",
        "Fig3f_SHAP特征重要性.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "g",
        "SHAP 蜂群图",
        "05_machine_learning/fig3g_shap_beeswarm.png",
        "Fig3g_SHAP蜂群图.png",
    ),
    PanelAlias(
        "Figure3_机器学习模型构建与SHAP核心特征",
        "h",
        "核心基因与纤维化分期关联",
        "05_machine_learning/fig3h_gse49541_core_genes_candidate_boxplots.png",
        "Fig3h_核心基因与纤维化分期.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "a",
        "小鼠 UMAP",
        "06_single_cell_mouse/fig4a_umap_cell_types.png",
        "Fig4a_小鼠UMAP.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "b",
        "细胞类型 marker 气泡图",
        "06_single_cell_mouse/fig4b_cell_type_markers.png",
        "Fig4b_细胞类型marker气泡图.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "c",
        "核心基因小提琴图",
        "06_single_cell_mouse/fig4c_core_gene_violin.png",
        "Fig4c_核心基因小提琴图.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "d",
        "核心基因 UMAP",
        "06_single_cell_mouse/fig4d_core_gene_umap.png",
        "Fig4d_核心基因UMAP.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "e",
        "细胞组成",
        "06_single_cell_mouse/fig4e_cell_composition.png",
        "Fig4e_细胞组成.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "f",
        "细胞通讯",
        "06_single_cell_mouse/fig4f_cellchat_network.png",
        "Fig4f_细胞通讯.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除 Top 10 变化数据",
        "06_single_cell_mouse/fig4g_virtual_knockout_top10.png",
        "Fig4g_虚拟敲除Top10.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除调控网络",
        "06_single_cell_mouse/fig4g_virtual_knockout_network.png",
        "Fig4g_虚拟敲除调控网络.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除 UMAP 偏移",
        "06_single_cell_mouse/fig4g_virtual_knockout_shift.png",
        "Fig4g_虚拟敲除UMAP偏移.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "h",
        "虚拟敲除富集结果",
        "06_single_cell_mouse/fig4h_virtual_knockout_enrichment.png",
        "Fig4h_虚拟敲除富集.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "i",
        "人类 snRNA UMAP",
        "07_single_cell_human/fig4i_umap_human_cell_types.png",
        "Fig4i_人类UMAP.png",
    ),
    PanelAlias(
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "j",
        "人类核心基因细胞类型验证",
        "07_single_cell_human/fig4j_human_core_gene_by_cell_type.png",
        "Fig4j_人类核心基因验证.png",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "a",
        "6PPD-Q 与核心靶点 3D 对接构象",
        "08_docking/fig5a_docking_pose_3d.png",
        "Fig5a_3D对接构象.png",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "b",
        "2D 相互作用平面图",
        "08_docking/fig5b_interaction_schematic.png",
        "Fig5b_2D相互作用平面图.png",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "c",
        "对接结合能热图",
        "08_docking/fig5c_docking_affinity_heatmap.png",
        "Fig5c_对接结合能热图.png",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "d",
        "蛋白主链 RMSD",
        "09_md_mmpbsa/fig5d_rmsd_NOT_RUN.png",
        "Fig5d_蛋白主链RMSD_未运行.png",
        "prepared_not_run",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "e",
        "配体 RMSD",
        "09_md_mmpbsa/fig5e_ligand_rmsd_NOT_RUN.png",
        "Fig5e_配体RMSD_未运行.png",
        "prepared_not_run",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "f",
        "关键残基 RMSF",
        "09_md_mmpbsa/fig5f_rmsf_NOT_RUN.png",
        "Fig5f_关键残基RMSF_未运行.png",
        "prepared_not_run",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "g",
        "回旋半径 Rg",
        "09_md_mmpbsa/fig5g_rg_NOT_RUN.png",
        "Fig5g_回旋半径_未运行.png",
        "prepared_not_run",
    ),
    PanelAlias(
        "Figure5_分子对接与分子动力学模拟",
        "h",
        "MM-PBSA 结合自由能分解",
        "09_md_mmpbsa/fig5h_mmpbsa_NOT_RUN.png",
        "Fig5h_MMPBSA_未运行.png",
        "prepared_not_run",
    ),
)


def classify_experiment_plan_results(
    output_root: Path,
    *,
    rebuild: bool = True,
) -> dict[str, Any]:
    """Create a readable Figure1-Figure5 mirror using hard links."""
    output_root = output_root.resolve()
    target = (output_root / CLASSIFIED_DIR).resolve()
    if target.parent != output_root:
        raise ValueError(f"unsafe classification target: {target}")
    if rebuild and target.exists():
        shutil.rmtree(target)
    ensure_dir(target)

    manifest: list[dict[str, Any]] = []
    for category, source_values in PLAN_GROUPS.items():
        category_dir = ensure_dir(target / category)
        ensure_dir(category_dir / "完整来源结果")
        for source_value in source_values:
            source = (output_root / source_value).resolve()
            if output_root not in source.parents and source != output_root:
                raise ValueError(f"source escapes output root: {source}")
            if not source.exists():
                manifest.append(
                    {
                        "分类": category,
                        "方案Panel": "整体来源",
                        "方案内容": "计划输出来源",
                        "分类文件": "",
                        "原文件": str(source),
                        "状态": "missing",
                        "说明": "源结果未生成",
                    }
                )
                continue
            destination_root = category_dir / "完整来源结果" / source.relative_to(
                output_root
            )
            if source.is_dir():
                for child in sorted(source.rglob("*")):
                    if not child.is_file():
                        continue
                    relative = child.relative_to(source)
                    destination = destination_root / relative
                    _link_or_copy(child, destination)
                    manifest.append(
                        _manifest_row(
                            category,
                            "整体来源",
                            "按方案归档的完整结果",
                            destination,
                            child,
                            target,
                        )
                    )
            else:
                _link_or_copy(source, destination_root)
                manifest.append(
                    _manifest_row(
                        category,
                        "整体来源",
                        "按方案归档的完整结果",
                        destination_root,
                        source,
                        target,
                    )
                )

    appendix = ensure_dir(target / "附录_数据_方法_质量状态")
    for source_value in ("00_data", "logs", "10_reports", "RESULTS_SUMMARY.md"):
        source = (output_root / source_value).resolve()
        if not source.exists():
            continue
        destination = appendix / source.relative_to(output_root)
        if source.is_dir():
            for child in sorted(source.rglob("*")):
                if child.is_file():
                    destination_file = destination / child.relative_to(source)
                    _link_or_copy(child, destination_file)
                    manifest.append(
                        _manifest_row(
                            appendix.name,
                            "附录",
                            "原始数据、运行日志与质量报告",
                            destination_file,
                            child,
                            target,
                        )
                    )
        else:
            _link_or_copy(source, destination)
            manifest.append(
                _manifest_row(
                    appendix.name,
                    "附录",
                    "原始数据、运行日志与质量报告",
                    destination,
                    source,
                    target,
                )
            )

    alias_records: list[dict[str, Any]] = []
    for alias in PANEL_ALIASES:
        source = (output_root / alias.source).resolve()
        destination = target / alias.figure / alias.destination
        destinations: list[Path] = []
        if source.exists():
            destinations.extend(_link_alias_variants(source, destination))
            actual_status = alias.status
            note = ""
        else:
            actual_status = "missing"
            note = "方案要求的文件名未找到，请检查分析日志"
        record = {
            "分类": alias.figure,
            "方案Panel": alias.panel,
            "方案内容": alias.description,
            "分类文件": ";".join(
                str(path.relative_to(target)) for path in destinations
            ),
            "原文件": str(source),
            "状态": actual_status,
            "说明": note,
        }
        alias_records.append(record)
        manifest.append(record)

    manifest_frame = pd.DataFrame(manifest)
    manifest_frame.to_csv(target / "分类清单.csv", index=False, encoding="utf-8-sig")
    _write_readme(target, alias_records)
    summary = {
        "classified_root": str(target),
        "categories": list(PLAN_GROUPS),
        "panel_entries": len(alias_records),
        "available_panels": sum(
            record["状态"] == "available" for record in alias_records
        ),
        "prepared_not_run_panels": sum(
            record["状态"] == "prepared_not_run" for record in alias_records
        ),
        "missing_panels": sum(
            record["状态"] == "missing" for record in alias_records
        ),
        "manifest": str(target / "分类清单.csv"),
        "readme": str(target / "分类说明.md"),
    }
    write_json(target / "分类汇总.json", summary)
    return summary


def _link_or_copy(source: Path, destination: Path) -> None:
    ensure_dir(destination.parent)
    if destination.exists():
        destination.unlink()
    try:
        os.link(source, destination)
    except OSError:
        shutil.copy2(source, destination)


def _link_alias_variants(source: Path, destination: Path) -> list[Path]:
    """Expose PNG, vector PDF when available, and an LZW TIFF derivative."""
    outputs: list[Path] = []
    _link_or_copy(source, destination)
    outputs.append(destination)
    for suffix in (".pdf", ".svg"):
        vector_source = source.with_suffix(suffix)
        if not vector_source.exists():
            continue
        vector_destination = destination.with_suffix(suffix)
        _link_or_copy(vector_source, vector_destination)
        outputs.append(vector_destination)
    if source.suffix.lower() == ".png":
        tiff_destination = destination.with_suffix(".tiff")
        ensure_dir(tiff_destination.parent)
        if tiff_destination.exists():
            tiff_destination.unlink()
        try:
            with Image.open(source) as image:
                source_dpi = image.info.get("dpi", (600.0, 600.0))
                dpi = (
                    float(source_dpi[0]) if len(source_dpi) >= 2 else 600.0,
                    float(source_dpi[1]) if len(source_dpi) >= 2 else 600.0,
                )
                image.convert("RGB").save(
                    tiff_destination,
                    format="TIFF",
                    compression="tiff_lzw",
                    dpi=dpi,
                )
            outputs.append(tiff_destination)
        except (UnidentifiedImageError, OSError):
            tiff_destination.unlink(missing_ok=True)
    return outputs


def _manifest_row(
    category: str,
    panel: str,
    description: str,
    destination: Path,
    source: Path,
    target: Path,
) -> dict[str, Any]:
    return {
        "分类": category,
        "方案Panel": panel,
        "方案内容": description,
        "分类文件": str(destination.relative_to(target)),
        "原文件": str(source),
        "状态": "available",
        "说明": "硬链接；原文件同步更新后内容通常保持一致",
    }


def _write_readme(target: Path, records: Iterable[dict[str, Any]]) -> None:
    lines = [
        "# 按实验方案分类的结果",
        "",
        "本目录把分析结果按方案中的 Figure 1-5 和 Panel 重新组织。",
        "每个 Figure 目录包含：",
        "",
        "- 以 `Fig...` 开头的方案 Panel 文件。",
        "- `完整来源结果/`：对应分析阶段的完整原始输出。",
        "- 原始结果仍在上一级目录保留，本目录使用硬链接，避免重复占用磁盘空间。",
        "",
        "## Panel 对照",
        "",
        "| Figure | Panel | 内容 | 状态 | 分类文件 |",
        "|---|---|---|---|---|",
    ]
    for record in records:
        lines.append(
            f"| {record['分类']} | {record['方案Panel']} | "
            f"{record['方案内容']} | {record['状态']} | "
            f"`{record['分类文件']}` |"
        )
    lines += [
        "",
        "## 重要说明",
        "",
        "- Figure 1d、1e 中 GeneCards、OMIM、TTD 无本地授权数据时，"
        "使用来源状态图而不是伪造 Venn 交集。",
        "- Figure 4f 是 CellChat 风格配体受体评分，不是完整 CellChat "
        "置换检验的替代声明。",
        "- Figure 5d-5h 标记为 `prepared_not_run`，表示 100 ns GROMACS "
        "生产轨迹未实际运行，输入文件和未运行说明均已归档。",
        "",
        "完整文件映射见 `分类清单.csv`，机器可读汇总见 `分类汇总.json`。",
        "",
    ]
    (target / "分类说明.md").write_text("\n".join(lines), encoding="utf-8")
