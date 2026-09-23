"""Automated quality screening and readable contact sheets for plan figures."""

from __future__ import annotations

import html
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

from .classify import CLASSIFIED_DIR, PANEL_ALIASES, PanelAlias
from .common import ensure_dir, read_json, write_json


@dataclass(frozen=True)
class FigureReview:
    scientific: str
    aesthetic: str
    verdict: str
    notes: str


REVIEWS: dict[tuple[str, str], FigureReview] = {
    ("Figure1_化合物表征_靶点预测与通路富集", "a"): FigureReview(
        "合理",
        "中等",
        "需重排",
        "流程完整且已压缩到7.2英寸画布，但框内文字仍偏小；正式主图建议精简为发现、建模、模拟三行，完整流程放扩展数据。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "b"): FigureReview(
        "合理",
        "良好",
        "可用",
        "二维结构已按600 dpi输出并提供SVG矢量版本。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "c"): FigureReview(
        "合理",
        "中等",
        "需修饰",
        "3D构象和理化性质均存在，但两者是独立文件，正式Figure需要合并为同一panel。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "d"): FigureReview(
        "部分合理",
        "中等",
        "不可替代Venn",
        "仅SwissTargetPrediction产生有效靶点，ChEMBL/Swiss与STITCH未形成可用交集，当前用来源覆盖图，未伪造Venn。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "e"): FigureReview(
        "部分合理",
        "中等",
        "不可替代Venn",
        "GeneCards、OMIM、TTD缺少授权数据，当前为来源状态图，不能当作三数据库交集结果。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "f"): FigureReview(
        "合理",
        "良好",
        "可用",
        "交集方向与数据一致，建议补充交集基因列表。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "g"): FigureReview(
        "合理",
        "良好",
        "需修饰",
        "Top10 KEGG已更新为600 dpi；长通路名在主图中仍偏小，建议缩写或缩短为机制相关条目。",
    ),
    ("Figure1_化合物表征_靶点预测与通路富集", "h"): FigureReview(
        "合理",
        "良好",
        "需修饰",
        "GO条目已更新为600 dpi；术语较多，主图建议保留Top5或拆分BP/CC/MF并缩短标签。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "a"): FigureReview(
        "合理",
        "中等",
        "需修饰",
        "已用STRING add_nodes=50扩展到60个节点、389条边，并区分种子基因与上下文节点；主图仍需减少标签，完整网络更适合扩展数据。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "b"): FigureReview(
        "合理",
        "良好",
        "可用",
        "按原方案运行R MCL；图注应记录inflation、上下文节点和聚类版本。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "c"): FigureReview(
        "合理",
        "良好",
        "可用",
        "Degree排序和Top5强调正确。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "d"): FigureReview(
        "合理",
        "良好",
        "可用",
        "Centrality排名正确，但当前网络小，betweenness区分度有限。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "e"): FigureReview(
        "合理",
        "良好",
        "可用",
        "MCC和Degree交集较小，图中数字与网络规模一致。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "f"): FigureReview(
        "合理",
        "良好",
        "可用",
        "已压缩到7.2英寸并增加分组色带；63个样本不再逐条显示编号，避免出版尺寸下标签重叠。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "g"): FigureReview(
        "合理",
        "良好",
        "可用",
        "纤维化分期箱线图合理，建议在图注报告FDR和效应量。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "h"): FigureReview(
        "合理",
        "良好",
        "可用",
        "GSE164441比较的是HCC肿瘤与癌旁组织，不能作为健康/NAFLD验证的替代终点。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "a"): FigureReview(
        "合理",
        "良好",
        "可用",
        "模型和特征组合覆盖充分，数值标注清晰。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "b"): FigureReview(
        "合理",
        "良好",
        "可用",
        "ROC由整个模型/特征选择流程的外折OOF预测绘制，而不是普通CV预测。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "c"): FigureReview(
        "部分合理",
        "良好",
        "需限定解释",
        "GSE49541检测的是纤维化分期，不是同一健康/NAFLD终点；仅作为跨终点探索，不设AUC通过阈值。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "d"): FigureReview(
        "部分合理",
        "良好",
        "需限定解释",
        "GSE135251是同一疾病类别候选但既往分析暴露需登记；GSE164441为HCC配对癌旁扩展终点。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "e"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "校准参数、Brier和H-L检验需联合解释；H-L p>0.05 不单独证明校准合格。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "f"): FigureReview(
        "合理",
        "良好",
        "可用",
        "SHAP值来自明确标识的未校准解释模型；需报告重复筛选稳定性。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "g"): FigureReview(
        "合理",
        "良好",
        "可用",
        "SHAP蜂群图方向与表达值一致。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "h"): FigureReview(
        "合理",
        "良好",
        "可用",
        "纤维化关联图可读，建议统一展示FDR。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "a"): FigureReview(
        "合理",
        "良好",
        "可用",
        "UMAP分区和细胞类型图例清楚。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "b"): FigureReview(
        "合理",
        "良好",
        "可用",
        "7.2英寸宽度下marker基因标签已放大，仍需在图注中说明marker与细胞类型对应关系。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "c"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "小提琴图仅描述细胞分布；GSE270583每个条件生物学重复不足，不输出组间显著性。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "d"): FigureReview(
        "合理",
        "良好",
        "可用",
        "核心基因UMAP可读，适合展示表达分布。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "e"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "组成比例先按独立文库/样本计算；单NCD对照不足以进行群体组成推断。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "按原方案运行R CellChat并保存版本、数据库和参数；细胞级概率不替代独立生物单位统计。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "g"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "按原方案运行R scTenifoldKnk并记录网络数、细胞数和随机种子；预测扰动仍不是湿实验敲除。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "h"): FigureReview(
        "合理",
        "中等",
        "阴性结果",
        "未发现显著GO/KEGG条目时保留为明确的阴性结果，不伪造富集通路。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "i"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "人类UMAP按当前SOFT解析展示；多文库患者需合并，未知donor不得参与患者级推断。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "j"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "热图为donor均值×细胞类型×疾病阶段探索性展示，并明确未知患者排除和重复文库合并。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "a"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "已重绘为结合口袋原子、蛋白主链轨迹、配体骨架和口袋残基标注；Discovery Studio/PyMOL 可作为正交复核，不是唯一合法实现。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "b"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "输出距离/几何候选接触并标注证据级别；未验证角度和环平面时不得写成已证实氢键、盐桥或π堆积。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "c"): FigureReview(
        "合理",
        "良好",
        "需限定解释",
        "热图保留预指定靶点顺序；Vina分数是口袋特异预测值，不跨蛋白排序，也不能替代实验亲和力。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "d"): FigureReview(
        "未运行",
        "不可评价",
        "不满足方案",
        "100 ns蛋白主链RMSD未运行。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "e"): FigureReview(
        "未运行",
        "不可评价",
        "不满足方案",
        "100 ns配体RMSD未运行。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "f"): FigureReview(
        "未运行",
        "不可评价",
        "不满足方案",
        "关键残基RMSF未运行。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "g"): FigureReview(
        "未运行",
        "不可评价",
        "不满足方案",
        "回旋半径Rg未运行。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "h"): FigureReview(
        "未运行",
        "不可评价",
        "不满足方案",
        "需要真实残基分解、TOTAL列、单位与误差；仅总能量不能满足本面板。",
    ),
}

LABEL_OVERLAP_REVIEWS: dict[tuple[str, str], tuple[str, str]] = {
    ("Figure1_化合物表征_靶点预测与通路富集", "f"): (
        "无",
        "两集合区域数字已分别放置，原交集中心的重叠问题已修复。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "a"): (
        "无",
        "仅保留交集基因标签，STRING上下文节点不再逐个标注，已消除中心标签拥挤。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "b"): (
        "无",
        "模块网络仅标注交集基因，模块归属由图例颜色表达，已消除中心标签拥挤。",
    ),
    ("Figure2_PPI网络与枢纽基因初步筛选", "e"): (
        "无",
        "两集合区域数字已分别放置，原交集中心的重叠问题已修复。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f"): (
        "轻微",
        "细胞和通路名称位于节点上，白色描边已提高可读性；最终组图可进一步放大。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "i"): (
        "无",
        "图例已增加边框并与散点区域分离。",
    ),
}

LABEL_OVERLAP_ALIAS_REVIEWS: dict[tuple[str, str, str], tuple[str, str]] = {
    (
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除 Top 10 变化数据",
    ): (
        "无",
        "表格文字在原分辨率下未发现重叠。",
    ),
    (
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除调控网络",
    ): (
        "轻微",
        "外部渲染的网络节点文字较满，建议仅用于补充材料。",
    ),
    (
        "Figure4_单细胞图谱_细胞通讯与虚拟扰动",
        "g",
        "虚拟敲除 UMAP 偏移",
    ): (
        "无",
        "原分辨率下未发现文字或标签重叠。",
    ),
}


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _format_number(value: Any, digits: int = 3) -> str:
    number = _number(value)
    return f"{number:.{digits}f}" if number is not None else "NA"


def _dynamic_result_reviews(
    output_root: Path,
) -> dict[tuple[str, str], FigureReview]:
    """Use fresh result metrics for panels whose interpretation is data-dependent."""
    summary = read_json(
        output_root / "05_machine_learning" / "ml_summary.json",
        {},
    )
    external = {
        str(row.get("dataset")): row
        for row in (summary.get("external_validation") or [])
        if isinstance(row, dict)
    }
    calibration = (
        summary.get("calibration")
        if isinstance(summary.get("calibration"), dict)
        else {}
    )
    figure = "Figure3_机器学习模型构建与SHAP核心特征"
    reviews: dict[tuple[str, str], FigureReview] = {}

    ppi_summary = read_json(
        output_root / "03_intersection_ppi" / "ppi_summary.json",
        {},
    )
    if isinstance(ppi_summary, dict) and ppi_summary:
        module_method = str(ppi_summary.get("module_method") or "")
        ppi_status = str(ppi_summary.get("status") or "unknown")
        reviews[("Figure2_PPI网络与枢纽基因初步筛选", "b")] = FigureReview(
            "合理",
            "良好",
            (
                "可用"
                if ppi_status == "completed" and module_method.upper() == "MCL"
                else "阴性结果"
                if ppi_status == "valid_negative"
                else "需限定解释"
            ),
            (
                f"status={ppi_status}; "
                f"module_method={module_method or 'unknown'}; "
                f"inflation={ppi_summary.get('mcl_inflation', 'NA')}; "
                "MCL为原方案方法，Louvain只能作为登记后的替代方法。"
            ),
        )

    mouse_dir = output_root / "06_single_cell_mouse"
    cellchat_summary = read_json(
        mouse_dir / "cellchat_permutation_summary.json",
        {},
    )
    if isinstance(cellchat_summary, dict) and cellchat_summary:
        communication_status = str(
            cellchat_summary.get("status") or "unknown"
        )
        method = str(cellchat_summary.get("method") or "unknown")
        reviews[("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f")] = FigureReview(
            "合理",
            "良好",
            (
                "阴性结果"
                if communication_status == "valid_negative"
                else "需限定解释"
            ),
            (
                f"method={method}; status={communication_status}; "
                f"interactions={cellchat_summary.get('n_interactions', 'NA')}; "
                "细胞级通讯概率不替代独立生物单位统计。"
            ),
        )

    insilico_candidates = sorted(mouse_dir.rglob("insilico_summary.json"))
    insilico_summary = (
        read_json(insilico_candidates[0], {})
        if insilico_candidates
        else {}
    )
    if isinstance(insilico_summary, dict) and insilico_summary:
        insilico_status = str(
            insilico_summary.get("status") or "unknown"
        )
        engine = str(insilico_summary.get("engine") or "unknown")
        scientific_role = str(
            insilico_summary.get("scientific_role")
            or "predicted_perturbation_response"
        )
        insilico_review = FigureReview(
            "合理",
            "良好",
            "可用" if insilico_status == "completed" else "需限定解释",
            (
                f"engine={engine}; status={insilico_status}; "
                f"role={scientific_role}; "
                "结果仍是预测扰动响应，不是湿实验敲除。"
            ),
        )
        reviews[
            ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "g")
        ] = insilico_review
        reviews[
            ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "h")
        ] = insilico_review

    interaction_status = read_json(
        output_root / "08_docking" / "fig5b_status.json",
        {},
    )
    if isinstance(interaction_status, dict) and interaction_status:
        validated = bool(interaction_status.get("validated"))
        interaction_state = str(
            interaction_status.get("status") or "unknown"
        )
        reviews[
            ("Figure5_分子对接与分子动力学模拟", "b")
        ] = FigureReview(
            "合理",
            "良好",
            (
                "可用"
                if validated
                else "阴性结果"
                if interaction_state == "valid_negative"
                else "需限定解释"
            ),
            (
                f"status={interaction_state}; "
                f"method={interaction_status.get('method') or 'geometric candidates'}; "
                f"validated={validated}; "
                "PLIP验证只针对所选docking pose，不等于实验结合证据。"
            ),
        )

    fibrosis = external.get("GSE49541_fibrosis")
    if fibrosis:
        reviews[(figure, "c")] = FigureReview(
            "合理",
            "良好",
            "需限定解释",
            (
                "GSE49541检测的是纤维化分期，不是同一健康/NAFLD终点；"
                f"AUC {_format_number(fibrosis.get('auc'))}，"
                "该终点不套用NAFLD分类AUC阈值，也不以跨终点性能判定成功。"
            ),
        )

    external_validation = [
        row
        for key in ("GSE164441_tumor", "GSE135251_NAFLD")
        if (row := external.get(key))
    ]
    if external_validation:
        details = []
        for row in external_validation:
            dataset = str(row.get("dataset") or "")
            suffix = (
                " (same endpoint category; prior analysis exposure recorded)"
                if dataset == "GSE135251_NAFLD"
                else " (HCC paired-adjacent endpoint; no NAFLD target)"
            )
            details.append(
                f"{dataset} AUC {_format_number(row.get('auc'))}{suffix}"
            )
        reviews[(figure, "d")] = FigureReview(
            "合理",
            "良好",
            "需限定解释",
            (
                "; ".join(details)
                + "。GSE135251为同一疾病类别候选但已登记既往分析暴露；"
                "GSE164441为HCC肿瘤与配对癌旁终点，不能称为同终点验证。"
            ),
        )

    if calibration:
        reviews[(figure, "e")] = FigureReview(
            "合理",
            "良好",
            "需限定解释",
            (
                "嵌套外折预测校准："
                f"Hosmer-Lemeshow p={_format_number(calibration.get('hosmer_lemeshow_p'), 3)}, "
                f"Brier={_format_number(calibration.get('brier'))}, "
                f"slope={_format_number(calibration.get('calibration_slope'))}, "
                f"intercept={_format_number(calibration.get('calibration_intercept'))}; "
                "这些指标需联合解释，H-L p>0.05 单独不证明校准合格。"
            ),
        )
    return reviews


def audit_figures(
    output_root: Path,
    *,
    rebuild: bool = True,
) -> dict[str, Any]:
    """Audit classified figures and create contact sheets for visual review."""
    output_root = output_root.resolve()
    classified = output_root / CLASSIFIED_DIR
    if not classified.exists():
        raise FileNotFoundError(
            f"classified result folder not found: {classified}; run classify first"
        )
    report_dir = ensure_dir(output_root / "10_reports" / "figure_quality_audit")
    contact_dir = ensure_dir(report_dir / "contact_sheets")
    aliases_by_figure: dict[str, list[PanelAlias]] = {}
    for alias in PANEL_ALIASES:
        aliases_by_figure.setdefault(alias.figure, []).append(alias)
    md_figure_status = read_json(
        output_root / "09_md_mmpbsa" / "md_plan_figures.json",
        {},
    )
    md_panel_status = md_figure_status.get("panels") or {}
    compound_venn = (
        output_root
        / "01_compound_characterization"
        / "fig1d_compound_target_source_venn.png"
    ).exists()
    disease_venn = (
        output_root
        / "02_disease_targets"
        / "fig1e_disease_target_venn.png"
    ).exists()
    dynamic_reviews = _dynamic_result_reviews(output_root)

    rows: list[dict[str, Any]] = []
    for figure, aliases in aliases_by_figure.items():
        for alias in aliases:
            source = (output_root / alias.source).resolve()
            review = dynamic_reviews.get(
                (figure, alias.panel),
                REVIEWS.get(
                    (figure, alias.panel),
                    FigureReview("不确定", "不确定", "人工复核", ""),
                ),
            )
            if (
                figure == "Figure5_分子对接与分子动力学模拟"
                and alias.panel in md_panel_status
                and bool(md_panel_status.get(alias.panel))
            ):
                review = FigureReview(
                    "合理",
                    "良好",
                    "可用",
                    (
                        "Generated from provenance-bound GROMACS/MM-PBSA "
                        "outputs; "
                        f"run={md_figure_status.get('run_manifest', {}).get('run_id', '')}, "
                        f"production={md_figure_status.get('production_ns', '')} ns."
                    ),
                )
            if (
                figure == "Figure1_化合物表征_靶点预测与通路富集"
                and alias.panel == "c"
                and "合并" in alias.description
            ):
                review = FigureReview(
                    "合理",
                    "良好",
                    "可用",
                    "3D 构象与理化性质已合并为单一 panel。",
                )
            if (
                figure == "Figure1_化合物表征_靶点预测与通路富集"
                and alias.panel == "d"
                and compound_venn
            ):
                review = FigureReview(
                    "合理",
                    "良好",
                    "可用",
                    "已基于至少两个可审计化合物靶点来源生成真实交集。",
                )
            if (
                figure == "Figure1_化合物表征_靶点预测与通路富集"
                and alias.panel == "e"
                and disease_venn
            ):
                review = FigureReview(
                    "合理",
                    "良好",
                    "可用",
                    "已基于至少两个可审计疾病靶点来源生成真实交集。",
                )
            overlap_severity, overlap_notes = LABEL_OVERLAP_ALIAS_REVIEWS.get(
                (figure, alias.panel, alias.description),
                LABEL_OVERLAP_REVIEWS.get(
                    (figure, alias.panel),
                    ("无", "原分辨率视觉复核未发现文字重叠。"),
                ),
            )
            metrics = _image_metrics(source) if source.exists() else {}
            row = {
                "figure": figure,
                "panel": alias.panel,
                "content": alias.description,
                "status": alias.status if source.exists() else "missing",
                "file": str(source),
                "scientific_reasonableness": review.scientific,
                "aesthetic_level": review.aesthetic,
                "audit_verdict": review.verdict,
                "recommended_output": _recommended_output(review.verdict),
                "notes": review.notes,
                "label_overlap_severity": overlap_severity,
                "label_overlap_notes": overlap_notes,
                **metrics,
            }
            row["automatic_quality_score"] = _automatic_score(row)
            rows.append(row)
        _contact_sheet(
            aliases,
            output_root,
            contact_dir / f"{figure}_contact_sheet.png",
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(report_dir / "figure_quality_audit.csv", index=False, encoding="utf-8-sig")
    summary = _summary(audit)
    write_json(report_dir / "figure_quality_audit.json", summary)
    (report_dir / "figure_quality_audit.md").write_text(
        _render_audit_markdown(audit, summary, classified),
        encoding="utf-8",
    )
    return summary


def _image_metrics(path: Path) -> dict[str, Any]:
    with Image.open(path) as image:
        image = image.convert("RGB")
        width, height = image.size
        dpi = image.info.get("dpi") or (72.0, 72.0)
        rgb = np.asarray(image.resize((max(1, width // 8), max(1, height // 8))))
        near_white = np.all(rgb >= 245, axis=2)
        non_white_fraction = float(1.0 - near_white.mean())
        svg_metrics = _svg_text_metrics(path.with_suffix(".svg"))
        return {
            "width_px": int(width),
            "height_px": int(height),
            "dpi_x": round(float(dpi[0]), 1),
            "dpi_y": round(float(dpi[1]), 1),
            "physical_width_in": round(width / max(float(dpi[0]), 1.0), 3),
            "physical_height_in": round(height / max(float(dpi[1]), 1.0), 3),
            "mode": image.mode,
            "format": image.format or path.suffix.lstrip(".").upper(),
            "file_size_mb": round(path.stat().st_size / (1024 * 1024), 3),
            "aspect_ratio": round(width / height, 3) if height else 0.0,
            "non_white_fraction": round(non_white_fraction, 4),
            "has_pdf": path.with_suffix(".pdf").exists(),
            "has_svg": path.with_suffix(".svg").exists(),
            "has_tiff": path.with_suffix(".tiff").exists(),
            **svg_metrics,
        }


def _svg_text_metrics(svg_path: Path) -> dict[str, Any]:
    if not svg_path.exists():
        return {
            "duplicate_svg_text_anchor_groups": None,
            "duplicate_svg_text_examples": "",
        }
    text = svg_path.read_text(encoding="utf-8", errors="replace")
    matches = re.findall(
        r'<text\b[^>]*\bx="([^"]+)"[^>]*\by="([^"]+)"[^>]*>(.*?)</text>',
        text,
        flags=re.DOTALL,
    )
    anchors: dict[tuple[float, float], list[str]] = {}
    for x_value, y_value, label in matches:
        value = re.sub(r"<[^>]+>", "", label)
        value = html.unescape(value).strip()
        if not value:
            continue
        key = (round(float(x_value), 3), round(float(y_value), 3))
        anchors.setdefault(key, []).append(value)
    duplicate_groups = {
        key: values
        for key, values in anchors.items()
        if len(values) > 1
    }
    examples = [
        f"({x:g},{y:g}): {' / '.join(values[:4])}"
        for (x, y), values in list(duplicate_groups.items())[:5]
    ]
    return {
        "duplicate_svg_text_anchor_groups": len(duplicate_groups),
        "duplicate_svg_text_examples": "; ".join(examples),
    }


def _automatic_score(row: dict[str, Any]) -> int | None:
    if row.get("status") == "prepared_not_run":
        return None
    if not row.get("width_px"):
        return 0
    score = 100
    dpi_x = float(row.get("dpi_x") or 0)
    width = int(row.get("width_px") or 0)
    height = int(row.get("height_px") or 0)
    aspect = float(row.get("aspect_ratio") or 0)
    non_white = float(row.get("non_white_fraction") or 0)
    size_mb = float(row.get("file_size_mb") or 0)
    physical_width = float(row.get("physical_width_in") or 0)
    if dpi_x < 300:
        score -= 30
    elif dpi_x < 450:
        score -= 8
    if min(width, height) < 600:
        score -= 15
    if aspect > 4.5 or aspect < 0.22:
        score -= 10
    if non_white < 0.01:
        score -= 45
    if size_mb > 10:
        score -= 10
    if physical_width > 8.5:
        score -= 12
    elif physical_width > 7.5:
        score -= 5
    if not row.get("has_pdf") and not row.get("has_svg"):
        score -= 12
    return max(0, min(100, score))


def _recommended_output(verdict: str) -> str:
    if verdict in {"不可替代Venn", "结果未达标", "不合理"}:
        return "不可作为主图"
    if verdict in {"结果未全达标", "需重绘"}:
        return "主图需重绘"
    if verdict in {"不满足方案"}:
        return "分析未完成"
    if verdict in {"需重排", "需修饰", "需补统计"}:
        return "主图需修饰"
    if verdict in {"需限定解释", "补充材料", "阴性结果"}:
        return "仅补充材料"
    if verdict == "可用":
        return "主图候选"
    return "人工复核"


def _summary(audit: pd.DataFrame) -> dict[str, Any]:
    available = audit[audit["status"] != "missing"]
    not_run = audit[audit["status"] == "prepared_not_run"]
    all_panel_keys = set(map(tuple, audit[["figure", "panel"]].to_numpy()))
    available_keys = set(
        map(tuple, available[["figure", "panel"]].to_numpy())
    )
    not_run_keys = set(map(tuple, not_run[["figure", "panel"]].to_numpy()))
    missing_keys = all_panel_keys - available_keys
    scores = pd.to_numeric(
        available["automatic_quality_score"],
        errors="coerce",
    ).dropna()
    oversize = available[
        pd.to_numeric(available["physical_width_in"], errors="coerce") > 7.5
    ]
    resolution_issues = available[
        (
            pd.to_numeric(available["dpi_x"], errors="coerce") < 300
        )
        | (
            pd.to_numeric(available["physical_width_in"], errors="coerce") > 7.5
        )
        | (
            ~available["has_pdf"].astype(bool)
            & ~available["has_svg"].astype(bool)
        )
    ]
    oversize_keys = set(
        map(tuple, oversize[["figure", "panel"]].to_numpy())
    )
    resolution_issue_keys = set(
        map(tuple, resolution_issues[["figure", "panel"]].to_numpy())
    )
    output_tiers = (
        available["recommended_output"]
        .value_counts()
        .to_dict()
    )
    severe_overlaps = available[
        available["label_overlap_severity"] == "严重"
    ]
    minor_overlaps = available[
        available["label_overlap_severity"] == "轻微"
    ]
    duplicate_anchor_panels = available[
        pd.to_numeric(
            available["duplicate_svg_text_anchor_groups"],
            errors="coerce",
        ).fillna(0)
        > 0
    ]
    major_issues = audit[
        audit["audit_verdict"].isin(
            {"需重绘", "不可替代Venn", "结果未达标", "结果未全达标", "不满足方案"}
        )
    ].drop_duplicates(subset=["figure", "panel"], keep="first")
    severe_overlap_keys = set(
        map(tuple, severe_overlaps[["figure", "panel"]].to_numpy())
    )
    minor_overlap_keys = set(
        map(tuple, minor_overlaps[["figure", "panel"]].to_numpy())
    )
    minor_overlap_keys -= severe_overlap_keys
    duplicate_anchor_keys = set(
        map(tuple, duplicate_anchor_panels[["figure", "panel"]].to_numpy())
    )
    return {
        "standard": "SCI IF>10 provisional figure screen",
        "target_journal": "not specified; exact journal requirements remain unverified",
        "overall_verdict": (
            "not_ready_for_if10"
            if (
                len(not_run)
                or len(major_issues)
                or len(resolution_issue_keys)
                or len(severe_overlap_keys)
                or len(duplicate_anchor_keys)
            )
            else "conditional_ready"
        ),
        "panel_entries": len(all_panel_keys),
        "image_entries": len(audit),
        "available_panels": len(available_keys),
        "prepared_not_run_panels": len(not_run_keys),
        "missing_panels": len(missing_keys),
        "mean_automatic_quality_score": float(scores.mean()) if len(scores) else None,
        "median_dpi_x": float(available["dpi_x"].median()) if len(available) else None,
        "panels_wider_than_7_5_inches": len(oversize_keys),
        "resolution_or_vector_issue_count": len(resolution_issue_keys),
        "recommended_output_counts": {
            str(key): int(value)
            for key, value in output_tiers.items()
        },
        "label_overlap_counts": {
            "严重": len(severe_overlap_keys),
            "轻微": len(minor_overlap_keys),
            "无": len(available_keys - severe_overlap_keys - minor_overlap_keys),
        },
        "duplicate_text_anchor_panel_count": len(duplicate_anchor_keys),
        "severe_label_overlap_panels": severe_overlaps[
            ["figure", "panel", "content", "label_overlap_notes"]
        ].to_dict(orient="records"),
        "minor_label_overlap_panels": minor_overlaps[
            ["figure", "panel", "content", "label_overlap_notes"]
        ].to_dict(orient="records"),
        "major_issue_count": len(major_issues),
        "major_issues": major_issues[
            ["figure", "panel", "content", "audit_verdict", "notes"]
        ].to_dict(orient="records"),
        "contact_sheets": {
            figure: str(
                Path("10_reports")
                / "figure_quality_audit"
                / "contact_sheets"
                / f"{figure}_contact_sheet.png"
            )
            for figure in sorted(audit["figure"].unique())
        },
        "csv": "10_reports/figure_quality_audit/figure_quality_audit.csv",
        "markdown": "10_reports/figure_quality_audit/figure_quality_audit.md",
    }


def _contact_sheet(
    aliases: list[PanelAlias],
    output_root: Path,
    destination: Path,
) -> Path:
    grouped: dict[str, list[PanelAlias]] = {}
    for alias in aliases:
        grouped.setdefault(alias.panel, []).append(alias)
    panels = list(grouped.items())
    columns = 3
    rows = math.ceil(len(panels) / columns)
    cell_width = 1100
    cell_height = 820
    gap = 34
    margin = 46
    canvas = Image.new(
        "RGB",
        (
            columns * cell_width + (columns + 1) * gap + 2 * margin,
            rows * cell_height + (rows + 1) * gap + 2 * margin,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    label_font = ImageFont.truetype(_font_path(bold=True), 54)
    small_font = ImageFont.truetype(_cjk_font_path(), 24)
    for index, (panel, panel_aliases) in enumerate(panels):
        row, column = divmod(index, columns)
        x = margin + gap + column * (cell_width + gap)
        y = margin + gap + row * (cell_height + gap)
        draw.text((x, y), panel, fill="#111111", font=label_font)
        image_top = y + 72
        image_height = cell_height - 96
        available = [
            (alias, (output_root / alias.source).resolve())
            for alias in panel_aliases
            if (output_root / alias.source).exists()
        ]
        if not available:
            draw.text(
                (x + 20, image_top + 100),
                "MISSING",
                fill="#b22222",
                font=label_font,
            )
            continue
        if len(available) == 1:
            sub_width, sub_height = cell_width, image_height
            positions = [(x, image_top)]
        else:
            sub_width = cell_width // len(available)
            sub_height = image_height
            positions = [
                (x + item * sub_width, image_top)
                for item in range(len(available))
            ]
        for (alias, source), (sub_x, sub_y) in zip(available, positions):
            with Image.open(source) as image:
                image = image.convert("RGB")
                image = ImageOps.contain(
                    image,
                    (max(1, sub_width - 12), max(1, sub_height - 42)),
                    Image.Resampling.LANCZOS,
                )
                paste_x = sub_x + (sub_width - image.width) // 2
                paste_y = sub_y + (sub_height - 36 - image.height) // 2
                canvas.paste(image, (paste_x, paste_y))
            draw.text(
                (sub_x + 4, sub_y + sub_height - 30),
                alias.description[:42],
                fill="#555555",
                font=small_font,
            )
    ensure_dir(destination.parent)
    canvas.save(destination, dpi=(300, 300), optimize=True)
    return destination


def _font_path(bold: bool = False) -> str:
    from matplotlib import font_manager

    candidate = "DejaVu Sans"
    return font_manager.findfont(
        font_manager.FontProperties(family=candidate, weight="bold" if bold else "normal")
    )


def _cjk_font_path() -> str:
    candidates = [
        Path("C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/simhei.ttf"),
        Path("C:/Windows/Fonts/NotoSansSC-Regular.otf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    return _font_path()


def _render_audit_markdown(
    audit: pd.DataFrame,
    summary: dict[str, Any],
    classified_root: Path,
) -> str:
    lines = [
        "# Figure quality audit for SCI IF>10 screening",
        "",
        f"- Classified root: `{classified_root}`",
        f"- Overall verdict: **{summary['overall_verdict']}**",
        f"- Panel entries: {summary['panel_entries']}",
        f"- Image aliases reviewed: {summary.get('image_entries', summary['panel_entries'])}",
        f"- Available panels: {summary['available_panels']}",
        f"- Prepared but not run: {summary['prepared_not_run_panels']}",
        f"- Major issues: {summary['major_issue_count']}",
        (
            f"- Panels wider than 7.5 inches at native scale: "
            f"{summary['panels_wider_than_7_5_inches']}"
        ),
        (
            f"- Resolution/vector issues: "
            f"{summary['resolution_or_vector_issue_count']}"
        ),
        (
            "- Recommended output tiers: "
            + ", ".join(
                f"{tier}={count}"
                for tier, count in summary["recommended_output_counts"].items()
            )
        ),
        (
            "- Label overlap review: "
            f"severe={summary['label_overlap_counts']['严重']}, "
            f"minor={summary['label_overlap_counts']['轻微']}, "
            f"none={summary['label_overlap_counts']['无']}"
        ),
        (
            "- Duplicate SVG text anchors: "
            f"{summary['duplicate_text_anchor_panel_count']}"
        ),
        (
            f"- Mean automated raster score: "
            f"{summary['mean_automatic_quality_score']:.1f}"
            if summary.get("mean_automatic_quality_score") is not None
            else "- Mean automated raster score: NA"
        ),
        "",
        "The automated score checks raster resolution, size, aspect ratio, blankness "
        "and file size. Scientific reasonableness and aesthetics are assigned by "
        "panel-specific review and do not represent journal acceptance.",
        "",
        "## Panel Review",
        "",
        "| Figure | Panel | Content | Scientific | Aesthetics | Verdict | Recommended output | Notes |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in audit.itertuples(index=False):
        lines.append(
            f"| {row.figure} | {row.panel} | {row.content} | "
            f"{row.scientific_reasonableness} | {row.aesthetic_level} | "
            f"{row.audit_verdict} | {row.recommended_output} | {row.notes} |"
        )
    lines += [
        "",
        "## Required Improvements",
        "",
        "- Run the 100 ns GROMACS production simulation and replace all `NOT_RUN` "
        "panels with measured trajectory and MM-PBSA outputs.",
            "- Report calibration intercept, slope, Brier and H-L together; do not "
            "treat p>0.05 as a standalone pass/fail test.",
        "- Rebuild Figure 1d/e only if licensed GeneCards, OMIM, TTD exports and "
        "independent ChEMBL/STITCH evidence become available.",
            "- Keep the R CellChat version, species database and cell-level "
            "probability limitations in the figure legend.",
            "- Re-render predicted perturbation-response panels from source data "
            "and state the R scTenifoldKnk parameters in the figure legend.",
        "- Assemble final multi-panel figures, add lowercase panel letters, export "
        "vector PDF/SVG for line art, and provide source data tables.",
        "- Add confidence intervals, effect sizes and exact FDR to key statistical panels.",
        "",
        "## Label Overlap Review",
        "",
        "| Figure | Panel | Content | Severity | Notes |",
        "|---|---|---|---|---|",
    ]
    overlap_rows = audit[
        (audit["label_overlap_severity"] != "无")
        | (
            pd.to_numeric(
                audit["duplicate_svg_text_anchor_groups"],
                errors="coerce",
            ).fillna(0)
            > 0
        )
    ]
    if overlap_rows.empty:
        lines.append("| - | - | - | 无 | 原分辨率视觉检查未发现标签重叠。 |")
    else:
        for row in overlap_rows.itertuples(index=False):
            examples = str(row.duplicate_svg_text_examples or "")
            note = row.label_overlap_notes
            if examples:
                note = f"{note} SVG重复锚点: {examples}"
            lines.append(
                f"| {row.figure} | {row.panel} | {row.content} | "
                f"{row.label_overlap_severity} | {note} |"
            )
    lines += [
        "",
        "## Contact Sheets",
        "",
    ]
    for figure, path in summary["contact_sheets"].items():
        lines.append(f"- `{figure}`: `{path}`")
    return "\n".join(lines) + "\n"
