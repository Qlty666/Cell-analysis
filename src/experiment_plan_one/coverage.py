"""Weighted implementation and result-coverage audit for experiment plan one."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .common import ensure_dir, write_json


@dataclass(frozen=True)
class PanelRequirement:
    figure: str
    panel: str
    label: str
    weight: float
    implementation: float
    dependency: str = ""
    note: str = ""


PANEL_REQUIREMENTS: tuple[PanelRequirement, ...] = (
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "a", "研究全局流程图", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "b", "6PPD-Q 2D 化学结构", 2, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "c", "3D 结构 + 理化性质", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "d", "化合物靶点数据库 Venn", 3, 0.80, "external_data", "需要至少两个可用化合物靶点来源；本地预测表可补齐。"),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "e", "疾病靶点数据库 Venn", 4, 0.80, "licensed_data", "GeneCards/OMIM/TTD 需要本地授权导出；无授权时可使用开放证据源替代并明确标注。"),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "f", "化合物靶点 ∩ 疾病靶点", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "g", "KEGG Top10 富集", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "h", "GO BP/CC/MF 富集", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "a", "STRING PPI 网络", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "b", "MCL/模块网络", 2, 0.95, "algorithm_equivalent", "当前实现为可复现社区检测；MCL 为可选增强。"),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "c", "Degree Top20", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "d", "Betweenness 排名", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "e", "MCC ∩ Degree", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "f", "GSE89632 候选基因热图", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "g", "GSE49541 纤维化箱线图", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "h", "GSE164441 表达验证", 2, 1.0, "endpoint_mismatch", "GSE164441 为肝癌与癌旁对照，必须按不同终点解释。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "a", "11 模型 AUC 热图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "b", "训练集 ROC", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "c", "GSE49541 外部 ROC", 4, 0.90, "outcome_dependent", "临床终点为纤维化分期；是否达到 AUC 目标取决于数据和模型表现。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "d", "外部 NAFLD ROC", 4, 0.90, "outcome_dependent", "外部泛化是否达到 AUC 目标不能由代码保证。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "e", "校准曲线", 4, 0.90, "outcome_dependent", "代码支持 sigmoid/isotonic 校准选择，但统计校准仍取决于样本量。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "f", "SHAP 条形图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "g", "SHAP 蜂群图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "h", "核心基因与纤维化关联", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "a", "小鼠 UMAP", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "b", "细胞类型 marker 气泡图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "c", "核心基因各亚群表达", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "d", "核心基因 UMAP 特征图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "e", "细胞组成图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f", "细胞通讯网络", 4, 0.90, "method_equivalent", "默认使用带置换检验和 FDR 的配体-受体评分；完整 CellChat 为可选外部工具。"),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "g", "虚拟敲除 Top10/网络/UMAP", 3, 0.95, "external_tool", "scTenifoldKnk 可用时执行；否则保留明确的网络模拟结果。"),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "h", "虚拟敲除后富集", 2, 0.95),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "i", "人类疾病谱 UMAP", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "j", "人类核心基因细胞类型验证", 3, 1.0),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "a", "3D 对接构象", 3, 0.90, "visualization_equivalent", "使用可复现的 3D 构象与结合口袋可视化；Discovery Studio 为可选复核工具。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "b", "2D 相互作用图", 3, 0.90, "visualization_equivalent", "按氢键、疏水、π-π、盐桥分类；PLIP/Discovery Studio 为可选复核工具。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "c", "对接能量热图", 2, 1.0),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "d", "蛋白 RMSD", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "e", "配体 RMSD", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "f", "关键残基 RMSF", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "g", "回旋半径 Rg", 3, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "h", "MM-PBSA 自由能分解", 5, 0.75, "gmx_mmpbsa", "需要可用的 gmx_MMPBSA 或等效命令；未安装时明确标记不可用。"),
)


GOOD_VERDICTS = {
    "可用",
    "主图候选",
    "需修饰",
    "需重排",
    "需补统计",
    "需限定解释",
    "阴性结果",
    "结果未达标",
    "结果未全达标",
}


def _read_audit_reviews(output_root: Path) -> dict[tuple[str, str], dict[str, Any]]:
    path = output_root / "10_reports" / "figure_quality_audit" / "figure_quality_audit.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    records: dict[tuple[str, str], dict[str, Any]] = {}
    for row in frame.to_dict("records"):
        records[(str(row.get("figure")), str(row.get("panel")))] = row
    return records


def _status(
    requirement: PanelRequirement,
    review: dict[str, Any] | None,
) -> str:
    if review is None:
        return "missing_audit"
    raw_status = str(review.get("status") or "")
    verdict = str(review.get("audit_verdict") or "")
    if raw_status == "missing":
        return "missing"
    if raw_status == "prepared_not_run":
        return "prepared_not_run"
    if verdict in {"不可替代Venn", "不合理", "需重绘", "不满足方案"}:
        return "needs_revision"
    return "available"


def audit_plan_coverage(output_root: Path) -> dict[str, Any]:
    """Score code support and actual result coverage against the 42 plan panels."""
    output_root = output_root.resolve()
    reviews = _read_audit_reviews(output_root)
    rows: list[dict[str, Any]] = []
    for requirement in PANEL_REQUIREMENTS:
        review = reviews.get((requirement.figure, requirement.panel))
        status = _status(requirement, review)
        verdict = str((review or {}).get("audit_verdict") or "")
        if status == "available" and verdict in GOOD_VERDICTS:
            current_credit = 1.0
        elif status == "available" and verdict == "补充材料":
            current_credit = 0.5
        else:
            current_credit = 0.0
        performance_target_met = True
        if (
            requirement.figure
            == "Figure3_机器学习模型构建与SHAP核心特征"
            and requirement.panel in {"c", "d", "e"}
        ):
            performance_target_met = verdict not in {
                "结果未达标",
                "结果未全达标",
                "不合理",
            }
        projected_credit = 1.0 if requirement.implementation >= 0.90 else requirement.implementation
        rows.append(
            {
                "figure": requirement.figure,
                "panel": requirement.panel,
                "label": requirement.label,
                "weight": requirement.weight,
                "implementation_support": requirement.implementation,
                "current_status": status,
                "audit_verdict": verdict,
                "dependency": requirement.dependency,
                "current_credit": current_credit,
                "projected_credit": projected_credit,
                "performance_target_met": performance_target_met,
                "note": requirement.note,
            }
        )
    frame = pd.DataFrame(rows)
    total_weight = float(frame["weight"].sum())
    weighted = lambda column: float(  # noqa: E731
        (frame[column] * frame["weight"]).sum() / total_weight
    )
    implementation_percent = weighted("implementation_support") * 100.0
    current_percent = weighted("current_credit") * 100.0
    projected_percent = weighted("projected_credit") * 100.0
    blockers = (
        frame[frame["current_credit"] < 1.0]
        .groupby(["dependency", "current_status"], as_index=False)
        .agg(
            panels=("panel", lambda values: ", ".join(sorted(set(map(str, values))))),
            count=("panel", "size"),
        )
        .to_dict("records")
    )
    summary = {
        "standard": "experiment-plan-one weighted coverage",
        "panels": int(len(frame)),
        "weight_total": total_weight,
        "implementation_completion_percent": round(implementation_percent, 2),
        "current_result_completion_percent": round(current_percent, 2),
        "projected_completion_percent_with_prerequisites": round(projected_percent, 2),
        "above_90_implementation": implementation_percent >= 90.0,
        "above_90_projected_with_prerequisites": projected_percent >= 90.0,
        "current_90pct_ready": current_percent >= 90.0,
        "performance_targets_met": int(
            frame["performance_target_met"].sum()
        ),
        "performance_targets_not_met": frame[
            ~frame["performance_target_met"]
        ][
            ["figure", "panel", "label", "audit_verdict", "note"]
        ].to_dict("records"),
        "blocking_groups": blockers,
        "missing_or_unmet_panels": frame[frame["current_credit"] < 1.0][
            [
                "figure",
                "panel",
                "label",
                "current_status",
                "audit_verdict",
                "dependency",
                "note",
            ]
        ].to_dict("records"),
        "csv": "10_reports/plan_coverage/plan_coverage.csv",
        "json": "10_reports/plan_coverage/plan_coverage.json",
        "markdown": "10_reports/plan_coverage/plan_coverage.md",
    }
    out_dir = ensure_dir(output_root / "10_reports" / "plan_coverage")
    frame.to_csv(out_dir / "plan_coverage.csv", index=False, encoding="utf-8-sig")
    write_json(out_dir / "plan_coverage.json", summary)
    (out_dir / "plan_coverage.md").write_text(
        _render_markdown(frame, summary),
        encoding="utf-8",
    )
    return summary


def _render_markdown(rows: pd.DataFrame, summary: dict[str, Any]) -> str:
    lines = [
        "# Experiment plan one weighted coverage",
        "",
        f"- Panels: {summary['panels']}",
        f"- Code/implementation coverage: {summary['implementation_completion_percent']:.2f}%",
        f"- Current executed-result coverage: {summary['current_result_completion_percent']:.2f}%",
        (
            "- Projected coverage with required external data/tools: "
            f"{summary['projected_completion_percent_with_prerequisites']:.2f}%"
        ),
        (
            "- Performance targets met: "
            f"{summary['performance_targets_met']}/{summary['panels']}"
        ),
        "",
        "## Panel Matrix",
        "",
        "| Figure | Panel | Content | Weight | Implementation | Current | Dependency |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for row in rows.itertuples(index=False):
        lines.append(
            f"| {row.figure} | {row.panel} | {row.label} | "
            f"{row.weight:g} | {row.implementation_support:.2f} | "
            f"{row.current_credit:.2f} | {row.dependency or '-'} |"
        )
    lines += ["", "## Required Actions", ""]
    if summary["missing_or_unmet_panels"]:
        for item in summary["missing_or_unmet_panels"]:
            lines.append(
                f"- {item['figure']} {item['panel']} ({item['label']}): "
                f"{item['current_status']} / {item['audit_verdict'] or 'no verdict'}; "
                f"{item['note']}"
            )
    else:
        lines.append("- All weighted panels are currently satisfied.")
    if summary["performance_targets_not_met"]:
        lines += ["", "## Performance Targets Not Met", ""]
        for item in summary["performance_targets_not_met"]:
            lines.append(
                f"- {item['figure']} {item['panel']} ({item['label']}): "
                f"{item['audit_verdict']}; {item['note']}"
            )
    return "\n".join(lines) + "\n"
