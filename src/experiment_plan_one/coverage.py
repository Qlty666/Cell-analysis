"""Weighted implementation and result-coverage audit for experiment plan one."""

from __future__ import annotations

import importlib.util
import re
import shutil
import subprocess
import sys
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
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "d", "化合物靶点数据库 Venn", 3, 1.00, "external_data", "内置 SwissTargetPrediction、ChEMBL、STITCH 与 SEA；任一联网源可用或使用缓存/本地预测表即可形成可审计多来源。"),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "e", "疾病靶点数据库 Venn", 4, 0.80, "licensed_data", "GeneCards/OMIM/TTD 需要本地授权导出；无授权时可使用开放证据源替代并明确标注。"),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "f", "化合物靶点 ∩ 疾病靶点", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "g", "KEGG Top10 富集", 3, 1.0),
    PanelRequirement("Figure1_化合物表征_靶点预测与通路富集", "h", "GO BP/CC/MF 富集", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "a", "STRING PPI 网络", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "b", "模块网络", 2, 0.95, "algorithm_equivalent", "按原方案运行 R MCL 并记录 inflation；若明确改用 Louvain，须单独登记方法替代。"),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "c", "Degree Top20", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "d", "Betweenness 排名", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "e", "MCC ∩ Degree", 3, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "f", "GSE89632 候选基因热图", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "g", "GSE49541 纤维化箱线图", 2, 1.0),
    PanelRequirement("Figure2_PPI网络与枢纽基因初步筛选", "h", "GSE164441 表达验证", 2, 1.0, "endpoint_mismatch", "GSE164441 为肝癌与癌旁对照，必须按不同终点解释。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "a", "11 模型 AUC 热图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "b", "训练集 ROC", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "c", "GSE49541 跨终点探索 ROC", 4, 0.90, "endpoint_mismatch", "纤维化分期与健康/NAFLD不是同一终点；报告实际AUC但不设成功阈值。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "d", "同终点候选与HCC扩展 ROC", 4, 0.90, "endpoint_mismatch", "GSE135251已登记既往分析暴露；GSE164441为HCC配对癌旁终点。AUC阈值不作为验收门禁。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "e", "校准曲线", 4, 0.90, "outcome_dependent", "联合报告校准截距、斜率、Brier和H-L结果；H-L p>0.05 单独不证明校准合格。"),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "f", "SHAP 条形图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "g", "SHAP 蜂群图", 3, 1.0),
    PanelRequirement("Figure3_机器学习模型构建与SHAP核心特征", "h", "核心基因与纤维化关联", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "a", "小鼠 UMAP", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "b", "细胞类型 marker 气泡图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "c", "核心基因各亚群表达", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "d", "核心基因 UMAP 特征图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "e", "细胞组成图", 2, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f", "细胞通讯网络", 4, 0.90, "method_equivalent", "按原方案运行 R CellChat；细胞级概率只作描述，组间统计仍需独立生物单位。"),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "g", "预测扰动响应", 3, 0.95, "external_tool", "按原方案运行 R scTenifoldKnk；工具缺失时标记 blocked，不静默替换为局部GRN。"),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "h", "虚拟敲除后富集", 2, 0.95),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "i", "人类疾病谱 UMAP", 3, 1.0),
    PanelRequirement("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "j", "人类核心基因细胞类型验证", 3, 1.0),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "a", "3D 对接构象", 3, 0.90, "visualization_equivalent", "使用可复现的 3D 构象与结合口袋可视化；Discovery Studio 为可选复核工具。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "b", "2D 相互作用图", 3, 0.90, "visualization_equivalent", "优先运行 PLIP 并保存机器可读相互作用表；几何候选不得写成经验证的氢键、盐桥或π堆积。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "c", "对接能量热图", 2, 1.0),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "d", "蛋白 RMSD", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "e", "配体 RMSD", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "f", "关键残基 RMSF", 4, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "g", "回旋半径 Rg", 3, 1.0, "gromacs", "需要 md.run=true 或外部 GROMACS 完成 100 ns 轨迹。"),
    PanelRequirement("Figure5_分子对接与分子动力学模拟", "h", "MM-PBSA 自由能分解", 5, 0.75, "gmx_mmpbsa", "需要残基分解、TOTAL列、单位和误差；仅总能量不能交付本面板。"),
)


GOOD_VERDICTS = {
    "可用",
    "主图候选",
    "需修饰",
    "需重排",
    "阴性结果",
    "结果未达标",
    "结果未全达标",
    "需限定解释",
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
        if requirement.dependency in {
            "external_data",
            "licensed_data",
            "external_tool",
            "gromacs",
            "gmx_mmpbsa",
            "method_equivalent",
            "algorithm_equivalent",
            "visualization_equivalent",
        }:
            return "blocked"
        return "missing_audit"
    raw_status = str(review.get("status") or "")
    verdict = str(review.get("audit_verdict") or "")
    if raw_status == "missing":
        return (
            "blocked"
            if requirement.dependency
            in {
                "external_data",
                "licensed_data",
                "external_tool",
                "gromacs",
                "gmx_mmpbsa",
                "method_equivalent",
                "algorithm_equivalent",
                "visualization_equivalent",
            }
            else "not_run"
        )
    if raw_status == "prepared_not_run":
        return "not_run"
    if verdict in {"不满足方案", "未运行"} or raw_status == "not_run":
        return "not_run"
    if verdict in {"不可替代Venn", "不合理", "需重绘", "不满足方案"}:
        return "needs_revision"
    if verdict == "需限定解释":
        return "exploratory"
    if verdict == "阴性结果":
        return "valid_negative"
    return "available"


def _command_output(command: list[str], timeout: int = 20) -> str:
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
    except Exception:  # noqa: BLE001
        return ""
    return (result.stdout or "") + "\n" + (result.stderr or "")


def audit_plan_environment() -> dict[str, Any]:
    """Record exact external-tool versions required by the plan."""
    root = Path(__file__).resolve().parents[2]
    vina_path = root / "dock" / "tools" / "vina.exe"
    if not vina_path.exists():
        discovered = shutil.which("vina") or shutil.which("vina.exe")
        vina_path = Path(discovered) if discovered else None
    vina_text = (
        _command_output([str(vina_path), "--version"])
        if vina_path
        else ""
    )
    vina_match = re.search(r"v?(\d+\.\d+(?:\.\d+)?)", vina_text)

    gmx_path = shutil.which("gmx") or shutil.which("gmx.exe")
    gmx_text = (
        _command_output([str(gmx_path), "--version"])
        if gmx_path
        else ""
    )
    gmx_match = re.search(r"GROMACS version:\s*([^\s]+)", gmx_text)

    gmx_mmpbsa = (
        shutil.which("gmx_MMPBSA")
        or shutil.which("gmx_MMPBSA.py")
        or shutil.which("gmx_MMPBSA.exe")
    )
    gmx_mmpbsa_text = (
        _command_output([str(gmx_mmpbsa), "--version"])
        if gmx_mmpbsa
        else ""
    )
    rscript = shutil.which("Rscript") or shutil.which("Rscript.exe")
    r_package_text = (
        _command_output(
            [
                str(rscript),
                "--vanilla",
                "-e",
                (
                    'cat(paste(c("MCL","CellChat","scTenifoldKnk"), '
                    'vapply(c("MCL","CellChat","scTenifoldKnk"), '
                    "requireNamespace, logical(1), quietly=TRUE), sep=':'))"
                ),
            ],
            timeout=30,
        )
        if rscript
        else ""
    )
    r_packages = {}
    for token in re.findall(r"([A-Za-z0-9.]+):(TRUE|FALSE)", r_package_text):
        r_packages[token[0]] = token[1] == "TRUE"
    plip_candidates = [
        shutil.which("plip"),
        shutil.which("plip.exe"),
        str(Path(sys.executable).parent / (
            "plip.exe" if sys.platform.startswith("win") else "plip"
        )),
        str(Path(sys.executable).parent / "Scripts" / (
            "plip.exe" if sys.platform.startswith("win") else "plip"
        )),
        str(Path(sys.executable).parent / "bin" / (
            "plip.exe" if sys.platform.startswith("win") else "plip"
        )),
    ]
    plip_path = next(
        (
            Path(candidate)
            for candidate in plip_candidates
            if candidate and Path(candidate).exists()
        ),
        None,
    )
    try:
        from openbabel import openbabel as ob

        openbabel_inchikey = bool(ob.OBConversion().FindFormat("inchikey"))
    except Exception:  # noqa: BLE001
        openbabel_inchikey = False
    mmpbsa_candidates = [
        shutil.which("gmx_MMPBSA"),
        shutil.which("gmx_MMPBSA.py"),
        shutil.which("gmx_MMPBSA.exe"),
        str(Path(sys.executable).parent / (
            "gmx_MMPBSA.exe"
            if sys.platform.startswith("win")
            else "gmx_MMPBSA"
        )),
        str(Path(sys.executable).parent / "Scripts" / (
            "gmx_MMPBSA.exe"
            if sys.platform.startswith("win")
            else "gmx_MMPBSA"
        )),
        str(Path(sys.executable).parent / "bin" / (
            "gmx_MMPBSA.exe"
            if sys.platform.startswith("win")
            else "gmx_MMPBSA"
        )),
    ]
    gmx_mmpbsa = next(
        (
            candidate
            for candidate in mmpbsa_candidates
            if candidate and Path(candidate).exists()
        ),
        None,
    )
    gmx_mmpbsa_text = (
        _command_output([str(gmx_mmpbsa), "--version"])
        if gmx_mmpbsa
        else ""
    )
    return {
        "vina": {
            "path": str(vina_path) if vina_path else "",
            "version": vina_match.group(1) if vina_match else "",
            "expected": "1.2.3",
            "match": bool(
                vina_match and vina_match.group(1).startswith("1.2.3")
            ),
        },
        "gromacs": {
            "path": str(gmx_path) if gmx_path else "",
            "version": gmx_match.group(1) if gmx_match else "",
            "expected": "2022",
            "match": bool(
                gmx_match and gmx_match.group(1).startswith("2022")
            ),
        },
        "gmx_mmpbsa": {
            "path": str(gmx_mmpbsa) if gmx_mmpbsa else "",
            "available": bool(gmx_mmpbsa),
            "version": (
                gmx_mmpbsa_text.strip().splitlines()[0]
                if gmx_mmpbsa_text.strip()
                else ""
            ),
        },
        "cellchat": {
            "rscript": str(rscript) if rscript else "",
            "available": bool(r_packages.get("CellChat")),
        },
        "mcl": {
            "rscript": str(rscript) if rscript else "",
            "available": bool(r_packages.get("MCL")),
        },
        "r_sctenifoldknk": {
            "rscript": str(rscript) if rscript else "",
            "available": bool(r_packages.get("scTenifoldKnk")),
        },
        "plip": {
            "path": str(plip_path) if plip_path else "",
            "available": bool(plip_path),
            "runtime_backend": "python_api_in_process",
            "openbabel_inchikey_format": openbabel_inchikey,
            "inchikey_writer_patch_required": bool(
                plip_path and not openbabel_inchikey
            ),
        },
        "r_packages": r_packages,
        "python_scripts": {
            "plip": str(plip_path) if plip_path else "",
            "gmx_mmpbsa": str(gmx_mmpbsa) if gmx_mmpbsa else "",
        },
        "python_packages": {
            name: bool(importlib.util.find_spec(name))
            for name in (
                "rdkit",
                "scanpy",
                "shap",
                "plip",
                "scTenifold",
                "torch",
                "posebusters",
                "meeko",
            )
        },
    }


def audit_plan_coverage(output_root: Path) -> dict[str, Any]:
    """Score code support and actual result coverage against the 42 plan panels."""
    output_root = output_root.resolve()
    reviews = _read_audit_reviews(output_root)
    rows: list[dict[str, Any]] = []
    for requirement in PANEL_REQUIREMENTS:
        review = reviews.get((requirement.figure, requirement.panel))
        status = _status(requirement, review)
        verdict = str((review or {}).get("audit_verdict") or "")
        if status in {"available", "valid_negative"} and verdict in GOOD_VERDICTS:
            current_credit = 1.0
        elif status == "exploratory" or (
            status == "available" and verdict == "补充材料"
        ):
            current_credit = 0.5
        else:
            current_credit = 0.0
        performance_target_met: bool | None = None
        performance_target_applicable = False
        performance_target_status = "not_applicable"
        if requirement.figure == "Figure3_机器学习模型构建与SHAP核心特征":
            if requirement.panel == "e":
                performance_target_applicable = True
                # A missing/ambiguous audit is never evidence of calibration.
                performance_target_status = "unknown"
                if status == "available" and verdict in {"可用", "主图候选"}:
                    performance_target_met = True
                    performance_target_status = "met"
                elif status == "available" and verdict in {"结果未达标", "结果未全达标", "不合理"}:
                    performance_target_met = False
                    performance_target_status = "not_met"
            elif requirement.panel in {"c", "d"}:
                # The available validation cohorts use fibrosis or HCC
                # endpoints rather than the planned healthy-versus-NAFLD one.
                performance_target_met = None
                performance_target_status = "not_evaluated"
        projected_credit = 1.0 if requirement.implementation >= 0.90 else requirement.implementation
        method_valid = bool(
            status in {"available", "valid_negative", "exploratory"}
            and verdict in GOOD_VERDICTS
        )
        evidence_sufficient = bool(
            status in {"available", "valid_negative"}
            and verdict
            in {"可用", "主图候选", "阴性结果"}
        )
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
                "file_generated": status
                in {"available", "prepared_not_run", "needs_revision"},
                "method_valid": method_valid,
                "evidence_sufficient": evidence_sufficient,
                "scientific_outcome_status": (
                    "valid_negative"
                    if verdict == "阴性结果"
                    else "reported_with_limitations"
                    if verdict == "需限定解释"
                    else "available_result"
                    if method_valid
                    else "unknown"
                ),
                "performance_target_met": performance_target_met,
                "performance_target_applicable": performance_target_applicable,
                "performance_target_status": performance_target_status,
                "note": requirement.note,
            }
        )
    frame = pd.DataFrame(rows)
    total_weight = float(frame["weight"].sum())
    weighted = lambda column: float(
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
        "standard": "experiment-plan-one weighted implementation/output coverage (not publication readiness)",
        "publication_readiness": "not_assessed",
        "evidence_state_vocabulary": [
            "blocked",
            "not_run",
            "failed",
            "exploratory",
            "valid_negative",
            "valid_positive",
            "unknown",
        ],
        "blocked_panels": int(frame["current_status"].eq("blocked").sum()),
        "not_run_panels": int(frame["current_status"].eq("not_run").sum()),
        "method_valid_panels": int(frame["method_valid"].sum()),
        "evidence_sufficient_panels": int(frame["evidence_sufficient"].sum()),
        "file_generated_panels": int(frame["file_generated"].sum()),
        "performance_targets_unknown": int(frame["performance_target_status"].eq("unknown").sum()),
        "panels": len(frame),
        "weight_total": total_weight,
        "implementation_completion_percent": round(implementation_percent, 2),
        "current_result_completion_percent": round(current_percent, 2),
        "projected_completion_percent_with_prerequisites": round(projected_percent, 2),
        "above_90_implementation": implementation_percent >= 90.0,
        "above_90_projected_with_prerequisites": projected_percent >= 90.0,
        "current_90pct_ready": current_percent >= 90.0,
        "performance_targets_met": int(
            frame["performance_target_status"].eq("met").sum()
        ),
        "performance_targets_evaluable": int(
            frame["performance_target_status"].isin({"met", "not_met"}).sum()
        ),
        "performance_targets_not_evaluated": int(
            frame["performance_target_status"].eq("not_evaluated").sum()
        ),
        "performance_targets_not_met": frame[
            frame["performance_target_status"].eq("not_met")
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
        "environment": audit_plan_environment(),
        "csv": "10_reports/plan_coverage/plan_coverage.csv",
        "json": "10_reports/plan_coverage/plan_coverage.json",
        "markdown": "10_reports/plan_coverage/plan_coverage.md",
    }
    out_dir = ensure_dir(output_root / "10_reports" / "plan_coverage")
    frame.to_csv(out_dir / "plan_coverage.csv", index=False, encoding="utf-8-sig")
    write_json(out_dir / "plan_coverage.json", summary)
    write_json(
        out_dir / "environment_audit.json",
        summary["environment"],
    )
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
        f"- Panels with generated files: {summary['file_generated_panels']}/{summary['panels']}",
        f"- Panels with method-valid audits: {summary['method_valid_panels']}/{summary['panels']}",
        f"- Panels with sufficient evidence audits: {summary['evidence_sufficient_panels']}/{summary['panels']}",
        f"- Blocked panels: {summary['blocked_panels']}",
        f"- Not-run panels: {summary['not_run_panels']}",
        (
            "- Projected coverage with required external data/tools: "
            f"{summary['projected_completion_percent_with_prerequisites']:.2f}%"
        ),
        (
            "- Performance targets met: "
            f"{summary['performance_targets_met']}/"
            f"{summary['performance_targets_evaluable']} evaluable"
        ),
        (
            "- Performance targets not evaluated because the available "
            "cohort measures a different endpoint: "
            f"{summary['performance_targets_not_evaluated']}"
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
