"""Automated quality screening and readable contact sheets for plan figures."""

from __future__ import annotations

import html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

from .classify import CLASSIFIED_DIR, PANEL_ALIASES, PanelAlias
from .common import ensure_dir, write_json


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
        "已改用可复现的Louvain社区检测，得到4个主要模块；应把方法和种子/上下文节点定义写入图注。",
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
        "交叉验证ROC绘制方式正确。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "c"): FigureReview(
        "部分合理",
        "良好",
        "结果未达标",
        "GSE49541检测的是纤维化分期，不是同一健康/NAFLD终点，AUC 0.459未达方案目标。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "d"): FigureReview(
        "部分合理",
        "良好",
        "结果未全达标",
        "GSE164441达到0.880；补充NAFLD队列GSE135251仅0.742，外部NAFLD泛化性不足。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "e"): FigureReview(
        "不合理",
        "中等",
        "结果未达标",
        "Hosmer-Lemeshow p=4.14e-11，校准明显失败，不能作为合格临床预测模型图。",
    ),
    ("Figure3_机器学习模型构建与SHAP核心特征", "f"): FigureReview(
        "合理",
        "良好",
        "可用",
        "SHAP排名清晰，但仅63例训练，应报告稳定性/重采样波动。",
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
        "可用",
        "GPAT3在主要细胞中表达稀疏，统计零膨胀明显，应与dropout一起解释。",
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
        "可用",
        "组成图清楚，图例已移至下方；建议增加每样本比例统计和β回归/FDR。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "f"): FigureReview(
        "部分合理",
        "中等",
        "需限定解释",
        "已改为显式配体-受体评分网络并减少边数，但仍不是完整CellChat置换检验；主图标题和图注必须保留CellChat-like限定。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "g"): FigureReview(
        "合理",
        "偏低",
        "补充材料",
        "虚拟敲除表格、网络和UMAP均存在，但当前外部包裹器输出为160 dpi且信息过密，主图应改用重绘版本或移至补充材料；正文需说明这是网络模拟而非湿实验敲除。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "h"): FigureReview(
        "合理",
        "偏低",
        "补充材料",
        "未发现显著GO/KEGG条目，纯文字面板不宜进入主图，应作为阴性结果写入正文或补充材料。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "i"): FigureReview(
        "合理",
        "良好",
        "可用",
        "人类疾病谱UMAP清楚，建议标注疾病分期比例。",
    ),
    ("Figure4_单细胞图谱_细胞通讯与虚拟扰动", "j"): FigureReview(
        "合理",
        "良好",
        "可用",
        "已重构为紧凑的基因×细胞类型×疾病阶段热图，显示相对Healthy的中位表达差并标注Kruskal-Wallis FDR。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "a"): FigureReview(
        "部分合理",
        "偏低",
        "需重绘",
        "当前为Cα散点与配体坐标，不是高质量蛋白-配体3D构象图；需要PyMOL/ChimeraX cartoon、结合口袋、相互作用残基和比例提示。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "b"): FigureReview(
        "部分合理",
        "偏低",
        "需重绘",
        "当前为近邻残基放射图，不等同于Discovery Studio/PLIP二维相互作用图；需要按氢键、疏水、π-π等相互作用类型区分。",
    ),
    ("Figure5_分子对接与分子动力学模拟", "c"): FigureReview(
        "合理",
        "良好",
        "可用",
        "五个靶点结合能热图清晰，但Vina分数不能替代实验亲和力。",
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
        "gmx_MMPBSA结合自由能分解未运行。",
    ),
}


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

    rows: list[dict[str, Any]] = []
    for figure, aliases in aliases_by_figure.items():
        for alias in aliases:
            source = (output_root / alias.source).resolve()
            review = REVIEWS.get((figure, alias.panel), FigureReview("不确定", "不确定", "人工复核", ""))
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
    output_tiers = (
        available["recommended_output"]
        .value_counts()
        .to_dict()
    )
    major_issues = audit[
        audit["audit_verdict"].isin(
            {"需重绘", "不可替代Venn", "结果未达标", "结果未全达标", "不满足方案"}
        )
    ]
    return {
        "standard": "SCI IF>10 provisional figure screen",
        "target_journal": "not specified; exact journal requirements remain unverified",
        "overall_verdict": (
            "not_ready_for_if10"
            if len(not_run) or len(major_issues) or len(resolution_issues)
            else "conditional_ready"
        ),
        "panel_entries": int(len(audit)),
        "available_panels": int(len(available)),
        "prepared_not_run_panels": int(len(not_run)),
        "missing_panels": int((audit["status"] == "missing").sum()),
        "mean_automatic_quality_score": float(scores.mean()) if len(scores) else None,
        "median_dpi_x": float(available["dpi_x"].median()) if len(available) else None,
        "panels_wider_than_7_5_inches": int(len(oversize)),
        "resolution_or_vector_issue_count": int(len(resolution_issues)),
        "recommended_output_counts": {
            str(key): int(value)
            for key, value in output_tiers.items()
        },
        "major_issue_count": int(len(major_issues)),
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
        "- Re-run Calibration after a calibration method is selected; the current "
        "Hosmer-Lemeshow test fails.",
        "- Rebuild Figure 1d/e only if licensed GeneCards, OMIM, TTD exports and "
        "independent ChEMBL/STITCH evidence become available.",
        "- Keep the CellChat-like label on the ligand-receptor network, or replace it "
        "with a full CellChat permutation analysis and validated interaction diagram.",
        "- Re-render the 160 dpi virtual-knockout wrapper panels from source data "
        "before considering them for a main figure.",
        "- Assemble final multi-panel figures, add lowercase panel letters, export "
        "vector PDF/SVG for line art, and provide source data tables.",
        "- Add confidence intervals, effect sizes and exact FDR to key statistical panels.",
        "",
        "## Contact Sheets",
        "",
    ]
    for figure, path in summary["contact_sheets"].items():
        lines.append(f"- `{figure}`: `{path}`")
    return "\n".join(lines) + "\n"
