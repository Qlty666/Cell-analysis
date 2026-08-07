#!/usr/bin/env python3
"""Generate an HTML summary report for a virtual screening run."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

from .config import ResolvedConfig
from .utils import DockingError


def generate_report(cfg: ResolvedConfig, log) -> Path:
    reports = cfg.reports_dir()
    ranked = reports / "ranked_results.csv"
    if not ranked.exists():
        raise DockingError(f"ranked results not found: {ranked}")
    top_n = int(cfg.get("report", "top_n", 20))
    frame = pd.read_csv(ranked, dtype={"id": str}).head(top_n)

    summary = _read_json(reports / "summary.json")
    ml_summary = _read_json(reports / "ml_predict_summary.json")
    figures = sorted(p.name for p in reports.glob("*.png"))
    rows_html = "".join(
        "<tr>"
        f"<td>{_esc(row.get('rank', ''))}</td>"
        f"<td>{_esc(row.get('id', ''))}</td>"
        f"<td>{_esc(row.get('affinity', ''))}</td>"
        f"<td>{_esc(row.get('smiles', ''))}</td>"
        "</tr>"
        for _, row in frame.iterrows()
    )
    figures_html = "".join(
        f'<figure><img src="{html.escape(name)}" alt="{html.escape(name)}">'
        f"<figcaption>{html.escape(name)}</figcaption></figure>"
        for name in figures
    )

    html_text = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>虚拟筛选报告</title>
<style>
body {{ font-family: "Segoe UI", Arial, sans-serif; margin: 24px; color: #1f2933; background: #f5f7fa; }}
h1 {{ font-size: 24px; }}
.card {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 18px; margin-bottom: 16px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th, td {{ border: 1px solid #e5e7eb; padding: 7px 8px; text-align: left; }}
th {{ background: #eef2f7; }}
.gallery {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }}
.gallery img {{ width: 100%; border: 1px solid #e5e7eb; border-radius: 6px; }}
.muted {{ color: #6b7280; }}
</style>
</head>
<body>
<h1>虚拟筛选报告</h1>
<div class="card">
  <p class="muted">项目：{_esc(cfg.data.get('name', 'virtual_screening'))}</p>
  <p>工作目录：{_esc(str(cfg.workdir))}</p>
  <p>受体：{_esc(str(cfg.receptor_output()))}</p>
  <p>配体库：{_esc(str(cfg.ligand_input()))}</p>
  <p>对接盒中心：{_esc(str(cfg.receptor_center()))}</p>
  <p>对接盒尺寸：{_esc(str(cfg.receptor_size()))}</p>
  <p>命中阈值：{_esc(str(cfg.get('analysis', 'cutoff', -7.0)))} kcal/mol</p>
</div>
<div class="card">
  <h2>统计</h2>
  <p>成功对接：{_esc(str(summary.get('total_docked', '')))}</p>
  <p>命中数：{_esc(str(summary.get('hits', '')))}</p>
  <p>最佳亲和力：{_esc(str(summary.get('best_affinity', '')))}</p>
  <p>ML 重打分：{_esc(str(ml_summary.get('scored', '未运行')))} 个配体</p>
</div>
<div class="card">
  <h2>Top {top_n} 结果</h2>
  <table>
    <thead><tr><th>排名</th><th>ID</th><th>亲和力</th><th>SMILES</th></tr></thead>
    <tbody>{rows_html}</tbody>
  </table>
</div>
<div class="card">
  <h2>结果图</h2>
  <div class="gallery">{figures_html or '<p class="muted">无结果图。</p>'}</div>
</div>
</body>
</html>
"""
    out_path = reports / "docking_report.html"
    out_path.write_text(html_text, encoding="utf-8")
    log.info("HTML report generated: %s", out_path)
    return out_path


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _esc(value) -> str:
    return html.escape(str(value if value is not None else ""))
