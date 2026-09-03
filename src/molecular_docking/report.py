"""HTML report generator for the standalone molecular docking board."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

from docking.config import ResolvedConfig
from docking.utils import DockingError


def generate_report(cfg: ResolvedConfig, log) -> Path:
    """Write a docking-only HTML report under the molecular docking outputs."""
    reports = cfg.reports_dir()
    ranked = reports / "01_analysis" / "data" / "fig_46_47_ranked_results.csv"
    if not ranked.exists():
        raise DockingError(f"ranked results not found: {ranked}")
    top_n = int(cfg.get("report", "top_n", 20))
    frame = pd.read_csv(ranked, dtype={"id": str}).head(top_n)

    summary = _read_json(reports / "01_analysis" / "summary.json")
    redock_summary = _read_json(reports / "02_redock" / "summary.json")
    figures = sorted(
        p.relative_to(reports).as_posix()
        for p in reports.rglob("*.png")
    )
    pose_files = sorted(
        p.name
        for p in (cfg.output_dir / "docked").glob("*.pdbqt")
        if p.is_file()
    )

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
    poses_html = "".join(
        f"<li><code>{html.escape(name)}</code></li>" for name in pose_files
    )

    html_text = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>分子对接报告</title>
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
<h1>分子对接报告</h1>
<div class="card">
  <p class="muted">独立分子对接板块，与虚拟筛选流程分开运行。</p>
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
  <p>重对接完成：{_esc(str(redock_summary.get('ok', '未运行')))}</p>
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
<div class="card">
  <h2>对接构象</h2>
  <ul>{poses_html or '<li class="muted">暂无 PDBQT 构象。</li>'}</ul>
</div>
</body>
</html>
"""
    out_path = reports / "molecular_docking_report.html"
    out_path.write_text(html_text, encoding="utf-8")
    log.info("molecular docking report generated: %s", out_path)
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
