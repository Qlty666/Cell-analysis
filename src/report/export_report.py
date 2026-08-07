#!/usr/bin/env python3
"""Export the pipeline report as DOCX and optionally PDF."""

import json
import sys
import textwrap
import zipfile
from pathlib import Path


def escape_xml(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def write_docx(path: Path, title: str, lines: list[str], rows: list[list[str]]) -> None:
    paras = "".join(
        f"<w:p><w:r><w:t>{escape_xml(line)}</w:t></w:r></w:p>"
        for line in lines
    )
    table_rows = ""
    for row in rows:
        cells = "".join(
            f"<w:tc><w:p><w:r><w:t>{escape_xml(cell)}</w:t></w:r></w:p></w:tc>"
            for cell in row
        )
        table_rows += f"<w:tr>{cells}</w:tr>"
    table = (
        f"<w:tbl><w:tblPr><w:tblBorders>"
        f"<w:top w:val='single'/><w:left w:val='single'/>"
        f"<w:bottom w:val='single'/><w:right w:val='single'/>"
        f"</w:tblBorders></w:tblPr>{table_rows}</w:tbl>"
    )
    document = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        f"<w:body>{paras}{table}</w:body></w:document>"
    )
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
            "</Types>",
        )
        zf.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
            "</Relationships>",
        )
        zf.writestr(
            "word/_rels/document.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>',
        )
        zf.writestr("word/document.xml", document)


def write_pdf(path: Path, title: str, lines: list[str]) -> None:
    try:
        from fpdf import FPDF
    except Exception:
        return
    pdf = FPDF()
    pdf.add_page()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.set_font("Helvetica", size=14)
    pdf.cell(0, 10, title)
    pdf.ln(12)
    pdf.set_font("Helvetica", size=10)
    width = pdf.w - pdf.l_margin - pdf.r_margin
    for line in lines:
        for wrapped in textwrap.wrap(line, width=90) or [line]:
            pdf.multi_cell(width, 6, wrapped)
    pdf.output(str(path))


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    summary_path = root / "results" / "summary.json"
    if not summary_path.exists():
        return 1
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    lines = [
        f"Dataset: {summary.get('dataset', '')}",
        f"Title: {summary.get('title', '')}",
        f"Cells raw: {summary.get('n_cells_raw', '')}",
        f"Cells after QC: {summary.get('n_cells_after_qc', '')}",
        f"Cells after doublet removal: {summary.get('n_cells_after_doublet_removal', '')}",
        f"Genes: {summary.get('n_genes', '')}",
        f"Clusters: {summary.get('n_clusters', '')}",
        f"Up DEGs: {summary.get('deg_up', '')}",
        f"Down DEGs: {summary.get('deg_down', '')}",
    ]
    top = summary.get("top_degs", [])
    if isinstance(top, dict):
        top = [top]
    rows = [["Gene", "Log2FC", "Padj"]]
    for item in top[:20]:
        rows.append([str(item.get("gene", "")), str(item.get("avg_log2FC", "")), str(item.get("p_val_adj", ""))])

    write_docx(root / "results" / "result_report.docx", "Single-Cell Analysis Report", lines, rows)
    write_pdf(root / "results" / "result_report.pdf", "Single-Cell Analysis Report", lines)
    (root / "results" / "export_status.txt").write_text(
        "DOCX/PDF export completed",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
