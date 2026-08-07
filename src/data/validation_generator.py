"""Generate small synthetic single-cell datasets for pipeline validation."""

import csv
import gzip
import json
import random
from pathlib import Path

CANONICAL_MARKERS = [
    "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "CD4",
    "CD79A", "MS4A1", "CD19", "IGHG1",
    "LYZ", "CD68", "C1QA", "C1QB", "FCGR3A",
    "ALB", "APOA1", "APOA2", "FGB", "SERPINA1",
    "KRT19", "KRT7", "EPCAM", "SOX9", "CFTR",
    "COL1A1", "COL1A2", "ACTA2", "RGS5", "PDGFRB",
    "PECAM1", "VWF", "CLDN5", "PLVAP",
    "KRT8", "KRT18", "AFP", "GPC3",
    "TPSAB1", "CPA3", "MS4A2",
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB",
]

CELLTYPES = [
    "T_NK", "B", "Myeloid", "Hepatocyte",
    "Cholangiocyte", "Hepatic_stellate", "Endothelial", "Fibroblast",
]

CELLTYPE_MARKERS = {
    "T_NK": ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "CD4"],
    "B": ["CD79A", "MS4A1", "CD19", "IGHG1"],
    "Myeloid": ["LYZ", "CD68", "C1QA", "C1QB", "FCGR3A"],
    "Hepatocyte": ["ALB", "APOA1", "APOA2", "FGB", "SERPINA1"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"],
    "Hepatic_stellate": ["COL1A1", "COL1A2", "ACTA2", "RGS5", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "PLVAP"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
}


def build_genes(n_genes: int = 3000) -> list[str]:
    genes = list(CANONICAL_MARKERS)
    for i in range(n_genes):
        name = f"Gene{i:05d}"
        if name not in genes:
            genes.append(name)
    return genes[:n_genes]


def build_cells(
    n_cells: int,
    n_samples: int,
    seed: int,
) -> list[dict]:
    rng = random.Random(seed)
    samples = []
    for i in range(n_samples):
        samples.append({
            "name": f"S{i + 1:02d}",
            "condition": "GroupA" if i < n_samples // 2 else "GroupB",
        })
    cells = []
    for i in range(n_cells):
        sample = samples[i % len(samples)]
        celltype = CELLTYPES[i % len(CELLTYPES)]
        cells.append({
            "barcode": f"ACGT{i:06d}",
            "sample": sample["name"],
            "condition": sample["condition"],
            "celltype": celltype,
        })
    return cells


def simulate_counts(
    cells: list[dict],
    genes: list[str],
    seed: int,
) -> tuple[list[list[int]], dict]:
    rng = random.Random(seed)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    marker_indices = {
        ct: [gene_index[g] for g in CELLTYPE_MARKERS[ct] if g in gene_index]
        for ct in CELLTYPES
    }
    mt_indices = [
        gene_index[g]
        for g in ["MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP6", "MT-CYB"]
        if g in gene_index
    ]
    diff_a = [gene_index[g] for g in genes if g.startswith("Gene0")][:80]
    diff_b = [gene_index[g] for g in genes if g.startswith("Gene1")][:80]
    matrix = [[0] * len(cells) for _ in genes]
    for j, cell in enumerate(cells):
        for idx in marker_indices[cell["celltype"]]:
            matrix[idx][j] += 6 + rng.randint(0, 4)
        for idx in mt_indices:
            matrix[idx][j] += 1 + rng.randint(0, 2)
        for idx in diff_a if cell["condition"] == "GroupA" else diff_b:
            matrix[idx][j] += rng.randint(1, 3)
        for _ in range(250):
            idx = rng.randrange(len(genes))
            matrix[idx][j] += 1
    return matrix, marker_indices


def write_csv_gzip(path: Path, rows) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)


def write_metadata(path: Path, cells: list[dict]) -> None:
    rows = [["Cell.Barcode", "Sample", "Condition", "CellType"]]
    rows.extend(
        [c["barcode"], c["sample"], c["condition"], c["celltype"]]
        for c in cells
    )
    write_csv_gzip(path, rows)


def write_mtx(
    raw_dir: Path,
    prefix: str,
    matrix: list[list[int]],
    genes: list[str],
    cells: list[dict],
) -> None:
    barcode_path = raw_dir / f"{prefix}_barcodes.tsv.gz"
    gene_path = raw_dir / f"{prefix}_genes.tsv.gz"
    mtx_path = raw_dir / f"{prefix}_matrix.mtx.gz"

    with gzip.open(barcode_path, "wt", encoding="utf-8") as fh:
        for cell in cells:
            fh.write(cell["barcode"] + "\n")
    with gzip.open(gene_path, "wt", encoding="utf-8") as fh:
        for idx, gene in enumerate(genes):
            fh.write(f"ENSG{idx:011d}\t{gene}\n")

    triplets = []
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if value:
                triplets.append((i + 1, j + 1, value))
    with gzip.open(mtx_path, "wt", encoding="utf-8") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{len(genes)} {len(cells)} {len(triplets)}\n")
        for i, j, value in triplets:
            fh.write(f"{i} {j} {value}\n")


def write_barcode_csv(
    path: Path,
    matrix: list[list[int]],
    genes: list[str],
    cells: list[dict],
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([c["barcode"] for c in cells])
        for gene, row in zip(genes, matrix):
            writer.writerow([gene] + row)


def write_gene_csv(
    path: Path,
    matrix: list[list[int]],
    genes: list[str],
    cells: list[dict],
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["gene"] + [c["barcode"] for c in cells])
        for gene, row in zip(genes, matrix):
            writer.writerow([gene] + row)


def write_manifest(data_dir: Path, accession: str, files: dict) -> None:
    manifest = {
        "accession": accession,
        "mode": "generic",
        "files": files,
    }
    (data_dir / f"{accession}_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )


def generate_dataset(
    accession: str,
    output_root: Path,
    kind: str,
    n_cells: int = 1200,
    n_genes: int = 3000,
    n_samples: int = 4,
    seed: int = 1,
) -> Path:
    out = output_root / accession
    raw_dir = out / "data" / "raw" / accession
    data_dir = out / "data"
    raw_dir.mkdir(parents=True, exist_ok=True)

    genes = build_genes(n_genes)
    cells = build_cells(n_cells, n_samples, seed)
    matrix, _ = simulate_counts(cells, genes, seed)
    write_metadata(raw_dir / "metadata.csv.gz", cells)

    files = {
        "matrix": [],
        "barcodes": [],
        "genes": [],
        "metadata": ["metadata.csv.gz"],
        "series_matrices": [],
    }

    if kind == "mtx_single":
        write_mtx(raw_dir, "counts", matrix, genes, cells)
        files["matrix"] = ["counts_matrix.mtx.gz"]
        files["barcodes"] = ["counts_barcodes.tsv.gz"]
        files["genes"] = ["counts_genes.tsv.gz"]
    elif kind == "barcode_csv_single":
        write_barcode_csv(raw_dir / "counts.csv.gz", matrix, genes, cells)
        files["matrix"] = ["counts.csv.gz"]
    elif kind == "gene_csv_single":
        write_gene_csv(raw_dir / "counts.csv.gz", matrix, genes, cells)
        files["matrix"] = ["counts.csv.gz"]
    elif kind == "mtx_multi":
        per = max(1, n_cells // n_samples)
        for i in range(n_samples):
            sub_cells = cells[i * per:(i + 1) * per]
            sub_matrix = [row[i * per:(i + 1) * per] for row in matrix]
            write_mtx(raw_dir, f"sample{i + 1}", sub_matrix, genes, sub_cells)
            files["matrix"].append(f"sample{i + 1}_matrix.mtx.gz")
            files["barcodes"].append(f"sample{i + 1}_barcodes.tsv.gz")
            files["genes"].append(f"sample{i + 1}_genes.tsv.gz")
    elif kind == "barcode_csv_multi":
        per = max(1, n_cells // n_samples)
        for i in range(n_samples):
            sub_cells = cells[i * per:(i + 1) * per]
            sub_matrix = [row[i * per:(i + 1) * per] for row in matrix]
            write_barcode_csv(
                raw_dir / f"sample{i + 1}.csv.gz",
                sub_matrix,
                genes,
                sub_cells,
            )
            files["matrix"].append(f"sample{i + 1}.csv.gz")
    elif kind == "gene_csv_multi":
        per = max(1, n_cells // n_samples)
        for i in range(n_samples):
            sub_cells = cells[i * per:(i + 1) * per]
            sub_matrix = [row[i * per:(i + 1) * per] for row in matrix]
            write_gene_csv(
                raw_dir / f"sample{i + 1}.csv.gz",
                sub_matrix,
                genes,
                sub_cells,
            )
            files["matrix"].append(f"sample{i + 1}.csv.gz")
    elif kind == "mtx_no_sample":
        write_mtx(raw_dir, "counts", matrix, genes, cells)
        files["matrix"] = ["counts_matrix.mtx.gz"]
        files["barcodes"] = ["counts_barcodes.tsv.gz"]
        files["genes"] = ["counts_genes.tsv.gz"]
        # remove Sample from metadata to exercise sample inference fallback
        rows = [["Cell.Barcode", "Condition", "CellType"]]
        rows.extend(
            [c["barcode"], c["condition"], c["celltype"]]
            for c in cells
        )
        write_csv_gzip(raw_dir / "metadata.csv.gz", rows)
    elif kind == "mixed_gene_mtx":
        half = n_cells // 2
        write_gene_csv(
            raw_dir / "part1.csv.gz",
            [row[:half] for row in matrix],
            genes,
            cells[:half],
        )
        write_mtx(
            raw_dir,
            "part2",
            [row[half:] for row in matrix],
            genes,
            cells[half:],
        )
        files["matrix"] = ["part1.csv.gz", "part2_matrix.mtx.gz"]
        files["barcodes"] = ["", "part2_barcodes.tsv.gz"]
        files["genes"] = ["", "part2_genes.tsv.gz"]

    write_manifest(data_dir, accession, files)
    return out
