#!/usr/bin/env Rscript
# Export a compact Seurat expression matrix + cell metadata for the
# single-cell in-silico knockout module.
# Usage: Rscript export_single_cell_insilico.R <rds> <out_dir> <max_cells> <max_genes> [seed]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript export_single_cell_insilico.R <rds> <out_dir> <max_cells> <max_genes> [seed]")
}
rds_path <- args[1]
out_dir <- args[2]
max_cells <- as.integer(args[3])
max_genes <- as.integer(args[4])
seed <- if (length(args) >= 5) as.integer(args[5]) else 123L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressWarnings(suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
}))

obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) {
  stop(rds_path, " is not a Seurat object")
}

assay <- DefaultAssay(obj)
counts <- tryCatch(
  GetAssayData(obj, assay = assay, layer = "counts"),
  error = function(e) NULL
)
if (is.null(counts) || ncol(counts) == 0 || nrow(counts) == 0 ||
    Matrix::nnzero(counts) == 0) {
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
}
if (is.null(counts)) {
  stop("no counts or data layer found for assay ", assay)
}
counts <- as(counts, "CsparseMatrix")

meta <- obj@meta.data
celltype_col <- NULL
for (candidate in c(
  "celltype_annot", "cell_type", "celltype", "cell_annotation",
  "singleR_labels", "singleR_annotation", "predicted.celltype.l2",
  "annotation", "cluster_label", "seurat_clusters"
)) {
  if (candidate %in% colnames(meta)) {
    celltype_col <- candidate
    break
  }
}
sample_col <- NULL
for (candidate in c("sample", "orig.ident", "condition", "group")) {
  if (candidate %in% colnames(meta)) {
    sample_col <- candidate
    break
  }
}

celltype <- if (!is.null(celltype_col)) {
  as.character(meta[[celltype_col]])
} else {
  as.character(Idents(obj))
}
celltype[is.na(celltype) | !nzchar(celltype)] <- "Unannotated"
sample_id <- if (!is.null(sample_col)) {
  as.character(meta[[sample_col]])
} else {
  rep("Sample", ncol(obj))
}

set.seed(seed)
cells_all <- colnames(obj)
if (length(cells_all) > max_cells) {
  by_type <- split(seq_along(cells_all), factor(celltype))
  picked <- c()
  for (idx in by_type) {
    take <- max(1L, round(length(idx) * max_cells / length(cells_all)))
    picked <- c(picked, sample(idx, min(length(idx), take)))
  }
  picked <- sort(unique(picked))
  if (length(picked) > max_cells) {
    picked <- sort(sample(picked, max_cells))
  }
  cells_all <- cells_all[picked]
  counts <- counts[, cells_all, drop = FALSE]
  celltype <- celltype[picked]
  sample_id <- sample_id[picked]
}

present <- Matrix::rowSums(counts > 0)
means <- Matrix::rowMeans(counts)
row_var <- function(x) {
  x2 <- x * x
  Matrix::rowMeans(x2) - Matrix::rowMeans(x)^2
}
vars <- row_var(counts)
keep <- which(present >= max(2, ncol(counts) * 0.01) & vars > 0)
keep <- keep[order(-means[keep] * vars[keep])]
if (length(keep) > max_genes) keep <- keep[seq_len(max_genes)]
mat <- counts[keep, , drop = FALSE]
if (length(keep) > 0) {
  mat <- as.matrix(mat)
  mat <- cbind(gene = rownames(mat), as.data.frame(mat))
  con <- gzfile(file.path(out_dir, "expression.csv.gz"), "w")
  write.csv(mat, con, row.names = FALSE)
  close(con)
}

umap <- NULL
if ("umap" %in% Reductions(obj)) {
  umap <- as.data.frame(Embeddings(obj, reduction = "umap")[cells_all, , drop = FALSE])
}
md <- data.frame(
  cell = cells_all,
  cell_type = celltype,
  sample = sample_id,
  stringsAsFactors = FALSE
)
if (!is.null(umap)) {
  md <- cbind(md, umap[, c(1, 2), drop = FALSE])
  colnames(md)[(ncol(md) - 1):ncol(md)] <- c("umap_1", "umap_2")
}
write.csv(md, file.path(out_dir, "metadata.csv"), row.names = FALSE)

message("exported ", length(cells_all), " cells, ", length(keep), " genes")
quit(status = 0)
