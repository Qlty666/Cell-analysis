#!/usr/bin/env Rscript
# Aggregate a Seurat object to a sample-level pseudobulk expression matrix.
# Usage: Rscript export_pseudobulk.R <single_cell_root> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: Rscript export_pseudobulk.R <single_cell_root> <out_dir>")
}
root <- args[1]
out_dir <- args[2]
suppressPackageStartupMessages(library(Seurat))

data_dir <- file.path(root, "results", "data")
rds_files <- list.files(
  data_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)
if (length(rds_files) == 0) {
  stop("no Seurat RDS object found under ", data_dir)
}
obj <- readRDS(rds_files[1])
if (!inherits(obj, "Seurat")) {
  stop(rds_files[1], " is not a Seurat object")
}
meta <- obj@meta.data
if (!"sample" %in% colnames(meta) || !"condition" %in% colnames(meta)) {
  stop("Seurat object has no 'sample' metadata column")
}

sample_orig <- as.character(meta[["sample"]])
condition <- as.character(meta[["condition"]])
if (any(is.na(sample_orig)) || any(is.na(condition))) {
  stop("Seurat metadata has missing sample/condition values")
}
if (length(unique(condition)) < 2) {
  stop("Seurat metadata must contain at least two conditions for pseudobulk export")
}

group_label <- paste(sample_orig, make.names(condition), sep = "__")
group_factor <- factor(group_label, levels = unique(group_label))

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
counts <- as(counts, "CsparseMatrix")
bulk <- vapply(
  levels(group_factor),
  function(grp) {
    cells <- which(group_factor == grp)
    as.numeric(Matrix::rowSums(counts[, cells, drop = FALSE]))
  },
  numeric(nrow(counts))
)
if (is.null(dim(bulk))) {
  bulk <- matrix(bulk, ncol = 1)
}
rownames(bulk) <- rownames(counts)
colnames(bulk) <- levels(group_factor)
counts_per_cell <- as.numeric(table(group_factor))
bulk <- sweep(bulk, 2, counts_per_cell, "/")

expr <- as.data.frame(bulk)
expr <- cbind(gene = rownames(expr), expr)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  expr,
  file.path(out_dir, "pseudobulk_expression.csv"),
  row.names = FALSE
)

meta_out <- unique(data.frame(
  sample = levels(group_factor),
  condition = condition[match(levels(group_factor), group_label)],
  original_sample = sample_orig[match(levels(group_factor), group_label)],
  stringsAsFactors = FALSE
))
if ("celltype_annot" %in% colnames(meta)) {
  celltype <- as.character(meta[["celltype_annot"]])
  type_tab <- table(group_label, celltype)
  dominant <- colnames(type_tab)[max.col(type_tab, ties.method = "first")]
  meta_out$cell_type <- dominant[match(meta_out$sample, rownames(type_tab))]
}
meta_out$n_cells <- counts_per_cell[match(meta_out$sample, levels(group_factor))]
write.csv(
  meta_out,
  file.path(out_dir, "pseudobulk_metadata.csv"),
  row.names = FALSE
)

message(
  "pseudobulk exported: ",
  nrow(expr),
  " genes x ",
  ncol(expr) - 1,
  " samples -> ",
  out_dir
)
