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
if (!"sample" %in% colnames(obj@meta.data)) {
  stop("Seurat object has no 'sample' metadata column")
}

bulk <- AggregateExpression(
  obj,
  group.by = "sample",
  assays = "RNA",
  return.seurat = FALSE
)$RNA
bulk <- as.matrix(bulk)
if (nrow(bulk) == 0 || ncol(bulk) == 0) {
  stop("pseudobulk aggregation returned an empty matrix")
}

# Seurat's AggregateExpression turns underscores into dashes; restore the
# original sample identifiers before matching the metadata table.
colnames(bulk) <- gsub("-", "_", colnames(bulk))
n_cells <- table(obj$sample)
sample_names <- colnames(bulk)
counts_per_cell <- as.numeric(n_cells[sample_names])
bulk <- sweep(bulk, 2, counts_per_cell, "/")

expr <- as.data.frame(bulk)
expr <- cbind(gene = rownames(expr), expr)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  expr,
  file.path(out_dir, "pseudobulk_expression.csv"),
  row.names = FALSE
)

meta <- unique(data.frame(
  sample = as.character(obj$sample),
  condition = as.character(obj$condition),
  stringsAsFactors = FALSE
))
if ("celltype_annot" %in% colnames(obj@meta.data)) {
  type_tab <- table(obj$sample, obj$celltype_annot)
  dominant <- colnames(type_tab)[max.col(type_tab, ties.method = "first")]
  meta$cell_type <- dominant[match(meta$sample, rownames(type_tab))]
}
meta$n_cells <- as.numeric(n_cells[meta$sample])
write.csv(
  meta,
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
