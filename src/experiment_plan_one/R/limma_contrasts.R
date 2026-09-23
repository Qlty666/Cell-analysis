#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(7, 9)) {
  stop(
    paste(
      "usage: limma_contrasts.R <expression.csv|gz> <metadata.csv>",
      "<output.csv> <condition_column> <group1_csv> <group2_csv> <comparison>"
    )
  )
}

suppressPackageStartupMessages({
  library(limma)
})

expression_path <- args[[1]]
metadata_path <- args[[2]]
output_path <- args[[3]]
condition_column <- args[[4]]
group1 <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
group2 <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
comparison <- args[[7]]
input_type <- if (length(args) >= 8) args[[8]] else "log_expression"
pair_column <- if (length(args) >= 9) args[[9]] else ""
if (!input_type %in% c("log_expression", "counts")) stop("invalid input_type")

read_table <- function(path, row_names = FALSE) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    read.csv(
      con,
      check.names = FALSE,
      row.names = if (row_names) 1 else NULL,
      stringsAsFactors = FALSE
    )
  } else {
    read.csv(
      path,
      check.names = FALSE,
      row.names = if (row_names) 1 else NULL,
      stringsAsFactors = FALSE
    )
  }
}

expression <- read_table(expression_path, row_names = TRUE)
metadata <- read_table(metadata_path, row_names = TRUE)
if (!condition_column %in% colnames(metadata)) {
  stop(sprintf("condition column not found: %s", condition_column))
}
sample_ids <- rownames(metadata)
if (is.null(sample_ids) || !length(sample_ids)) {
  stop("metadata sample IDs are missing")
}
if (anyDuplicated(sample_ids)) stop("duplicate metadata sample IDs")
if (length(intersect(group1, group2))) stop("comparison groups overlap")

labels <- as.character(metadata[[condition_column]])
names(labels) <- rownames(metadata)
wanted <- c(group1, group2)
selected <- names(labels)[labels %in% wanted]
selected <- intersect(selected, colnames(expression))
if (length(selected) < 4) {
  stop("fewer than four samples remain after metadata/expression matching")
}

expression <- expression[, selected, drop = FALSE]
group <- factor(ifelse(labels[selected] %in% group1, "reference", "case"),
                levels = c("reference", "case"))
if (any(table(group) < 2)) {
  stop("each comparison group must contain at least two samples")
}
if (nzchar(pair_column)) {
  if (!pair_column %in% names(metadata)) stop("pair column missing")
  ids <- as.character(metadata[selected, pair_column])
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) stop("missing patient IDs")
  patient <- factor(ids)
  if (any(table(patient, group) != 1L)) stop("expected one sample per patient per group")
  design <- model.matrix(~ 0 + group + patient)
} else {
  design <- model.matrix(~ 0 + group)
}
if (qr(design)$rank != ncol(design)) stop("design is not full rank")
if (nrow(design) <= ncol(design)) stop("no residual degrees of freedom")
contrast_matrix <- makeContrasts(groupcase - groupreference, levels = design)
expression <- as.matrix(expression)
if (any(!is.finite(expression))) stop("expression contains non-finite values")
if (input_type == "counts") {
  if (any(expression < 0) || any(colSums(expression) <= 0)) stop("invalid counts")
  suppressPackageStartupMessages(library(edgeR))
  counts <- DGEList(counts = expression)
  keep <- filterByExpr(counts, design = design)
  if (sum(keep) < 2) stop("too few expressed genes")
  counts <- calcNormFactors(counts[keep, , keep.lib.sizes = FALSE])
  modeled <- voom(counts, design, plot = FALSE)
  expression <- modeled$E
  fit <- lmFit(modeled, design)
} else {
  fit <- lmFit(expression, design)
}
fit <- eBayes(contrasts.fit(fit, contrast_matrix))
result <- topTable(fit, number = Inf, sort.by = "P", adjust.method = "BH")
result$gene <- rownames(result)
group1_mean <- rowMeans(expression[, group == "reference", drop = FALSE])
group2_mean <- rowMeans(expression[, group == "case", drop = FALSE])
result$group1_mean <- group1_mean[rownames(result)]
result$group2_mean <- group2_mean[rownames(result)]
result$delta <- result$group2_mean - result$group1_mean
result$comparison <- comparison
result$group1 <- paste(group1, collapse = "+")
result$group2 <- paste(group2, collapse = "+")
result <- result[
  ,
  c(
    "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
    "group1_mean", "group2_mean", "delta", "comparison", "group1", "group2"
  ),
  drop = FALSE
]
result <- result[order(result$P.Value, result$logFC, decreasing = c(FALSE, TRUE)), , drop = FALSE]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_path, row.names = FALSE)
write.csv(data.frame(sample_id = selected, group = group, design, check.names = FALSE),
          paste0(output_path, ".design.csv"), row.names = FALSE)
