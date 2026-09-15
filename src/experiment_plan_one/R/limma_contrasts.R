#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
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
metadata <- metadata[!duplicated(sample_ids), , drop = FALSE]

labels <- as.character(metadata[[condition_column]])
names(labels) <- rownames(metadata)
wanted <- c(group1, group2)
selected <- names(labels)[labels %in% wanted]
selected <- intersect(selected, colnames(expression))
if (length(selected) < 4) {
  stop("fewer than four samples remain after metadata/expression matching")
}

expression <- expression[, selected, drop = FALSE]
group <- factor(labels[selected], levels = wanted)
if (any(table(group) < 2)) {
  stop("each comparison group must contain at least two samples")
}
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast_text <- paste0(make.names(tail(wanted, 1)), "-", make.names(head(wanted, 1)))
contrast_matrix <- makeContrasts(contrasts = contrast_text, levels = design)
fit <- lmFit(expression, design)
fit <- eBayes(contrasts.fit(fit, contrast_matrix))
result <- topTable(fit, number = Inf, sort.by = "P", adjust.method = "BH")
result$gene <- rownames(result)
group1_mean <- rowMeans(expression[, group == head(wanted, 1), drop = FALSE])
group2_mean <- rowMeans(expression[, group == tail(wanted, 1), drop = FALSE])
result$group1_mean <- group1_mean
result$group2_mean <- group2_mean
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
