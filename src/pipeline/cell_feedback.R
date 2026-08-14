#!/usr/bin/env Rscript
# Re-analyze a Seurat object using virtual knockout and docking results.
# Usage: Rscript cell_feedback.R <single_cell_root> <out_dir> <top_n> <max_features>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: Rscript cell_feedback.R <single_cell_root> <out_dir> <top_n> <max_features>")
}
root <- args[1]
out_dir <- args[2]
top_n <- as.integer(args[3])
max_features <- as.integer(args[4])

suppressWarnings(suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
}))

data_dir <- file.path(out_dir, "data")
fig_dir <- file.path(out_dir, "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(file, plot, width, height, dpi = 150, plot_margin = NULL) {
  if (is.null(plot_margin)) {
    plot_margin <- ggplot2::margin(26, 26, 24, 18, "pt")
  }
  if (inherits(plot, "patchwork")) {
    plot <- plot & theme(plot.margin = plot_margin)
  } else {
    plot <- plot + theme(plot.margin = plot_margin)
  }
  ggsave(file, plot, width = width, height = height, dpi = dpi, bg = "white")
}

write_summary <- function(x) {
  writeLines(
    jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE),
    file.path(out_dir, "cell_feedback_summary.json")
  )
}

manifest_path <- file.path(out_dir, "manifest.csv")
if (!file.exists(manifest_path)) {
  write_summary(list(status = "skipped", reason = "manifest.csv not found"))
  quit(status = 0)
}
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (nrow(manifest) == 0) {
  write_summary(list(status = "skipped", reason = "feedback manifest is empty"))
  quit(status = 0)
}

find_seurat_rds <- function(root_path) {
  ckpt <- file.path(root_path, "results", "checkpoints")
  preferred <- c(
    "seurat_annotated.rds",
    "seurat_clustered.rds",
    "seurat_singlet.rds",
    "seurat_raw.rds"
  )
  for (name in preferred) {
    path <- file.path(ckpt, name)
    if (file.exists(path)) return(path)
  }
  if (dir.exists(ckpt)) {
    found <- list.files(ckpt, pattern = "\\.rds$", full.names = TRUE)
    if (length(found) > 0) return(found[1])
  }
  res_data <- file.path(root_path, "results", "data")
  if (dir.exists(res_data)) {
    found <- list.files(res_data, pattern = "\\.rds$", full.names = TRUE)
    if (length(found) > 0) return(found[1])
  }
  NULL
}

rds_path <- find_seurat_rds(root)
if (is.null(rds_path)) {
  write_summary(list(
    status = "skipped",
    reason = "no Seurat RDS object found under results/checkpoints or results/data"
  ))
  quit(status = 0)
}

obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) {
  write_summary(list(
    status = "skipped",
    reason = paste(rds_path, "is not a Seurat object")
  ))
  quit(status = 0)
}

available <- intersect(manifest$gene, rownames(obj))
if (length(available) == 0) {
  write_summary(list(
    status = "skipped",
    reason = "no feedback genes were found in the Seurat object"
  ))
  quit(status = 0)
}

manifest <- manifest[match(available, manifest$gene), , drop = FALSE]
manifest$feedback_score <- as.numeric(manifest$feedback_score)
manifest$feedback_score[is.na(manifest$feedback_score)] <- 0.5

assay_name <- DefaultAssay(obj)
data_mat <- GetAssayData(obj, assay = assay_name, layer = "data")
expr_mat <- as.matrix(data_mat[available, , drop = FALSE])
expr_df <- as.data.frame(t(expr_mat))
colnames(expr_df) <- paste0("expr_", make.names(available))

if (length(available) >= 3) {
  ctrl_n <- min(50, max(2, floor(ncol(obj) * 0.1)))
  module_added <- tryCatch({
    obj <- AddModuleScore(
      obj,
      features = list(available),
      name = "Feedback_",
      ctrl = ctrl_n
    )
    TRUE
  }, error = function(e) FALSE)
  if (module_added && "Feedback_1" %in% colnames(obj@meta.data)) {
    module_col <- "Feedback_1"
  } else {
    scaled_expr <- t(scale(t(expr_mat)))
    scaled_expr[!is.finite(scaled_expr)] <- 0
    obj$Feedback_1 <- as.numeric(Matrix::colMeans(scaled_expr))
    module_col <- "Feedback_1"
  }
} else {
  module_score <- as.numeric(Matrix::colMeans(expr_mat))
  obj$Feedback_1 <- module_score
  module_col <- "Feedback_1"
}

celltype_col <- if ("celltype_annot" %in% colnames(obj@meta.data)) {
  "celltype_annot"
} else if ("cluster_label" %in% colnames(obj@meta.data)) {
  "cluster_label"
} else {
  "seurat_clusters"
}
celltype_vec <- as.character(obj@meta.data[[celltype_col]])
condition_vec <- as.character(obj$condition)
sample_vec <- as.character(obj$sample)
if (any(is.na(celltype_vec))) celltype_vec[is.na(celltype_vec)] <- "Unannotated"
if (any(is.na(condition_vec))) condition_vec[is.na(condition_vec)] <- "Unknown"

score_df <- data.frame(
  cell = colnames(obj),
  condition = condition_vec,
  sample = sample_vec,
  celltype = celltype_vec,
  module_score = as.numeric(obj@meta.data[[module_col]]),
  stringsAsFactors = FALSE
)
umap <- tryCatch(
  Embeddings(obj, reduction = "umap"),
  error = function(e) NULL
)
if (!is.null(umap) && ncol(umap) >= 2) {
  score_df$UMAP_1 <- umap[, 1]
  score_df$UMAP_2 <- umap[, 2]
}
score_df <- cbind(score_df, expr_df)
write.csv(
  score_df,
  file.path(data_dir, "cell_scores.csv"),
  row.names = FALSE
)

celltypes <- unique(celltype_vec)
ct_summary <- do.call(rbind, lapply(celltypes, function(ct) {
  idx <- which(celltype_vec == ct)
  means <- colMeans(expr_df[idx, , drop = FALSE])
  row <- data.frame(
    celltype = ct,
    n_cells = length(idx),
    module_mean = mean(score_df$module_score[idx]),
    module_median = median(score_df$module_score[idx]),
    stringsAsFactors = FALSE
  )
  for (col in colnames(expr_df)) row[[paste0("mean_", col)]] <- means[[col]]
  row
}))
ct_summary <- ct_summary[order(-ct_summary$module_mean), , drop = FALSE]
write.csv(
  ct_summary,
  file.path(data_dir, "celltype_summary.csv"),
  row.names = FALSE
)

enrichment <- data.frame()
if (length(celltypes) >= 2) {
  enrichment <- do.call(rbind, lapply(celltypes, function(ct) {
    idx <- celltype_vec == ct
    p_value <- tryCatch(
      suppressWarnings(
        wilcox.test(
          score_df$module_score[idx],
          score_df$module_score[!idx]
        )$p.value
      ),
      error = function(e) NA_real_
    )
    data.frame(
      celltype = ct,
      n_cells = sum(idx),
      module_mean = mean(score_df$module_score[idx]),
      other_module_mean = mean(score_df$module_score[!idx]),
      module_diff = mean(score_df$module_score[idx]) -
        mean(score_df$module_score[!idx]),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }))
  enrichment$p_adjust <- stats::p.adjust(enrichment$p_value, method = "BH")
  enrichment <- enrichment[order(enrichment$p_adjust), , drop = FALSE]
  write.csv(
    enrichment,
    file.path(data_dir, "celltype_enrichment.csv"),
    row.names = FALSE
  )
}

condition_groups <- unique(data.frame(
  condition = condition_vec,
  celltype = celltype_vec,
  stringsAsFactors = FALSE
))
condition_summary <- do.call(rbind, lapply(seq_len(nrow(condition_groups)), function(i) {
  cond <- condition_groups$condition[i]
  ct <- condition_groups$celltype[i]
  idx <- condition_vec == cond & celltype_vec == ct
  means <- colMeans(expr_df[idx, , drop = FALSE])
  row <- data.frame(
    condition = cond,
    celltype = ct,
    n_cells = sum(idx),
    module_mean = mean(score_df$module_score[idx]),
    stringsAsFactors = FALSE
  )
  for (col in colnames(expr_df)) row[[paste0("mean_", col)]] <- means[[col]]
  row
}))
write.csv(
  condition_summary,
  file.path(data_dir, "condition_summary.csv"),
  row.names = FALSE
)

detection_rate <- vapply(available, function(g) {
  mean(expr_mat[g, ] > 0)
}, numeric(1))
ct_means_by_gene <- do.call(rbind, lapply(celltypes, function(ct) {
  idx <- celltype_vec == ct
  rowMeans(expr_mat[, idx, drop = FALSE])
}))
rownames(ct_means_by_gene) <- celltypes
colnames(ct_means_by_gene) <- available
if (length(celltypes) >= 2) {
  prop <- sweep(
    ct_means_by_gene,
    2,
    colSums(ct_means_by_gene),
    "/"
  )
  prop[is.na(prop)] <- 0
  prop[prop <= 0] <- 1e-12
  entropy <- -colSums(prop * log2(prop))
  specificity <- 1 - entropy / log2(length(celltypes))
  specificity <- pmax(0, pmin(1, specificity))
} else {
  specificity <- rep(0.5, length(available))
  names(specificity) <- available
}
top_celltype <- vapply(available, function(g) {
  celltypes[which.max(ct_means_by_gene[, g])]
}, character(1))

targets <- data.frame(
  gene = available,
  source = manifest$source,
  feedback_score = manifest$feedback_score,
  target_score = manifest$target_score,
  knockout_score = manifest$knockout_score,
  docking_status = manifest$docking_status,
  docking_hits = manifest$docking_hits,
  best_affinity = manifest$best_affinity,
  cell_detection_rate = detection_rate,
  celltype_specificity = unname(specificity),
  top_celltype = top_celltype,
  stringsAsFactors = FALSE
)
targets$cell_support_score <- pmin(
  1,
  0.7 * targets$feedback_score + 0.3 * targets$celltype_specificity
)
targets <- targets[order(-targets$cell_support_score), , drop = FALSE]
write.csv(
  targets,
  file.path(data_dir, "feedback_targets.csv"),
  row.names = FALSE
)

figures <- character()
if (!is.null(umap) && module_col %in% colnames(obj@meta.data)) {
  p_module <- FeaturePlot(
    obj,
    features = module_col,
    cols = c("grey90", "#B31B1B")
  ) + ggtitle("Screening feedback module score")
  fig_module <- file.path(fig_dir, "fig_54_feedback_module_umap.png")
  save_fig(fig_module, p_module, width = 8, height = 7)
  figures <- c(figures, basename(fig_module))
}

feature_genes <- head(available, max_features)
if (length(feature_genes) > 0 && !is.null(umap)) {
  p_features <- FeaturePlot(
    obj,
    features = feature_genes,
    ncol = 3,
    cols = c("grey90", "#2E86AB")
  )
  fig_features <- file.path(fig_dir, "fig_55_feedback_target_expression_umap.png")
  save_fig(
    fig_features,
    p_features,
    width = 12,
    height = 4 * ceiling(length(feature_genes) / 3)
  )
  figures <- c(figures, basename(fig_features))
}

if (length(feature_genes) > 0) {
  p_dot <- DotPlot(obj, features = feature_genes, group.by = celltype_col) +
    RotatedAxis() +
    ggtitle("Screening target expression by cell type")
  fig_dot <- file.path(fig_dir, "fig_56_feedback_celltype_dotplot.png")
  save_fig(fig_dot, p_dot, width = 11, height = 7)
  figures <- c(figures, basename(fig_dot))
}

if (length(celltypes) >= 2) {
  p_box <- ggplot(
    score_df,
    aes(x = celltype, y = module_score, fill = condition)
  ) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Cell type",
      y = "Feedback module score",
      title = "Screening feedback module by cell type"
    )
  fig_box <- file.path(fig_dir, "fig_57_feedback_celltype_boxplot.png")
  save_fig(fig_box, p_box, width = 10, height = 7)
  figures <- c(figures, basename(fig_box))
}

if (length(available) > 1 && requireNamespace("pheatmap", quietly = TRUE)) {
  heat_mat <- t(ct_means_by_gene)
  heat_mat <- heat_mat[, order(colSums(heat_mat), decreasing = TRUE), drop = FALSE]
  fig_heat <- file.path(fig_dir, "fig_58_feedback_celltype_heatmap.png")
  pheatmap::pheatmap(
    heat_mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("#f7fbff", "#d6e7ff", "#2E86AB", "#7B241C"))(100),
    main = "Mean screening target expression by cell type",
    filename = fig_heat,
    width = 8,
    height = 7
  )
  figures <- c(figures, basename(fig_heat))
}

write_summary(list(
  status = "completed",
  seurat_rds = rds_path,
  genes_requested = nrow(manifest),
  genes_matched = length(available),
  module_score_column = module_col,
  celltype_column = celltype_col,
  n_cells = ncol(obj),
  n_celltypes = length(celltypes),
  top_celltypes = head(
    ct_summary$celltype[order(-ct_summary$module_mean)],
    10
  ),
  figures = figures,
  outputs = list(
    cell_scores_csv = file.path(data_dir, "cell_scores.csv"),
    celltype_summary_csv = file.path(data_dir, "celltype_summary.csv"),
    celltype_enrichment_csv = file.path(data_dir, "celltype_enrichment.csv"),
    condition_summary_csv = file.path(data_dir, "condition_summary.csv"),
    feedback_targets_csv = file.path(data_dir, "feedback_targets.csv")
  )
))
