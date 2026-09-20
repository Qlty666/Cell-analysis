# Biological samples are the replicate unit; cell types are fitted separately.
liver_celltype_pseudobulk <- function(counts, metadata, covariates = character()) {
  required <- c("sample", "condition", "celltype_annot", covariates)
  if (!all(required %in% colnames(metadata))) stop("Missing pseudobulk metadata columns")
  if (!identical(colnames(counts), rownames(metadata))) stop("Counts and metadata are not aligned")
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("DESeq2 is required")
  records <- list()
  statuses <- list()
  for (cell_type in sort(unique(as.character(metadata$celltype_annot)))) {
    result <- tryCatch({
      keep <- which(metadata$celltype_annot == cell_type)
      meta <- metadata[keep, required, drop = FALSE]
      if (anyNA(meta)) stop("Missing sample/condition/covariate values")
      smeta <- unique(meta[, c("sample", "condition", covariates), drop = FALSE])
      if (anyDuplicated(smeta$sample)) stop("Sample has conflicting metadata")
      if (length(unique(smeta$condition)) != 2 || any(table(smeta$condition) < 2)) {
        stop("At least two biological samples per condition are required")
      }
      rownames(smeta) <- as.character(smeta$sample)
      assignment <- Matrix::sparseMatrix(i = seq_along(keep), j = match(meta$sample, smeta$sample), x = 1,
                                         dims = c(length(keep), nrow(smeta)))
      bulk <- counts[, keep, drop = FALSE] %*% assignment
      colnames(bulk) <- smeta$sample
      values <- if (inherits(bulk, "sparseMatrix")) bulk@x else as.numeric(bulk)
      if (any(!is.finite(values)) || any(values < 0) || any(abs(values - round(values)) > 1e-6)) {
        stop("Cell-type pseudobulk requires raw counts")
      }
      for (name in c(covariates, "condition")) smeta[[name]] <- factor(smeta[[name]])
      design <- reformulate(c(covariates, "condition"))
      model <- model.matrix(design, smeta)
      if (qr(model)$rank < ncol(model) || nrow(model) <= ncol(model)) stop("Confounded or saturated design")
      dds <- DESeq2::DESeqDataSetFromMatrix(round(as.matrix(bulk)), smeta, design)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      levels <- levels(smeta$condition)
      res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", levels[1], levels[2])))
      res$gene <- rownames(res)
      res$cell_type <- cell_type
      res$inference_unit <- "biological_sample"
      res$design <- paste(deparse(design), collapse = "")
      res
    }, error = function(e) e)
    if (inherits(result, "error")) {
      statuses[[cell_type]] <- list(status = "not_estimable", reason = conditionMessage(result))
    } else {
      records[[cell_type]] <- result
      statuses[[cell_type]] <- list(status = "completed", genes = nrow(result))
    }
  }
  combined <- if (length(records)) do.call(rbind, records) else data.frame()
  if (nrow(combined)) combined$padj_across_celltypes <- p.adjust(combined$pvalue, method = "BH")
  list(results = combined, status = statuses)
}
