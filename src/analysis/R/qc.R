# qc.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

qc_percentage <- function(object, features, label = "QC") {
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  features <- intersect(features, rownames(counts))
  if (length(features) == 0 || ncol(counts) == 0) {
    # No matching features means the metric is undefined, not zero; returning 0
    # would let every cell pass the QC filter for this metric.
    log_msg(
      "WARNING: no features matched for ", label,
      "; setting ", label, " to NA for all cells"
    )
    return(rep(NA_real_, ncol(object)))
  }
  total_counts <- Matrix::colSums(counts)
  100 * Matrix::colSums(counts[features, , drop = FALSE]) /
    pmax(total_counts, 1)
}

qc_plot_features <- function(object, cols) {
  # All-NA QC metrics (no matching features) cannot be drawn by VlnPlot, which
  # errors on `all(x == x[1])`; drop them from the figure and say so.
  keep <- character(0)
  dropped <- character(0)
  for (col in cols) {
    values <- object[[]][[col]]
    if (length(values) == 0 || all(is.na(values))) {
      dropped <- c(dropped, col)
    } else {
      keep <- c(keep, col)
    }
  }
  if (length(dropped) > 0) {
    log_msg(
      "WARNING: skipping all-NA QC metric(s) in violin plots: ",
      paste(dropped, collapse = ", ")
    )
  }
  keep
}

ribo_features <- function(object) {
  grep("^(RP[SL]|Rp[ls])", rownames(object), value = TRUE)
}

hemoglobin_features <- function(object) {
  genes <- rownames(object)
  if (species == "mm") {
    patterns <- c(
      "^Hba-a[12]$", "^Hbb-[a-z0-9]+$", "^Hbd$", "^Hbe1?$",
      "^Hbg[12]$", "^Hbm$", "^Hbq1[ab]?$", "^Hbz$"
    )
  } else {
    patterns <- c(
      "^HBA1$", "^HBA2$", "^HBB$", "^HBD$", "^HBE1$",
      "^HBG1$", "^HBG2$", "^HBM$", "^HBQ1$", "^HBZ$"
    )
  }
  unique(unlist(lapply(patterns, grep, x = genes, value = TRUE)))
}

mt_features <- function(object) {
  genes <- rownames(object)
  hits <- grep(mt_pattern, genes, value = TRUE)
  if (length(hits) > 0) {
    return(hits)
  }
  # Ensembl-ID datasets (e.g. "ENSG00000198888") have no "MT-" symbols, so fall
  # back to the organism annotation package and keep the CHR == "MT" genes.
  ensembl <- genes[grepl("^(ENSG|ENSMUSG)", genes)]
  if (length(ensembl) == 0) {
    log_msg(
      "WARNING: no mitochondrial features found (no MT- symbols, no Ensembl IDs)"
    )
    return(character(0))
  }
  org_db <- if (species == "mm") {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    } else {
      NULL
    }
  } else if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
  } else {
    NULL
  }
  if (is.null(org_db)) {
    log_msg(
      "WARNING: no MT- features and organism annotation package unavailable; ",
      "percent.mt will be NA"
    )
    return(character(0))
  }
  mapped <- tryCatch(
    AnnotationDbi::select(
      org_db,
      keys = ensembl,
      columns = "CHR",
      keytype = "ENSEMBL"
    ),
    error = function(e) {
      log_msg(
        "WARNING: Ensembl MT annotation lookup failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  if (is.null(mapped) || !"CHR" %in% colnames(mapped)) {
    log_msg("WARNING: Ensembl MT annotation lookup returned no usable rows")
    return(character(0))
  }
  mt_ids <- unique(as.character(mapped$ENSEMBL[
    !is.na(mapped$CHR) & as.character(mapped$CHR) == "MT"
  ]))
  hits <- intersect(mt_ids, genes)
  if (length(hits) == 0) {
    log_msg(
      "WARNING: no mitochondrial features detected; percent.mt will be NA"
    )
  } else {
    log_msg("mitochondrial features resolved via Ensembl CHR=MT: ", length(hits))
  }
  hits
}

umi_feature_correlation_stats <- function(qc_frame, stage_label) {
  x <- log1p(as.numeric(qc_frame$nFeature_RNA))
  y <- log1p(as.numeric(qc_frame$nCount_RNA))
  n <- length(x)
  if (n < 3) {
    return(data.frame(
      stage = stage_label,
      n_cells = n,
      loglog_slope = NA_real_,
      loglog_intercept = NA_real_,
      loglog_correlation = NA_real_,
      residual_sd = NA_real_,
      residual_mad = NA_real_,
      review_threshold = NA_real_,
      n_review = 0L,
      review_pct = 0,
      stringsAsFactors = FALSE
    ))
  }
  fit <- lm(y ~ x)
  residuals_fit <- residuals(fit)
  residual_sd <- stats::sd(residuals_fit)
  residual_mad <- suppressWarnings(stats::mad(residuals_fit))
  if (!is.finite(residual_mad) || residual_mad <= .Machine$double.eps) {
    residual_mad <- residual_sd
  }
  review_threshold <- if (
    is.finite(residual_mad) && residual_mad > .Machine$double.eps
  ) {
    4 * residual_mad
  } else {
    max(abs(residuals_fit), na.rm = TRUE)
  }
  n_review <- sum(abs(residuals_fit) > review_threshold, na.rm = TRUE)
  data.frame(
    stage = stage_label,
    n_cells = n,
    loglog_slope = unname(coef(fit)[2]),
    loglog_intercept = unname(coef(fit)[1]),
    loglog_correlation = tryCatch(
      suppressWarnings(cor(x, y, use = "complete.obs")),
      error = function(e) NA_real_
    ),
    residual_sd = residual_sd,
    residual_mad = residual_mad,
    review_threshold = review_threshold,
    n_review = n_review,
    review_pct = 100 * n_review / max(1, n),
    stringsAsFactors = FALSE
  )
}

