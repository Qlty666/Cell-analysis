#!/usr/bin/env Rscript
# Advanced bulk analyses: ComBat, limma, WGCNA, immune deconvolution and
# survival. Optional packages are checked independently so the Python
# ML/target-ranking layer can still run when an R package is unavailable.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: advanced_bulk.R <config.json> <output_dir>")
}
config_path <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(jsonlite)
})

config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
base_dir <- dirname(config_path)
log_lines <- character()
`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) fallback else value
}
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(...))
  log_lines <<- c(log_lines, line)
  message(line)
}
write_log <- function() {
  writeLines(log_lines, file.path(out_dir, "advanced_r.log"), useBytes = TRUE)
}
on.exit(write_log(), add = TRUE)

resolve_path <- function(value) {
  if (is.null(value) || !nzchar(as.character(value))) return(NULL)
  path <- as.character(value)
  if (!grepl("^([A-Za-z]:[/\\\\]|/)", path)) {
    path <- file.path(base_dir, path)
  }
  normalizePath(path, mustWork = FALSE)
}

read_table <- function(path) {
  if (!file.exists(path)) stop("table not found: ", path)
  sep <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

first_existing <- function(values, candidates) {
  lower <- tolower(trimws(values))
  for (candidate in candidates) {
    idx <- match(tolower(candidate), lower)
    if (!is.na(idx)) return(values[[idx]])
  }
  NULL
}

read_expression <- function(path, gene_column = NULL) {
  frame <- read_table(path)
  gene_col <- gene_column
  if (is.null(gene_col)) {
    gene_col <- first_existing(
      colnames(frame),
      c("gene", "symbol", "gene_symbol", "hgnc", "hugo", "feature", "id")
    )
  }
  if (!is.null(gene_col)) {
    if (!gene_col %in% colnames(frame)) {
      stop("gene column not found: ", gene_col)
    }
    rownames(frame) <- toupper(trimws(as.character(frame[[gene_col]])))
    frame[[gene_col]] <- NULL
  } else {
    rownames(frame) <- toupper(trimws(as.character(frame[[1]])))
    frame[[1]] <- NULL
  }
  frame <- frame[!is.na(rownames(frame)) & nzchar(rownames(frame)), , drop = FALSE]
  frame <- frame[!duplicated(rownames(frame)), , drop = FALSE]
  frame[] <- lapply(frame, function(x) suppressWarnings(as.numeric(x)))
  frame
}

read_metadata <- function(path, sample_column = NULL) {
  frame <- read_table(path)
  sample_col <- sample_column
  if (is.null(sample_col)) {
    sample_col <- first_existing(
      colnames(frame),
      c("sample", "sample_id", "sampleid", "donor", "patient", "id")
    )
  }
  if (is.null(sample_col)) {
    rownames(frame) <- as.character(frame[[1]])
    frame[[1]] <- NULL
  } else {
    rownames(frame) <- as.character(frame[[sample_col]])
    frame[[sample_col]] <- NULL
  }
  frame
}

infer_case_control <- function(labels) {
  values <- sort(unique(as.character(labels[!is.na(labels)])))
  if (length(values) != 2) {
    stop("classification requires exactly two labels: ", paste(values, collapse = ", "))
  }
  case_tokens <- c("case", "tumor", "tumour", "disease", "treated", "exposed", "hcc")
  control_tokens <- c("control", "normal", "healthy", "vehicle", "untreated")
  find_token <- function(tokens) {
    hit <- values[vapply(values, function(value) {
      any(vapply(tokens, function(token) grepl(token, value, ignore.case = TRUE), logical(1)))
    }, logical(1))]
    if (length(hit)) hit[[1]] else NA_character_
  }
  case_label <- find_token(case_tokens)
  control_label <- find_token(control_tokens)
  if (is.na(case_label)) case_label <- values[[1]]
  if (is.na(control_label)) control_label <- setdiff(values, case_label)[[1]]
  list(case = case_label, control = control_label)
}

load_cohort <- function(spec) {
  expression_path <- resolve_path(spec$expression)
  metadata_path <- resolve_path(spec$metadata)
  if (is.null(expression_path) || is.null(metadata_path)) {
    stop("cohort requires expression and metadata")
  }
  expression <- read_expression(expression_path, spec$gene_column)
  metadata <- read_metadata(metadata_path, spec$sample_column)
  condition_col <- spec$condition_column
  if (is.null(condition_col) || !condition_col %in% colnames(metadata)) {
    condition_col <- first_existing(
      colnames(metadata),
      c("condition", "group", "disease", "status", "tissue", "treatment")
    )
  }
  if (is.null(condition_col) || !condition_col %in% colnames(metadata)) {
    stop("cohort has no condition column")
  }
  common <- intersect(colnames(expression), rownames(metadata))
  if (length(common) < 4) {
    norm_expr <- gsub("\\.", "-", toupper(colnames(expression)))
    norm_meta <- gsub("\\.", "-", toupper(rownames(metadata)))
    idx <- match(norm_meta, norm_expr)
    common <- colnames(expression)[idx[!is.na(idx)]]
  }
  if (length(common) < 4) stop("cohort has fewer than four matched samples")
  expression <- expression[, common, drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  labels <- as.character(metadata[[condition_col]])
  case_label <- spec$case_label
  control_label <- spec$control_label
  if (is.null(case_label) || is.null(control_label)) {
    inferred <- infer_case_control(labels)
    if (is.null(case_label)) case_label <- inferred$case
    if (is.null(control_label)) control_label <- inferred$control
  }
  keep <- labels %in% c(case_label, control_label)
  if (sum(keep) < 4) stop("fewer than four samples in requested labels")
  list(
    name = if (is.null(spec$name)) "cohort" else as.character(spec$name),
    expression = expression[, keep, drop = FALSE],
    metadata = metadata[keep, , drop = FALSE],
    labels = factor(labels[keep], levels = c(control_label, case_label)),
    condition_column = condition_col,
    case_label = as.character(case_label),
    control_label = as.character(control_label)
  )
}

looks_like_counts <- function(frame) {
  values <- as.matrix(frame)
  values <- values[is.finite(values)]
  if (!length(values) || min(values) < 0) return(FALSE)
  all(abs(values - round(values)) < 1e-6) &&
    stats::median(colSums(as.matrix(frame), na.rm = TRUE)) >= 100
}

normalize_expression <- function(frame) {
  if (looks_like_counts(frame)) {
    library_size <- colSums(frame, na.rm = TRUE)
    library_size[library_size <= 0] <- NA_real_
    normalized <- log2(sweep(frame, 2, library_size, "/") * 1e6 + 1)
    return(list(matrix = normalized, kind = "counts"))
  }
  normalized <- as.matrix(frame)
  if (min(normalized, na.rm = TRUE) < 0) {
    normalized <- normalized - min(normalized, na.rm = TRUE) + 1
  }
  list(matrix = normalized, kind = "normalized_or_microarray")
}

write_matrix <- function(matrix, path) {
  out <- data.frame(gene = rownames(matrix), matrix, check.names = FALSE)
  utils::write.csv(out, path, row.names = FALSE, quote = FALSE)
}

primary <- load_cohort(config$primary)
primary_expr <- normalize_expression(primary$expression)
primary_matrix <- primary_expr$matrix
validation <- list()
if (!is.null(config$validation) && length(config$validation) > 0) {
  validation <- lapply(config$validation, load_cohort)
}

r_summary <- list(
  status = "completed",
  expression_kind = primary_expr$kind,
  n_discovery_samples = ncol(primary_matrix),
  n_genes = nrow(primary_matrix),
  wgcna = list(status = "skipped", reason = "not run"),
  immune = list(status = "skipped", reason = "not run"),
  survival = list(status = "skipped", reason = "not run"),
  de = list(status = "skipped", reason = "limma unavailable")
)

# ---------------------------------------------------------------------------
# Discovery/validation batch correction and limma differential expression.
# ---------------------------------------------------------------------------
if (requireNamespace("limma", quietly = TRUE)) {
  corrected <- primary_matrix
  batch_cfg <- config$batch
  if (
    !is.null(batch_cfg) &&
    isTRUE(batch_cfg$enabled) &&
    requireNamespace("sva", quietly = TRUE) &&
    length(validation) > 0
  ) {
    common_genes <- Reduce(
      intersect,
      c(list(rownames(primary_matrix)), lapply(validation, function(x) rownames(x$expression)))
    )
    combined <- cbind(
      primary_matrix[common_genes, , drop = FALSE],
      do.call(cbind, lapply(validation, function(x) {
        normalized <- normalize_expression(x$expression)$matrix
        normalized[common_genes, , drop = FALSE]
      }))
    )
    batch <- c(
      rep(primary$name, ncol(primary_matrix)),
      unlist(lapply(validation, function(x) rep(x$name, ncol(x$expression))))
    )
    condition <- c(
      as.character(primary$labels),
      unlist(lapply(validation, function(x) as.character(x$labels)))
    )
    n_primary <- ncol(primary_matrix)
    keep <- stats::complete.cases(t(combined))
    combined <- combined[, keep, drop = FALSE]
    batch <- batch[keep]
    condition <- condition[keep]
    if (length(unique(batch)) > 1 && length(unique(condition)) > 1) {
      mod <- stats::model.matrix(~ condition)
      corrected_all <- sva::ComBat(
        dat = as.matrix(combined),
        batch = batch,
        mod = mod,
        par.prior = TRUE,
        prior.plots = FALSE
      )
      primary_keep <- keep[seq_len(n_primary)]
      if (!any(primary_keep)) {
        stop("ComBat removed every discovery sample; check missing values")
      }
      corrected <- corrected_all[
        ,
        which(primary_keep),
        drop = FALSE
      ]
      colnames(corrected) <- colnames(primary_matrix)[primary_keep]
      primary_matrix <- corrected
      primary$expression <- primary$expression[
        ,
        primary_keep,
        drop = FALSE
      ]
      primary$metadata <- primary$metadata[
        primary_keep,
        ,
        drop = FALSE
      ]
      primary$labels <- primary$labels[primary_keep]
      write_matrix(
        primary_matrix,
        file.path(out_dir, "expression_primary_corrected.csv")
      )
      r_summary$batch_correction <- list(
        status = "completed",
        method = "ComBat",
        cohorts = unique(batch),
        discovery_samples_removed = sum(!primary_keep)
      )
    }
  }
  design <- stats::model.matrix(~ primary$labels)
  fit <- limma::lmFit(corrected, design)
  fit <- limma::eBayes(fit)
  de <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  de$gene <- rownames(de)
  utils::write.csv(de, file.path(out_dir, "deg_primary_limma.csv"), row.names = FALSE)
  r_summary$de <- list(
    status = "completed",
    method = "limma",
    genes = nrow(de),
    significant = sum(de$adj.P.Val <= 0.05, na.rm = TRUE)
  )
} else {
  log_msg("limma is not installed; Python Welch test fallback will be used")
}

# ---------------------------------------------------------------------------
# WGCNA.
# ---------------------------------------------------------------------------
if (
  !is.null(config$wgcna) &&
  isTRUE(config$wgcna$enabled) &&
  requireNamespace("WGCNA", quietly = TRUE)
) {
  wgcna_cfg <- config$wgcna
  dat_expr <- t(primary_matrix)
  dat_expr <- dat_expr[stats::complete.cases(dat_expr), , drop = FALSE]
  vars <- apply(dat_expr, 2, stats::var, na.rm = TRUE)
  keep_genes <- names(sort(vars, decreasing = TRUE))[
    seq_len(min(length(vars), as.integer(wgcna_cfg$top_genes %||% 8000)))
  ]
  dat_expr <- dat_expr[, keep_genes, drop = FALSE]
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  if (!gsg$allOK) {
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  powers <- c(1:20)
  sft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = "unsigned",
    verbose = 0
  )
  fit <- sft$fitIndices
  valid <- which(fit[, "SFT.R.sq"] >= 0.8)
  power <- if (length(valid)) fit[valid[[1]], "Power"] else fit[which.max(fit[, "SFT.R.sq"]), "Power"]
  net <- WGCNA::blockwiseModules(
    dat_expr,
    power = power,
    TOMType = "unsigned",
    minModuleSize = as.integer(wgcna_cfg$min_module_size %||% 30),
    mergeCutHeight = as.numeric(wgcna_cfg$merge_cut_height %||% 0.25),
    numericLabels = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    verbose = 0
  )
  traits <- data.frame(
    condition = as.numeric(primary$labels == primary$case_label),
    row.names = rownames(dat_expr)
  )
  module_trait_cor <- stats::cor(net$MEs, traits, use = "p")
  module_trait_p <- WGCNA::corPvalueStudent(module_trait_cor, nrow(dat_expr))
  module_table <- data.frame(
    module = rownames(module_trait_cor),
    correlation = module_trait_cor[, 1],
    pvalue = module_trait_p[, 1],
    row.names = NULL
  )
  utils::write.csv(
    module_table,
    file.path(out_dir, "wgcna_module_trait.csv"),
    row.names = FALSE
  )
  kme <- as.data.frame(
    stats::cor(dat_expr, net$MEs, use = "p"),
    check.names = FALSE
  )
  hub_rows <- list()
  for (module in module_table$module[module_table$pvalue < 0.05]) {
    if (!module %in% colnames(kme)) next
    membership <- kme[[module]]
    selected <- names(membership)[
      abs(membership) >= as.numeric(wgcna_cfg$kme_cutoff %||% 0.8)
    ]
    if (length(selected)) {
      hub_rows[[length(hub_rows) + 1]] <- data.frame(
        gene = selected,
        module = module,
        kME = membership[selected],
        module_correlation = module_table$correlation[
          match(module, module_table$module)
        ],
        module_pvalue = module_table$pvalue[
          match(module, module_table$module)
        ],
        row.names = NULL
      )
    }
  }
  if (length(hub_rows)) {
    hubs <- do.call(rbind, hub_rows)
    utils::write.csv(
      hubs,
      file.path(out_dir, "wgcna_hubs.csv"),
      row.names = FALSE
    )
  } else {
    utils::write.csv(
      data.frame(gene = character(), module = character(), kME = numeric()),
      file.path(out_dir, "wgcna_hubs.csv"),
      row.names = FALSE
    )
  }
  utils::write.csv(
    data.frame(power = fit$Power, fit = fit$SFT.R.sq, mean_k = fit$mean.k.),
    file.path(out_dir, "wgcna_soft_threshold.csv"),
    row.names = FALSE
  )
  r_summary$wgcna <- list(
    status = "completed",
    power = power,
    modules = length(unique(net$colors)),
    hub_genes = if (length(hub_rows)) nrow(do.call(rbind, hub_rows)) else 0
  )
} else if (!is.null(config$wgcna) && isTRUE(config$wgcna$enabled)) {
  r_summary$wgcna <- list(status = "skipped", reason = "WGCNA not installed")
}

# ---------------------------------------------------------------------------
# Immune deconvolution. immunedeconv is preferred; GSVA custom signatures
# and a clear skipped status are supported for reproducible offline runs.
# ---------------------------------------------------------------------------
if (
  !is.null(config$immune) &&
  isTRUE(config$immune$enabled)
) {
  immune_cfg <- config$immune
  method <- if (is.null(immune_cfg$method)) "ssgsea" else tolower(as.character(immune_cfg$method))
  if (requireNamespace("immunedeconv", quietly = TRUE) && method != "custom") {
    immune <- tryCatch(
      immunedeconv::deconvolute(
        as.matrix(primary_matrix),
        method = method,
        tumor = FALSE
      ),
      error = function(e) e
    )
    if (!inherits(immune, "error")) {
      utils::write.csv(
        immune,
        file.path(out_dir, "immune_deconvolution.csv"),
        row.names = FALSE
      )
      r_summary$immune <- list(status = "completed", method = method)
    } else {
      r_summary$immune <- list(status = "failed", method = method, reason = conditionMessage(immune))
    }
  } else if (
    requireNamespace("GSVA", quietly = TRUE) &&
    !is.null(config$immune$gene_sets)
  ) {
    gene_sets <- lapply(config$immune$gene_sets, function(x) unique(as.character(x)))
    scores <- GSVA::gsva(
      as.matrix(primary_matrix),
      gene_sets,
      method = "ssgsea",
      verbose = FALSE
    )
    utils::write.csv(
      data.frame(cell_type = rownames(scores), scores, check.names = FALSE),
      file.path(out_dir, "immune_deconvolution.csv"),
      row.names = FALSE
    )
    r_summary$immune <- list(status = "completed", method = "GSVA ssGSEA")
  } else {
    r_summary$immune <- list(
      status = "skipped",
      reason = "immunedeconv/GSVA or custom immune gene sets unavailable"
    )
  }
}

# ---------------------------------------------------------------------------
# Survival analysis. This is cohort-specific and only runs when the metadata
# provides time and event columns.
# ---------------------------------------------------------------------------
if (
  !is.null(config$survival) &&
  isTRUE(config$survival$enabled) &&
  requireNamespace("survival", quietly = TRUE)
) {
  survival_cfg <- config$survival
  time_col <- survival_cfg$time_column
  event_col <- survival_cfg$event_column
  if (
    !is.null(time_col) &&
    !is.null(event_col) &&
    all(c(time_col, event_col) %in% colnames(primary$metadata))
  ) {
    time <- suppressWarnings(as.numeric(primary$metadata[[time_col]]))
    event <- suppressWarnings(as.numeric(primary$metadata[[event_col]]))
    valid <- is.finite(time) & is.finite(event) & time > 0
    survival_rows <- list()
    for (gene in rownames(primary_matrix)) {
      values <- suppressWarnings(as.numeric(primary_matrix[gene, ]))
      keep <- valid & is.finite(values)
      if (sum(keep) < 10 || length(unique(values[keep])) < 2) next
      frame <- data.frame(
        time = time[keep],
        event = event[keep],
        expression = values[keep]
      )
      fit <- tryCatch(
        survival::coxph(
          survival::Surv(time, event) ~ expression,
          data = frame
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      summary_fit <- summary(fit)
      survival_rows[[length(survival_rows) + 1]] <- data.frame(
        gene = gene,
        hazard_ratio = exp(stats::coef(fit))[[1]],
        lower95 = exp(stats::confint(fit))[[1]],
        upper95 = exp(stats::confint(fit))[[2]],
        pvalue = summary_fit$coefficients[, "Pr(>|z|)"][[1]],
        row.names = NULL
      )
    }
    if (length(survival_rows)) {
      survival_frame <- do.call(rbind, survival_rows)
      survival_frame$padj <- stats::p.adjust(survival_frame$pvalue, method = "BH")
      survival_frame <- survival_frame[
        order(survival_frame$padj, survival_frame$pvalue),
        ,
        drop = FALSE
      ]
      utils::write.csv(
        survival_frame,
        file.path(out_dir, "survival_cox.csv"),
        row.names = FALSE
      )
      r_summary$survival <- list(
        status = "completed",
        genes = nrow(survival_frame),
        significant = sum(survival_frame$padj <= 0.05, na.rm = TRUE)
      )
    } else {
      r_summary$survival <- list(status = "skipped", reason = "no valid survival fits")
    }
  } else {
    r_summary$survival <- list(
      status = "skipped",
      reason = "time/event columns not configured"
    )
  }
} else if (!is.null(config$survival) && isTRUE(config$survival$enabled)) {
  r_summary$survival <- list(status = "skipped", reason = "survival not installed")
}

jsonlite::write_json(
  r_summary,
  file.path(out_dir, "advanced_r_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
