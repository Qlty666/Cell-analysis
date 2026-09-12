#!/usr/bin/env Rscript
# Optional R backend for Mendelian randomisation and colocalisation.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: mr_coloc.R <harmonised.csv> <config.json> <output_dir>")
}
harmonised_path <- normalizePath(args[[1]], mustWork = TRUE)
config_path <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(jsonlite)
})

config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])) fallback else value
}
data <- utils::read.csv(harmonised_path, check.names = FALSE)
required <- c(
  "snp",
  "beta_exposure",
  "se_exposure",
  "beta_outcome",
  "se_outcome"
)
if (!all(required %in% colnames(data))) {
  stop("harmonised table is missing columns: ",
       paste(setdiff(required, colnames(data)), collapse = ", "))
}

summary <- list(status = "completed", methods = list(), coloc = list(status = "skipped"))
methods <- list()

if (requireNamespace("MendelianRandomization", quietly = TRUE)) {
  input <- MendelianRandomization::mr_input(
    bx = data$beta_exposure,
    bxse = data$se_exposure,
    by = data$beta_outcome,
    byse = data$se_outcome,
    snps = data$snp
  )
  fit_method <- function(name, fit) {
    data.frame(
      method = name,
      estimate = as.numeric(fit@Estimate),
      se = as.numeric(fit@StdError),
      pvalue = as.numeric(fit@Pvalue),
      nsnp = length(data$snp)
    )
  }
  methods[[length(methods) + 1]] <- fit_method(
    "inverse_variance_weighted",
    MendelianRandomization::mr_ivw(input)
  )
  methods[[length(methods) + 1]] <- fit_method(
    "weighted_median",
    MendelianRandomization::mr_median(input)
  )
  if (nrow(data) >= 3) {
    methods[[length(methods) + 1]] <- fit_method(
      "mr_egger",
      MendelianRandomization::mr_egger(input)
    )
  }
  methods_frame <- do.call(rbind, methods)
  utils::write.csv(
    methods_frame,
    file.path(out_dir, "mr_r_methods.csv"),
    row.names = FALSE
  )
  summary$methods <- methods_frame
} else {
  summary$methods <- list(status = "skipped", reason = "MendelianRandomization not installed")
}

if (
  isTRUE(config$coloc$enabled) &&
  requireNamespace("coloc", quietly = TRUE)
) {
  exposure_n <- config$coloc$exposure_n
  outcome_n <- config$coloc$outcome_n
  if (is.null(exposure_n) || is.null(outcome_n)) {
    summary$coloc <- list(
      status = "skipped",
      reason = "coloc.exposure_n and coloc.outcome_n are required"
    )
  } else {
    dataset1 <- list(
      beta = data$beta_exposure,
      varbeta = data$se_exposure^2,
      snp = data$snp,
      type = "quant",
      N = as.numeric(exposure_n)
    )
    dataset2 <- list(
      beta = data$beta_outcome,
      varbeta = data$se_outcome^2,
      snp = data$snp,
      type = "quant",
      N = as.numeric(outcome_n)
    )
    result <- tryCatch(
      coloc::coloc.abf(
        dataset1 = dataset1,
        dataset2 = dataset2,
        p12 = as.numeric(config$coloc$p12 %||% 1e-5)
      ),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      summary$coloc <- list(status = "failed", reason = conditionMessage(result))
    } else {
      posterior <- as.data.frame(t(result$summary))
      posterior$term <- rownames(posterior)
      utils::write.csv(
        posterior,
        file.path(out_dir, "coloc_posterior.csv"),
        row.names = FALSE
      )
      summary$coloc <- list(
        status = "completed",
        PP.H4 = as.numeric(result$summary["PP.H4.abf"]),
        output = file.path(out_dir, "coloc_posterior.csv")
      )
    }
  }
} else if (isTRUE(config$coloc$enabled)) {
  summary$coloc <- list(status = "skipped", reason = "coloc not installed")
}

jsonlite::write_json(
  summary,
  file.path(out_dir, "mr_r_summary.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
