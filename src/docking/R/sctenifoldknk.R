#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, required = TRUE) {
  prefix <- paste0("--", name, "=")
  value <- args[startsWith(args, prefix)]
  if (!length(value)) {
    if (required) stop(sprintf("missing required argument --%s", name))
    return(NULL)
  }
  sub(prefix, "", value[[1]], fixed = TRUE)
}

matrix_path <- get_arg("matrix")
gene <- get_arg("gene")
output_dir <- get_arg("output")
get_arg_or <- function(name, default) {
  value <- get_arg(name, required = FALSE)
  if (is.null(value)) default else value
}
n_networks <- as.integer(get_arg_or("n-networks", "5"))
n_cells <- as.integer(get_arg_or("n-cells", "500"))
n_comp <- as.integer(get_arg_or("n-comp", "3"))
q_value <- as.numeric(get_arg_or("q", "0.95"))
k_value <- as.integer(get_arg_or("k", "3"))
ma_dim <- as.integer(get_arg_or("ma-dim", "2"))
min_lib_size <- as.numeric(get_arg_or("min-lib-size", "0"))
min_percent <- as.numeric(get_arg_or("min-percent", "0"))
min_exp_sum <- as.numeric(get_arg_or("min-exp-sum", "0"))
max_mito_ratio <- as.numeric(get_arg_or("max-mito-ratio", "1"))
remove_outliers <- tolower(
  get_arg_or("remove-outliers", "false")
) %in% c("true", "1", "yes")
seed <- as.integer(get_arg_or("seed", "123"))
n_cores <- as.integer(get_arg_or("n-cores", "1"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(jsonlite)
  library(scTenifoldKnk)
})

counts_frame <- read.csv(
  matrix_path,
  check.names = FALSE,
  row.names = 1
)
counts <- as.matrix(counts_frame)
storage.mode(counts) <- "double"
if (!gene %in% rownames(counts)) {
  stop(sprintf("knockout gene '%s' is absent from the count matrix", gene))
}
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("scTenifoldKnk requires finite non-negative raw counts")
}

set.seed(seed)
result <- scTenifoldKnk::scTenifoldKnk(
  countMatrix = counts,
  gKO = gene,
  transcriptomeWide = FALSE,
  qc = min_lib_size > 0 || min_percent > 0 || min_exp_sum > 0,
  qc_minLibSize = min_lib_size,
  qc_removeOutlierCells = remove_outliers,
  qc_minPCT = min_percent,
  qc_maxMTratio = max_mito_ratio,
  nc_nNet = n_networks,
  nc_nCells = min(n_cells, ncol(counts)),
  nc_nComp = n_comp,
  nc_q = q_value,
  td_K = k_value,
  ma_nDim = ma_dim,
  nCores = n_cores
)

diff <- result$diffRegulation
if (is.null(diff) || !nrow(diff)) {
  stop("scTenifoldKnk returned no differential regulation result")
}
write.csv(
  diff,
  file.path(output_dir, "sctenifoldknk_results.csv"),
  row.names = FALSE,
  na = ""
)
writeLines(
  jsonlite::toJSON(
    list(
      engine = "scTenifoldKnk (R package)",
      package_version = as.character(utils::packageVersion("scTenifoldKnk")),
      gene = gene,
      genes_scored = nrow(diff),
      n_networks = n_networks,
      n_cells_per_network = min(n_cells, ncol(counts)),
      n_comp = n_comp,
      q = q_value,
      k = k_value,
      ma_dim = ma_dim,
      qc = min_lib_size > 0 || min_percent > 0 || min_exp_sum > 0
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  ),
  file.path(output_dir, "sctenifoldknk_summary.json")
)
