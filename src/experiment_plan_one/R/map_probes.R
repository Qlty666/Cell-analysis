#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: map_probes.R <GPL14951|GPL570> <probe_ids.txt> <output.csv>")
}

platform <- args[[1]]
ids <- unique(trimws(readLines(args[[2]], warn = FALSE)))
ids <- ids[nzchar(ids)]

if (platform == "GPL14951") {
  suppressPackageStartupMessages(library(illuminaHumanv4.db))
  package <- "illuminaHumanv4.db"
} else if (platform == "GPL570") {
  suppressPackageStartupMessages(library(hgu133plus2.db))
  package <- "hgu133plus2.db"
} else {
  stop(sprintf("unsupported platform: %s", platform))
}

suppressPackageStartupMessages(library(AnnotationDbi))
symbols <- suppressMessages(
  AnnotationDbi::mapIds(
    get(package),
    keys = ids,
    keytype = "PROBEID",
    column = "SYMBOL",
    multiVals = "first"
  )
)
output <- data.frame(
  probe_id = ids,
  symbol = toupper(ifelse(is.na(symbols), "", symbols)),
  stringsAsFactors = FALSE
)
output <- output[nzchar(output$symbol), , drop = FALSE]
dir.create(dirname(args[[3]]), recursive = TRUE, showWarnings = FALSE)
write.csv(output, args[[3]], row.names = FALSE)
