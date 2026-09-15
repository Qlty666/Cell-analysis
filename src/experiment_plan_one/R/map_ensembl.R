#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: map_ensembl.R <input_ids.txt> <output.csv>")
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

ids <- unique(trimws(readLines(args[[1]], warn = FALSE)))
ids <- ids[nzchar(ids)]
clean <- sub("\\..*$", "", ids)
symbols <- suppressMessages(
  AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = clean,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
)
output <- data.frame(
  ensembl_id = clean,
  symbol = toupper(ifelse(is.na(symbols), "", symbols)),
  stringsAsFactors = FALSE
)
output <- output[nzchar(output$symbol), , drop = FALSE]
dir.create(dirname(args[[2]]), recursive = TRUE, showWarnings = FALSE)
write.csv(output, args[[2]], row.names = FALSE)
