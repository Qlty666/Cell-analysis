#!/usr/bin/env Rscript
# GO/KEGG enrichment for in-silico knockout target genes.
# Usage: Rscript insilico_enrichment.R <gene_csv> <out_dir> <species>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: Rscript insilico_enrichment.R <gene_csv> <out_dir> <species>")
}
gene_csv <- args[1]
out_dir <- args[2]
species <- tolower(args[3])
if (!species %in% c("hs", "mm")) species <- "hs"

suppressWarnings(suppressPackageStartupMessages({
  library(clusterProfiler)
  library(httr)
  library(R.utils)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
}))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

kegg_cache_dir <- file.path(getwd(), "data_cache", "kegg")
dir.create(kegg_cache_dir, recursive = TRUE, showWarnings = FALSE)

kegg_fetch_text <- function(url, timeout = 600, attempts = 3) {
  last_error <- NULL
  for (i in seq_len(attempts)) {
    resp <- tryCatch(
      httr::GET(url, httr::timeout(timeout)),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) {
      return(httr::content(resp, as = "text", encoding = "UTF-8"))
    }
    if (i < attempts) Sys.sleep(3 * i)
  }
  NULL
}

kegg_rest_cached <- function(rest_url) {
  cache_file <- file.path(
    kegg_cache_dir,
    paste0(gsub("[^A-Za-z0-9]", "_", rest_url), ".rds")
  )
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  content <- kegg_fetch_text(rest_url)
  if (is.null(content)) {
    stop("KEGG REST request failed for ", rest_url)
  }
  lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("KEGG REST returned empty content")
  }
  mat <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
  res <- data.frame(from = mat[, 1], to = mat[, 2], stringsAsFactors = FALSE)
  saveRDS(res, cache_file)
  res
}

suppressWarnings(
  tryCatch(
    utils::assignInNamespace(
      "kegg_rest",
      kegg_rest_cached,
      ns = "clusterProfiler"
    ),
    error = function(e) NULL
  )
)

genes <- read.csv(gene_csv, stringsAsFactors = FALSE)$gene
genes <- unique(as.character(genes))
genes <- genes[nzchar(genes) & !is.na(genes)]
if (length(genes) < 3) quit(status = 0)

org_pkg <- if (species == "mm") "org.Mm.eg.db" else "org.Hs.eg.db"
org_db <- getExportedValue(org_pkg, org_pkg)
kegg_org <- ifelse(species == "mm", "mmu", "hsa")

eg <- tryCatch(
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db),
  error = function(e) NULL
)
if (is.null(eg) || nrow(eg) == 0) {
  write.csv(data.frame(note = "no gene ID mapping"),
            file.path(out_dir, "insilico_go_enrichment.csv"), row.names = FALSE)
  write.csv(data.frame(note = "no gene ID mapping"),
            file.path(out_dir, "insilico_kegg_enrichment.csv"), row.names = FALSE)
  quit(status = 0)
}

go_parts <- lapply(c("BP", "CC", "MF"), function(ont) {
  res <- tryCatch(
    enrichGO(
      gene = unique(eg$ENTREZID),
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.3,
      readable = TRUE
    ),
    error = function(e) NULL
  )
  if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
    df <- as.data.frame(res)
    df$ONTOLOGY <- ont
    return(df)
  }
  NULL
})
go_df <- do.call(rbind, go_parts)
if (!is.null(go_df) && nrow(go_df) > 0) {
  write.csv(go_df, file.path(out_dir, "insilico_go_enrichment.csv"),
            row.names = FALSE)
} else {
  write.csv(data.frame(note = "no significant GO terms"),
            file.path(out_dir, "insilico_go_enrichment.csv"), row.names = FALSE)
}

kegg <- tryCatch(
  enrichKEGG(gene = unique(eg$ENTREZID), organism = kegg_org,
             pvalueCutoff = 0.2),
  error = function(e) NULL
)
if (!is.null(kegg)) {
  kegg <- tryCatch(
    setReadable(kegg, OrgDb = org_db, keyType = "ENTREZID"),
    error = function(e) kegg
  )
}
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  write.csv(as.data.frame(kegg), file.path(out_dir, "insilico_kegg_enrichment.csv"),
            row.names = FALSE)
} else {
  write.csv(data.frame(note = "no significant KEGG terms"),
            file.path(out_dir, "insilico_kegg_enrichment.csv"), row.names = FALSE)
}
quit(status = 0)
