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
universe_ids <- NULL
if (length(args) >= 4 && nzchar(args[4])) {
  background <- read.csv(args[4], stringsAsFactors = FALSE)$gene
  unique_background <- unique(as.character(background))
  unique_background <- unique_background[
    !is.na(unique_background) & nzchar(unique_background)
  ]
  mapped_background <- suppressWarnings(bitr(
    unique_background,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  ))
  background_mapping <- data.frame(
    gene = unique_background,
    entrez_id = mapped_background$ENTREZID[
      match(unique_background, mapped_background$SYMBOL)
    ],
    mapping_status = ifelse(
      unique_background %in% mapped_background$SYMBOL,
      "mapped",
      "unmapped"
    ),
    stringsAsFactors = FALSE
  )
  universe_ids <- unique(mapped_background$ENTREZID)
  if (length(universe_ids) < 3) stop("fewer than three background genes mapped")
  write.csv(background_mapping, file.path(out_dir, "background_id_mapping.csv"), row.names = FALSE)
}

eg <- tryCatch(
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db),
  error = function(e) NULL
)
query_mapping <- data.frame(
  gene = genes,
  entrez_id = NA_character_,
  mapping_status = "unmapped",
  stringsAsFactors = FALSE
)
if (!is.null(eg) && nrow(eg) > 0) {
  query_mapping$entrez_id <- eg$ENTREZID[match(genes, eg$SYMBOL)]
  query_mapping$mapping_status <- ifelse(
    genes %in% eg$SYMBOL,
    "mapped",
    "unmapped"
  )
}
write.csv(query_mapping, file.path(out_dir, "query_id_mapping.csv"), row.names = FALSE)
if (is.null(eg) || nrow(eg) == 0) {
  write.csv(data.frame(note = "no gene ID mapping"),
            file.path(out_dir, "insilico_go_enrichment.csv"), row.names = FALSE)
  write.csv(data.frame(note = "no gene ID mapping"),
            file.path(out_dir, "insilico_kegg_enrichment.csv"), row.names = FALSE)
  quit(status = 0)
}
if (!is.null(universe_ids) && any(!eg$ENTREZID %in% universe_ids)) {
  stop("query genes are outside the specified testable background")
}

package_versions <- data.frame(
  package = c(
    "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "KEGG_REST"
  ),
  version = c(
    as.character(utils::packageVersion("clusterProfiler")),
    as.character(utils::packageVersion("org.Hs.eg.db")),
    as.character(utils::packageVersion("org.Mm.eg.db")),
    paste0("query_date=", Sys.Date())
  ),
  stringsAsFactors = FALSE
)
write.csv(
  package_versions,
  file.path(out_dir, "enrichment_database_versions.csv"),
  row.names = FALSE
)

go_parts <- lapply(c("BP", "CC", "MF"), function(ont) {
  res <- tryCatch(
    enrichGO(
      gene = unique(eg$ENTREZID),
      OrgDb = org_db,
      keyType = "ENTREZID",
      ont = ont,
      universe = universe_ids,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
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
             universe = universe_ids, pAdjustMethod = "BH", pvalueCutoff = 1,
             qvalueCutoff = 1),
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
