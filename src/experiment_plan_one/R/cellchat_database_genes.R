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

output_path <- get_arg("output")
summary_path <- get_arg("summary")
species_arg <- get_arg("species", required = FALSE)
species <- tolower(if (is.null(species_arg)) "mm" else species_arg)

suppressPackageStartupMessages({
  library(CellChat)
  library(jsonlite)
})

database <- if (species %in% c("hs", "human", "homo_sapiens")) {
  CellChatDB.human
} else if (species %in% c("mm", "mouse", "mus_musculus")) {
  CellChatDB.mouse
} else {
  stop(sprintf("unsupported CellChat species: %s", species))
}

gene_info <- database$geneInfo
symbol_column <- if ("Symbol" %in% names(gene_info)) {
  "Symbol"
} else if ("symbol" %in% names(gene_info)) {
  "symbol"
} else {
  names(gene_info)[[1]]
}
all_symbols <- unique(as.character(gene_info[[symbol_column]]))
all_symbols <- all_symbols[!is.na(all_symbols) & nzchar(all_symbols)]

collect_tokens <- function(x) {
  text <- as.character(unlist(x, use.names = FALSE))
  text <- text[!is.na(text) & nzchar(text)]
  tokens <- unlist(
    strsplit(text, "[^A-Za-z0-9_.-]+", perl = TRUE),
    use.names = FALSE
  )
  unique(tokens[nzchar(tokens)])
}

tokens <- collect_tokens(
  list(database$interaction, database$complex, database$cofactor)
)
symbols <- all_symbols[toupper(all_symbols) %in% toupper(tokens)]
if (length(symbols) < 10) symbols <- all_symbols

writeLines(symbols, output_path, useBytes = TRUE)
writeLines(
  jsonlite::toJSON(
    list(
      method = "R CellChat",
      cellchat_version = as.character(utils::packageVersion("CellChat")),
      database = if (species %in% c("hs", "human", "homo_sapiens")) {
        "CellChatDB.human"
      } else {
        "CellChatDB.mouse"
      },
      symbol_column = symbol_column,
      interaction_database_genes = length(symbols)
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  ),
  summary_path
)
