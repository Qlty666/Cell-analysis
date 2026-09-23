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
genes_path <- get_arg("genes")
cells_path <- get_arg("cells")
metadata_path <- get_arg("metadata")
output_dir <- get_arg("output")
species_arg <- get_arg("species", required = FALSE)
group_arg <- get_arg("group-column", required = FALSE)
celltype_arg <- get_arg("celltype-column", required = FALSE)
min_cells_arg <- get_arg("min-cells", required = FALSE)
nboot_arg <- get_arg("nboot", required = FALSE)
groups_arg <- get_arg("groups", required = FALSE)
species <- tolower(if (is.null(species_arg)) "mm" else species_arg)
group_column <- if (is.null(group_arg)) "condition" else group_arg
celltype_column <- if (is.null(celltype_arg)) "cell_type" else celltype_arg
min_cells <- as.integer(if (is.null(min_cells_arg)) "10" else min_cells_arg)
nboot <- as.integer(if (is.null(nboot_arg)) "100" else nboot_arg)
requested_groups <- if (is.null(groups_arg)) {
  character()
} else {
  trimws(strsplit(groups_arg, ",", fixed = TRUE)[[1]])
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(jsonlite)
})

counts <- as(Matrix::readMM(matrix_path), "dgCMatrix")
genes <- readLines(genes_path, warn = FALSE)
cells <- readLines(cells_path, warn = FALSE)
if (nrow(counts) != length(genes) || ncol(counts) != length(cells)) {
  stop("matrix dimensions do not match gene and cell identifiers")
}
rownames(counts) <- genes
colnames(counts) <- cells

metadata <- read.delim(
  metadata_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!"cell" %in% names(metadata)) stop("metadata must contain a cell column")
if (!group_column %in% names(metadata)) {
  stop(sprintf("metadata is missing group column '%s'", group_column))
}
if (!celltype_column %in% names(metadata)) {
  stop(sprintf("metadata is missing cell-type column '%s'", celltype_column))
}
rownames(metadata) <- as.character(metadata$cell)
metadata <- metadata[cells, , drop = FALSE]
if (anyNA(rownames(metadata))) stop("metadata cells do not match matrix cells")

database <- if (species %in% c("hs", "human", "homo_sapiens")) {
  CellChatDB.human
} else if (species %in% c("mm", "mouse", "mus_musculus")) {
  CellChatDB.mouse
} else {
  stop(sprintf("unsupported CellChat species: %s", species))
}

input_gene_count <- nrow(counts)
gene_info <- database$geneInfo
if (!is.data.frame(gene_info) || !nrow(gene_info)) {
  stop("CellChat species database has no geneInfo table")
}
symbol_column <- if ("Symbol" %in% names(gene_info)) {
  "Symbol"
} else if ("symbol" %in% names(gene_info)) {
  "symbol"
} else {
  names(gene_info)[[1]]
}
all_database_symbols <- unique(as.character(gene_info[[symbol_column]]))
all_database_symbols <- all_database_symbols[
  !is.na(all_database_symbols) & nzchar(all_database_symbols)
]
collect_tokens <- function(x) {
  text <- as.character(unlist(x, use.names = FALSE))
  text <- text[!is.na(text) & nzchar(text)]
  tokens <- unlist(
    strsplit(text, "[^A-Za-z0-9_.-]+", perl = TRUE),
    use.names = FALSE
  )
  unique(tokens[nzchar(tokens)])
}
database_tokens <- collect_tokens(
  list(database$interaction, database$complex, database$cofactor)
)
database_symbols <- all_database_symbols[
  toupper(all_database_symbols) %in% toupper(database_tokens)
]
if (length(database_symbols) < 10) {
  database_symbols <- all_database_symbols
}
database_symbols <- database_symbols[
  !is.na(database_symbols) & nzchar(database_symbols)
]
match_index <- match(
  toupper(database_symbols),
  toupper(rownames(counts))
)
matched <- !is.na(match_index)
canonical_symbols <- database_symbols[matched]
matched_rows <- match_index[matched]
duplicated_symbol <- duplicated(canonical_symbols)
canonical_symbols <- canonical_symbols[!duplicated_symbol]
matched_rows <- matched_rows[!duplicated_symbol]
counts <- counts[matched_rows, , drop = FALSE]
rownames(counts) <- canonical_symbols
mapped_gene_count <- nrow(counts)
gene_mapping <- list(
  input_genes = as.integer(input_gene_count),
  mapped_genes = as.integer(mapped_gene_count),
  interaction_database_genes = as.integer(length(database_symbols)),
  mapping_rate = if (input_gene_count) {
    mapped_gene_count / input_gene_count
  } else {
    0
  },
  database_gene_coverage = if (length(database_symbols)) {
    mapped_gene_count / length(database_symbols)
  } else {
    0
  },
  case_insensitive = TRUE,
  database_symbol_column = symbol_column
)
if (mapped_gene_count == 0) {
  write.csv(
    data.frame(),
    file.path(output_dir, "cellchat_r_interactions.csv"),
    row.names = FALSE,
    na = ""
  )
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
        group_column = group_column,
        celltype_column = celltype_column,
        requested_groups = requested_groups,
        excluded_groups = character(),
        gene_mapping = gene_mapping,
        min_cells = min_cells,
        n_bootstrap = nboot,
        inference_status = "cell_level_probability_only",
        group_statistics = FALSE,
        status = "valid_negative",
        reason = "no input genes matched the CellChat species database",
        groups = list()
      ),
      auto_unbox = TRUE,
      pretty = TRUE,
      null = "null"
    ),
    file.path(output_dir, "cellchat_r_summary.json")
  )
  quit(status = 0)
}

groups <- unique(as.character(metadata[[group_column]]))
groups <- groups[nzchar(groups)]
excluded_groups <- character()
if (length(requested_groups)) {
  excluded_groups <- setdiff(groups, requested_groups)
  groups <- intersect(groups, requested_groups)
}
records <- list()
group_summary <- list()
for (group in groups) {
  group_mask <- as.character(metadata[[group_column]]) == group
  group_meta <- metadata[group_mask, , drop = FALSE]
  cell_types <- as.character(group_meta[[celltype_column]])
  type_counts <- table(cell_types)
  retained_types <- names(type_counts)[type_counts >= min_cells]
  keep <- group_mask & cell_types %in% retained_types
  if (sum(keep) < min_cells || length(retained_types) < 2) {
    group_summary[[group]] <- list(
      status = "insufficient_cells",
      n_cells = as.integer(sum(keep)),
      n_cell_types = as.integer(length(retained_types))
    )
    next
  }

  group_counts <- counts[, keep, drop = FALSE]
  group_meta <- metadata[keep, , drop = FALSE]
  cellchat <- CellChat::createCellChat(
    object = group_counts,
    meta = group_meta,
    group.by = celltype_column
  )
  cellchat@DB <- database
  cellchat@data <- CellChat::normalizeData(cellchat@data)
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  interaction_table <- database$interaction
  required_genes <- strsplit(
    paste(interaction_table$ligand, interaction_table$receptor),
    "[^A-Za-z0-9_.-]+"
  )
  keep_lr <- vapply(
    required_genes,
    function(genes) {
      genes <- genes[nzchar(genes)]
      length(genes) > 0 && all(genes %in% rownames(group_counts))
    },
    logical(1)
  )
  lr_use <- interaction_table[keep_lr, , drop = FALSE]
  if (!nrow(lr_use)) {
    group_summary[[group]] <- list(
      status = "valid_negative",
      n_cells = as.integer(sum(keep)),
      n_cell_types = as.integer(length(retained_types)),
      n_interactions = 0L,
      reason = "no CellChat ligand-receptor pair genes were available"
    )
    next
  }
  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "triMean",
    population.size = FALSE,
    LR.use = lr_use,
    nboot = nboot
  )
  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = min_cells
  )
  interactions <- CellChat::subsetCommunication(cellchat)
  if (is.null(interactions) || !nrow(interactions)) {
    group_summary[[group]] <- list(
      status = "valid_negative",
      n_cells = as.integer(sum(keep)),
      n_cell_types = as.integer(length(retained_types)),
      n_interactions = 0L
    )
    next
  }
  pathways <- tryCatch(
    CellChat::subsetCommunication(cellchat, slot.name = "netP"),
    error = function(error) data.frame()
  )
  interactions$group <- group
  interactions$inference_status <- "cell_level_probability_only"
  interactions$n_cells <- as.integer(sum(keep))
  records[[group]] <- interactions
  group_summary[[group]] <- list(
    status = "completed",
    n_cells = as.integer(sum(keep)),
    n_cell_types = as.integer(length(retained_types)),
    n_interactions = as.integer(nrow(interactions)),
    n_pathways = as.integer(nrow(pathways))
  )
}

if (!length(records)) {
  output <- data.frame()
} else {
  columns <- unique(unlist(lapply(records, names)))
  output <- do.call(
    rbind,
    lapply(records, function(frame) {
      missing <- setdiff(columns, names(frame))
      for (column in missing) frame[[column]] <- NA
      frame[columns]
    })
  )
}
write.csv(
  output,
  file.path(output_dir, "cellchat_r_interactions.csv"),
  row.names = FALSE,
  na = ""
)
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
      group_column = group_column,
      celltype_column = celltype_column,
      requested_groups = requested_groups,
      excluded_groups = excluded_groups,
      gene_mapping = gene_mapping,
      min_cells = min_cells,
      n_bootstrap = nboot,
      inference_status = "cell_level_probability_only",
      group_statistics = FALSE,
      groups = group_summary
    ),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  ),
  file.path(output_dir, "cellchat_r_summary.json")
)
