# read_10x.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

read_10x_set <- function(set_name) {
  prefix <- file.path(raw_dir, paste0("GSE125449_Set", set_name))
  barcodes <- as.character(
    fread(paste0(prefix, "_barcodes.tsv.gz"), header = FALSE)$V1
  )
  genes <- fread(paste0(prefix, "_genes.tsv.gz"), header = FALSE)
  samples_tbl <- fread(paste0(prefix, "_samples.txt.gz"), header = TRUE)

  con <- gzfile(paste0(prefix, "_matrix.mtx.gz"), "rt")
  mtx <- readMM(con)
  close(con)
  mtx <- as(mtx, "CsparseMatrix")
  rownames(mtx) <- make.unique(genes$V2)
  colnames(mtx) <- barcodes

  meta <- data.frame(
    row.names = barcodes,
    sample = samples_tbl$Sample[match(barcodes, samples_tbl[["Cell Barcode"]])],
    published_type = samples_tbl$Type[match(barcodes, samples_tbl[["Cell Barcode"]])],
    set = set_name,
    stringsAsFactors = FALSE
  )
  list(counts = mtx, meta = meta)
}

load_manifest <- function() {
  path <- file.path(root, "data", paste0(accession, "_manifest.json"))
  if (!file.exists(path)) {
    stop("Dataset manifest not found: ", path)
  }
  jsonlite::fromJSON(path)
}

read_mtx_fallback <- function(path) {
  con <- gzfile(path, "rt")
  head_lines <- readLines(con, n = 10, warn = FALSE)
  close(con)
  first_data <- which(!startsWith(head_lines, "%"))[1]
  dims <- as.integer(strsplit(trimws(head_lines[first_data]), "\\s+")[[1]])
  if (length(dims) >= 2) {
    dims <- dims[1:2]
  }
  tab <- fread(path, header = FALSE, skip = first_data, fill = TRUE)
  tab <- tab[!is.na(V1) & !is.na(V2)]
  if (length(dims) != 2 || any(is.na(dims))) {
    dims <- c(max(tab$V1), max(tab$V2))
  }
  m <- sparseMatrix(
    i = tab$V1,
    j = tab$V2,
    x = tab$V3,
    dims = dims
  )
  as(m, "CsparseMatrix")
}

normalize_ensembl_ids <- function(ids) {
  ids <- as.character(ids)
  hit <- grepl("^(ENSG|ENSMUSG|ENST|ENSMUST)", ids)
  ids[hit] <- sub("\\.[0-9]+$", "", ids[hit])
  ids
}

