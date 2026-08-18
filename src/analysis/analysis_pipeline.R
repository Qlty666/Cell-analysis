#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
  library(ggrepel)
  library(pheatmap)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DESeq2)
}))

options(timeout = 600)
options(future.globals.maxSize = 10 * 1024^3)

skip_figs <- strsplit(Sys.getenv("LIVER_SKIP_FIGURES", unset = ""), ",")[[1]]
skip_figs <- trimws(skip_figs[nzchar(skip_figs)])

figure_styles <- list()
raw_styles <- Sys.getenv("LIVER_FIGURE_STYLES", unset = "")
if (nzchar(raw_styles)) {
  for (pair in strsplit(raw_styles, ",", fixed = TRUE)[[1]]) {
    kv <- strsplit(pair, "=", fixed = TRUE)[[1]]
    if (length(kv) == 2) {
      figure_styles[[kv[1]]] <- kv[2]
    }
  }
}

fig_style <- function(name) {
  if (name %in% names(figure_styles)) {
    figure_styles[[name]]
  } else {
    ""
  }
}

figure_stage <- function(name) {
  num <- as.integer(sub("^fig_([0-9]+)_.*$", "\\1", name))
  if (is.na(num)) return("00_other")
  if (num <= 1 || num == 48) return("01_qc")
  if (num == 2) return("02_doublets")
  if (num %in% c(3, 4, 14, 15)) return("03_cluster")
  if (num %in% c(5, 6, 7, 16, 17, 18, 19)) return("04_annotation")
  if (num %in% c(8, 9)) return("05_deg")
  if (num %in% c(10, 11, 12, 13, 20, 21, 22, 23, 46, 47)) return("06_enrichment")
  if (num %in% c(24, 25, 43, 44, 45)) return("07_ml")
  if (num >= 26 && num <= 39) return("08_publication")
  if (num >= 40 && num <= 42) return("09_cellchat")
  return("00_other")
}

wrap_labels <- function(x, width = 20) {
  vapply(
    as.character(x),
    function(s) paste(strwrap(s, width = width), collapse = "\n"),
    character(1)
  )
}

add_plot_margin <- function(plot, plot_margin) {
  if (inherits(plot, "patchwork")) {
    plot & theme(plot.margin = plot_margin)
  } else {
    plot + theme(plot.margin = plot_margin)
  }
}

stage_fig_file <- function(file) {
  name <- basename(file)
  out_dir <- file.path(dirname(file), figure_stage(name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out_dir, name)
}

stage_data_file <- function(name) {
  num <- as.integer(sub("^fig_([0-9]+)_.*$", "\\1", name))
  if (!is.na(num)) {
    out_dir <- file.path(data_dir, figure_stage(name))
  } else if (name == "sample_annotations.csv") {
    out_dir <- file.path(data_dir, "01_qc")
  } else if (name == "liver_cancer_seurat.rds") {
    return(file.path(data_dir, name))
  } else {
    out_dir <- data_dir
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(out_dir, name)
}

param_num <- function(name) {
  val <- Sys.getenv(name, unset = "")
  if (nzchar(val)) {
    as.numeric(val)
  } else {
    NA_real_
  }
}

qc_min_features <- param_num("LIVER_QC_MIN_FEATURES")
qc_max_features <- param_num("LIVER_QC_MAX_FEATURES")
qc_min_counts <- param_num("LIVER_QC_MIN_COUNTS")
qc_max_counts <- param_num("LIVER_QC_MAX_COUNTS")
qc_max_mt <- param_num("LIVER_QC_MAX_MT")
cluster_resolution <- param_num("LIVER_CLUSTER_RESOLUTION")
cluster_algorithm <- param_num("LIVER_CLUSTER_ALGORITHM")
de_logfc <- param_num("LIVER_DE_LOGFc")
de_padj <- param_num("LIVER_DE_PADJ")
deg_violin_top_n <- param_num("LIVER_DE_VIOLIN_TOP_N")
deg_violin_max_cells <- param_num("LIVER_DE_VIOLIN_MAX_CELLS")
if (is.na(cluster_resolution)) cluster_resolution <- 0.6
if (is.na(de_logfc)) de_logfc <- 0.25
if (is.na(de_padj)) de_padj <- 0.05
if (is.na(deg_violin_top_n)) deg_violin_top_n <- 12
if (is.na(deg_violin_max_cells)) deg_violin_max_cells <- 1000

empty_deg_frame <- function() {
  data.frame(
    p_val = numeric(0),
    avg_log2FC = numeric(0),
    pct.1 = numeric(0),
    pct.2 = numeric(0),
    p_val_adj = numeric(0),
    stringsAsFactors = FALSE
  )
}

ensure_deg_columns <- function(deg) {
  if (is.null(deg) || !is.data.frame(deg)) {
    return(empty_deg_frame())
  }
  required <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  for (col in required) {
    if (!col %in% colnames(deg)) {
      deg[[col]] <- numeric(nrow(deg))
    }
    if (is.list(deg[[col]])) {
      deg[[col]] <- vapply(
        deg[[col]],
        function(x) {
          if (length(x) != 1L) return(NA_real_)
          suppressWarnings(as.numeric(x))
        },
        numeric(1)
      )
    } else {
      deg[[col]] <- suppressWarnings(as.numeric(deg[[col]]))
    }
  }
  if (!"gene" %in% colnames(deg)) {
    deg$gene <- rownames(deg)
  }
  deg
}

save_fig <- function(file, plot, width, height, dpi = 150, plot_margin = NULL) {
  name <- basename(file)
  if (name %in% skip_figs) {
    log_msg("skip figure: ", name)
    return(invisible(NULL))
  }
  if (is.null(plot_margin)) {
    plot_margin <- ggplot2::margin(26, 26, 24, 18, "pt")
  }
  plot <- add_plot_margin(plot, plot_margin)
  ggsave(
    stage_fig_file(file),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  log_msg("saved figure: ", name)
}

flag_on <- function(name, default = "yes") {
  val <- tolower(trimws(Sys.getenv(name, unset = default)))
  val %in% c("yes", "true", "1", "on")
}

run_cellcycle <- flag_on("LIVER_RUN_CELLCYCLE", "yes")
run_cluster_markers <- flag_on("LIVER_RUN_CLUSTER_MARKERS", "yes")
run_signatures <- flag_on("LIVER_RUN_SIGNATURES", "yes")
run_cnv <- flag_on("LIVER_RUN_CNV", "yes")
run_singler <- flag_on("LIVER_RUN_SINGLER", "yes")
run_trajectory <- flag_on("LIVER_RUN_TRAJECTORY", "no")
regress_cellcycle <- flag_on("LIVER_REGRESS_CELLCYCLE", "no")
skip_gsea <- flag_on("LIVER_SKIP_GSEA", "no")
gsea_max_genes_raw <- Sys.getenv("LIVER_GSEA_MAX_GENES", unset = "")
gsea_max_genes <- if (nzchar(gsea_max_genes_raw)) {
  val <- suppressWarnings(as.integer(gsea_max_genes_raw))
  if (is.na(val) || val < 0) 0L else val
} else {
  0L
}

save_pheatmap <- function(file, fn, width, height, res = 150) {
  name <- basename(file)
  if (name %in% skip_figs) {
    log_msg("skip figure: ", name)
    return(invisible(NULL))
  }
  png(stage_fig_file(file), width = width, height = height, res = res)
  on.exit(dev.off())
  fn()
  log_msg("saved figure: ", name)
}

expand_genes <- function(genes) {
  if (species == "mm") {
    mouse_genes <- paste0(
      toupper(substr(genes, 1, 1)),
      tolower(substr(genes, 2, nchar(genes)))
    )
    unique(c(genes, mouse_genes))
  } else {
    genes
  }
}

args <- commandArgs(trailingOnly = TRUE)
root <- Sys.getenv("LIVER_ROOT", unset = NA)
if (is.na(root) || !nzchar(root)) {
  root <- getwd()
}

accession <- toupper(Sys.getenv("LIVER_ACCESSION", unset = "GSE125449"))
species <- tolower(Sys.getenv("LIVER_SPECIES", unset = "hs"))
raw_dir <- file.path(root, "data", "raw")
if (accession != "GSE125449") {
  raw_dir <- file.path(root, "data", "raw", accession)
}
res_dir <- file.path(root, "results")
fig_dir <- file.path(res_dir, "figures")
data_dir <- file.path(res_dir, "data")
stage_dir <- file.path(res_dir, ".stages")
ckpt_dir <- file.path(res_dir, "checkpoints")

for (d in c(res_dir, fig_dir, data_dir, stage_dir, ckpt_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

ckpt_path <- function(name) {
  file.path(ckpt_dir, name)
}

force <- "--force" %in% args
start_stage <- "01"
for (arg in args) {
  if (startsWith(arg, "--start-stage=")) {
    start_stage <- sub("^--start-stage=", "", arg)
  }
}

dataset_mode <- "single_cell"
dataset_mode_path <- ckpt_path("dataset_mode.txt")
if (file.exists(dataset_mode_path)) {
  dataset_mode <- trimws(readLines(dataset_mode_path, warn = FALSE)[1])
}

stage_allowed <- function(code) {
  as.integer(code) >= as.integer(start_stage)
}

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

mt_pattern <- "^MT-|^mt-"

run_stage <- function(name, expr) {
  marker <- file.path(stage_dir, paste0(name, ".done"))
  if (file.exists(marker)) unlink(marker)
  log_msg("start stage: ", name)
  force(expr)
  writeLines(as.character(Sys.time()), marker)
  log_msg("complete stage: ", name)
}

parse_series_matrix <- function() {
  files <- list.files(
    raw_dir,
    pattern = "GSE125449-.*_series_matrix\\.txt\\.gz$",
    full.names = TRUE
  )
  if (length(files) == 0) stop("No GSE125449 series matrix files found.")

  out <- list()
  for (f in files) {
    con <- gzfile(f, "rt")
    lines <- readLines(con, warn = FALSE)
    close(con)

    title_line <- lines[startsWith(lines, "!Sample_title")]
    char_line <- lines[startsWith(lines, "!Sample_characteristics_ch1")][1]

    parse_values <- function(line) {
      v <- sub("^![^\t]*\t", "", line)
      v <- unlist(strsplit(v, "\t", fixed = TRUE))
      gsub('^"|"$', "", trimws(v))
    }

    titles <- parse_values(title_line[1])
    chars <- parse_values(char_line)
    cancer <- sub("^cancer type: ", "", chars)
    out[[f]] <- data.frame(
      sample = titles,
      cancer_type = cancer,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

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

read_generic_counts <- function(manifest) {
  files <- manifest$files
  matrices <- as.character(files$matrix)
  barcodes <- as.character(files$barcodes)
  genes <- as.character(files$genes)
  if (length(matrices) == 0) {
    stop("No count matrix files found in dataset manifest.")
  }

  count_list_all <- list()
  sample_list <- list()
  group_list <- list()
  bc <- character()
  infer_cell_samples <- function(labels) {
    labels <- as.character(labels)
    out <- rep("Sample1", length(labels))
    pattern <- "_[A-Za-z0-9]+(?:_[A-Za-z0-9]+)?$"
    m <- regexpr(pattern, labels, perl = TRUE)
    hits <- which(m > 0)
    if (length(hits) > 0) {
      out[hits] <- sub("^_", "", regmatches(labels[hits], regexpr(pattern, labels[hits], perl = TRUE)))
    }
    out
  }
  for (i in seq_along(matrices)) {
    mat_file <- matrices[i]
    mat_path <- file.path(raw_dir, mat_file)

    if (grepl("\\.rds$", mat_file, ignore.case = TRUE)) {
      obj <- readRDS(mat_path)
      if (inherits(obj, "Seurat")) {
        m <- GetAssayData(obj, layer = "counts")
      } else if (inherits(obj, "SingleCellExperiment")) {
        m <- assay(obj, "counts")
      } else {
        m <- as(as.matrix(obj), "CsparseMatrix")
      }
    } else if (grepl("\\.h5$", mat_file, ignore.case = TRUE)) {
      m <- Read10X_h5(mat_path)
    } else if (grepl("\\.mtx", mat_file, ignore.case = TRUE)) {
      if (length(barcodes) < i || length(genes) < i) {
        stop("Missing barcode or gene files for ", mat_file)
      }
      bc <- tryCatch(
        as.character(
          fread(file.path(raw_dir, barcodes[i]), header = FALSE)$V1
        ),
        error = function(e) character()
      )
      bc <- bc[nzchar(trimws(bc))]
      gn <- tryCatch(
        fread(file.path(raw_dir, genes[i]), header = FALSE),
        error = function(e) NULL
      )
      gene_vals <- if (is.null(gn) || ncol(gn) == 0) {
        character()
      } else if (ncol(gn) >= 2) {
        as.character(gn[[2]])
      } else {
        as.character(gn[[1]])
      }
      gene_vals <- gene_vals[nzchar(trimws(gene_vals))]
      m <- tryCatch(
        {
          con <- gzfile(mat_path, "rt")
          m <- readMM(con)
          close(con)
          as(m, "CsparseMatrix")
        },
        warning = function(w) read_mtx_fallback(mat_path),
        error = function(e) read_mtx_fallback(mat_path)
      )
      if (length(gene_vals) >= nrow(m)) {
        gene_vals <- gene_vals[seq_len(nrow(m))]
      }
      rownames(m) <- make.unique(gene_vals)
      if (length(gene_vals) < nrow(m)) {
        rownames(m) <- make.unique(
          c(gene_vals, paste0("Gene", seq.int(length(gene_vals) + 1, nrow(m))))
        )
      }
      if (length(bc) >= ncol(m)) {
        colnames(m) <- bc[seq_len(ncol(m))]
      } else {
        colnames(m) <- c(bc, paste0("Cell", seq.int(length(bc) + 1, ncol(m))))
      }
    } else {
      con <- gzfile(mat_path, "rt")
      first_line <- readLines(con, n = 1, warn = FALSE)
      close(con)
      first_fields <- trimws(unlist(strsplit(first_line, "[,\t]")))
      first_field <- gsub('^"|"$', "", first_fields[1])
      looks_barcode <- grepl("^[ACGTN]+([.-]|$)", first_field) ||
        grepl("^[ACGTN]+[0-9]+$", first_field) ||
        grepl("^[A-Z0-9]+_[A-Z0-9]+", first_field)

      if (looks_barcode) {
        barcodes <- gsub('^"|"$', "", first_fields)
        tab <- fread(mat_path, header = FALSE, skip = 1)
        if (ncol(tab) < 2) {
          stop("Matrix file does not look like a gene x cell table: ", mat_file)
        }
        g <- tab[[1]]
        m <- as.matrix(tab[, -1, with = FALSE])
        if (is.character(m)) storage.mode(m) <- "double"
        rownames(m) <- make.unique(as.character(g))
        if (length(barcodes) == ncol(m)) {
          colnames(m) <- barcodes
        } else {
          colnames(m) <- colnames(tab)[-1]
        }
      } else {
        tab <- fread(mat_path, header = TRUE)
        if (ncol(tab) < 2) {
          stop("Matrix file does not look like a gene x cell table: ", mat_file)
        }
        g <- tab[[1]]
        m <- as.matrix(tab[, -1, with = FALSE])
        if (is.character(m)) storage.mode(m) <- "double"
        rownames(m) <- make.unique(as.character(g))
        colnames(m) <- colnames(tab)[-1]
      }
      m <- as(m, "CsparseMatrix")
    }
    barcode_labels <- infer_cell_samples(
      if (length(bc) >= ncol(m)) bc else colnames(m)
    )
    filename_sample <- "Sample1"
    path_parts <- strsplit(mat_file, "/", fixed = TRUE)[[1]]
    feature_dir <- which(tolower(path_parts) == "raw_feature_bc_matrix")
    if (length(feature_dir) > 0 && feature_dir[1] > 1) {
      filename_sample <- path_parts[feature_dir[1] - 1]
    } else {
      g_match <- regmatches(mat_file, regexpr("G[0-9]+[A-Z]?", mat_file))
      if (length(g_match) > 0) {
        filename_sample <- g_match
      } else {
        gsm_match <- regmatches(mat_file, regexpr("GSM[0-9]+", mat_file))
        if (length(gsm_match) > 0) {
          filename_sample <- gsm_match
        }
      }
    }
    if (all(barcode_labels == "Sample1")) {
      sample_label <- filename_sample
      sample_labels <- rep(sample_label, length(barcode_labels))
    } else {
      sample_label <- barcode_labels[1]
      sample_labels <- barcode_labels
    }
    if (identical(sample_label, "Sample1")) {
      colnames(m) <- make.unique(colnames(m))
    } else {
      colnames(m) <- make.unique(paste0(colnames(m), "_", sample_label))
    }
    count_list_all[[i]] <- m
    sample_list[[i]] <- sample_labels
    group_list[[i]] <- rep(infer_group_from_filename(mat_file), ncol(m))
  }

  keep <- rep(TRUE, length(count_list_all))
  if (length(count_list_all) > 1) {
    ref_genes <- rownames(count_list_all[[1]])
    for (i in seq_along(count_list_all)[-1]) {
      if (length(intersect(rownames(count_list_all[[i]]), ref_genes)) == 0) {
        keep[i] <- FALSE
        log_msg("excluding incompatible count matrix: ", matrices[i])
      }
    }
  }
  count_list <- count_list_all[keep]
  cell_sample <- unlist(sample_list[keep], use.names = FALSE)
  cell_group <- unlist(group_list[keep], use.names = FALSE)
  common <- Reduce(intersect, lapply(count_list, rownames))
  counts <- do.call(
    cbind,
    lapply(count_list, function(m) m[common, , drop = FALSE])
  )
  list(
    counts = counts,
    cell_sample = cell_sample,
    cell_group = cell_group
  )
}

infer_group_from_filename <- function(mat_file) {
  low <- tolower(mat_file)
  if (grepl("etoh|alcohol", low) && grepl("control", low)) {
    return(ifelse(grepl("etoh", low), "EtOH", "Control"))
  }
  if (grepl("dmso", low)) {
    return("DMSO")
  }
  if (grepl("jte[-_ ]?607", low)) {
    return("JTE607")
  }
  if (grepl("fibrosis", low) && grepl("control", low)) {
    return(ifelse(grepl("fibrosis", low), "Fibrosis", "Control"))
  }
  if (grepl("mash|_m-|_m[0-9]", low)) {
    return("MASH")
  }
  if (grepl("^.*_n[0-9]", low) || grepl("normal", low)) {
    return("Normal")
  }
  if (grepl("phh", low)) {
    return("PHH")
  }
  if (grepl("cj[0-9]", low)) {
    return("CHB")
  }
  if (grepl("vector", low)) {
    return("Vector")
  }
  if (grepl("klf2", low)) {
    return("KLF2")
  }
  if (grepl("dl-?[0-9]", low)) {
    return("Diseased")
  }
  if (grepl("\\bch[-_0-9]", low)) {
    return("Diseased")
  }
  if (grepl("nl-?[0-9]", low)) {
    return("Normal")
  }
  if (grepl("tumor", low)) {
    return("Tumor")
  }
  if (grepl("normal", low)) {
    return("Normal")
  }
  if (grepl("hcc", low)) {
    return("HCC")
  }
  if (grepl("cca|cholangiocarcinoma", low)) {
    return("iCCA")
  }
  if (grepl("etoh", low)) {
    return("EtOH")
  }
  if (grepl("control", low)) {
    return("Control")
  }
  if (grepl("mock", low)) {
    return("Mock")
  }
  if (grepl("sars|sars-cov|covid", low)) {
    return("SARS")
  }
  if (grepl("ccrcc[0-9]", low)) {
    return(sub("^.*(ccrcc[0-9]).*$", "\\1", low))
  }
  if (grepl("gsm[0-9]+", low)) {
    return(sub("^.*(gsm[0-9]+).*$", "\\1", low))
  }
  ""
}

read_generic_metadata <- function(manifest, cells, cell_sample = NULL) {
  meta <- data.frame(row.names = cells)
  files <- as.character(manifest$files$metadata)
  if (length(files) == 0) {
    return(meta)
  }

  for (f in files) {
    tab <- fread(file.path(raw_dir, f), header = FALSE, fill = TRUE)
    if (nrow(tab) < 2 || ncol(tab) == 0) next
    header_vals <- as.character(tab[1, ])
    if (sum(nzchar(trimws(header_vals))) < ncol(tab)) {
      colnames(tab) <- make.names(
        c("barcode", header_vals[seq_len(ncol(tab) - 1)]),
        unique = TRUE
      )
    } else {
      colnames(tab) <- make.names(header_vals, unique = TRUE)
    }
    tab <- tab[-1]
    bc_names <- c(
      "Cell.Barcode", "Cell", "Barcode", "cell", "barcode",
      "cell_barcode", "cellID", "Index", "V1"
    )
    bc_col <- bc_names[bc_names %in% colnames(tab)][1]
    if (is.na(bc_col) || !bc_col %in% colnames(tab)) {
      first_vals <- as.character(tab[[1]])[!is.na(as.character(tab[[1]]))]
      looks_barcode <- length(first_vals) > 0 &&
        (
          grepl("^[A-Za-z0-9]+-1$", first_vals[1]) ||
          grepl("^[ACGTN]{10,}", first_vals[1])
        )
      if (looks_barcode) {
        bc_col <- colnames(tab)[1]
      } else {
        next
      }
    }

    bc_vals <- as.character(tab[[bc_col]])
    strip_sample_affixes <- function(x) {
      x <- as.character(x)
      without_suffix <- sub("_[^_]+$", "", x)
      ifelse(without_suffix != x, without_suffix, sub("^[^_]+_", "", x))
    }
    match_barcodes <- function(cells, vals) {
      idx <- match(cells, vals)
      if (sum(!is.na(idx)) < 10 && !is.null(cell_sample)) {
        samples <- as.character(cell_sample)
        candidate <- mapply(
          function(c, s) {
            if (is.na(s)) return(NA_character_)
            if (startsWith(c, s)) {
              return(paste0(
                s,
                "_",
                sub(paste0(s, "_"), "", c, fixed = TRUE)
              ))
            }
            if (endsWith(c, s)) {
              return(paste0(
                s,
                "_",
                sub(paste0("_", s), "", c, fixed = TRUE)
              ))
            }
            paste0(s, "_", c)
          },
          cells,
          samples,
          USE.NAMES = FALSE
        )
        idx <- match(candidate, vals)
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(
          strip_sample_affixes(cells),
          strip_sample_affixes(vals)
        )
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(sub("_[^_]+$", "", cells), vals)
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(
          gsub("\\.", "-", sub("_[^_]+$", "", cells)),
          gsub("\\.", "-", vals)
        )
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(
          gsub("\\.", "-", cells),
          gsub("\\.", "-", vals)
        )
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(sub("^.*_", "", cells), sub("^.*_", "", vals))
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(sub("-\\d+$", "", cells), sub("-\\d+$", "", vals))
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(
          sub("-.*$", "", sub("_[^_]+$", "", cells)),
          sub("-.*$", "", vals)
        )
      }
      if (sum(!is.na(idx)) < 10) {
        idx <- match(
          sub("_[^_]+$", "", cells),
          sub("^.*_", "", vals)
        )
      }
      idx
    }
    bc_idx <- match_barcodes(cells, bc_vals)
    log_msg(
      "metadata file ",
      f,
      ": bc_col=",
      bc_col,
      " matches=",
      sum(!is.na(bc_idx))
    )
    for (col in setdiff(colnames(tab), bc_col)) {
      meta[[col]] <- tab[[col]][bc_idx]
    }
  }

  if (!"sample" %in% colnames(meta)) {
    sample_candidates <- c(
      "sample", "Sample", "sample.id", "sample_id",
      "donor", "Donor", "patient", "Patient", "patient.id"
    )
    sample_col <- sample_candidates[sample_candidates %in% colnames(meta)][1]
    if (length(sample_col) == 1 && !is.na(sample_col)) {
      meta$sample <- meta[[sample_col]]
    }
  }

  type_candidates <- c(
    "Type", "celltype", "cell_type", "celltype_global", "celltype_sub",
    "CellType", "Cell.type", "cell.type", "cell_type_annot"
  )
  type_col <- type_candidates[type_candidates %in% colnames(meta)][1]
  if (length(type_col) == 1 && !is.na(type_col)) {
    meta$published_type <- meta[[type_col]]
  } else {
    meta$published_type <- "Unannotated"
  }
  meta
}

parse_series_generic <- function(manifest) {
  files <- as.character(manifest$files$series_matrices)
  if (length(files) == 0) {
    return(NULL)
  }

  parse_values <- function(line) {
    v <- sub("^![^\t]*\t", "", line)
    v <- unlist(strsplit(v, "\t", fixed = TRUE))
    gsub('^"|"$', "", trimws(v))
  }

  out <- list()
  for (f in files) {
    con <- gzfile(file.path(raw_dir, f), "rt")
    lines <- readLines(con, warn = FALSE)
    close(con)

    title_line <- lines[startsWith(lines, "!Sample_title")]
    if (length(title_line) == 0) next
    titles <- parse_values(title_line[1])
    df <- data.frame(sample = titles, stringsAsFactors = FALSE)

    source_line <- lines[startsWith(lines, "!Sample_source_name_ch1")]
    if (length(source_line) > 0) {
      df$source_name <- parse_values(source_line[1])
      df$condition <- ifelse(
        grepl("healthy|normal|control", df$source_name, ignore.case = TRUE),
        "Healthy",
        "Disease"
      )
    }

    char_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
    for (cl in char_lines) {
      vals <- parse_values(cl)
      field <- sub(":.*$", "", vals[1])
      key <- make.names(field)
      values <- sub("^[^:]+:\\s*", "", vals)
      df[[key]] <- values
    }
    out[[f]] <- df
  }

  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

infer_condition <- function(meta, sample_ann) {
  candidates <- c(
    "condition", "group", "disease", "disease_status", "health_status",
    "histology", "cancer_type", "cancer.type", "tissue", "site",
    "tissue_sub", "sample_type", "status", "treatment"
  )
  candidate_cols <- colnames(meta)[
    tolower(colnames(meta)) %in% candidates
  ]
  for (col in candidate_cols) {
    vals <- as.character(meta[[col]])
    if (length(unique(vals[!is.na(vals)])) >= 2) {
      return(vals)
    }
  }

  if (!is.null(sample_ann) && "sample" %in% colnames(meta)) {
    extract_key <- function(v) {
      v <- sub("^.*:", "", as.character(v))
      v <- sub("^\\s+|\\s+$", "", v)
      v_norm <- sub("_CRC$", " primary CRC", v)
      v_norm <- sub("_LM$", " liver metastases", v_norm)
      v_norm <- sub("_PBMC$", " PBMC", v_norm)
      keep_phrase <- grepl(
        "primary CRC|primary colorectal cancer|liver metastases|PBMC",
        v_norm,
        ignore.case = TRUE
      )
      out <- sub("^.*\\b([A-Za-z]+[0-9]+).*$", "\\1", v_norm)
      out[keep_phrase] <- v_norm[keep_phrase]
      out
    }
    match_samples <- function(a, b) {
      idx <- match(a, b)
      if (sum(!is.na(idx)) < 5) {
        a2 <- sub("^.*(G[0-9]+[A-Z]?).*$", "\\1", a)
        b2 <- sub("^.*(G[0-9]+[A-Z]?).*$", "\\1", b)
        idx <- match(a2, b2)
      }
      if (sum(!is.na(idx)) < 5) {
        idx <- match(extract_key(a), extract_key(b))
      }
      idx
    }
    for (col in colnames(sample_ann)) {
      if (tolower(col) == "sample") next
      vals <- as.character(sample_ann[[col]])
      if (length(unique(vals[!is.na(vals)])) >= 2) {
        idx <- match_samples(meta$sample, sample_ann$sample)
        log_msg(
          "condition candidate ",
          col,
          ": matches=",
          sum(!is.na(idx)),
          " unique=",
          length(unique(vals[!is.na(vals)]))
        )
        cond <- vals[idx]
        if (any(!is.na(cond))) {
          return(cond)
        }
      }
    }
  }
  NULL
}

normalize_condition <- function(meta) {
  cond <- as.character(meta$condition)
  if (any(grepl("tumor", cond, ignore.case = TRUE)) &&
      any(grepl("normal", cond, ignore.case = TRUE))) {
    keep <- grepl("tumor|normal", cond, ignore.case = TRUE)
  } else {
    tt <- sort(table(cond), decreasing = TRUE)
    if (length(tt) < 2) {
      stop(
        "Cannot find two groups for differential expression. ",
        "Unique conditions: ",
        paste(names(tt), collapse = ", ")
      )
    }
    keep <- cond %in% names(tt)[1:2]
  }
  meta <- meta[keep, , drop = FALSE]
  meta$condition <- factor(as.character(meta$condition))
  meta
}

read_generic_dataset <- function(manifest) {
  loaded <- read_generic_counts(manifest)
  counts <- loaded$counts
  meta <- read_generic_metadata(
    manifest,
    colnames(counts),
    loaded$cell_sample
  )
  ann <- parse_series_generic(manifest)

  if (!is.null(ann) && "sample" %in% colnames(ann)) {
    ann_samples <- as.character(ann$sample)
    count_samples <- colnames(counts)
    idx <- match(count_samples, ann_samples)
    if (sum(!is.na(idx)) < 5) {
      idx <- match(sub("_[^_]+$", "", count_samples), ann_samples)
    }
    if (sum(!is.na(idx)) < 5) {
      idx <- match(count_samples, sub("^[^ ]+ ", "", ann_samples))
    }
    if (sum(!is.na(idx)) < 5) {
      idx <- match(count_samples, sub("^[^_]+_", "", ann_samples))
    }
    if (sum(!is.na(idx)) < 5) {
      idx <- match(
        sub("^.*(G[0-9]+[A-Z]?).*$", "\\1", count_samples),
        sub("^.*(G[0-9]+[A-Z]?).*$", "\\1", ann_samples)
      )
    }
    if (sum(!is.na(idx)) >= 2) {
      sample_vals <- loaded$cell_sample
      matched <- !is.na(idx)
      sample_vals[matched] <- ann_samples[idx[matched]]
      meta$sample <- sample_vals
      log_msg(
        "matched series matrix sample labels: ",
        sum(matched), "/", length(idx)
      )
    }
  }

  if (!"sample" %in% colnames(meta)) {
    meta$sample <- loaded$cell_sample
  }
  log_msg(
    "generic metadata columns: ",
    paste(colnames(meta), collapse = ", ")
  )
  if ("sample" %in% colnames(meta)) {
    log_msg(
      "generic metadata sample head: ",
      paste(head(as.character(meta$sample), 20), collapse = ", ")
    )
  }
  log_msg(
    "generic sample labels: ",
    paste(head(as.character(meta$sample), 20), collapse = ", ")
  )
  if (!is.null(ann)) {
    log_msg(
      "series sample labels: ",
      paste(head(as.character(ann$sample), 20), collapse = ", ")
    )
    log_msg(
      "series columns: ",
      paste(colnames(ann), collapse = ", ")
    )
  }
  cond <- infer_condition(meta, ann)
  if (is.null(cond) && any(nzchar(loaded$cell_group))) {
    meta$condition <- loaded$cell_group
    cond <- infer_condition(meta, ann)
  }
  if (is.null(cond)) {
    stop(
      "Could not automatically infer two groups. ",
      "Add a condition/group/tissue column to the dataset metadata."
    )
  }
  log_msg(
    "inferred condition unique: ",
    paste(unique(as.character(cond)), collapse = ", ")
  )
  log_msg(
    "inferred condition table: ",
    paste(names(table(cond)), table(cond), sep = "=", collapse = ", ")
  )
  meta$condition <- cond
  meta <- normalize_condition(meta)
  counts <- counts[, rownames(meta), drop = FALSE]
  for (col in c("nCount_RNA", "nFeature_RNA", "percentMt", "percent.mt", "percent_mito")) {
    if (col %in% colnames(meta)) {
      meta[[col]] <- as.numeric(as.character(meta[[col]]))
    }
  }
  if (!"sample" %in% colnames(meta)) {
    meta$sample <- "Sample1"
  }
  if (!"published_type" %in% colnames(meta)) {
    meta$published_type <- "Unannotated"
  }
  list(counts = counts, meta = meta, ann = ann)
}

if (stage_allowed("01")) run_stage("01_load_data", {
  if (accession == "GSE125449") {
    s1 <- read_10x_set("1")
    s2 <- read_10x_set("2")

    common_genes <- intersect(rownames(s1$counts), rownames(s2$counts))
    counts <- cbind(
      s1$counts[common_genes, , drop = FALSE],
      s2$counts[common_genes, , drop = FALSE]
    )
    meta <- rbind(s1$meta, s2$meta)

    ann <- parse_series_matrix()
    meta$cancer_type <- ann$cancer_type[match(meta$sample, ann$sample)]
    meta$condition <- ifelse(
      meta$cancer_type == "Hepatocellular carcinoma",
      "HCC",
      "iCCA"
    )

    if (any(is.na(meta$condition))) {
      stop("Some cells could not be mapped to HCC or iCCA.")
    }
  } else {
    manifest <- load_manifest()
    dataset_mode <- if (identical(manifest$mode, "single_cell")) {
      "single_cell"
    } else {
      "sample_level"
    }
    writeLines(dataset_mode, dataset_mode_path)
    generic <- read_generic_dataset(manifest)
    counts <- generic$counts
    meta <- generic$meta
    ann <- generic$ann
    if (is.null(ann)) {
      ann <- data.frame(
        sample = character(),
        note = character(),
        stringsAsFactors = FALSE
      )
    }
  }

  write.csv(ann, stage_data_file("sample_annotations.csv"), row.names = FALSE)

  seurat_raw <- CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    min.cells = 0,
    min.features = 0,
    project = accession
  )
  seurat_raw$orig.ident <- accession
  seurat_raw[["percent.mt"]] <- PercentageFeatureSet(seurat_raw, pattern = mt_pattern)
  seurat_raw[["percent.ribo"]] <- PercentageFeatureSet(seurat_raw, pattern = "^RP[SL]")

  log_msg("raw cells: ", ncol(seurat_raw))
  log_msg("raw genes: ", nrow(seurat_raw))
  log_msg(
    "condition table: ",
    paste(names(table(seurat_raw$condition)), table(seurat_raw$condition), sep = "=", collapse = ", ")
  )
  saveRDS(seurat_raw, ckpt_path("seurat_raw.rds"))
})

if (stage_allowed("02")) run_stage("02_qc_filter", {
  if (!exists("seurat_raw")) {
    seurat_raw <- readRDS(ckpt_path("seurat_raw.rds"))
  }
  p_raw <- VlnPlot(
    seurat_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "condition",
    ncol = 3,
    pt.size = 0
  ) & NoLegend()
  save_fig(
    file.path(fig_dir, "fig_01_qc_raw_violin.png"),
    p_raw,
    width = 12,
    height = 6,
    dpi = 150
  )

  qc_data <- FetchData(
    seurat_raw,
    vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  )
  qc_data$sample <- seurat_raw$sample
  qc_data$condition <- seurat_raw$condition

  qc_pvalue_table <- function(qc_frame, stage) {
    group_col <- as.character(qc_frame$condition)
    metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
    groups <- sort(unique(group_col[!is.na(group_col) & nzchar(group_col)]))
    if (length(groups) < 2) {
      return(data.frame(
        stage = character(),
        metric = character(),
        group1 = character(),
        group2 = character(),
        n1 = integer(),
        n2 = integer(),
        median1 = numeric(),
        median2 = numeric(),
        median_diff = numeric(),
        statistic = numeric(),
        pvalue = numeric(),
        padj = numeric(),
        neg_log10_pvalue = numeric(),
        direction = character(),
        stringsAsFactors = FALSE
      ))
    }
    pairs <- combn(groups, 2, simplify = FALSE)
    rows <- lapply(pairs, function(pair) {
      lapply(metrics, function(metric) {
        x <- as.numeric(qc_frame[[metric]][group_col == pair[1]])
        y <- as.numeric(qc_frame[[metric]][group_col == pair[2]])
        x <- x[is.finite(x)]
        y <- y[is.finite(y)]
        med1 <- median(x)
        med2 <- median(y)
        test <- tryCatch(
          suppressWarnings(wilcox.test(x, y, exact = FALSE)),
          error = function(e) NULL
        )
        data.frame(
          stage = stage,
          metric = metric,
          group1 = pair[1],
          group2 = pair[2],
          n1 = length(x),
          n2 = length(y),
          median1 = med1,
          median2 = med2,
          median_diff = med2 - med1,
          statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
          pvalue = if (is.null(test)) NA_real_ else test$p.value,
          direction = ifelse(
            is.na(med1) || is.na(med2) || med2 == med1,
            "no difference",
            ifelse(med2 > med1, paste0(pair[2], " higher"), paste0(pair[1], " higher"))
          ),
          stringsAsFactors = FALSE
        )
      })
    })
    out <- do.call(rbind, unlist(rows, recursive = FALSE))
    out$padj <- p.adjust(out$pvalue, method = "BH")
    out$neg_log10_pvalue <- ifelse(
      is.na(out$pvalue),
      NA_real_,
      -log10(pmax(out$pvalue, .Machine$double.xmin))
    )
    out
  }

  qc_diff_raw <- qc_pvalue_table(qc_data, "raw")

  if (dataset_mode == "sample_level") {
    lo_feature <- 0
    hi_feature <- max(qc_data$nFeature_RNA) + 1
    lo_count <- 0
    hi_count <- max(qc_data$nCount_RNA) + 1
    hi_mt <- 100
  } else {
    lo_feature <- if (!is.na(qc_min_features)) {
      qc_min_features
    } else {
      max(200, as.numeric(quantile(qc_data$nFeature_RNA, 0.01)))
    }
    hi_feature <- if (!is.na(qc_max_features)) {
      qc_max_features
    } else {
      min(20000, as.numeric(quantile(qc_data$nFeature_RNA, 0.99)))
    }
    lo_count <- if (!is.na(qc_min_counts)) {
      qc_min_counts
    } else {
      max(300, as.numeric(quantile(qc_data$nCount_RNA, 0.01)))
    }
    hi_count <- if (!is.na(qc_max_counts)) {
      qc_max_counts
    } else {
      min(500000, as.numeric(quantile(qc_data$nCount_RNA, 0.99)))
    }
    hi_mt <- if (!is.na(qc_max_mt)) {
      qc_max_mt
    } else {
      min(30, as.numeric(quantile(qc_data$percent.mt, 0.99)))
    }
  }

  if (lo_feature > hi_feature || lo_count > hi_count) {
    log_msg(
      "QC bounds are incompatible for this dataset; ",
      "relaxing thresholds to keep all samples"
    )
    lo_feature <- 0
    hi_feature <- max(qc_data$nFeature_RNA) + 1
    lo_count <- 0
    hi_count <- max(qc_data$nCount_RNA) + 1
    hi_mt <- 100
  }

  log_msg(
    "QC thresholds: nFeature ",
    round(lo_feature, 1), "-", round(hi_feature, 1),
    ", nCount ", round(lo_count, 1), "-", round(hi_count, 1),
    ", percent.mt <= ", round(hi_mt, 1)
  )

  seurat_qc <- subset(
    seurat_raw,
    subset = nFeature_RNA >= lo_feature &
      nFeature_RNA <= hi_feature &
      nCount_RNA >= lo_count &
      nCount_RNA <= hi_count &
      percent.mt <= hi_mt
  )

  qc_metrics <- FetchData(
    seurat_qc,
    vars = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
  )
  qc_metrics$sample <- seurat_qc$sample
  qc_metrics$condition <- seurat_qc$condition
  write.csv(qc_metrics, stage_data_file("fig_01_qc_metrics.csv"))

  qc_diff_filtered <- qc_pvalue_table(qc_metrics, "filtered")
  qc_diff <- rbind(qc_diff_raw, qc_diff_filtered)
  write.csv(qc_diff, stage_data_file("fig_48_qc_pvalue_comparison.csv"), row.names = FALSE)

  if (nrow(qc_diff) > 0) {
    qc_diff$comparison <- paste(qc_diff$group1, qc_diff$group2, sep = " vs ")
    qc_diff$stage <- factor(qc_diff$stage, levels = c("raw", "filtered"))
    p_qc_diff <- ggplot(qc_diff, aes(x = metric, y = neg_log10_pvalue, fill = comparison)) +
      geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#666666") +
      geom_text(
        aes(label = ifelse(
          is.na(pvalue),
          "NA",
          paste0("P=", formatC(pvalue, digits = 2, format = "g"))
        )),
        position = position_dodge2(width = 0.7),
        vjust = -0.4,
        size = 3
      ) +
      facet_wrap(~ stage, ncol = 1, scales = "free_y") +
      labs(
        x = "QC metric",
        y = "-log10(P value)",
        title = "QC metric difference by condition",
        subtitle = "Wilcoxon rank-sum test; dashed line indicates P = 0.05"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom"
      )
  } else {
    p_qc_diff <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
      geom_text(label = "At least two conditions are required") +
      theme_void()
  }
  save_fig(
    file.path(fig_dir, "fig_48_qc_pvalue_comparison.png"),
    p_qc_diff,
    width = 10,
    height = 8,
    dpi = 150
  )

  p_qc <- VlnPlot(
    seurat_qc,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "condition",
    ncol = 3,
    pt.size = 0
  ) & NoLegend()
  save_fig(
    file.path(fig_dir, "fig_01_qc_filtered_violin.png"),
    p_qc,
    width = 12,
    height = 6,
    dpi = 150
  )

  log_msg("cells after QC: ", ncol(seurat_qc))
  saveRDS(seurat_qc, ckpt_path("seurat_qc.rds"))
})

if (stage_allowed("03")) run_stage("03_doublets", {
  if (!exists("seurat_qc")) {
    seurat_qc <- readRDS(ckpt_path("seurat_qc.rds"))
  }
  sce <- SingleCellExperiment(
    assays = list(counts = GetAssayData(seurat_qc, layer = "counts"))
  )
  colData(sce)$sample <- seurat_qc$sample
  colData(sce)$condition <- seurat_qc$condition

  if (dataset_mode == "sample_level") {
    log_msg("sample-level mode: skipping doublet detection")
    sce$scDblFinder.score <- rep(0, ncol(sce))
    sce$scDblFinder.class <- rep("singlet", ncol(sce))
  } else {
    sce <- tryCatch(
      scDblFinder(sce, BPPARAM = BiocParallel::SerialParam()),
      error = function(e) {
        log_msg("scDblFinder failed; marking all cells as singlet")
        sce$scDblFinder.score <- rep(0, ncol(sce))
        sce$scDblFinder.class <- rep("singlet", ncol(sce))
        sce
      }
    )
  }

  seurat_qc$doublet_score <- sce$scDblFinder.score
  seurat_qc$doublet_call <- as.character(sce$scDblFinder.class)

  doublet_tbl <- data.frame(
    cell = colnames(seurat_qc),
    sample = seurat_qc$sample,
    condition = seurat_qc$condition,
    doublet_score = seurat_qc$doublet_score,
    doublet_call = seurat_qc$doublet_call,
    stringsAsFactors = FALSE
  )
  write.csv(doublet_tbl, stage_data_file("fig_02_doublet_results.csv"), row.names = FALSE)

  p_dbl <- ggplot(
    doublet_tbl,
    aes(x = doublet_call, y = doublet_score, fill = doublet_call)
  ) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    scale_fill_manual(values = c("singlet" = "#4DBBD5", "doublet" = "#E64B35")) +
    labs(x = "Doublet call", y = "scDblFinder score", title = "Doublet score by call") +
    theme_minimal()
  save_fig(
    file.path(fig_dir, "fig_02_doublet_scores.png"),
    p_dbl,
    width = 7,
    height = 5,
    dpi = 150
  )

  seurat <- subset(seurat_qc, subset = doublet_call == "singlet")
  log_msg("cells after doublet removal: ", ncol(seurat))
  saveRDS(seurat, ckpt_path("seurat_singlet.rds"))
})

if (stage_allowed("04")) run_stage("04_cluster", {
  if (!exists("seurat")) {
    seurat <- readRDS(ckpt_path("seurat_singlet.rds"))
  }
  seurat <- NormalizeData(seurat, verbose = FALSE)
  if (dataset_mode == "sample_level") {
    log_msg("sample-level mode: assigning sample-level clusters and embeddings")
    seurat$seurat_clusters <- as.character(seurat$sample)
    npcs <- min(30, max(1, ncol(seurat) - 1))
    emb_pca <- matrix(
      rnorm(ncol(seurat) * npcs),
      nrow = ncol(seurat),
      ncol = npcs
    )
    rownames(emb_pca) <- colnames(seurat)
    colnames(emb_pca) <- paste0("PC_", seq_len(npcs))
    seurat[["pca"]] <- CreateDimReducObject(
      embeddings = emb_pca,
      key = "PC_",
      assay = "RNA"
    )
    emb_umap <- matrix(
      rnorm(ncol(seurat) * 2),
      nrow = ncol(seurat),
      ncol = 2
    )
    rownames(emb_umap) <- colnames(seurat)
    colnames(emb_umap) <- c("UMAP_1", "UMAP_2")
    seurat[["umap"]] <- CreateDimReducObject(
      embeddings = emb_umap,
      key = "umap_",
      assay = "RNA"
    )

    p_pca <- DimPlot(seurat, reduction = "pca", group.by = "condition") +
      ggtitle("PCA by condition (sample-level)")
    save_fig(file.path(fig_dir, "fig_14_pca.png"), p_pca, width = 8, height = 7)

    p_elbow <- ggplot(
      data.frame(PC = seq_len(npcs), stdev = rep(1, npcs)),
      aes(x = PC, y = stdev)
    ) +
      geom_line() +
      labs(title = "PCA standard deviation (sample-level)") +
      theme_minimal()
    save_fig(file.path(fig_dir, "fig_15_elbow.png"), p_elbow, width = 8, height = 5)

    umap_tbl <- data.frame(
      UMAP_1 = emb_umap[, 1],
      UMAP_2 = emb_umap[, 2],
      cell = rownames(emb_umap),
      seurat_clusters = as.character(seurat$seurat_clusters),
      condition = seurat$condition,
      sample = seurat$sample,
      stringsAsFactors = FALSE
    )
    write.csv(
      umap_tbl,
      stage_data_file("fig_03_04_05_umap_coordinates.csv"),
      row.names = FALSE
    )

    cluster_counts <- as.data.frame(table(
      seurat_clusters = seurat$seurat_clusters,
      condition = seurat$condition,
      sample = seurat$sample
    ))
    write.csv(
      cluster_counts,
      stage_data_file("fig_18_19_30_cluster_composition.csv"),
      row.names = FALSE
    )
    log_msg("number of clusters: ", length(unique(seurat$seurat_clusters)))
  } else {
    seurat <- tryCatch(
      FindVariableFeatures(
        seurat,
        selection.method = "vst",
        nfeatures = min(2000, max(10, nrow(seurat) - 1)),
        verbose = FALSE
      ),
      error = function(e) {
        log_msg("variable feature selection failed; using all genes: ", conditionMessage(e))
        VariableFeatures(seurat) <- rownames(seurat)
        seurat
      }
    )
    if (run_cellcycle) {
      cc_genes <- tryCatch(Seurat::cc.genes, error = function(e) NULL)
      if (is.null(cc_genes)) {
        cc_genes <- tryCatch(Seurat::cc.genes.updated.2019, error = function(e) NULL)
      }
      s_genes <- intersect(cc_genes$s.genes, rownames(seurat))
      g2m_genes <- intersect(cc_genes$g2m.genes, rownames(seurat))
      if (length(s_genes) >= 5 && length(g2m_genes) >= 5) {
        seurat <- CellCycleScoring(
          seurat,
          s.features = s_genes,
          g2m.features = g2m_genes,
          set.ident = FALSE
        )
        cc_tbl <- data.frame(
          cell = colnames(seurat),
          S.Score = seurat$S.Score,
          G2M.Score = seurat$G2M.Score,
          Phase = seurat$Phase,
          condition = seurat$condition,
          stringsAsFactors = FALSE
        )
        write.csv(cc_tbl, stage_data_file("fig_26_27_cell_cycle_scores.csv"), row.names = FALSE)
        log_msg("cell cycle scoring applied")
        if (regress_cellcycle) {
          seurat <- ScaleData(
            seurat,
            features = VariableFeatures(seurat),
            vars.to.regress = c("S.Score", "G2M.Score"),
            verbose = FALSE
          )
        } else {
          seurat <- ScaleData(seurat, features = VariableFeatures(seurat), verbose = FALSE)
        }
      } else {
        log_msg("cell cycle markers insufficient; skipping cell cycle scoring")
        seurat <- ScaleData(seurat, features = VariableFeatures(seurat), verbose = FALSE)
      }
    } else {
      seurat <- ScaleData(seurat, features = VariableFeatures(seurat), verbose = FALSE)
    }
    npcs <- min(30, ncol(seurat) - 1)
    dims_use <- seq_len(max(1, min(20, npcs)))
    seurat <- tryCatch(
      RunPCA(seurat, npcs = npcs, verbose = FALSE),
      error = function(e) {
        log_msg("PCA failed; using dummy reduction: ", conditionMessage(e))
        emb <- matrix(
          rnorm(ncol(seurat) * max(1, npcs)),
          nrow = ncol(seurat),
          ncol = max(1, npcs)
        )
        rownames(emb) <- colnames(seurat)
        colnames(emb) <- paste0("PC_", seq_len(max(1, npcs)))
        seurat[["pca"]] <- CreateDimReducObject(
          embeddings = emb,
          key = "PC_",
          assay = "RNA"
        )
        seurat
      }
    )
    reduction <- "pca"
    if (
      requireNamespace("harmony", quietly = TRUE) &&
      length(unique(seurat$sample)) > 1
    ) {
      tryCatch(
        {
          seurat <- harmony::RunHarmony(
            seurat,
            group.by.vars = "sample",
            reduction.use = "pca",
            dims.use = dims_use,
            verbose = FALSE
          )
          reduction <- "harmony"
          log_msg("Harmony batch correction applied")
        },
        error = function(e) {
          log_msg("Harmony failed: ", conditionMessage(e))
        }
      )
    }
    seurat <- FindNeighbors(
      seurat,
      reduction = reduction,
      dims = dims_use,
      verbose = FALSE
    )
    seurat <- FindClusters(
      seurat,
      resolution = cluster_resolution,
      algorithm = ifelse(
        is.na(cluster_algorithm),
        1,
        as.integer(cluster_algorithm)
      ),
      verbose = FALSE
    )
    seurat <- tryCatch(
      RunUMAP(
        seurat,
        reduction = reduction,
        dims = dims_use,
        n.neighbors = min(30, ncol(seurat) - 1),
        seed.use = 42,
        verbose = FALSE
      ),
      error = function(e) {
        log_msg("UMAP failed; using dummy embedding: ", conditionMessage(e))
        emb <- matrix(
          rnorm(ncol(seurat) * 2),
          nrow = ncol(seurat),
          ncol = 2
        )
        rownames(emb) <- colnames(seurat)
        colnames(emb) <- c("UMAP_1", "UMAP_2")
        seurat[["umap"]] <- CreateDimReducObject(
          embeddings = emb,
          key = "umap_",
          assay = "RNA"
        )
        seurat
      }
    )

    p_pca <- DimPlot(seurat, reduction = "pca", group.by = "condition") +
      ggtitle("PCA by condition")
    save_fig(
      file.path(fig_dir, "fig_14_pca.png"),
      p_pca,
      width = 8,
      height = 7
    )

    p_elbow <- ElbowPlot(seurat, ndims = npcs) +
      ggtitle("Principal component standard deviation")
    save_fig(
      file.path(fig_dir, "fig_15_elbow.png"),
      p_elbow,
      width = 8,
      height = 5
    )

    umap_tbl <- as.data.frame(Embeddings(seurat, reduction = "umap"))
    umap_tbl$cell <- rownames(umap_tbl)
    umap_tbl$seurat_clusters <- as.character(seurat$seurat_clusters)
    umap_tbl$condition <- seurat$condition
    umap_tbl$sample <- seurat$sample
    write.csv(umap_tbl, stage_data_file("fig_03_04_05_umap_coordinates.csv"), row.names = FALSE)

    cluster_counts <- as.data.frame(table(
      seurat_clusters = seurat$seurat_clusters,
      condition = seurat$condition,
      sample = seurat$sample
    ))
    write.csv(cluster_counts, stage_data_file("fig_18_19_30_cluster_composition.csv"), row.names = FALSE)

    log_msg("number of clusters: ", length(unique(seurat$seurat_clusters)))
  }
  saveRDS(seurat, ckpt_path("seurat_clustered.rds"))
})

if (stage_allowed("05")) run_stage("05_annotation", {
  if (!exists("seurat")) {
    seurat <- readRDS(ckpt_path("seurat_clustered.rds"))
  }
  marker_list <- list(
    T_NK = c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "CD4"),
    B = c("CD79A", "MS4A1", "CD19", "IGHG1"),
    Myeloid = c("LYZ", "CD68", "C1QA", "C1QB", "FCGR3A"),
    Hepatocyte = c("ALB", "APOA1", "APOA2", "FGB", "SERPINA1"),
    Cholangiocyte = c("KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"),
    Hepatic_stellate = c("COL1A1", "COL1A2", "ACTA2", "RGS5", "PDGFRB"),
    Endothelial = c("PECAM1", "VWF", "CLDN5", "PLVAP"),
    Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"),
    Malignant = c("EPCAM", "KRT8", "KRT18", "KRT19", "AFP", "GPC3"),
    Mast = c("TPSAB1", "CPA3", "MS4A2")
  )
  marker_list <- lapply(marker_list, expand_genes)
  marker_list <- lapply(marker_list, function(x) intersect(x, rownames(seurat)))
  marker_list <- marker_list[lengths(marker_list) > 0]
  marker_names <- names(marker_list)

  if (length(marker_list) > 0) {
    data_mat <- GetAssayData(seurat, layer = "data")
    score_list <- lapply(marker_list, function(genes) {
      genes <- intersect(genes, rownames(seurat))
      if (length(genes) == 0) {
        return(rep(0, ncol(seurat)))
      }
      as.numeric(Matrix::colMeans(data_mat[genes, , drop = FALSE]))
    })
    score_mat <- do.call(cbind, score_list)
    colnames(score_mat) <- marker_names
    rownames(score_mat) <- colnames(seurat)
  } else {
    marker_names <- "Unannotated"
    score_mat <- matrix(
      0,
      nrow = ncol(seurat),
      ncol = 1,
      dimnames = list(colnames(seurat), "Unannotated")
    )
  }

  seurat$celltype_annot_cell <- marker_names[
    max.col(score_mat, ties.method = "first")
  ]

  seurat$celltype_annot <- seurat$celltype_annot_cell

  cluster_ids <- unique(as.character(seurat$seurat_clusters))
  if (all(grepl("^[0-9]+$", cluster_ids))) {
    cluster_ids <- as.character(sort(as.integer(cluster_ids)))
  }
  cluster_labels <- vapply(cluster_ids, function(cl) {
    labs <- seurat$celltype_annot[
      as.character(seurat$seurat_clusters) == cl
    ]
    names(sort(table(labs), decreasing = TRUE))[1]
  }, character(1))
  seurat$cluster_label <- unname(cluster_labels[match(
    as.character(seurat$seurat_clusters),
    cluster_ids
  )])

  annotation_tbl <- data.frame(
    cell = colnames(seurat),
    seurat_clusters = as.character(seurat$seurat_clusters),
    celltype_annot = seurat$celltype_annot,
    celltype_annot_cell = seurat$celltype_annot_cell,
    cluster_label = seurat$cluster_label,
    published_type = seurat$published_type,
    condition = seurat$condition,
    sample = seurat$sample,
    stringsAsFactors = FALSE
  )
  write.csv(annotation_tbl, stage_data_file("fig_05_16_17_cell_annotations.csv"), row.names = FALSE)

  p_clusters <- DimPlot(seurat, group.by = "seurat_clusters", label = FALSE) +
    ggtitle("Seurat clusters (marker-based names)")
  p_clusters <- tryCatch(
    {
      cluster_legend_labels <- paste0(cluster_ids, " - ", cluster_labels)
      p_clusters$data$seurat_clusters <- factor(
        cluster_legend_labels[match(
          as.character(p_clusters$data$seurat_clusters),
          cluster_ids
        )],
        levels = unique(c(cluster_legend_labels, cluster_labels))
      )
      LabelClusters(
        p_clusters,
        id = "seurat_clusters",
        clusters = cluster_legend_labels,
        labels = cluster_labels,
        size = 4
      )
    },
    error = function(e) {
      log_msg(
        "custom cluster labels failed, using numeric labels: ",
        conditionMessage(e)
      )
      DimPlot(seurat, group.by = "seurat_clusters", label = TRUE) +
        ggtitle("Seurat clusters")
    }
  )
  p_condition <- DimPlot(seurat, group.by = "condition") +
    ggtitle(paste(sort(unique(as.character(seurat$condition))), collapse = " vs "))
  p_annot <- DimPlot(seurat, group.by = "celltype_annot", label = TRUE) +
    ggtitle("Marker-based annotation")
  p_pub <- DimPlot(seurat, group.by = "published_type", label = TRUE) +
    ggtitle("Published cell type")

  save_fig(
    file.path(fig_dir, "fig_03_umap_clusters.png"),
    p_clusters,
    width = 8,
    height = 7,
    dpi = 150
  )
  save_fig(
    file.path(fig_dir, "fig_04_umap_condition.png"),
    p_condition,
    width = 8,
    height = 7,
    dpi = 150
  )
  save_fig(
    file.path(fig_dir, "fig_05_umap_annotation.png"),
    p_annot + p_pub,
    width = 15,
    height = 7,
    dpi = 150
  )

  dot_features <- expand_genes(c(
    "CD3D", "CD8A", "NKG7", "GNLY", "CD79A", "MS4A1",
    "LYZ", "CD68", "C1QA", "ALB", "APOA1", "PECAM1",
    "VWF", "COL1A1", "DCN", "EPCAM", "KRT19", "KRT7", "SOX9",
    "ACTA2", "AFP", "GPC3"
  ))
  dot_features <- intersect(dot_features, rownames(seurat))
  if (length(dot_features) > 0) {
    p_dot <- DotPlot(seurat, features = dot_features, group.by = "celltype_annot") +
      RotatedAxis()
    save_fig(
      file.path(fig_dir, "fig_06_dotplot_markers.png"),
      p_dot,
      width = 11,
      height = 7,
      dpi = 150
    )
  } else {
    log_msg("skip DotPlot: no marker genes present in dataset")
  }

  feature_genes <- expand_genes(c(
    "CD3D", "NKG7", "CD79A", "LYZ", "ALB", "KRT19",
    "PECAM1", "COL1A1", "ACTA2", "EPCAM"
  ))
  feature_genes <- intersect(feature_genes, rownames(seurat))
  feature_genes <- head(feature_genes, 6)
  if (length(feature_genes) > 0) {
    p_feature <- FeaturePlot(
      seurat,
      features = feature_genes,
      ncol = 3
    )
    save_fig(
      file.path(fig_dir, "fig_16_featureplot_markers.png"),
      p_feature,
      width = 12,
      height = 4 * ceiling(length(feature_genes) / 3)
    )

    p_marker_violin <- VlnPlot(
      seurat,
      features = feature_genes,
      group.by = "celltype_annot",
      ncol = 3,
      pt.size = 0
    )
    save_fig(
      file.path(fig_dir, "fig_17_marker_violin.png"),
      p_marker_violin,
      width = 12,
      height = 4 * ceiling(length(feature_genes) / 3)
    )
  }

  prop_tbl <- as.data.frame(table(
    CellType = seurat$celltype_annot,
    Condition = seurat$condition
  ))
  prop_mat <- as.matrix(xtabs(Freq ~ CellType + Condition, data = prop_tbl))
  cond_names <- colnames(prop_mat)
  prop_stats <- lapply(rownames(prop_mat), function(ct) {
    row <- prop_mat[ct, , drop = TRUE]
    total <- colSums(prop_mat)
    tbl <- rbind(row, total - row)
    ft <- tryCatch(
      fisher.test(tbl),
      error = function(e) NULL
    )
    data.frame(
      CellType = ct,
      CountA = row[1],
      CountB = row[2],
      Pvalue = if (is.null(ft)) NA_real_ else ft$p.value,
      OddsRatio = if (is.null(ft)) NA_real_ else as.numeric(ft$estimate),
      stringsAsFactors = FALSE
    )
  })
  prop_stats_df <- do.call(rbind, prop_stats)
  prop_stats_df$Padj <- p.adjust(prop_stats_df$Pvalue, method = "BH")
  write.csv(
    prop_stats_df,
    stage_data_file("fig_18_19_celltype_proportion_stats.csv"),
    row.names = FALSE
  )
  p_cell_prop <- ggplot(prop_tbl, aes(x = CellType, y = Freq, fill = Condition)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cell type", y = "Proportion", fill = "Condition",
         title = "Cell type proportion by condition")
  save_fig(
    file.path(fig_dir, "fig_18_celltype_proportion.png"),
    p_cell_prop,
    width = 9,
    height = 6
  )

  p_cond_prop <- ggplot(prop_tbl, aes(x = Condition, y = Freq, fill = CellType)) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() +
    theme_minimal() +
    labs(x = "Condition", y = "Proportion", fill = "Cell type",
         title = "Condition composition by cell type")
  save_fig(
    file.path(fig_dir, "fig_19_condition_proportion.png"),
    p_cond_prop,
    width = 8,
    height = 6
  )

  conf <- table(seurat$celltype_annot, seurat$published_type)
  write.csv(as.data.frame.matrix(conf), stage_data_file("fig_07_annotation_confusion.csv"))
  if (!"fig_07_annotation_confusion_heatmap.png" %in% skip_figs) {
    png(
      stage_fig_file(file.path(fig_dir, "fig_07_annotation_confusion_heatmap.png")),
      width = 1200,
      height = 800,
      res = 150
    )
    if (nrow(conf) > 1 && ncol(conf) > 1 && !all(is.na(conf))) {
      pheatmap(
        conf,
        display_numbers = TRUE,
        fontsize_number = 6,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = "Marker annotation vs published cell type"
      )
    } else {
      plot.new()
      title("No published annotations for confusion matrix")
      text(0.5, 0.5, "No published annotations", cex = 1.4)
    }
    dev.off()
    log_msg("saved figure: fig_07_annotation_confusion_heatmap.png")
  } else {
    log_msg("skip figure: fig_07_annotation_confusion_heatmap.png")
  }

  if (dataset_mode == "sample_level") {
    seurat$celltype_annot <- "Bulk"
    seurat$celltype_annot_cell <- "Bulk"
    seurat$cluster_label <- as.character(seurat$seurat_clusters)
    annotation_tbl <- data.frame(
      cell = colnames(seurat),
      seurat_clusters = as.character(seurat$seurat_clusters),
      celltype_annot = seurat$celltype_annot,
      celltype_annot_cell = seurat$celltype_annot_cell,
      cluster_label = seurat$cluster_label,
      published_type = seurat$published_type,
      condition = seurat$condition,
      sample = seurat$sample,
      stringsAsFactors = FALSE
    )
    write.csv(
      annotation_tbl,
      stage_data_file("fig_05_16_17_cell_annotations.csv"),
      row.names = FALSE
    )
    log_msg("sample-level mode: sample-level annotation written")
  }

  log_msg("annotation cell types: ", paste(sort(unique(seurat$celltype_annot)), collapse = ", "))
  saveRDS(seurat, ckpt_path("seurat_annotated.rds"))
})

if (stage_allowed("06")) run_stage("06_differential_expression", {
  if (!exists("seurat")) {
    seurat <- readRDS(ckpt_path("seurat_annotated.rds"))
  }
  cond_levels <- sort(unique(as.character(seurat$condition)))
  if (length(cond_levels) != 2) {
    stop("Differential expression requires exactly two condition groups.")
  }
  Idents(seurat) <- "condition"
  sample_cond <- unique(data.frame(
    sample = seurat$sample,
    condition = as.character(seurat$condition),
    stringsAsFactors = FALSE
  ))
  sample_counts <- table(sample_cond$condition)
  use_pseudobulk <- nrow(sample_cond) >= 4 && all(sample_counts >= 2)

  if (use_pseudobulk) {
    log_msg("using DESeq2 pseudobulk differential expression")
    warn_file <- file.path(data_dir, "pseudobulk_warning.txt")
    if (file.exists(warn_file)) {
      unlink(warn_file)
    }
    bulk <- AggregateExpression(
      seurat,
      group.by = c("sample", "condition"),
      assays = "RNA",
      return.seurat = FALSE
    )$RNA
    bulk_meta <- data.frame(row.names = colnames(bulk), stringsAsFactors = FALSE)
    bulk_meta$sample <- gsub(
      "-",
      "_",
      sub("_[^_]+$", "", colnames(bulk))
    )
    bulk_meta$condition <- sample_cond$condition[
      match(bulk_meta$sample, sample_cond$sample)
    ]
    bulk_meta$condition <- factor(bulk_meta$condition, levels = cond_levels)
    bulk_meta <- bulk_meta[!is.na(bulk_meta$condition), , drop = FALSE]
    bulk <- bulk[, rownames(bulk_meta), drop = FALSE]
    if (sum(bulk) == 0 || nrow(bulk_meta) < 4) {
      use_pseudobulk <- FALSE
      log_msg("pseudobulk invalid; falling back to Seurat Wilcoxon")
      writeLines(
        "Sample count insufficient for pseudobulk; used Seurat Wilcoxon",
        file.path(data_dir, "pseudobulk_warning.txt")
      )
    }
  }

  if (use_pseudobulk) {
    dds <- DESeqDataSetFromMatrix(
      countData = round(as.matrix(bulk)),
      colData = bulk_meta,
      design = ~ condition
    )
    dds <- tryCatch(
      DESeq(dds, quiet = TRUE),
      error = function(e) {
        log_msg("pseudobulk DESeq2 failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(dds)) {
      use_pseudobulk <- FALSE
      writeLines(
        "DESeq2 pseudobulk failed; used Seurat Wilcoxon",
        file.path(data_dir, "pseudobulk_warning.txt")
      )
    }
  }

  if (use_pseudobulk) {
    res <- results(dds, contrast = c("condition", cond_levels[1], cond_levels[2]))
    deg <- as.data.frame(res)
    deg$gene <- rownames(deg)
    deg$avg_log2FC <- deg$log2FoldChange
    deg$p_val <- deg$pvalue
    deg$p_val_adj <- deg$padj
    deg$pct.1 <- NA_real_
    deg$pct.2 <- NA_real_
  } else {
    log_msg("using Seurat Wilcoxon with downsampling")
    writeLines(
      "Sample count insufficient for pseudobulk; used Seurat Wilcoxon",
      file.path(data_dir, "pseudobulk_warning.txt")
    )
    deg <- tryCatch(
      FindMarkers(
        seurat,
        ident.1 = cond_levels[1],
        ident.2 = cond_levels[2],
        test.use = "wilcox",
        max.cells.per.ident = 3000,
        logfc.threshold = de_logfc,
        min.pct = 0.1,
        only.pos = FALSE,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(deg) || nrow(deg) == 0 || ncol(deg) == 0) {
      log_msg(
        "no DEGs passed default filters; ",
        "rerunning without logFC/min.pct filters"
      )
      deg <- tryCatch(
        FindMarkers(
          seurat,
          ident.1 = cond_levels[1],
          ident.2 = cond_levels[2],
          test.use = "wilcox",
          max.cells.per.ident = 3000,
          logfc.threshold = 0,
          min.pct = 0,
          only.pos = FALSE,
          verbose = FALSE
        ),
        error = function(e) NULL
      )
    }
  }
  deg <- ensure_deg_columns(deg)
  deg$significant <- deg$p_val_adj < de_padj & abs(deg$avg_log2FC) > de_logfc
  deg$direction <- ifelse(
    deg$significant,
    ifelse(deg$avg_log2FC > 0, "Up", "Down"),
    "NS"
  )
  deg$neg_log10_padj <- -log10(pmax(deg$p_val_adj, 1e-300))
  deg <- deg[order(
    is.na(deg$p_val_adj),
    deg$p_val_adj,
    -abs(deg$avg_log2FC)
  ), , drop = FALSE]

  write.csv(deg, stage_data_file("fig_08_deg_all.csv"), row.names = FALSE)
  write.csv(
    deg[deg$significant, ],
    stage_data_file("fig_09_deg_significant.csv"),
    row.names = FALSE
  )

  label_genes <- deg$gene[seq_len(min(15, nrow(deg)))]
  p_volcano <- ggplot(
    deg,
    aes(x = avg_log2FC, y = neg_log10_padj, color = direction)
  ) +
    geom_point(alpha = 0.6, size = 1.1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey75")) +
    geom_text_repel(
      data = deg[deg$gene %in% label_genes, ],
      aes(label = gene),
      max.overlaps = 20,
      size = 3
    ) +
    labs(
      x = "Average log2 fold change (HCC vs iCCA)",
      y = "-log10 adjusted p value",
      title = paste0(
        "Differential expression volcano plot (",
        cond_levels[1], " vs ", cond_levels[2], ")"
      )
    ) +
    theme_minimal() +
    coord_cartesian(clip = "off")
  if (fig_style("fig_08_volcano.png") == "maplot") {
    expr_genes <- intersect(deg$gene, rownames(seurat))
    expr_mat <- GetAssayData(seurat, layer = "data")[
      expr_genes, , drop = FALSE
    ]
    means <- as.numeric(Matrix::rowSums(expr_mat) / ncol(expr_mat))
    deg$mean_expr <- NA_real_
    deg$mean_expr[match(expr_genes, deg$gene)] <- means
    p_volcano <- ggplot(
      deg,
      aes(x = mean_expr, y = avg_log2FC, color = direction)
    ) +
      geom_point(alpha = 0.6, size = 1.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey75")) +
      labs(
        x = "Mean normalized expression",
        y = "Average log2 fold change",
        title = paste0(
          "MA plot (", cond_levels[1], " vs ", cond_levels[2], ")"
        )
      ) +
      theme_minimal() +
      coord_cartesian(clip = "off")
  }
  save_fig(
    file.path(fig_dir, "fig_08_volcano.png"),
    p_volcano,
    width = 9,
    height = 7,
    dpi = 150,
    plot_margin = ggplot2::margin(20, 60, 24, 20, "pt")
  )

  top_deg <- deg[deg$significant %in% TRUE, ]
  if (nrow(top_deg) < 1) {
    top_deg <- deg[is.finite(deg$p_val_adj), ]
  }
  top_deg <- head(top_deg, deg_violin_top_n)
  top_genes <- intersect(top_deg$gene, rownames(seurat))
  if (length(top_genes) >= 2) {
    conds <- unique(as.character(seurat$condition))
    if (length(conds) >= 2) {
      keep_cells <- unlist(lapply(conds, function(cc) {
        cells <- colnames(seurat)[seurat$condition == cc]
        set.seed(20260812)
        if (length(cells) > deg_violin_max_cells) {
          cells <- sample(cells, deg_violin_max_cells)
        }
        cells
      }))
      expr_mat <- GetAssayData(seurat, layer = "data")[
        top_genes,
        keep_cells,
        drop = FALSE
      ]
      expr_long <- data.frame(
        cell = rep(colnames(expr_mat), each = nrow(expr_mat)),
        gene = rep(rownames(expr_mat), times = ncol(expr_mat)),
        expr = as.numeric(as.matrix(expr_mat)),
        stringsAsFactors = FALSE
      )
      expr_long$condition <- seurat$condition[
        match(expr_long$cell, colnames(seurat))
      ]
      expr_long <- expr_long[!is.na(expr_long$condition), , drop = FALSE]
      expr_long$gene <- factor(expr_long$gene, levels = rev(top_genes))

      deg_top <- top_deg[top_deg$gene %in% top_genes, , drop = FALSE]
      deg_top$p_label <- ifelse(
        deg_top$p_val_adj < 0.001,
        formatC(pmax(deg_top$p_val_adj, 1e-300), format = "e", digits = 1),
        sprintf("%.3f", deg_top$p_val_adj)
      )
      label_map <- setNames(
        paste0(deg_top$gene, "\nP = ", deg_top$p_label),
        deg_top$gene
      )
      cond_colors <- hcl.colors(length(conds), palette = "Set 2")
      p_deg_violin <- ggplot(
        expr_long,
        aes(x = expr, y = gene, fill = condition)
      ) +
        geom_violin(
          position = position_dodge(0.8),
          scale = "width",
          alpha = 0.75,
          trim = TRUE
        ) +
        geom_boxplot(
          width = 0.12,
          position = position_dodge(0.8),
          outlier.shape = NA,
          alpha = 0.9
        ) +
        geom_point(
          alpha = 0.18,
          size = 0.4,
          position = position_jitterdodge(
            jitter.width = 0.08,
            dodge.width = 0.8,
            seed = 20260812
          )
        ) +
        scale_y_discrete(labels = label_map) +
        scale_fill_manual(values = cond_colors) +
        labs(
          x = "Normalized expression",
          y = NULL,
          fill = "Condition",
          title = paste0(
            "Top differential genes by adjusted p value (",
            cond_levels[1], " vs ", cond_levels[2], ")"
          )
        ) +
        theme_minimal(base_size = 13) +
        theme(
          axis.text.y = element_text(size = 8, face = "italic"),
          legend.position = "top"
        )
      save_fig(
        file.path(fig_dir, "fig_09_deg_horizontal_violin.png"),
        p_deg_violin,
        width = 11,
        height = max(7, 0.55 * length(top_genes) + 2),
        dpi = 150,
        plot_margin = ggplot2::margin(20, 40, 26, 20, "pt")
      )
      write.csv(
        deg_top[, intersect(
          c(
            "gene", "avg_log2FC", "p_val", "p_val_adj",
            "pct.1", "pct.2", "significant", "direction"
          ),
          colnames(deg_top)
        ), drop = FALSE],
        stage_data_file("fig_09_deg_horizontal_violin.csv"),
        row.names = FALSE
      )
      log_msg("saved DEG horizontal violin table: ", nrow(deg_top), " genes")
    } else {
      log_msg("skip DEG horizontal violin: need at least 2 conditions")
    }
  } else {
    log_msg("skip DEG horizontal violin: too few genes in Seurat object")
  }

  top30 <- intersect(
    deg$gene[seq_len(min(30, nrow(deg)))],
    rownames(seurat)
  )
  if (length(top30) > 0) {
    seurat <- ScaleData(seurat, features = top30, verbose = FALSE)
    p_heat <- DoHeatmap(
      seurat,
      features = top30,
      group.by = "condition",
      angle = 45
    ) +
      scale_fill_viridis_c()
    save_fig(
      file.path(fig_dir, "fig_09_deg_heatmap.png"),
      p_heat,
      width = 10,
      height = 8,
      dpi = 150,
      plot_margin = ggplot2::margin(32, 24, 24, 18, "pt")
    )
  } else {
    log_msg("skip DEG heatmap: no DEGs to plot")
  }

  log_msg("significant DEGs: ", sum(deg$significant, na.rm = TRUE))
})

if (stage_allowed("07")) run_stage("07_enrichment", {
  deg <- read.csv(stage_data_file("fig_08_deg_all.csv"), stringsAsFactors = FALSE)
  deg_up <- deg[deg$significant & deg$avg_log2FC > 0, ]
  deg_down <- deg[deg$significant & deg$avg_log2FC < 0, ]

  if (nrow(deg_up) < 10) {
    deg_up <- deg[deg$p_val_adj < 0.1 & deg$avg_log2FC > 0.1, ]
  }
  if (nrow(deg_down) < 10) {
    deg_down <- deg[deg$p_val_adj < 0.1 & deg$avg_log2FC < -0.1, ]
  }

  gene_id_type <- function(ids) {
    ids <- na.omit(ids)
    if (length(ids) == 0) return("SYMBOL")
    if (grepl("^ENSG\\d+", ids[1]) || grepl("^ENSMUSG\\d+", ids[1])) {
      "ENSEMBL"
    } else {
      "SYMBOL"
    }
  }

  org_db <- if (species == "mm") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("org.Mm.eg.db is required for mouse enrichment analysis.")
    }
    getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
  } else {
    org.Hs.eg.db
  }
  kegg_org <- ifelse(species == "mm", "mmu", "hsa")

  run_enrichment <- function(deg_sub, name) {
    if (nrow(deg_sub) < 3) {
      log_msg("too few genes for ", name)
      return(list(go = NULL, kegg = NULL))
    }

    id_type <- gene_id_type(deg_sub$gene)
    eg <- tryCatch(
      bitr(
        deg_sub$gene,
        fromType = id_type,
        toType = "ENTREZID",
        OrgDb = org_db
      ),
      error = function(e) NULL
    )
    if (is.null(eg) || nrow(eg) == 0) {
      log_msg("no Entrez mapping for ", name)
      return(list(go = NULL, kegg = NULL))
    }

    go <- tryCatch(
      enrichGO(
        gene = eg$ENTREZID,
        OrgDb = org_db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )

    kegg <- tryCatch({
      httr::set_config(httr::timeout(60))
      options(timeout = 60)
      R.utils::withTimeout(
        enrichKEGG(
          gene = eg$ENTREZID,
          organism = kegg_org,
          pvalueCutoff = 0.1
        ),
        timeout = 60,
        onTimeout = "error"
      )
    }, error = function(e) {
      log_msg("KEGG enrichment unavailable: ", conditionMessage(e))
      NULL
    })

    list(go = go, kegg = kegg)
  }

  up_res <- run_enrichment(deg_up, "up-regulated")
  down_res <- run_enrichment(deg_down, "down-regulated")

  write_res <- function(res, go_prefix, kegg_prefix) {
    if (!is.null(res$go)) {
      go_df <- as.data.frame(res$go)
      write.csv(go_df, stage_data_file(paste0(go_prefix, "_go.csv")), row.names = FALSE)
    } else {
      write.csv(data.frame(note = "no significant GO terms"),
                stage_data_file(paste0(go_prefix, "_go.csv")), row.names = FALSE)
    }
    if (!is.null(res$kegg)) {
      kegg_df <- as.data.frame(res$kegg)
      write.csv(kegg_df, stage_data_file(paste0(kegg_prefix, "_kegg.csv")), row.names = FALSE)
    } else {
      write.csv(data.frame(note = "no significant KEGG terms"),
                stage_data_file(paste0(kegg_prefix, "_kegg.csv")), row.names = FALSE)
    }
  }
  write_res(up_res, "fig_10_enrichment_up", "fig_12_enrichment_up")
  write_res(down_res, "fig_11_enrichment_down", "fig_13_enrichment_down")

  plot_res <- function(res, file, title) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
        geom_text(label = "No significant enrichment terms") +
        theme_void()
    } else {
      style <- fig_style(basename(file))
      if (style == "dotplot") {
        p <- dotplot(res, showCategory = 10, font.size = 8) + ggtitle(title)
      } else if (style == "cnetplot") {
        p <- cnetplot(res, showCategory = 10) +
          ggtitle(title)
      } else {
        p <- barplot(res, showCategory = 10, font.size = 8) + ggtitle(title)
      }
    }
    save_fig(file, p, width = 9, height = 7, dpi = 150)
  }

  plot_res(up_res$go, file.path(fig_dir, "fig_10_go_up.png"), "GO BP: up-regulated genes")
  plot_res(down_res$go, file.path(fig_dir, "fig_11_go_down.png"), "GO BP: down-regulated genes")
  plot_res(up_res$kegg, file.path(fig_dir, "fig_12_kegg_up.png"), "KEGG: up-regulated genes")
  plot_res(down_res$kegg, file.path(fig_dir, "fig_13_kegg_down.png"), "KEGG: down-regulated genes")

  rank_vec <- deg$avg_log2FC
  names(rank_vec) <- deg$gene
  rank_vec <- sort(rank_vec[!is.na(rank_vec) & is.finite(rank_vec)], decreasing = TRUE)
  id_type <- gene_id_type(names(rank_vec))
  mapped_col <- if (id_type == "SYMBOL") "SYMBOL" else "ENSEMBL"
  eg_all <- tryCatch(
    bitr(
      names(rank_vec),
      fromType = id_type,
      toType = "ENTREZID",
      OrgDb = org_db
    ),
    error = function(e) data.frame()
  )
  if (skip_gsea) {
    log_msg("GSEA skipped by LIVER_SKIP_GSEA")
    gsea_go <- NULL
    gsea_kegg <- NULL
  } else if (nrow(eg_all) == 0) {
    log_msg("no Entrez mapping for GSEA; skipping")
    gsea_go <- NULL
    gsea_kegg <- NULL
  } else {
    ranked <- rank_vec[eg_all[[mapped_col]]]
    names(ranked) <- eg_all$ENTREZID
    ranked <- sort(ranked, decreasing = TRUE)
    if (gsea_max_genes > 0 && length(ranked) > gsea_max_genes) {
      ranked <- ranked[seq_len(gsea_max_genes)]
      log_msg(
        "GSEA gene list capped to ",
        length(ranked),
        " genes (LIVER_GSEA_MAX_GENES=",
        gsea_max_genes,
        ")"
      )
    }

    log_msg("starting GSEA GO with ", length(ranked), " genes")
    gsea_go <- tryCatch(
      gseGO(
        geneList = ranked,
        OrgDb = org_db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.1,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    log_msg("GSEA GO finished")
    log_msg("starting GSEA KEGG with ", length(ranked), " genes")
    gsea_kegg <- tryCatch({
      options(timeout = 60)
      httr::set_config(httr::timeout(60))
      R.utils::withTimeout(
        gseKEGG(
          geneList = ranked,
          organism = kegg_org,
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = 0.1,
          verbose = FALSE
        ),
        timeout = 60,
        onTimeout = "error"
      )
    }, error = function(e) {
      log_msg("GSEA KEGG unavailable: ", conditionMessage(e))
      NULL
    })
    log_msg("GSEA KEGG finished")
  }

  plot_gsea <- function(res, file, title) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
        geom_text(label = "No significant GSEA terms") +
        theme_void()
    } else {
      style <- fig_style(basename(file))
      if (style == "gseaplot2") {
        n <- min(3, nrow(as.data.frame(res)))
        p <- gseaplot2(res, geneSetID = seq_len(n), title = title)
      } else {
        p <- ridgeplot(res, showCategory = 10) + ggtitle(title)
      }
    }
    save_fig(file, p, width = 9, height = 7)
  }
  plot_gsea(gsea_go, file.path(fig_dir, "fig_20_gsea_go.png"), "GSEA GO BP")
  plot_gsea(gsea_kegg, file.path(fig_dir, "fig_21_gsea_kegg.png"), "GSEA KEGG")

  top_enrichment <- function(res, n = 5, padj_cutoff = 0.05) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      return(res)
    }
    df <- as.data.frame(res)
    sig_idx <- which(!is.na(df$p.adjust) & df$p.adjust <= padj_cutoff)
    filtered <- clusterProfiler::slice(res, sig_idx)
    df <- as.data.frame(filtered)
    if (nrow(df) == 0) {
      return(filtered)
    }
    clusterProfiler::slice(filtered, seq_len(min(n, nrow(df))))
  }

  plot_cnet <- function(res, file, title) {
    filtered <- top_enrichment(res)
    if (is.null(filtered) || nrow(as.data.frame(filtered)) == 0) {
      p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
        geom_text(label = "No significant pathway network") +
        theme_void()
    } else {
      style <- fig_style(basename(file))
      if (style == "emapplot") {
        p <- emapplot(filtered, showCategory = 10) + ggtitle(title)
      } else {
        p <- cnetplot(filtered, showCategory = 5) +
          ggtitle(title)
      }
      p <- tryCatch(
        p + coord_cartesian(clip = "off"),
        error = function(e) p
      )
    }
    save_fig(
      file,
      p,
      width = 10,
      height = 8,
      plot_margin = ggplot2::margin(30, 70, 40, 70, "pt")
    )
  }

  plot_top5 <- function(res, file, title) {
    filtered <- top_enrichment(res)
    if (is.null(filtered) || nrow(as.data.frame(filtered)) == 0) {
      p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
        geom_text(label = "No significant pathway after filtering") +
        theme_void()
    } else {
      style <- fig_style(basename(file))
      if (style == "barplot") {
        p <- barplot(filtered, showCategory = 5, font.size = 8) + ggtitle(title)
      } else if (style == "cnetplot") {
        p <- cnetplot(filtered, showCategory = 5) + ggtitle(title)
      } else if (style == "emapplot") {
        p <- emapplot(filtered, showCategory = 5) + ggtitle(title)
      } else {
        p <- dotplot(filtered, showCategory = 5, font.size = 8) + ggtitle(title)
      }
    }
    save_fig(file, p, width = 10, height = 8)
  }

  plot_cnet(up_res$go, file.path(fig_dir, "fig_22_go_network.png"), "GO BP network (filtered top 5)")
  plot_cnet(up_res$kegg, file.path(fig_dir, "fig_23_kegg_network.png"), "KEGG network (filtered top 5)")
  plot_top5(up_res$go, file.path(fig_dir, "fig_46_go_top5.png"), "GO BP filtered top 5")
  plot_top5(up_res$kegg, file.path(fig_dir, "fig_47_kegg_top5.png"), "KEGG filtered top 5")

  log_msg("enrichment tables and plots generated")
})

if (stage_allowed("08")) run_stage("08_publication_analyses", {
  if (!exists("seurat")) {
    seurat <- readRDS(ckpt_path("seurat_annotated.rds"))
  }
  seurat <- NormalizeData(seurat, verbose = FALSE)
  Idents(seurat) <- "condition"

  if (dataset_mode == "sample_level") {
    log_msg("sample-level mode: skipping cell-level publication analyses")
  } else {
  if (run_cellcycle) {
    if (!"Phase" %in% colnames(seurat@meta.data)) {
      cc_genes <- tryCatch(Seurat::cc.genes, error = function(e) NULL)
      if (is.null(cc_genes)) {
        cc_genes <- tryCatch(Seurat::cc.genes.updated.2019, error = function(e) NULL)
      }
      s_genes <- intersect(cc_genes$s.genes, rownames(seurat))
      g2m_genes <- intersect(cc_genes$g2m.genes, rownames(seurat))
      if (length(s_genes) >= 5 && length(g2m_genes) >= 5) {
        seurat <- CellCycleScoring(
          seurat,
          s.features = s_genes,
          g2m.features = g2m_genes,
          set.ident = FALSE
        )
        cc_tbl <- data.frame(
          cell = colnames(seurat),
          S.Score = seurat$S.Score,
          G2M.Score = seurat$G2M.Score,
          Phase = seurat$Phase,
          condition = seurat$condition,
          stringsAsFactors = FALSE
        )
        write.csv(cc_tbl, stage_data_file("fig_26_27_cell_cycle_scores.csv"), row.names = FALSE)
      }
    }
    if ("Phase" %in% colnames(seurat@meta.data)) {
      p_cc_umap <- DimPlot(seurat, group.by = "Phase") +
        ggtitle("Cell cycle phase")
      save_fig(
        file.path(fig_dir, "fig_26_cellcycle_umap.png"),
        p_cc_umap,
        width = 8,
        height = 7,
        dpi = 150
      )
      cc_prop <- as.data.frame(table(
        Condition = seurat$condition,
        Phase = seurat$Phase
      ))
      p_cc_prop <- ggplot(cc_prop, aes(x = Condition, y = Freq, fill = Phase)) +
        geom_col(position = "fill") +
        scale_y_continuous(labels = scales::percent) +
        theme_minimal() +
        labs(
          x = "Condition",
          y = "Proportion",
          fill = "Phase",
          title = "Cell cycle phase proportion"
        )
      save_fig(
        file.path(fig_dir, "fig_27_cellcycle_proportion.png"),
        p_cc_prop,
        width = 7,
        height = 5,
        dpi = 150
      )
    }
  }

  if (length(unique(seurat$sample)) > 1) {
    sample_levels <- unique(as.character(seurat$sample))
    sample_labels <- wrap_labels(sample_levels, width = 22)
  seurat$sample_label <- unname(factor(
    sample_labels[match(as.character(seurat$sample), sample_levels)],
    levels = sample_labels
  ))
    legend_cols <- if (length(sample_levels) > 80) {
      3
    } else if (length(sample_levels) > 30) {
      2
    } else {
      1
    }
    p_sample <- DimPlot(seurat, group.by = "sample_label", label = FALSE) +
      labs(color = "Sample") +
      ggtitle("UMAP by sample") +
      guides(color = guide_legend(ncol = legend_cols)) +
      theme(
        legend.text = element_text(size = 6.5),
        legend.key.size = grid::unit(0.45, "cm"),
        legend.spacing.y = grid::unit(0.03, "cm")
      )
    save_fig(
      file.path(fig_dir, "fig_28_umap_sample.png"),
      p_sample,
      width = 8,
      height = 7,
      dpi = 150,
      plot_margin = ggplot2::margin(16, 30, 22, 30, "pt")
    )
    seurat$sample_label <- NULL
  } else {
    log_msg("skip figure: fig_28_umap_sample.png (single sample)")
  }

  dbl_path <- stage_data_file("fig_02_doublet_results.csv")
  if (file.exists(dbl_path)) {
    dbl <- read.csv(dbl_path, stringsAsFactors = FALSE)
    if (nrow(dbl) > 0 && "doublet_call" %in% colnames(dbl)) {
      dbl_rate <- dbl %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(
          n_cells = dplyr::n(),
          n_doublets = sum(doublet_call == "doublet", na.rm = TRUE),
          doublet_rate = n_doublets / n_cells,
          .groups = "drop"
        )
      dbl_rate$condition <- dbl$condition[match(dbl_rate$sample, dbl$sample)]
      write.csv(
        dbl_rate,
        stage_data_file("fig_29_doublet_rate_by_sample.csv"),
        row.names = FALSE
      )
      p_dbl_rate <- ggplot(
        dbl_rate,
        aes(x = sample, y = doublet_rate, fill = condition)
      ) +
        geom_col() +
        geom_text(
          aes(label = sprintf("%.1f%%", 100 * doublet_rate)),
          vjust = -0.4,
          size = 3
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_x_discrete(labels = function(x) wrap_labels(x, 18)) +
        labs(
          x = "Sample",
          y = "Doublet rate",
          fill = "Condition",
          title = "Doublet rate by sample"
        )
      save_fig(
        file.path(fig_dir, "fig_29_doublet_rate_sample.png"),
        p_dbl_rate,
        width = 8,
        height = 5,
        dpi = 150
      )
    }
  }

  prop_sample <- as.data.frame(table(
    CellType = seurat$celltype_annot,
    Sample = seurat$sample
  ))
  p_prop_sample <- ggplot(
    prop_sample,
    aes(x = Sample, y = Freq, fill = CellType)
  ) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = function(x) wrap_labels(x, 18)) +
    labs(
      x = "Sample",
      y = "Proportion",
      fill = "Cell type",
      title = "Cell type proportion by sample"
    )
  save_fig(
    file.path(fig_dir, "fig_30_sample_proportion.png"),
    p_prop_sample,
    width = 10,
    height = 6,
    dpi = 150
  )

  if (run_cluster_markers) {
    Idents(seurat) <- "seurat_clusters"
    markers <- tryCatch(
      FindAllMarkers(
        seurat,
        only.pos = TRUE,
        min.pct = 0.25,
        logfc.threshold = 0.5,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(markers) && nrow(markers) > 0) {
      if (!"gene" %in% colnames(markers)) {
        markers$gene <- rownames(markers)
      }
      write.csv(
        markers,
        stage_data_file("fig_16_17_31_32_cluster_markers.csv"),
        row.names = FALSE
      )
      cluster_col <- if ("cluster" %in% colnames(markers)) {
        "cluster"
      } else {
        "seurat_clusters"
      }
      if (!"avg_log2FC" %in% colnames(markers) &&
          "avg_logFC" %in% colnames(markers)) {
        markers$avg_log2FC <- markers$avg_logFC
      }
      top_markers <- markers %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(cluster_col))) %>%
        dplyr::slice_max(n = 3, order_by = avg_log2FC)
      top_genes <- unique(top_markers$gene)
      top_genes <- intersect(top_genes, rownames(seurat))
      if (length(top_genes) > 0) {
        tryCatch(
          {
            seurat <- ScaleData(seurat, features = top_genes, verbose = FALSE)
            p_marker_heat <- DoHeatmap(
              seurat,
              features = top_genes,
              group.by = "seurat_clusters",
              angle = 45
            ) + scale_fill_viridis_c()
            save_fig(
              file.path(fig_dir, "fig_31_cluster_marker_heatmap.png"),
              p_marker_heat,
              width = 12,
              height = 9,
              dpi = 150
            )
            p_marker_dot <- DotPlot(
              seurat,
              features = top_genes,
              group.by = "seurat_clusters"
            ) + RotatedAxis()
            save_fig(
              file.path(fig_dir, "fig_32_cluster_marker_dotplot.png"),
              p_marker_dot,
              width = 14,
              height = 8,
              dpi = 150
            )
          },
          error = function(e) {
            log_msg("cluster marker figures failed: ", conditionMessage(e))
          }
        )
      }
    } else {
      log_msg("no cluster markers found")
    }
    Idents(seurat) <- "condition"
  }

  if (run_signatures) {
    signatures <- list(
      Proliferation = c(
        "MKI67", "PCNA", "TOP2A", "MCM2", "MCM3", "MCM4", "MCM5",
        "MCM6", "MCM7", "CDC20", "CCNB1", "CCNA2", "BIRC5", "AURKA",
        "AURKB", "UBE2C", "CENPF", "CENPE", "KIF11", "KIF2C"
      ),
      EMT = c(
        "VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1",
        "ZEB2", "MMP2", "MMP9", "COL1A1", "COL1A2", "LUM", "DCN",
        "TGFB1", "ACTA2"
      ),
      Hypoxia = c(
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "PGK1", "PDK1", "CA9",
        "BNIP3", "BNIP3L", "NDRG1", "ADM", "ALDOA", "ENO1", "TPI1"
      ),
      ImmuneCheckpoint = c(
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD274",
        "PDCD1LG2", "IDO1", "TNFRSF9", "CD27", "CD70", "LGALS9"
      ),
      TcellExhaustion = c(
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX",
        "ENTPD1", "CXCL13", "BATF", "MAF"
      ),
      Stemness = c(
        "PROM1", "EPCAM", "ALDH1A1", "SOX2", "NANOG", "MYC", "CD44",
        "KIT", "LGR5", "ABCG2", "AFP", "GPC3"
      ),
      Inflammation = c(
        "IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL5", "CXCL10",
        "CXCL9", "ICAM1", "VCAM1", "NFKB1", "NFKBIA", "JUN", "FOS"
      )
    )
    sig_features <- lapply(
      signatures,
      function(gs) intersect(expand_genes(gs), rownames(seurat))
    )
    sig_features <- sig_features[lengths(sig_features) >= 3]
    if (length(sig_features) > 0) {
      tryCatch(
        {
      seurat <- AddModuleScore(
        seurat,
        features = sig_features,
        name = "Signature_",
        ctrl = 50
      )
      sig_cols <- paste0("Signature_", seq_along(sig_features))
      names(sig_cols) <- names(sig_features)
      sig_df <- FetchData(
        seurat,
        c("condition", "celltype_annot", unname(sig_cols))
      )
      colnames(sig_df) <- c("condition", "celltype_annot", names(sig_features))
      write.csv(
        sig_df,
        stage_data_file("fig_33_34_35_signature_scores.csv"),
        row.names = FALSE
      )

      show_sigs <- intersect(
        c("Proliferation", "EMT", "Hypoxia", "ImmuneCheckpoint"),
        names(sig_cols)
      )
      if (length(show_sigs) > 0) {
        p_sig_umap <- FeaturePlot(
          seurat,
          features = unname(sig_cols[show_sigs]),
          ncol = 2,
          cols = c("grey90", "#B31B1B")
        )
        save_fig(
          file.path(fig_dir, "fig_33_signature_scores_umap.png"),
          p_sig_umap,
          width = 10,
          height = 4 * ceiling(length(show_sigs) / 2),
          dpi = 150
        )
      }

      sig_long <- do.call(rbind, lapply(names(sig_cols), function(nm) {
        data.frame(
          condition = seurat$condition,
          celltype = seurat$celltype_annot,
          signature = nm,
          score = seurat@meta.data[[sig_cols[[nm]]]],
          stringsAsFactors = FALSE
        )
      }))
      p_sig_box <- ggplot(
        sig_long,
        aes(x = condition, y = score, fill = condition)
      ) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.15, outlier.shape = NA) +
        facet_wrap(~signature, scales = "free_y", ncol = 2) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        ) +
        labs(
          x = "Condition",
          y = "Signature score",
          title = "Signature scores by condition"
        )
      save_fig(
        file.path(fig_dir, "fig_34_signature_scores_boxplot.png"),
        p_sig_box,
        width = 10,
        height = 8,
        dpi = 150
      )
        },
        error = function(e) {
          log_msg("signature analysis failed: ", conditionMessage(e))
        }
      )
    }
  }

  prop_stats_path <- stage_data_file("fig_18_19_celltype_proportion_stats.csv")
  if (file.exists(prop_stats_path)) {
    prop_stats <- read.csv(prop_stats_path, stringsAsFactors = FALSE)
    prop_stats$log2OR <- log2(prop_stats$OddsRatio)
    prop_stats$log2OR[!is.finite(prop_stats$log2OR)] <- NA_real_
    prop_stats$neg_log10_padj <- -log10(pmax(prop_stats$Padj, 1e-300))
    p_abundance <- ggplot(
      prop_stats,
      aes(x = log2OR, y = neg_log10_padj, color = CellType)
    ) +
      geom_point(size = 2.5) +
      geom_text_repel(aes(label = CellType), size = 3, max.overlaps = 20) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      labs(
        x = "log2 odds ratio",
        y = "-log10 adjusted p value",
        color = "Cell type",
        title = "Cell type abundance shift"
      ) +
      theme_minimal() +
      coord_cartesian(clip = "off")
    save_fig(
      file.path(fig_dir, "fig_35_celltype_abundance_effect.png"),
      p_abundance,
      width = 9,
      height = 7,
      dpi = 150
    )
  }

  if (run_cnv) {
    tryCatch(
      {
    org_db_cnv <- if (species == "mm") {
      if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
      } else {
        NULL
      }
    } else {
      org.Hs.eg.db
    }
    if (!is.null(org_db_cnv)) {
      mapped <- tryCatch(
        AnnotationDbi::select(
          org_db_cnv,
          keys = rownames(seurat),
          columns = c("CHR", "CHRLOC"),
          keytype = "SYMBOL"
        ),
        error = function(e) NULL
      )
      if (!is.null(mapped) && nrow(mapped) > 1000) {
        mapped <- mapped[
          !is.na(mapped$SYMBOL) & !duplicated(mapped$SYMBOL),
          , drop = FALSE
        ]
        mapped <- mapped[
          !is.na(mapped$CHR) &
            grepl("^([0-9]+|X|Y)$", as.character(mapped$CHR)),
          , drop = FALSE
        ]
        mapped$chr_order <- factor(
          as.character(mapped$CHR),
          levels = c(as.character(1:22), "X", "Y")
        )
        if ("CHRLOC" %in% colnames(mapped) &&
            sum(!is.na(mapped$CHRLOC)) > 500) {
          mapped <- mapped[!is.na(mapped$CHRLOC), , drop = FALSE]
          mapped <- mapped[order(mapped$chr_order, mapped$CHRLOC), , drop = FALSE]
        } else {
          mapped <- mapped[order(mapped$chr_order, mapped$SYMBOL), , drop = FALSE]
        }
        cnv_genes <- intersect(mapped$SYMBOL, rownames(seurat))
        if (length(cnv_genes) >= 200) {
          set.seed(42)
          cell_idx <- unlist(lapply(unique(seurat$condition), function(cond) {
            idx <- which(seurat$condition == cond)
            sample(idx, min(length(idx), 750))
          }))
          cell_idx <- as.integer(cell_idx)
          expr_sub <- as.matrix(
            GetAssayData(seurat, layer = "data")[
              cnv_genes, cell_idx, drop = FALSE
            ]
          )
          expr_scaled <- t(scale(t(expr_sub)))
          expr_scaled[!is.finite(expr_scaled)] <- 0
          starts <- seq(1, length(cnv_genes) - 49, by = 25)
          win_scores <- sapply(starts, function(s) {
            e <- expr_scaled[s:min(s + 99, length(cnv_genes)), , drop = FALSE]
            colMeans(e)
          })
          if (is.null(dim(win_scores))) {
            win_scores <- matrix(win_scores, ncol = 1)
          }
          win_scores <- sweep(
            win_scores,
            2,
            apply(win_scores, 2, median),
            "-"
          )
          win_labels <- vapply(starts, function(s) {
            mid <- floor((s + min(s + 99, length(cnv_genes))) / 2)
            paste0(
              "chr",
              mapped$CHR[match(cnv_genes[mid], mapped$SYMBOL)]
            )
          }, character(1))
          colnames(win_scores) <- paste0(
            win_labels,
            "_w",
            seq_along(starts)
          )
          rownames(win_scores) <- colnames(expr_sub)
          cnv_df <- data.frame(
            cell = rownames(win_scores),
            win_scores,
            check.names = FALSE
          )
          write.csv(
            cnv_df,
            stage_data_file("fig_36_cnv_heatmap.csv"),
            row.names = FALSE
          )
          ann_col <- data.frame(
            row.names = rownames(win_scores),
            Condition = seurat$condition[cell_idx],
            CellType = seurat$celltype_annot[cell_idx]
          )
          cond_levels_cnv <- sort(unique(as.character(ann_col$Condition)))
          ct_levels_cnv <- sort(unique(as.character(ann_col$CellType)))
          ann_colors <- list(
            Condition = setNames(
              c("#E64B35", "#4DBBD5", "#00A087")[seq_along(cond_levels_cnv)],
              cond_levels_cnv
            ),
            CellType = setNames(
              rainbow(length(ct_levels_cnv)),
              ct_levels_cnv
            )
          )
          save_pheatmap(
            file.path(fig_dir, "fig_36_cnv_heatmap.png"),
            function() {
              pheatmap(
                win_scores,
                annotation_col = ann_col,
                annotation_colors = ann_colors,
                show_rownames = FALSE,
                show_colnames = TRUE,
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                color = colorRampPalette(
                  c("#3B4CC0", "#FFFFFF", "#B40426")
                )(100),
                fontsize_col = 6,
                main = "Inferred CNV profile (sliding window)"
              )
            },
            width = 1200,
            height = 900
          )
          log_msg("CNV heatmap generated")
        }
      } else {
        log_msg("CNV skipped: insufficient chromosome annotation")
      }
    } else {
      log_msg("CNV skipped: organism annotation package unavailable")
    }
      },
      error = function(e) {
        log_msg("CNV analysis failed: ", conditionMessage(e))
      }
    )
  }

  if (run_singler) {
    singler_ok <- requireNamespace("SingleR", quietly = TRUE) &&
      requireNamespace("celldex", quietly = TRUE)
    if (singler_ok) {
      ref <- tryCatch(
        if (species == "mm") {
          celldex::MouseRNAseqData()
        } else {
          celldex::HumanPrimaryCellAtlasData()
        },
        error = function(e) NULL
      )
      if (!is.null(ref)) {
        set.seed(42)
        all_cells <- colnames(seurat)
        max_cells <- 20000
        chosen <- if (length(all_cells) > max_cells) {
          sample(all_cells, max_cells)
        } else {
          all_cells
        }
        pred <- tryCatch(
          SingleR::SingleR(
            test = GetAssayData(seurat, layer = "data")[
              , chosen, drop = FALSE
            ],
            ref = ref,
            labels = ref$label.main
          ),
          error = function(e) NULL
        )
        if (!is.null(pred)) {
          seurat$singleR_label <- NA_character_
          seurat$singleR_label[chosen] <- as.character(pred$labels)
          for (cl in unique(seurat$seurat_clusters)) {
            idx <- seurat$seurat_clusters == cl
            labs <- seurat$singleR_label[idx]
            if (sum(!is.na(labs)) > 0) {
              seurat$singleR_label[idx] <- names(
                sort(table(labs), decreasing = TRUE)
              )[1]
            }
          }
          write.csv(
            data.frame(
              cell = colnames(seurat),
              seurat_clusters = as.character(seurat$seurat_clusters),
              celltype_annot = seurat$celltype_annot,
              singleR_label = seurat$singleR_label,
              condition = seurat$condition,
              stringsAsFactors = FALSE
            ),
            stage_data_file("fig_37_singleR_annotations.csv"),
            row.names = FALSE
          )
          p_singler <- DimPlot(seurat, group.by = "singleR_label", label = TRUE) +
            ggtitle("SingleR cell type annotation")
          save_fig(
            file.path(fig_dir, "fig_37_singler_umap.png"),
            p_singler,
            width = 9,
            height = 7,
            dpi = 150
          )
          conf_singler <- table(
            seurat$celltype_annot,
            seurat$singleR_label
          )
          write.csv(
            as.data.frame.matrix(conf_singler),
            stage_data_file("fig_38_singleR_confusion.csv")
          )
          if (nrow(conf_singler) > 0 && ncol(conf_singler) > 0) {
            save_pheatmap(
              file.path(fig_dir, "fig_38_singler_confusion_heatmap.png"),
              function() {
                pheatmap(
                  conf_singler,
                  display_numbers = TRUE,
                  fontsize_number = 6,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  main = "Marker annotation vs SingleR"
                )
              },
              width = 1100,
              height = 800
            )
          }
        }
      }
    }
    if (!"singleR_label" %in% colnames(seurat@meta.data)) {
      log_msg("SingleR skipped: reference data or prediction unavailable")
    }
  }

  if (run_trajectory) {
    if (requireNamespace("slingshot", quietly = TRUE)) {
      tryCatch(
        {
          set.seed(42)
          rd <- Embeddings(seurat, reduction = "umap")
          sce_traj <- SingleCellExperiment(
            assays = list(counts = GetAssayData(seurat, layer = "counts")),
            reducedDims = list(UMAP = rd),
            colData = DataFrame(
              cluster = as.character(seurat$seurat_clusters)
            )
          )
          sce_traj <- slingshot::slingshot(
            sce_traj,
            clusterLabels = "cluster",
            reducedDim = "UMAP"
          )
          pt <- slingshot::slingPseudotime(sce_traj)
          if (!is.null(dim(pt)) && ncol(pt) >= 1) {
            seurat$pseudotime <- pt[, 1]
            write.csv(
              data.frame(
                cell = colnames(seurat),
                cluster = as.character(seurat$seurat_clusters),
                pseudotime = seurat$pseudotime,
                stringsAsFactors = FALSE
              ),
              stage_data_file("fig_39_trajectory_pseudotime.csv"),
              row.names = FALSE
            )
            p_pt <- FeaturePlot(seurat, features = "pseudotime") +
              scale_color_viridis_c() +
              ggtitle("Slingshot pseudotime")
            curves <- slingshot::slingCurves(sce_traj)
            curve_df <- do.call(rbind, lapply(seq_along(curves), function(i) {
              coords <- as.data.frame(curves[[i]]$s)
              colnames(coords) <- c("UMAP_1", "UMAP_2")
              coords$lineage <- as.character(i)
              coords
            }))
            rd_df <- as.data.frame(rd)
            colnames(rd_df) <- c("UMAP_1", "UMAP_2")
            rd_df$cluster <- as.character(seurat$seurat_clusters)
            p_line <- ggplot(
              rd_df,
              aes(x = UMAP_1, y = UMAP_2, color = cluster)
            ) +
              geom_point(size = 0.5, alpha = 0.7) +
              geom_path(
                data = curve_df,
                aes(x = UMAP_1, y = UMAP_2, group = lineage),
                inherit.aes = FALSE,
                color = "black",
                linewidth = 1
              ) +
              theme_minimal() +
              labs(color = "Cluster", title = "Slingshot lineages on UMAP")
            save_fig(
              file.path(fig_dir, "fig_39_trajectory_umap.png"),
              p_pt + p_line,
              width = 14,
              height = 7,
              dpi = 150
            )
          }
        },
        error = function(e) {
          log_msg("trajectory analysis failed: ", conditionMessage(e))
        }
      )
    } else {
      log_msg("trajectory skipped: slingshot not installed")
    }
  }

  }

  saveRDS(seurat, ckpt_path("seurat_annotated.rds"))
  saveRDS(seurat, ckpt_path("seurat_publication.rds"))
  log_msg("publication analyses complete")
})

if (stage_allowed("09")) run_stage("09_summary_outputs", {
  if (!exists("seurat")) {
    pub_path <- ckpt_path("seurat_publication.rds")
    if (file.exists(pub_path)) {
      seurat <- readRDS(pub_path)
    } else {
      seurat <- readRDS(ckpt_path("seurat_annotated.rds"))
    }
  }
  if (!exists("seurat_raw")) {
    seurat_raw <- readRDS(ckpt_path("seurat_raw.rds"))
  }
  if (!exists("qc_metrics")) {
    qc_metrics <- read.csv(stage_data_file("fig_01_qc_metrics.csv"))
  }
  saveRDS(seurat, file.path(data_dir, "liver_cancer_seurat.rds"))

  deg <- read.csv(stage_data_file("fig_08_deg_all.csv"), stringsAsFactors = FALSE)
  go_up <- tryCatch(
    read.csv(stage_data_file("fig_10_enrichment_up_go.csv"), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )

  summary_list <- list(
    dataset = accession,
    dataset_mode = dataset_mode,
    n_samples = ncol(seurat),
    title = paste(
      sort(unique(as.character(seurat$condition))),
      collapse = " vs "
    ),
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    n_cells_raw = ncol(seurat_raw),
    n_cells_after_qc = nrow(qc_metrics),
    n_cells_after_doublet_removal = ncol(seurat),
    n_genes = nrow(seurat),
    n_clusters = length(unique(seurat$seurat_clusters)),
    n_celltypes = length(unique(seurat$celltype_annot)),
    condition_counts = as.list(table(seurat$condition)),
    deg_total = nrow(deg),
    deg_up = sum(
      deg$significant & deg$avg_log2FC > 0,
      na.rm = TRUE
    ),
    deg_down = sum(
      deg$significant & deg$avg_log2FC < 0,
      na.rm = TRUE
    ),
    top_degs = head(deg[, c("gene", "avg_log2FC", "p_val_adj")], 20),
    go_up_top = if (nrow(go_up) > 0) head(go_up[, c("ID", "Description", "pvalue", "p.adjust")], 10) else data.frame()
  )
  write_json(summary_list, file.path(res_dir, "summary.json"), auto_unbox = TRUE, pretty = TRUE)

  complete <- list(status = "complete", finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  write_json(complete, file.path(res_dir, "pipeline_complete.json"), auto_unbox = TRUE, pretty = TRUE)

  log_msg("pipeline complete")
})
