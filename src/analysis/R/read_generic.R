# read_generic.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

read_generic_counts <- function(manifest) {
  files <- manifest$files
  matrices <- as.character(files$matrix)
  barcodes <- as.character(files$barcodes)
  genes <- as.character(files$genes)
  if (length(matrices) == 0) {
    stop("No count matrix files found in dataset manifest.")
  }
  if (identical(manifest$mode, "single_cell") && length(matrices) > 1) {
    h5_suffixes <- grepl(
      "\\.(h5ad(\\.gz)?|h5|loom|rds)$",
      matrices,
      ignore.case = TRUE
    )
    if (any(h5_suffixes) && any(!h5_suffixes)) {
      dropped <- matrices[!h5_suffixes]
      matrices <- matrices[h5_suffixes]
      barcodes <- character()
      genes <- character()
      log_msg(
        "single-cell mode: dropped sample-level matrix files ",
        paste(dropped, collapse = ", ")
      )
    }
  }

  count_list_all <- list()
  sample_list <- list()
  group_list <- list()
  meta_list <- list()
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
    embedded_meta <- NULL
    if (exists("sample_label", inherits = FALSE)) {
      rm(sample_label)
    }

    if (grepl("\\.h5ad(\\.gz)?$", mat_file, ignore.case = TRUE)) {
      h5_loaded <- read_h5ad_matrix(mat_path)
      m <- h5_loaded$counts
      embedded_meta <- h5_loaded$meta
      sample_label <- "Sample1"
      sample_labels <- rep("Sample1", ncol(m))
    } else if (grepl("\\.rds$", mat_file, ignore.case = TRUE)) {
      obj <- readRDS(mat_path)
      if (inherits(obj, "Seurat")) {
        m <- GetAssayData(obj, layer = "counts")
      } else if (inherits(obj, "SingleCellExperiment")) {
        m <- assay(obj, "counts")
      } else {
        m <- as(as.matrix(obj), "CsparseMatrix")
      }
      embedded_meta <- if (inherits(obj, "Seurat")) obj[[]] else NULL
    } else if (grepl("\\.h5$", mat_file, ignore.case = TRUE)) {
      m <- Read10X_h5(mat_path)
      if (is.list(m) && !inherits(m, "Matrix")) {
        pick_h5_gex <- function(x) {
          if (inherits(x, "Matrix") || is.matrix(x)) {
            return(x)
          }
          if (is.list(x) && "Gene Expression" %in% names(x)) {
            return(x[["Gene Expression"]])
          }
          NULL
        }
        picked <- pick_h5_gex(m)
        if (is.null(picked)) {
          for (entry in m) {
            picked <- pick_h5_gex(entry)
            if (!is.null(picked)) break
          }
        }
        if (is.null(picked)) {
          stop("H5 file has no Gene Expression matrix: ", mat_file)
        }
        log_msg("selecting Gene Expression matrix from multi-modal H5: ", mat_file)
        m <- picked
      }
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
        g <- normalize_ensembl_ids(as.character(tab[[1]]))
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
        g <- normalize_ensembl_ids(as.character(tab[[1]]))
        drop_cols <- 1
        if (
          ncol(tab) > 2 &&
          tolower(trimws(as.character(colnames(tab)[2]))) %in% c(
            "gene", "genename", "gene_name", "symbol", "genesymbol"
          )
        ) {
          drop_cols <- 2
        }
        m <- as.matrix(tab[, -seq_len(drop_cols), with = FALSE])
        if (is.character(m)) storage.mode(m) <- "double"
        rownames(m) <- make.unique(as.character(g))
        colnames(m) <- colnames(tab)[-seq_len(drop_cols)]
      }
      m <- as(m, "CsparseMatrix")
      if (ncol(m) > 0) {
        sample_labels <- make.unique(as.character(colnames(m)))
      }
      sample_label <- sample_labels[1]
    }
    bulk_sample_mode <- exists("sample_label")
    if (!bulk_sample_mode) {
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
        gsm_token <- regmatches(
          mat_file,
          regexpr(
            "GSM[0-9]+_[A-Za-z0-9-]+(?=-matrix\\.mtx)",
            mat_file,
            perl = TRUE
          )
        )
        if (length(gsm_token) > 0) {
          filename_sample <- gsub("-", "_", gsm_token)
        } else {
          filename_sample <- gsm_match
        }
      } else {
        token_match <- regmatches(
          mat_file,
          regexpr(
            "[A-Za-z0-9]+_[A-Za-z0-9]+(?=_[^_/]*\\.mtx)",
            mat_file,
            perl = TRUE
          )
        )
        if (length(token_match) > 0) {
          filename_sample <- token_match
        }
      }
    }
  }
    }
    if (bulk_sample_mode) {
      colnames(m) <- make.unique(colnames(m))
    } else {
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
    }
    count_list_all[[i]] <- m
    sample_list[[i]] <- sample_labels
    group_list[[i]] <- rep(infer_group_from_filename(mat_file), ncol(m))
    meta_list[[i]] <- embedded_meta
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
  kept_meta <- meta_list[keep]
  common <- Reduce(intersect, lapply(count_list, rownames))
  counts <- do.call(
    cbind,
    lapply(count_list, function(m) m[common, , drop = FALSE])
  )
  list(
    counts = counts,
    cell_sample = cell_sample,
    cell_group = cell_group,
    meta_list = kept_meta
  )
}

read_generic_dataset <- function(manifest) {
  loaded <- read_generic_counts(manifest)
  counts <- loaded$counts
  meta <- read_generic_metadata(
    manifest,
    colnames(counts),
    loaded$cell_sample,
    loaded$meta_list
  )
  ann <- parse_series_generic(manifest)

  if (!is.null(ann) && "sample" %in% colnames(ann)) {
    ann_samples <- as.character(ann$sample)
    count_samples <- colnames(counts)
    idx <- if ("geo_accession" %in% colnames(ann)) {
      count_gsm <- sub(
        "^(GSM[0-9]+).*$",
        "\\1",
        count_samples
      )
      ann_gsm <- as.character(ann$geo_accession)
      gsm_idx <- match(count_gsm, ann_gsm)
      if (sum(!is.na(gsm_idx)) >= 2) gsm_idx else match(count_samples, ann_samples)
    } else {
      match(count_samples, ann_samples)
    }
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
    if (sum(!is.na(idx)) < 5) {
      extract_patient_id <- function(x) {
        m <- regexpr(
          "(?:patient[ _-]?)?(?:P|Patient)[ _-]?[0-9]+",
          x,
          ignore.case = TRUE
        )
        out <- rep(NA_character_, length(x))
        hits <- m > 0
        if (any(hits)) {
          raw <- regmatches(x, m)
          out[hits] <- sub(
            "^[^0-9]*",
            "",
            gsub("[ _-]", "", raw[hits])
          )
        }
        toupper(out)
      }
      idx <- match(
        extract_patient_id(count_samples),
        extract_patient_id(ann_samples)
      )
    }
    if (sum(!is.na(idx)) >= 2) {
      sample_vals <- loaded$cell_sample
      matched <- !is.na(idx)
      sample_vals[matched] <- ann_samples[idx[matched]]
      meta$sample_annotation <- sample_vals
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
    cond <- infer_condition_from_sample_names(as.character(meta$sample))
    if (!is.null(cond)) {
      log_msg("inferred condition from sample name tokens")
      meta$condition <- cond
    }
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
  meta$condition <- decode_condition_labels(meta$condition, accession)
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

