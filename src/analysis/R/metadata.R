# metadata.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

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

read_generic_metadata <- function(
  manifest,
  cells,
  cell_sample = NULL,
  embedded_meta = NULL
) {
  meta <- data.frame(row.names = cells)
  if (!is.null(embedded_meta) && length(embedded_meta) > 0) {
    for (frame in embedded_meta) {
      if (is.null(frame) || nrow(frame) == 0 || ncol(frame) == 0) next
      frame_cells <- make.unique(as.character(rownames(frame)))
      idx <- match(cells, frame_cells)
      for (col in colnames(frame)) {
        meta[[col]] <- frame[[col]][idx]
      }
    }
  }
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
    "CellType", "Cell.type", "cell.type", "cell_type_annot",
    "major_cluster", "sub_cluster"
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

  out <- list()
  for (f in files) {
    con <- gzfile(file.path(raw_dir, f), "rt")
    lines <- readLines(con, warn = FALSE)
    close(con)

    title_line <- lines[startsWith(lines, "!Sample_title")]
    if (length(title_line) == 0) next
    titles <- parse_values(title_line[1])
    df <- data.frame(sample = titles, stringsAsFactors = FALSE)

    geo_line <- lines[startsWith(lines, "!Sample_geo_accession")]
    if (length(geo_line) > 0) {
      df$geo_accession <- parse_values(geo_line[1])
    }

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
  tryCatch(
    do.call(rbind, out),
    error = function(e) {
      log_msg("series matrix column mismatch; using rbind.fill: ", conditionMessage(e))
      if (requireNamespace("plyr", quietly = TRUE)) {
        plyr::rbind.fill(out)
      } else if (requireNamespace("dplyr", quietly = TRUE)) {
        dplyr::bind_rows(out)
      } else {
        stop("Cannot combine series matrix tables with mismatched columns")
      }
    }
  )
}

infer_condition_from_sample_names <- function(samples) {
  samples <- as.character(samples)
  if (length(samples) < 4) return(NULL)
  tokenize <- function(s) {
    parts <- unlist(strsplit(s, "[_.-]", perl = TRUE))
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    norm <- tolower(parts)
    norm <- sub("[0-9]+$", "", norm)
    norm <- sub("^[0-9]+", "", norm)
    parts[norm != ""]
  }
  tokens <- lapply(samples, tokenize)
  max_pos <- max(lengths(tokens))
  if (max_pos < 1) return(NULL)
  best_groups <- 0L
  best_pos <- 0L
  for (pos in seq_len(max_pos)) {
    vals <- vapply(tokens, function(t) {
      if (length(t) >= pos) {
        v <- tolower(t[[pos]])
        sub("[0-9]+$", "", sub("^[0-9]+", "", trimws(v)))
      } else {
        NA_character_
      }
    }, character(1))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    tt <- table(vals)
    tt <- tt[tt >= 2]
    if (length(tt) >= 2 &&
        sum(tt) >= length(samples) * 0.6 &&
        length(tt) >= best_groups) {
      best_groups <- length(tt)
      best_pos <- pos
    }
  }
  if (best_pos == 0L) return(NULL)
  vapply(tokens, function(t) {
    if (length(t) >= best_pos) {
      v <- tolower(t[[best_pos]])
      sub("[0-9]+$", "", sub("^[0-9]+", "", trimws(v)))
    } else {
      NA_character_
    }
  }, character(1))
}

reduce_condition_to_two_groups <- function(vals) {
  vals <- as.character(vals)
  uni <- unique(vals[!is.na(vals) & nzchar(trimws(vals))])
  if (length(uni) <= 2) {
    return(vals)
  }

  disease_prefixes <- c("ICC", "HCC")
  prefix_vals <- sub("^\\s*([A-Za-z]+).*$", "\\1", vals)
  if (all(prefix_vals %in% disease_prefixes) &&
      length(unique(prefix_vals)) >= 2) {
    return(prefix_vals)
  }

  tissue_vals <- vals
  tissue_vals <- ifelse(
    grepl("normal|adjacent|control", tissue_vals, ignore.case = TRUE),
    "Normal",
    tissue_vals
  )
  tissue_vals <- ifelse(
    grepl("tumor|cancer|metastasis|lymph node", tissue_vals, ignore.case = TRUE),
    "Tumor",
    tissue_vals
  )
  tissue_uni <- unique(
    tissue_vals[!is.na(tissue_vals) & nzchar(trimws(tissue_vals))]
  )
  if (length(tissue_uni) == 2 && all(c("Normal", "Tumor") %in% tissue_uni)) {
    return(tissue_vals)
  }

  tt <- sort(table(vals), decreasing = TRUE)
  if (length(tt) >= 2) {
    # Behaviour choice: minority labels are kept as "Other" instead of being
    # silently relabelled as the majority group, which would fabricate data.
    top2 <- names(tt)[1:2]
    out <- ifelse(vals %in% top2, vals, "Other")
    other_counts <- sort(table(out[out == "Other"]))
    log_msg(
      "WARNING: condition column has more than two groups; keeping ",
      paste(top2, collapse = ", "), " and labelling the remaining ",
      sum(out == "Other"), " samples as 'Other'"
    )
    warn_path <- file.path(data_dir, "condition_warning.txt")
    write(
      paste0(
        "condition_reduction: kept [", paste(top2, collapse = ", "),
        "]; relabelled ", sum(out == "Other"), " samples as 'Other'",
        if (length(other_counts) > 0) {
          paste0(
            " (original labels: ",
            paste(
              names(other_counts), as.integer(other_counts),
              sep = "=", collapse = ", "
            ),
            ")"
          )
        } else {
          ""
        }
      ),
      file = warn_path,
      append = file.exists(warn_path)
    )
    return(out)
  }
  vals
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

  match_meta_sample <- if (
    !is.null(sample_ann) && "sample_annotation" %in% colnames(meta)
  ) {
    as.character(meta$sample_annotation)
  } else if (!is.null(sample_ann) && "sample" %in% colnames(meta)) {
    as.character(meta$sample)
  } else {
    NULL
  }
  if (!is.null(match_meta_sample)) {
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
        idx <- match_samples(match_meta_sample, sample_ann$sample)
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
          reduced <- reduce_condition_to_two_groups(cond)
          if (length(unique(reduced[!is.na(reduced)])) >= 2) {
            return(reduced)
          }
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
    dropped <- sort(table(cond[!keep]))
    if (length(dropped) > 0) {
      # Behaviour choice: samples outside the two largest condition groups are
      # dropped explicitly (never relabelled into another group) and recorded.
      log_msg(
        "WARNING: dropping ", sum(dropped), " samples whose condition is not ",
        "one of the two largest groups (",
        paste(names(tt)[1:2], collapse = ", "), "): ",
        paste(names(dropped), as.integer(dropped), sep = "=", collapse = ", ")
      )
      warn_path <- file.path(data_dir, "condition_warning.txt")
      write(
        paste0(
          "condition_filter: dropped ", sum(dropped), " samples not in [",
          paste(names(tt)[1:2], collapse = ", "), "]; dropped counts: ",
          paste(names(dropped), as.integer(dropped), sep = "=", collapse = ", ")
        ),
        file = warn_path,
        append = file.exists(warn_path)
      )
    }
  }
  meta <- meta[keep, , drop = FALSE]
  meta$condition <- factor(as.character(meta$condition))
  meta
}

decode_condition_labels <- function(vals, accession) {
  vals <- as.character(vals)
  unique_vals <- unique(vals[!is.na(vals) & nzchar(trimws(vals))])
  if (length(unique_vals) == 0 || !all(unique_vals %in% c("P", "T", "N"))) {
    return(vals)
  }
  mapping <- c(
    P = if (identical(accession, "GSE235863")) "blood" else "PB",
    T = if (identical(accession, "GSE235863")) "liver tumor" else "Tumor",
    N = "Normal"
  )
  mapped <- unname(mapping[vals])
  mapped[is.na(mapped)] <- vals[is.na(mapped)]
  mapped
}

