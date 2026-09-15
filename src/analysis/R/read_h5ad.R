# read_h5ad.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

h5_vector <- function(item) {
  if (is.null(item)) return(NULL)
  value <- item$read()
  if (is.data.frame(value)) {
    value <- value[[1]]
  }
  if (is.character(value)) return(value)
  if (is.factor(value)) return(as.character(value))
  if (is.raw(value)) return(rawToChar(value))
  as.vector(value)
}

h5_obs_column <- function(obs_group, name) {
  item <- obs_group[[name]]
  if (inherits(item, "H5D")) {
    value <- h5_vector(item)
    if (is.integer(value) &&
        length(unique(value[!is.na(value)])) <= 2 &&
        all(unique(value[!is.na(value)]) %in% c(0L, 1L))) {
      value <- as.logical(value)
    }
    return(value)
  }
  if (inherits(item, "H5Group")) {
    categories <- h5_vector(item[["categories"]])
    codes <- as.integer(h5_vector(item[["codes"]]))
    if (length(categories) == 0 || length(codes) == 0) {
      return(rep(NA_character_, length(codes)))
    }
    out <- rep(NA_character_, length(codes))
    valid <- !is.na(codes) & codes >= 0 & codes < length(categories)
    out[valid] <- categories[codes[valid] + 1L]
    return(out)
  }
  NULL
}

h5_index_name <- function(group) {
  if (is.null(group)) return(NULL)
  if ("_index" %in% names(group)) return("_index")
  attr_names <- tryCatch(hdf5r::h5attr_names(group), error = function(e) character())
  if ("_index" %in% attr_names) {
    value <- tryCatch(group$attr_open("_index")$read(), error = function(e) NULL)
    if (!is.null(value) && length(value) == 1 && nzchar(as.character(value))) {
      return(as.character(value))
    }
  }
  NULL
}

read_h5ad_matrix <- function(path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("hdf5r is required to read .h5ad inputs; run install_deps.R")
  }
  plain_path <- path
  temp_plain <- NULL
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    if (!requireNamespace("R.utils", quietly = TRUE)) {
      stop("R.utils is required to decompress .h5ad.gz inputs")
    }
    temp_plain <- tempfile(fileext = ".h5ad")
    R.utils::gunzip(
      filename = path,
      destname = temp_plain,
      remove = FALSE,
      overwrite = TRUE
    )
    plain_path <- temp_plain
  }

  h5 <- hdf5r::h5file(plain_path, mode = "r")
  on.exit({
    try(h5$close_all(), silent = TRUE)
    if (!is.null(temp_plain) && file.exists(temp_plain)) {
      unlink(temp_plain)
    }
  }, add = TRUE)

  h5ad_read_group <- function(x_group, label) {
    shape <- as.integer(x_group$attr_open("shape")$read())
    if (length(shape) != 2) {
      stop("h5ad ", label, " does not have a 2D shape: ", path)
    }
    attr_names <- tryCatch(
      hdf5r::h5attr_names(x_group),
      error = function(e) character()
    )
    encoding <- if ("encoding-type" %in% attr_names) {
      as.character(x_group$attr_open("encoding-type")$read())
    } else {
      NULL
    }
    if ("data" %in% names(x_group)) {
      data <- h5_vector(x_group[["data"]])
      indices <- as.integer(h5_vector(x_group[["indices"]]))
      indptr <- as.integer(h5_vector(x_group[["indptr"]]))
      if (identical(encoding, "csc_matrix")) {
        matrix_c <- new(
          "dgCMatrix",
          p = indptr,
          i = indices,
          x = as.numeric(data),
          Dim = c(shape[1], shape[2])
        )
      } else {
        # AnnData defaults to CSR when encoding-type is absent/unrecognised.
        matrix_c <- new(
          "dgRMatrix",
          p = indptr,
          j = indices,
          x = as.numeric(data),
          Dim = c(shape[1], shape[2])
        )
      }
      rm(data, indices, indptr)
      gc(FALSE)
    } else {
      dense <- h5_vector(x_group)
      matrix_c <- Matrix::Matrix(
        dense,
        nrow = shape[1],
        ncol = shape[2],
        sparse = TRUE
      )
      rm(dense)
      gc(FALSE)
    }
    # h5ad stores cells x genes; the pipeline needs genes x cells.
    counts <- as(t(matrix_c), "CsparseMatrix")
    rm(matrix_c)
    gc(FALSE)
    counts
  }

  h5ad_is_integer <- function(m) {
    values <- as.numeric(m@x)
    all(is.finite(values)) &&
      all(values >= 0) &&
      all(abs(values - round(values)) < 1e-8)
  }

  h5ad_candidates <- list()
  layers <- tryCatch(h5[["layers"]], error = function(e) NULL)
  if (inherits(layers, "H5Group")) {
    layer_candidates <- c("counts", "raw_counts", "count")
    layer_name <- layer_candidates[layer_candidates %in% names(layers)][1]
    if (length(layer_name) == 1 && !is.na(layer_name)) {
      h5ad_candidates[[length(h5ad_candidates) + 1L]] <- list(
        group = layers[[layer_name]],
        label = paste0("layers/", layer_name),
        var_group = tryCatch(h5[["var"]], error = function(e) NULL)
      )
    }
  }
  raw_group <- tryCatch(h5[["raw"]], error = function(e) NULL)
  if (inherits(raw_group, "H5Group") && "X" %in% names(raw_group)) {
    h5ad_candidates[[length(h5ad_candidates) + 1L]] <- list(
      group = raw_group[["X"]],
      label = "raw/X",
      var_group = if ("var" %in% names(raw_group)) {
        raw_group[["var"]]
      } else {
        tryCatch(h5[["var"]], error = function(e) NULL)
      }
    )
  }
  x_main <- tryCatch(h5[["X"]], error = function(e) NULL)
  if (!is.null(x_main)) {
    h5ad_candidates[[length(h5ad_candidates) + 1L]] <- list(
      group = x_main,
      label = "X",
      var_group = tryCatch(h5[["var"]], error = function(e) NULL)
    )
  }

  # Precedence: layers["counts"] -> .raw/X -> X, but only matrices whose stored
  # values are non-negative integers are accepted as counts. X is frequently
  # log-normalized, so it must never be silently rounded into fake counts.
  counts <- NULL
  counts_label <- NULL
  counts_var <- NULL
  tried_labels <- character()
  non_integer_labels <- character()
  for (candidate in h5ad_candidates) {
    tried_labels <- c(tried_labels, candidate$label)
    candidate_counts <- h5ad_read_group(candidate$group, candidate$label)
    if (h5ad_is_integer(candidate_counts)) {
      counts <- candidate_counts
      counts_label <- candidate$label
      counts_var <- candidate$var_group
      break
    }
    non_integer_labels <- c(non_integer_labels, candidate$label)
    rm(candidate_counts)
    gc(FALSE)
  }
  if (is.null(counts)) {
    stop(
      "No integer count matrix found in h5ad ", basename(path), ". Checked: ",
      paste(tried_labels, collapse = ", "),
      ". Non-integer (likely log-normalized) matrices: ",
      paste(non_integer_labels, collapse = ", "),
      ". Provide raw counts in layers['counts'] or .raw/X."
    )
  }
  log_msg(
    "h5ad: using ", counts_label, " as integer counts from ", basename(path)
  )

  obs <- h5[["obs"]]
  index_name <- h5_index_name(obs)
  n_cells <- ncol(counts)
  cells <- if (is.null(index_name)) {
    rep("", n_cells)
  } else {
    h5_vector(obs[[index_name]])
  }
  if (length(cells) == 0 || all(!nzchar(cells))) {
    cells <- paste0("Cell", seq_len(n_cells))
  }
  cells <- make.unique(as.character(cells))

  meta <- data.frame(row.names = cells, stringsAsFactors = FALSE)
  if (!is.null(obs)) {
    for (name in names(obs)) {
      if (identical(name, index_name)) next
      value <- h5_obs_column(obs, name)
      if (is.null(value)) next
      if (length(value) == length(cells)) {
        meta[[name]] <- value
      }
    }
  }

  var <- if (!is.null(counts_var)) counts_var else h5[["var"]]
  var_index <- h5_index_name(var)
  genes <- if (is.null(var_index)) {
    h5_vector(var[["gene_ids"]])
  } else {
    h5_vector(var[[var_index]])
  }
  if (length(genes) != nrow(counts)) {
    genes <- h5_vector(var[["gene_ids"]])
  }
  if (length(genes) != nrow(counts)) {
    genes <- paste0("Gene", seq_len(nrow(counts)))
  }
  if (length(genes) != 0 && length(genes) == nrow(counts) && "exclude" %in% names(var)) {
    exclude <- as.logical(h5_vector(var[["exclude"]]))
    exclude[is.na(exclude)] <- FALSE
    if (any(exclude)) {
      counts <- counts[!exclude, , drop = FALSE]
      genes <- genes[!exclude]
      log_msg("h5ad: dropped ", sum(exclude), " genes flagged exclude in ", basename(path))
    }
  }
  genes <- normalize_ensembl_ids(genes)
  rownames(counts) <- make.unique(as.character(genes))
  colnames(counts) <- cells

  log_msg(
    "h5ad loaded: ", basename(path), " cells=", ncol(counts),
    " genes=", nrow(counts)
  )
  list(counts = counts, meta = meta)
}

