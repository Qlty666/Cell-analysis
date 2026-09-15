# zz_utils.R - extracted from analysis_pipeline.R
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

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
  if (num <= 1 || num %in% c(48, 49)) return("01_qc")
  if (num == 2) return("02_doublets")
  if (num %in% c(3, 4, 14, 15)) return("03_cluster")
  if (num %in% c(5, 6, 7, 16, 17, 18, 19, 50, 51)) return("04_annotation")
  if (num %in% c(8, 9)) return("05_deg")
  if (num %in% c(10, 11, 12, 13, 20, 21, 22, 23, 46, 47)) return("06_enrichment")
  if (num %in% c(24, 25, 43, 44, 45)) return("07_ml")
  if (num >= 26 && num <= 39) return("08_publication")
  if (num >= 40 && num <= 42) return("09_cellchat")
  return("00_other")
}

sample_short_label <- function(x) {
  sub(":.*$", "", trimws(as.character(x)))
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
      # Missing statistics must stay missing: filling p-values with 0 would
      # silently mark every gene as significant.
      deg[[col]] <- rep(NA_real_, nrow(deg))
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

ckpt_path <- function(name) {
  file.path(ckpt_dir, name)
}

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

parse_values <- function(line) {
  v <- sub("^![^\t]*\t", "", line)
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(v, tmp)
  tab <- tryCatch(
    data.table::fread(
      tmp,
      header = FALSE,
      sep = "\t",
      quote = '"',
      fill = TRUE,
      colClasses = "character"
    ),
    error = function(e) NULL
  )
  if (is.null(tab) || ncol(tab) == 0) {
    return(character())
  }
  vals <- unname(unlist(as.list(tab[1, ]), use.names = FALSE))
  gsub('^"|"$', "", trimws(vals))
}

