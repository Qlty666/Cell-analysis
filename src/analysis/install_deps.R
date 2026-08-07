#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

r_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(r_lib)) {
  dir.create(r_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(r_lib, .libPaths()))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

cran_pkgs <- c(
  "Seurat", "dplyr", "ggplot2", "patchwork", "Matrix", "data.table",
  "jsonlite", "ggrepel", "pheatmap", "RColorBrewer", "harmony"
)

bioc_pkgs <- c(
  "scDblFinder", "SingleCellExperiment", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "BiocParallel",
  "SingleR", "celldex", "DESeq2"
)

if (nzchar(r_lib)) {
  all_libs <- unique(c(r_lib, .libPaths()))
  installed <- unique(unlist(lapply(
    all_libs,
    function(lib) rownames(installed.packages(lib.loc = lib))
  )))
} else {
  installed <- rownames(installed.packages())
}

to_cran <- setdiff(cran_pkgs, installed)
if (length(to_cran) > 0) {
  cat("Installing missing CRAN packages:", to_cran, sep = "\n")
  if (nzchar(r_lib)) {
    install.packages(to_cran, lib = r_lib, quiet = TRUE)
  } else {
    install.packages(to_cran, quiet = TRUE)
  }
}

to_bioc <- setdiff(bioc_pkgs, installed)
if (length(to_bioc) > 0) {
  cat("Installing missing Bioconductor packages:", to_bioc, sep = "\n")
  if (nzchar(r_lib)) {
    BiocManager::install(
      to_bioc,
      lib = r_lib,
      update = FALSE,
      ask = FALSE
    )
  } else {
    BiocManager::install(to_bioc, update = FALSE, ask = FALSE)
  }
}

cat("All pipeline dependencies are available.\n")
