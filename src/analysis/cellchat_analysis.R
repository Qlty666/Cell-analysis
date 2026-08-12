#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
species <- if (length(args) > 1) args[2] else "hs"

if (!requireNamespace("CellChat", quietly = TRUE)) {
  status_dir <- file.path(root, "results", "data", "09_cellchat")
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("CellChat not installed", file.path(status_dir, "cellchat_status.txt"))
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
})

seurat <- readRDS(file.path(root, "results", "checkpoints", "seurat_annotated.rds"))
data_mat <- GetAssayData(seurat, layer = "data")
meta <- seurat@meta.data
cellchat <- createCellChat(
  object = data_mat,
  meta = meta,
  group.by = "celltype_annot"
)

CellChatDB <- if (species == "mm") {
  CellChatDB.mouse
} else {
  CellChatDB.human
}
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

data_dir <- file.path(root, "results", "data", "09_cellchat")
fig_dir <- file.path(root, "results", "figures", "09_cellchat")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_cc_fig <- function(name, fn, width = 1200, height = 900) {
  file <- file.path(fig_dir, name)
  png(file, width = width, height = height, res = 150)
  on.exit(dev.off())
  fn()
  cat("saved figure: ", name, "\n")
}

tryCatch(
  {
    net_count <- cellchat@net$count
    net_df <- as.data.frame(as.table(net_count))
    colnames(net_df) <- c("Source", "Target", "Count")
    write.csv(
      net_df,
      file.path(data_dir, "fig_40_cellchat_communication.csv"),
      row.names = FALSE
    )
    net_weight <- cellchat@net$weight
    weight_df <- as.data.frame(as.table(net_weight))
    colnames(weight_df) <- c("Source", "Target", "Weight")
    write.csv(
      weight_df,
      file.path(data_dir, "fig_40_cellchat_communication_weight.csv"),
      row.names = FALSE
    )
  },
  error = function(e) {
    cat("cellchat table export failed: ", conditionMessage(e), "\n")
  }
)

tryCatch(
  {
    path_df <- subsetCommunication(cellchat)
    write.csv(
      path_df,
      file.path(data_dir, "fig_42_cellchat_pathways.csv"),
      row.names = FALSE
    )
  },
  error = function(e) {
    cat("cellchat pathway export failed: ", conditionMessage(e), "\n")
  }
)

tryCatch(
  save_cc_fig(
    "fig_40_cellchat_network.png",
    function() {
      netVisual_circle(
        cellchat@net$count,
        vertex.weight = as.numeric(table(cellchat@idents)),
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = "Number of interactions"
      )
    },
    1200,
    1200
  ),
  error = function(e) {
    cat("cellchat network figure failed: ", conditionMessage(e), "\n")
  }
)

tryCatch(
  save_cc_fig(
    "fig_41_cellchat_heatmap.png",
    function() {
      netVisual_heatmap(cellchat, measure = "count")
    },
    1000,
    800
  ),
  error = function(e) {
    cat("cellchat heatmap figure failed: ", conditionMessage(e), "\n")
  }
)

tryCatch(
  save_cc_fig(
    "fig_42_cellchat_bubble.png",
    function() {
      netVisual_bubble(cellchat, remove.isolate = TRUE)
    },
    1100,
    900
  ),
  error = function(e) {
    cat("cellchat bubble figure failed: ", conditionMessage(e), "\n")
  }
)

saveRDS(cellchat, file.path(data_dir, "cellchat_object.rds"))
writeLines("CellChat completed", file.path(data_dir, "cellchat_status.txt"))
