#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
species <- if (length(args) > 1) args[2] else "hs"

if (!requireNamespace("CellChat", quietly = TRUE)) {
  dir.create(file.path(root, "results", "data"), recursive = TRUE, showWarnings = FALSE)
  writeLines("CellChat not installed", file.path(root, "results", "data", "cellchat_status.txt"))
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

saveRDS(cellchat, file.path(root, "results", "data", "cellchat_object.rds"))
writeLines("CellChat completed", file.path(root, "results", "data", "cellchat_status.txt"))
