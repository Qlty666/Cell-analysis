#!/usr/bin/env Rscript
# Cell-cell communication, including condition-specific and differential
# analyses when a condition column is available.

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
data_dir <- file.path(root, "results", "data", "09_cellchat")
fig_dir <- file.path(root, "results", "figures", "09_cellchat")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cellchat_db <- if (species == "mm") {
  CellChatDB.mouse
} else {
  CellChatDB.human
}

save_cc_fig <- function(name, fn, width = 1200, height = 900) {
  file <- file.path(fig_dir, name)
  png(file, width = width, height = height, res = 150)
  on.exit(dev.off(), add = TRUE)
  fn()
  cat("saved figure: ", name, "\n")
}

build_cellchat <- function(object, group.by = "celltype_annot") {
  data_mat <- GetAssayData(object, layer = "data")
  chat <- createCellChat(
    object = data_mat,
    meta = object@meta.data,
    group.by = group.by
  )
  chat@DB <- cellchat_db
  chat <- subsetData(chat)
  chat <- identifyOverExpressedGenes(chat)
  chat <- identifyOverExpressedInteractions(chat)
  chat <- computeCommunProb(chat)
  chat <- filterCommunication(chat, min.cells = 10)
  chat <- computeCommunProbPathway(chat)
  chat <- aggregateNet(chat)
  chat
}

export_network <- function(chat, prefix) {
  tryCatch(
    {
      net_count <- chat@net$count
      net_df <- as.data.frame(as.table(net_count))
      colnames(net_df) <- c("Source", "Target", "Count")
      write.csv(
        net_df,
        file.path(data_dir, paste0(prefix, "_communication.csv")),
        row.names = FALSE
      )
      net_weight <- chat@net$weight
      weight_df <- as.data.frame(as.table(net_weight))
      colnames(weight_df) <- c("Source", "Target", "Weight")
      write.csv(
        weight_df,
        file.path(data_dir, paste0(prefix, "_communication_weight.csv")),
        row.names = FALSE
      )
    },
    error = function(e) {
      cat("cellchat network export failed: ", conditionMessage(e), "\n")
    }
  )
  tryCatch(
    {
      path_df <- subsetCommunication(chat)
      write.csv(
        path_df,
        file.path(data_dir, paste0(prefix, "_pathways.csv")),
        row.names = FALSE
      )
      interactions <- subsetCommunication(chat)
      if (nrow(interactions) > 0) {
        write.csv(
          interactions,
          file.path(data_dir, paste0(prefix, "_ligand_receptor.csv")),
          row.names = FALSE
        )
      }
    },
    error = function(e) {
      cat("cellchat pathway export failed: ", conditionMessage(e), "\n")
    }
  )
}

global_chat <- build_cellchat(seurat)
export_network(global_chat, "fig_40_cellchat")

tryCatch(
  save_cc_fig(
    "fig_40_cellchat_network.png",
    function() {
      netVisual_circle(
        global_chat@net$count,
        vertex.weight = as.numeric(table(global_chat@idents)),
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = "Number of interactions"
      )
    },
    1200,
    1200
  ),
  error = function(e) cat("global network figure failed: ", conditionMessage(e), "\n")
)

tryCatch(
  save_cc_fig(
    "fig_41_cellchat_heatmap.png",
    function() netVisual_heatmap(global_chat, measure = "count"),
    1000,
    800
  ),
  error = function(e) cat("global heatmap figure failed: ", conditionMessage(e), "\n")
)

tryCatch(
  save_cc_fig(
    "fig_42_cellchat_bubble.png",
    function() netVisual_bubble(global_chat, remove.isolate = TRUE),
    1100,
    900
  ),
  error = function(e) cat("global bubble figure failed: ", conditionMessage(e), "\n")
)

saveRDS(global_chat, file.path(data_dir, "cellchat_object.rds"))

# Condition-specific and differential communication analysis.
if ("condition" %in% colnames(seurat@meta.data)) {
  conditions <- unique(as.character(seurat$condition))
  conditions <- conditions[!is.na(conditions) & nzchar(conditions)]
  if (length(conditions) >= 2) {
    chat_list <- list()
    for (condition in conditions) {
      cells <- colnames(seurat)[as.character(seurat$condition) == condition]
      if (length(cells) < 10) next
      subset_obj <- subset(seurat, cells = cells)
      tryCatch(
        {
          chat <- build_cellchat(subset_obj)
          chat_list[[condition]] <- chat
          export_network(chat, paste0("fig_40_condition_", make.names(condition)))
          saveRDS(
            chat,
            file.path(data_dir, paste0("cellchat_", make.names(condition), ".rds"))
          )
        },
        error = function(e) {
          cat(
            "condition-specific CellChat failed for ",
            condition,
            ": ",
            conditionMessage(e),
            "\n",
            sep = ""
          )
        }
      )
    }
    if (length(chat_list) >= 2) {
      tryCatch(
        {
          merged <- mergeCellChat(chat_list, add.names = names(chat_list))
          saveRDS(merged, file.path(data_dir, "cellchat_merged.rds"))
          if (exists("netVisual_diffInteraction", where = asNamespace("CellChat"))) {
            save_cc_fig(
              "fig_40_cellchat_diff_interactions.png",
              function() {
                netVisual_diffInteraction(
                  merged,
                  weight.scale = TRUE,
                  measure = "count"
                )
              },
              1200,
              1000
            )
          }
          if (exists("rankNet", where = asNamespace("CellChat"))) {
            rank_df <- rankNet(
              merged,
              mode = "comparison",
              stacked = TRUE,
              do.stat = TRUE
            )
            write.csv(
              rank_df,
              file.path(data_dir, "fig_40_cellchat_pathway_rank.csv"),
              row.names = FALSE
            )
          }
          if (exists("netAnalysis_computeCentrality", where = asNamespace("CellChat"))) {
            for (condition in names(chat_list)) {
              tryCatch(
                {
                  chat <- netAnalysis_computeCentrality(
                    chat_list[[condition]],
                    slot.name = "netP"
                  )
                  centrality_path <- file.path(
                    data_dir,
                    paste0("cellchat_centrality_", make.names(condition), ".csv")
                  )
                  centrality <- chat@netP$centr
                  if (is.list(centrality) && length(centrality) > 0) {
                    rows <- list()
                    for (measure in names(centrality)) {
                      value <- centrality[[measure]]
                      if (is.null(value)) next
                      rows[[measure]] <- data.frame(
                        cell_type = names(value),
                        measure = measure,
                        value = as.numeric(value),
                        condition = condition,
                        row.names = NULL
                      )
                    }
                    if (length(rows)) {
                      write.csv(
                        do.call(rbind, rows),
                        centrality_path,
                        row.names = FALSE
                      )
                    }
                  }
                },
                error = function(e) {
                  cat(
                    "centrality export failed for ",
                    condition,
                    ": ",
                    conditionMessage(e),
                    "\n",
                    sep = ""
                  )
                }
              )
            }
          }
        },
        error = function(e) {
          cat("merged CellChat analysis failed: ", conditionMessage(e), "\n")
        }
      )
    }
  }
}

writeLines("CellChat completed", file.path(data_dir, "cellchat_status.txt"))
