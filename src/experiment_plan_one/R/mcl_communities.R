#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, required = TRUE) {
  prefix <- paste0("--", name, "=")
  value <- args[startsWith(args, prefix)]
  if (!length(value)) {
    if (required) stop(sprintf("missing required argument --%s", name))
    return(NULL)
  }
  sub(prefix, "", value[[1]], fixed = TRUE)
}

edges_path <- get_arg("edges")
nodes_path <- get_arg("nodes", required = FALSE)
output_path <- get_arg("output")
inflation_arg <- get_arg("inflation", required = FALSE)
inflation <- if (is.null(inflation_arg)) 2 else as.numeric(inflation_arg)

read_nodes <- function(path) {
  if (is.null(path)) return(character())
  nodes <- readLines(path, warn = FALSE)
  unique(nodes[nzchar(nodes)])
}

edges <- read.delim(
  edges_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("node1", "node2") %in% names(edges))) {
  stop("edge table must contain node1 and node2 columns")
}
if (!"score" %in% names(edges)) edges$score <- 1
edges$node1 <- as.character(edges$node1)
edges$node2 <- as.character(edges$node2)
edges$score <- as.numeric(edges$score)

all_nodes <- unique(c(read_nodes(nodes_path), edges$node1, edges$node2))
if (!length(all_nodes)) {
  write.table(
    data.frame(node = character(), cluster = integer()),
    output_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  quit(status = 0)
}

adjacency <- matrix(
  0,
  nrow = length(all_nodes),
  ncol = length(all_nodes),
  dimnames = list(all_nodes, all_nodes)
)
for (index in seq_len(nrow(edges))) {
  source <- edges$node1[[index]]
  target <- edges$node2[[index]]
  weight <- edges$score[[index]]
  if (!is.finite(weight) || weight <= 0) next
  adjacency[source, target] <- adjacency[source, target] + weight
  adjacency[target, source] <- adjacency[target, source] + weight
}

if (!any(adjacency > 0)) {
  clusters <- seq_along(all_nodes)
} else {
  result <- MCL::mcl(
    adjacency,
    addLoops = TRUE,
    inflation = inflation,
    max.iter = 100
  )
  clusters <- as.integer(result$Cluster)
  if (length(clusters) != length(all_nodes)) {
    stop("MCL returned a cluster vector with unexpected length")
  }
  clusters[clusters < 1] <- 0L
}

sizes <- table(clusters)
major <- as.integer(names(sizes)[sizes >= 3])
ordering <- order(-as.integer(sizes)[match(major, as.integer(names(sizes)))], major)
remap <- seq_along(ordering)
names(remap) <- as.character(major[ordering])
for (index in seq_along(all_nodes)) {
  value <- clusters[[index]]
  if (as.character(value) %in% names(remap)) {
    clusters[[index]] <- unname(remap[[as.character(value)]])
  } else {
    clusters[[index]] <- 0L
  }
}

write.table(
  data.frame(node = all_nodes, cluster = as.integer(clusters)),
  output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
