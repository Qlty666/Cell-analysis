source("src/analysis/R/celltype_pseudobulk.R")
set.seed(42)
counts <- matrix(rnbinom(100 * 80, mu = 20, size = 3), nrow = 100)
rownames(counts) <- paste0("G", seq_len(100))
colnames(counts) <- paste0("C", seq_len(80))
meta <- data.frame(sample = rep(paste0("S", 1:8), each = 10),
                   condition = rep(c("Control", "Disease"), each = 40),
                   celltype_annot = rep(rep(c("T", "B"), each = 5), 8), row.names = colnames(counts))
counts[1:10, 41:80] <- counts[1:10, 41:80] * 4
result <- liver_celltype_pseudobulk(Matrix::Matrix(counts, sparse = TRUE), meta)
stopifnot(nrow(result$results) > 0, all(vapply(result$status, function(s) s$status == "completed", logical(1))))
meta$batch <- meta$condition
bad <- liver_celltype_pseudobulk(Matrix::Matrix(counts, sparse = TRUE), meta, "batch")
stopifnot(nrow(bad$results) == 0, all(vapply(bad$status, function(s) s$status == "not_estimable", logical(1))))
cat("CELLTYPE_PSEUDOBULK_OK\n")
