# enrichment.R - KEGG REST helpers for the enrichment stage.
# Function bodies are unchanged; see analysis_pipeline.R for the stages.

# Robust KEGG REST access with a local disk cache and retries. clusterProfiler
# downloads pathway annotations from rest.kegg.jp during enrichKEGG/gseKEGG,
# which can be slow or flaky; cache the parsed tables so reruns are offline.

kegg_fetch_text <- function(url, timeout = 600, attempts = 3) {
  last_error <- NULL
  for (i in seq_len(attempts)) {
    resp <- tryCatch(
      httr::GET(url, httr::timeout(timeout)),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) {
      return(httr::content(resp, as = "text", encoding = "UTF-8"))
    }
    if (i < attempts) Sys.sleep(3 * i)
  }
  stop(
    "KEGG REST request failed after ", attempts, " attempts: ",
    if (!is.null(last_error)) conditionMessage(last_error) else "non-200 response"
  )
}

kegg_rest_cached <- function(rest_url) {
  cache_file <- file.path(
    kegg_cache_dir,
    paste0(gsub("[^A-Za-z0-9]", "_", rest_url), ".rds")
  )
  if (file.exists(cache_file)) {
    return(readRDS(cache_file))
  }
  message("Reading KEGG annotation online: \"", rest_url, "\"...")
  content <- kegg_fetch_text(rest_url)
  lines <- strsplit(content, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  mat <- do.call(rbind, strsplit(lines, "\t", fixed = TRUE))
  res <- data.frame(from = mat[, 1], to = mat[, 2], stringsAsFactors = FALSE)
  saveRDS(res, cache_file)
  res
}

kegg_call_retry <- function(fn, attempts = 3, timeout = 600) {
  last_error <- NULL
  for (i in seq_len(attempts)) {
    res <- tryCatch(
      {
        options(timeout = timeout)
        R.utils::withTimeout(fn(), timeout = timeout, onTimeout = "error")
      },
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(res)) return(res)
    if (i < attempts) {
      message(
        "KEGG attempt ", i, " failed (", conditionMessage(last_error),
        "); retrying"
      )
      Sys.sleep(3 * i)
    }
  }
  if (!is.null(last_error)) {
    message("KEGG failed after ", attempts, " attempts: ", conditionMessage(last_error))
  }
  NULL
}

suppressWarnings(
  tryCatch(
    utils::assignInNamespace(
      "kegg_rest",
      kegg_rest_cached,
      ns = "clusterProfiler"
    ),
    error = function(e) {
      message("KEGG cache patch unavailable: ", conditionMessage(e))
    }
  )
)
