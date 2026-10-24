# Offline ranking of marker gene pairs for every cluster of a dataset
# (replaces the legacy find_markers_h5.R brute-force search with one matrix
# product per cluster). Load with source("script/load.R", chdir = TRUE).

find_pairs <- function(clusters, out = NULL, cut = 0.5, required = character(0),
                       min_n = 1, top = 20) {
  #' Best marker pairs per cluster: data.frame with cluster (ID),
  #' cluster_name (NA when unnamed), gene1, gene2, offtargets (other clusters
  #' in which both genes are ON in any split) and single_offtargets; rows
  #' sorted best first per cluster, at most `top` per cluster.
  #' cut: probability above which a gene is ON.
  #' required, min_n: the pair must be ON in every `required` split and in at
  #'   least `min_n` splits.
  #' out: also write the table there as CSV (empty cluster_name = unnamed).
  ds <- ds_load("find_pairs", list(label = "find_pairs", clusters = clusters))
  b <- ds_binarize(ds, cut)
  results <- lapply(ds$clusters, function(cl) {
    res <- pair_offtargets(b, cl, required, as.integer(min_n),
                           limit = as.integer(top), genes = ds$analysis_genes)
    if (nrow(res) == 0) return(NULL)
    cbind(cluster = cl, cluster_name = unname(ds$cluster_names[cl]), res)
  })
  pairs <- do.call(rbind, results)
  if (!is.null(out)) {
    dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
    write.csv(pairs, out, row.names = FALSE)
    message(sprintf("Wrote %d pairs for %d clusters to %s",
                    if (is.null(pairs)) 0 else nrow(pairs),
                    length(unique(pairs$cluster)), out))
  }
  if (is.null(out)) pairs else invisible(pairs)
}
