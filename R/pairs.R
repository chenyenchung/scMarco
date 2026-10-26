# Offline ranking of marker gene pairs for every cluster of a dataset
# (replaces the legacy find_markers_h5.R brute-force search with one matrix
# product per cluster).

#' Rank marker pairs for every cluster
#'
#' Runs [pair_offtargets()] for every cluster of a dataset, e.g. to screen a
#' whole atlas offline.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param out Optional path: also write the table there as CSV (empty
#'   `cluster_name` = unnamed).
#' @param cut Probability above which a gene is ON (the confidence margin of
#'   [ds_binarize()]).
#' @param required,min_n Both genes must be ON in every `required` split and
#'   in at least `min_n` splits.
#' @param top Maximum number of pairs per cluster.
#'
#' @return A data.frame with columns `cluster` (ID), `cluster_name` (`NA`
#'   when unnamed), `gene1`, `gene2`, `offtargets` (other clusters in which
#'   both genes are ON in any split) and `single_offtargets`, best first per
#'   cluster; invisibly when `out` is given.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' pairs <- find_pairs(ex, top = 3)
#' head(pairs)
#' @export
find_pairs <- function(clusters, out = NULL, cut = 0.5, required = character(0),
                       min_n = 1, top = 20) {
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
