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
#' @param band Ambiguity band of [ds_binarize()]; the default 0 calls a gene
#'   ON above P(On) = 0.5 with no ambiguous band.
#' @param required,min_n Both genes must be ON in every `required` split and
#'   in at least `min_n` splits.
#' @param top Maximum number of pairs per cluster.
#' @param cut Deprecated: the former confidence margin, `0.5 - band`.
#' @param cores Number of clusters ranked at once, in forked processes
#'   ([parallel::mclapply()]); one on Windows, which cannot fork. The result
#'   does not depend on it.
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
find_pairs <- function(clusters, out = NULL, band = 0, required = character(0),
                       min_n = 1, top = 20, cores = 1, cut = NULL) {
    if (!is.null(cut)) {
        warning("`cut` is deprecated; use `band = 0.5 - cut` (here ", 0.5 - cut, ")",
            call. = FALSE)
        band <- 0.5 - cut
    }
    ds <- ds_load("find_pairs", list(label = "find_pairs", clusters = clusters))
    b <- ds_binarize(ds, band)
    rank_cluster <- function(cl) {
        res <- pair_offtargets(b, cl, required, as.integer(min_n),
            limit = as.integer(top), genes = ds$analysis_genes)
        if (nrow(res) == 0) return(NULL)
        cbind(cluster = cl, cluster_name = unname(ds$cluster_names[cl]), res)
    }
    if (.Platform$OS.type == "windows") cores <- 1
    results <- if (cores > 1) {
        parallel::mclapply(ds$clusters, rank_cluster, mc.cores = cores)
    } else {
        lapply(ds$clusters, rank_cluster)
    }
    # mclapply() returns a worker's error instead of raising it
    failed <- vapply(results, inherits, logical(1), "try-error")
    if (any(failed)) {
        stop("Ranking pairs failed for cluster ", ds$clusters[which(failed)[1]], ": ",
            conditionMessage(attr(results[[which(failed)[1]]], "condition")), call. = FALSE)
    }
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
