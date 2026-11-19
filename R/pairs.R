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
#' @param band_width Width of the ambiguity band of [ds_binarize()]; the
#'   default 0 calls a gene ON above P(On) = 0.5 with no ambiguous band.
#' @param required,min_n Both genes must be ON in every `required` split and
#'   in at least `min_n` splits.
#' @param splits Optional splits to search, e.g. one stage: genes are ON and
#'   clusters are off-targets only by their calls there (see
#'   [restrict_splits()]); clusters absent from them are not ranked.
#' @param genes Optional genes to choose pairs from (e.g. genes with a driver
#'   line); default: every gene used for marker finding.
#' @param top Maximum number of pairs per cluster (`Inf` for all).
#' @param max_onon Optional number: per cluster, every pair with at most this
#'   many ON-ON off-targets, before `top` applies (see [pair_offtargets()]).
#' @param cut Deprecated: the former confidence margin, `band_width = 1 - 2 * cut`.
#' @param cores Number of clusters ranked at once, in forked processes
#'   ([parallel::mclapply()]); one on Windows, which cannot fork. The result
#'   does not depend on it.
#'
#' @return A data.frame with columns `cluster` (ID), `cluster_name` (`NA`
#'   when unnamed), `gene1`, `gene2`, `onon_offtargets` (other clusters in
#'   which both genes are ON), `ambiguous_offtargets` (other clusters in which
#'   neither is OFF and at least one is ambiguous), `offtargets` (their sum),
#'   `separation`,
#'   `target_fold`, `worst_leak`, `ashman_D_lower` (see [pair_offtargets()];
#'   `NA` for datasets without the binarizer 1.1 ranking data) and
#'   `single_offtargets`, best first per cluster; invisibly when `out` is
#'   given.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' pairs <- find_pairs(ex, top = 3)
#' head(pairs)
#' @export
find_pairs <- function(clusters, out = NULL, band_width = 0, required = character(0),
                       min_n = 1, top = 20, cores = 1, cut = NULL, splits = NULL,
                       genes = NULL, max_onon = NULL) {
    if (!is.null(cut)) {
        warning("`cut` is deprecated; use `band_width = 1 - 2 * cut` (here ", 1 - 2 * cut, ")",
            call. = FALSE)
        band_width <- 1 - 2 * cut
    }
    ds <- ds_load("find_pairs", list(label = "find_pairs", clusters = clusters))
    b <- ds_binarize(ds, band_width)
    if (!is.null(splits)) b <- restrict_splits(b, splits)
    use <- if (is.null(genes)) ds$analysis_genes else intersect(ds$analysis_genes, genes)
    limit <- if (is.finite(top)) as.integer(top) else Inf
    rank_cluster <- function(cl) {
        res <- pair_offtargets(b, cl, required, as.integer(min_n),
            limit = limit, genes = use, max_onon = max_onon)
        if (nrow(res) == 0) return(NULL)
        cbind(cluster = cl, cluster_name = unname(ds$cluster_names[cl]), res)
    }
    if (.Platform$OS.type == "windows") cores <- 1
    results <- if (cores > 1) {
        parallel::mclapply(b$clusters, rank_cluster, mc.cores = cores)
    } else {
        lapply(b$clusters, rank_cluster)
    }
    # mclapply() returns a worker's error instead of raising it
    failed <- vapply(results, inherits, logical(1), "try-error")
    if (any(failed)) {
        stop("Ranking pairs failed for cluster ", b$clusters[which(failed)[1]], ": ",
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
