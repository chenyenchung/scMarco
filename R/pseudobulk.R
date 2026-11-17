# Pseudobulk samples for the mixture model. make_pseudobulk() pools each
# cluster's cells into samples as Seurat::AggregateExpression() does (summed
# counts, then LogNormalize): random parts (sampling without replacement,
# which assumes the cluster is homogeneous) or the values of a cell metadata
# column.

with_seed <- function(seed, code) {
    # Evaluate `code` after set.seed(seed), leaving the caller's random state
    # as it was; seed NULL uses the current random state
    if (is.null(seed)) return(code)
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
    } else {
        rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
    code
}

#' Pool cells into pseudobulk samples for the mixture model
#'
#' Makes the input of [fit_mixture()] for one split, as
#' `Seurat::AggregateExpression(return.seurat = TRUE)` does: each cluster's
#' cells are pooled into up to `n_samples` samples, the raw counts of each
#' sample are summed, and each sample is then log-normalized like a cell,
#' `log1p(summed counts / sample total * 1e4)`.
#'
#' With `method = "random"` the cells of a cluster are shuffled and dealt
#' into samples of (nearly) equal size, without replacement. This assumes
#' the cluster is homogeneous, and gives the model several samples per
#' cluster when there are few clusters.
#'
#' With `method = "metadata"` each cluster is split by the values of
#' `group`, a cell metadata column (e.g. library, batch, genotype): each
#' value is one sample, and also its replicate, the same in every cluster
#' (in factor level order, else sorted).
#' To pool libraries into two replicates, as the mixture_modeling_smk
#' pipeline did, map them to two groups first.
#'
#' A cluster gets fewer random parts when it has fewer than `min_cells`
#' cells per part (one sample below `2 * min_cells` cells). With metadata, a
#' value with fewer than `min_cells` cells in a cluster is left out of it
#' (those cells are in no sample).
#'
#' @param counts Raw counts of one split's cells, genes x cells (a matrix or
#'   sparse `Matrix`, gene row names), e.g. the `counts` layer of a Seurat
#'   object.
#' @param cluster Cluster ID of each cell (column of `counts`).
#' @param n_samples Samples per cluster for `"random"`.
#' @param method `"random"` or `"metadata"` (see Details).
#' @param group For `"metadata"`, the value of each cell that splits the
#'   clusters (e.g. `obj$orig.ident`).
#' @param min_cells Fewest cells per sample (see Details).
#' @param seed Random seed for the shuffle; the caller's random
#'   state is restored afterwards. `NULL` uses the current random state.
#'
#' @return A list: `lognorm` (genes x samples), `cluster` and `replicate` of
#'   each sample (the `group` value, or the sample's number within its
#'   cluster), `n_cells` per sample, `sample` (the
#'   sample of each cell, `NA` for cells left out) and `method`. Pass `lognorm`, `cluster` and
#'   `replicate` to [fit_mixture()].
#' @seealso [fit_mixture()].
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(40 * 60, 2), 40, 60,
#'     dimnames = list(paste0("g", 1:40), paste0("c", 1:60)))
#' cluster <- rep(c("1", "2"), each = 30)
#' pb <- make_pseudobulk(counts, cluster, n_samples = 3, min_cells = 5)
#' pb$cluster
#' pb$n_cells
#' # one sample per library and cluster
#' pb <- make_pseudobulk(counts, cluster, method = "metadata",
#'     group = rep(c("L1", "L2"), 30))
#' pb$replicate
#' @export
make_pseudobulk <- function(counts, cluster, n_samples = 2,
                            method = c("random", "metadata"), group = NULL,
                            min_cells = 10, seed = 42) {
    method <- match.arg(method)
    counts <- as_count_matrix(counts)
    n <- ncol(counts)
    if (length(cluster) != n) fail("cluster has ", length(cluster), " entries for ", n, " cells")
    if (anyNA(cluster)) fail("cluster has missing values")
    if (n_samples < 1 || min_cells < 1) fail("n_samples and min_cells must be at least 1")
    if (method == "metadata") {
        if (is.null(group)) fail("method = \"metadata\" needs group, a value per cell")
        if (length(group) != n) fail("group has ", length(group), " entries for ", n, " cells")
        if (anyNA(group)) fail("group has missing values")
        # samples in factor level order, else sorted
        group_levels <- if (is.factor(group)) levels(group) else sort(unique(as.character(group)))
        group <- as.character(group)
    } else if (!is.null(group)) {
        fail("group is used by method = \"metadata\" only")
    }
    cluster <- as.character(cluster)
    clusters <- unique(cluster)

    # sample of each cell within its cluster (NA: in none)
    groups <- with_seed(seed, {
        lapply(clusters, function(cl) {
            idx <- which(cluster == cl)
            if (method == "metadata") {
                g <- group[idx]
                g[g %in% names(which(table(g) < min_cells))] <- NA
                return(g)
            }
            k <- min(n_samples, max(1, length(idx) %/% min_cells))
            if (k == 1) return(rep(1L, length(idx)))
            rep_len(seq_len(k), length(idx))[sample(length(idx))]
        })
    })

    cell_sample <- character(n)
    replicate <- character(0)
    sample_cluster <- character(0)
    for (i in seq_along(clusters)) {
        idx <- which(cluster == clusters[i])
        present <- if (method == "metadata") {
            intersect(group_levels, groups[[i]])
        } else {
            sort(unique(groups[[i]]))
        }
        cell_sample[idx] <- ifelse(is.na(groups[[i]]), NA, paste0(clusters[i], "_", groups[[i]]))
        sample_cluster <- c(sample_cluster, rep(clusters[i], length(present)))
        replicate <- c(replicate, as.character(present))
    }
    samples <- paste0(sample_cluster, "_", replicate)

    dropped <- sum(is.na(cell_sample))
    if (dropped > 0) {
        message(dropped, " cell(s) are in a group with fewer than ", min_cells,
            " cells in their cluster and are left out")
    }

    # summed counts of each sample, then LogNormalize per sample
    j <- match(cell_sample, samples)
    kept <- !is.na(j)
    pool <- Matrix::sparseMatrix(i = which(kept), j = j[kept], x = 1,
        dims = c(n, length(samples)))
    summed <- as.matrix(counts %*% pool)
    out <- log1p(t(t(summed) / pmax(colSums(summed), 1)) * 1e4)
    dimnames(out) <- list(rownames(counts), samples)

    list(lognorm = out, cluster = sample_cluster, replicate = replicate,
        n_cells = setNames(tabulate(j, length(samples)), samples), sample = cell_sample,
        method = method)
}
