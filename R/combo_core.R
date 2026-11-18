# Combination-finding core
#
# Pure functions (no Shiny, no I/O) shared by every tab. They operate on a
# "binarized dataset" produced by binarize():
#
#   on[gene, cluster, split]  logical  confidently ON: prob >= 0.5 + band
#                                      (FALSE when the cluster is absent from
#                                      the split or the gene is non-bimodal in
#                                      that split)
#   ambiguous[gene, cluster, split]
#                             logical  0.5 - band < prob < 0.5 + band: neither
#                                      confidently ON nor confidently OFF; or
#                                      a call that the binarizer 1.1 flags
#                                      (chains disagree, ON below the detection
#                                      floor, OFF above ambient RNA)
#   present[cluster, split]   logical  the cluster exists in the split
#   bimodal[gene, split]      logical  the mixture model called the gene bimodal
#
# Semantics shared by all functions (ambiguous calls behave like NA, so they
# never help a gene):
#   - A gene "passes" in a cluster when it is ON in every `required` split and
#     ON in at least `min_n` splits overall; ambiguous is not ON.
#   - Off-target expression is counted if the gene is ON or ambiguous in ANY
#     split of the other cluster (a conservative choice kept from both legacy
#     variants).

#' Binarize mixture-model probabilities
#'
#' Turns per-cluster probabilities of expression into ON, OFF and ambiguous
#' calls. Most users get this from [ds_binarize()], which reads the inputs
#' from a dataset and caches the result.
#'
#' A gene is ON in a cluster and split when its probability is at least
#' `0.5 + band`, OFF when it is at most `0.5 - band`, and ambiguous in
#' between (the ambiguity band).
#' Ambiguous calls never count as ON, but they do count as off-target
#' expression. A gene the mixture model did not call bimodal in a split is
#' neither ON nor ambiguous there.
#'
#' With `flags` (binarizer 1.1, from [call_mixture()]), calls the binarizer
#' does not trust become ambiguous whatever the band: any call whose chains
#' disagree (flag 1), an ON call below the detection floor (flag 2) and an
#' OFF call above ambient RNA (flag 4).
#'
#' @param prob Numeric matrix, one row per split x cluster pair that exists
#'   and one column per gene (column names are the genes). `NA` is OFF.
#' @param obs_split,obs_cluster Character vectors: the split and cluster of
#'   each row of `prob`.
#' @param splits,clusters The full, ordered sets of splits and clusters.
#' @param band Half-width of the ambiguity band around P(On) = 0.5, in
#'   \[0, 0.5). At 0 there is no ambiguous band and ON is `prob > 0.5`.
#'   The default elsewhere in scMarco is 0.3 (ON at 0.8, OFF at 0.2).
#' @param bimodal Optional logical matrix, genes x splits: whether the
#'   mixture model called the gene bimodal in that split.
#' @param flags Optional integer matrix shaped like `prob`: the flags of
#'   [call_mixture()] (the sum of 1 chains disagree, 2 below the detection
#'   floor, 4 above ambient RNA), as stored in `layers/call_flags`.
#' @param fold Optional numeric matrix shaped like `prob`: each entry's log2
#'   fold over ambient RNA (`layers/fold_ambient`), used by
#'   [pair_offtargets()] to rank pairs.
#' @param ashman_D_lower Optional numeric matrix, genes x splits: the lower
#'   bound of each gene's Ashman's D (`varm/call_ashman_D_lower`), used by
#'   [pair_offtargets()] to break ties.
#'
#' @return An object of class `scmarco_binarized`: a list with the logical
#'   arrays `on` and `ambiguous` (genes x clusters x splits), the logical
#'   matrices `present` (clusters x splits), `bimodal` (genes x splits) and
#'   `anyon` (genes x clusters: ON or ambiguous in any split, see [any_on()]),
#'   `genes`, `clusters`, `splits`, `band`, `flags_applied` (whether
#'   `flags` were given), and `fold` (genes x clusters x splits, `NULL`
#'   without `fold`) and `ashman_D_lower` as given.
#' @param margin Deprecated: the former confidence margin, `0.5 - band`.
#' @examples
#' prob <- matrix(c(0.95, 0.5, 0.02, 0.99), 2,
#'     dimnames = list(NULL, c("geneA", "geneB")))
#' b <- binarize(prob, obs_split = c("adult", "adult"),
#'     obs_cluster = c("c1", "c2"), splits = "adult",
#'     clusters = c("c1", "c2"), band = 0.4)
#' b$on[, , "adult"]
#' b$ambiguous[, , "adult"]
#' @export
binarize <- function(prob, obs_split, obs_cluster, splits, clusters,
                     band, bimodal = NULL, flags = NULL, fold = NULL,
                     ashman_D_lower = NULL, margin = NULL) {
    band <- band_from_margin(if (missing(band)) NULL else band, margin)
    stopifnot(
        nrow(prob) == length(obs_split),
        length(obs_split) == length(obs_cluster),
        band >= 0, band < 0.5
    )
    genes <- colnames(prob)
    present <- matrix(
        FALSE, length(clusters), length(splits),
        dimnames = list(clusters, splits)
    )

    ci <- match(obs_cluster, clusters)
    si <- match(obs_split, splits)
    if (anyNA(ci) || anyNA(si)) {
        stop("obs contains clusters or splits missing from the level sets")
    }
    present[cbind(ci, si)] <- TRUE

    # tolerance so slider steps such as 1 - 0.3 land on 0.7
    eps <- 1e-9
    off <- is.na(prob) | prob <= 0.5 - band + eps
    hit <- !off & prob >= 0.5 + band - eps
    unsure <- !off & !hit
    if (!is.null(flags)) {
        stopifnot(identical(dim(flags), dim(prob)))
        has <- function(bit) !is.na(flags) & (flags %/% bit) %% 2 == 1
        chains <- has(CALL_FLAG_CHAINS)
        # NA probabilities (no call) stay OFF
        to_unsure <- (hit & (chains | has(CALL_FLAG_DETECTION))) |
            (off & !is.na(prob) & (chains | has(CALL_FLAG_AMBIENT)))
        hit <- hit & !to_unsure
        off <- off & !to_unsure
        unsure <- unsure | to_unsure
    }

    # Fill genes x (clusters * splits) matrices, whose column (si - 1) *
    # n_clusters + ci is slice [, ci, si] of the final array, with one
    # assignment each
    n_cl <- length(clusters)
    col <- (si - 1L) * n_cl + ci
    on <- matrix(FALSE, length(genes), n_cl * length(splits))
    ambiguous <- on
    on[, col] <- t(hit)
    ambiguous[, col] <- t(unsure)

    if (is.null(bimodal)) {
        bimodal <- matrix(
            TRUE, length(genes), length(splits),
            dimnames = list(genes, splits)
        )
    } else {
        bimodal <- bimodal[genes, splits, drop = FALSE]
        # a gene non-bimodal in a split is neither ON nor ambiguous there
        by_col <- bimodal[, rep(seq_along(splits), each = n_cl), drop = FALSE]
        on <- on & by_col
        ambiguous <- ambiguous & by_col
    }
    dim(on) <- dim(ambiguous) <- c(length(genes), n_cl, length(splits))
    dimnames(on) <- dimnames(ambiguous) <- list(genes, clusters, splits)
    fold_arr <- NULL
    if (!is.null(fold)) {
        stopifnot(identical(dim(fold), dim(prob)))
        fold_arr <- matrix(NA_real_, length(genes), n_cl * length(splits))
        fold_arr[, col] <- t(fold)
        dim(fold_arr) <- c(length(genes), n_cl, length(splits))
        dimnames(fold_arr) <- list(genes, clusters, splits)
    }

    structure(
        list(
            on = on, ambiguous = ambiguous, present = present, bimodal = bimodal,
            anyon = anyon_matrix(on, ambiguous),
            genes = genes, clusters = clusters, splits = splits, band = band,
            flags_applied = !is.null(flags), fold = fold_arr,
            ashman_D_lower = ashman_D_lower
        ),
        class = "scmarco_binarized"
    )
}

band_from_margin <- function(band, margin) {
    # `margin` (the former confidence margin) is deprecated: band = 0.5 - margin.
    # band: NULL when the caller did not give it
    if (!is.null(margin)) {
        if (!is.null(band)) stop("give band or margin, not both", call. = FALSE)
        warning("`margin` is deprecated; use `band = 0.5 - margin` (here ", 0.5 - margin, ")",
            call. = FALSE)
        return(0.5 - margin)
    }
    band
}

#' Limit a binarized dataset to some clusters
#'
#' Every search on the result ignores the clusters left out, both as targets
#' and as off-targets.
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param keep Cluster IDs to keep; the dataset's order is kept.
#'
#' @return A binarized dataset with only the clusters in `keep`.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' sub <- restrict_clusters(b, b$clusters[1:20])
#' length(sub$clusters)
#' @export
restrict_clusters <- function(b, keep) {
    keep <- b$clusters[b$clusters %in% keep]
    if (identical(keep, b$clusters)) return(b)
    b$on <- b$on[, keep, , drop = FALSE]
    b$ambiguous <- b$ambiguous[, keep, , drop = FALSE]
    b$present <- b$present[keep, , drop = FALSE]
    if (!is.null(b$anyon)) b$anyon <- b$anyon[, keep, drop = FALSE]
    if (!is.null(b$fold)) b$fold <- b$fold[, keep, , drop = FALSE]
    b$clusters <- keep
    b
}

check_split_args <- function(b, required, min_n) {
    unknown <- setdiff(required, b$splits)
    if (length(unknown) > 0) {
        stop("Unknown split(s): ", paste(unknown, collapse = ", "))
    }
    if (min_n < 1 || min_n > length(b$splits)) {
        stop("min_n must be between 1 and the number of splits")
    }
}

# gene x cluster logical: gene passes the split criteria in that cluster
pass_matrix <- function(b, required, min_n, clusters = b$clusters,
                        genes = b$genes) {
    check_split_args(b, required, min_n)
    on <- b$on[genes, clusters, , drop = FALSE]
    n_on <- rowSums(on, dims = 2)
    pass <- n_on >= min_n
    if (length(required) > 0) {
        n_req <- rowSums(on[, , required, drop = FALSE], dims = 2)
        pass <- pass & n_req == length(required)
    }
    dim(pass) <- c(length(genes), length(clusters))
    dimnames(pass) <- list(genes, clusters)
    pass
}

# gene x cluster logical: gene is ON, or ambiguous, in any split
#' Which genes are ON or ambiguous in which clusters
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param genes Genes (default: every gene of `b`).
#' @param clusters Clusters (default: every cluster of `b`).
#'
#' @return Logical matrix, genes x clusters: `TRUE` where the gene is ON, or
#'   ambiguous, in any split of the cluster. This is what counts as
#'   off-target expression.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' any_on(b, genes = b$genes[1:3], clusters = b$clusters[1:4])
#' @export
any_on <- function(b, genes = b$genes, clusters = b$clusters) {
    # binarize() computes it once for every gene and cluster; objects made
    # without it (e.g. kept from an older scMarco) get it computed here
    if (!is.null(b$anyon)) {
        return(b$anyon[genes, clusters, drop = FALSE])
    }
    anyon_matrix(b$on[genes, clusters, , drop = FALSE],
        b$ambiguous[genes, clusters, , drop = FALSE])
}

anyon_matrix <- function(on, ambiguous) {
    # genes x clusters logical from genes x clusters x splits arrays: ON or
    # ambiguous in any split
    out <- rowSums(on | ambiguous, dims = 2) > 0
    dim(out) <- dim(on)[1:2]
    dimnames(out) <- dimnames(on)[1:2]
    out
}

target_genes <- function(b, ident, required, min_n, genes = b$genes) {
    # Genes that pass the split criteria in the single cluster `ident`
    if (length(ident) != 1 || !ident %in% b$clusters) {
        stop("ident must be a single known cluster")
    }
    genes <- intersect(genes, b$genes)
    pass <- pass_matrix(b, required, min_n, clusters = ident, genes = genes)
    genes[pass[, 1]]
}

offtarget_counts <- function(b, genes, exclude) {
    # Number of other clusters in which each gene is ON (or ambiguous) in any
    # split
    others <- setdiff(b$clusters, exclude)
    if (length(genes) == 0) {
        return(setNames(integer(0), character(0)))
    }
    counts <- rowSums(any_on(b, genes = genes, clusters = others))
    setNames(as.integer(counts), genes)
}

#' Candidate markers of a cluster
#'
#' The "From cluster" search: genes that pass the split criteria in one
#' cluster, ranked by the number of other clusters that also express them.
#'
#' A gene passes in a cluster when it is ON in every `required` split and in
#' at least `min_n` splits. A cluster counts as an off-target when the gene
#' is ON or ambiguous in any of its splits.
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param ident One cluster ID.
#' @param required Character vector of splits in which a gene must be ON
#'   (may be empty).
#' @param min_n Minimum number of splits in which a gene must be ON.
#' @param limit Maximum number of genes returned.
#' @param genes Genes to consider (default: every gene of `b`).
#'
#' @return Named integer vector: the number of off-target clusters of each
#'   candidate gene, fewest first (ties alphabetical).
#' @seealso [pair_offtargets()] to rank gene pairs, [distinct_within()] to
#'   tell the off-targets apart.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' from_cluster(b, ident = b$clusters[1], required = character(0), min_n = 1)
#' @export
from_cluster <- function(b, ident, required, min_n, limit = 10,
                         genes = b$genes) {
    cand <- target_genes(b, ident, required, min_n, genes = genes)
    counts <- offtarget_counts(b, cand, exclude = ident)
    counts <- counts[order(counts, names(counts))]
    head(counts, limit)
}

#' Clusters that express a set of genes
#'
#' The "From gene" search.
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param genes Gene names.
#' @param required Character vector of splits in which a gene must be ON
#'   (may be empty).
#' @param min_n Minimum number of splits in which a gene must be ON.
#' @param mode `"all"`: clusters in which every gene passes; `"any"`:
#'   clusters in which at least one does.
#'
#' @return Character vector of cluster IDs, in dataset order.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' clusters_expressing(b, b$genes[1:2], required = character(0), min_n = 1,
#'     mode = "any")
#' @export
clusters_expressing <- function(b, genes, required, min_n,
                                mode = c("all", "any")) {
    mode <- match.arg(mode)
    genes <- intersect(genes, b$genes)
    if (length(genes) == 0) {
        return(character(0))
    }
    pass <- pass_matrix(b, required, min_n, genes = genes)
    keep <- if (mode == "all") {
        colSums(pass) == length(genes)
    } else {
        colSums(pass) > 0
    }
    b$clusters[keep]
}

#' Genes that tell a set of clusters apart
#'
#' The "Find distinct genes within" search: genes that pass the split
#' criteria in exactly `k` of `idents` and are neither ON nor ambiguous (in
#' any split) in the other clusters of `idents`.
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param idents At least two cluster IDs.
#' @param required Character vector of splits in which a gene must be ON
#'   (may be empty).
#' @param min_n Minimum number of splits in which a gene must be ON.
#' @param k Number of clusters of `idents` a gene may mark, from 1 to
#'   `length(idents) - 1`.
#' @param genes Genes to consider (default: every gene of `b`).
#'
#' @return A data.frame with columns `cluster`, `gene` and `offtargets` (the
#'   number of clusters outside `idents` that also express the gene), best
#'   first within each cluster; zero rows when nothing is found.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' hits <- distinct_within(b, b$clusters[1:3], required = character(0),
#'     min_n = 1)
#' head(hits)
#' @export
distinct_within <- function(b, idents, required, min_n, k = 1,
                            genes = b$genes) {
    idents <- intersect(idents, b$clusters)
    if (length(idents) < 2) {
        stop("Need at least 2 clusters")
    }
    if (k < 1 || k >= length(idents)) {
        stop("k must be between 1 and the number of clusters - 1")
    }
    genes <- intersect(genes, b$genes)
    pass <- pass_matrix(b, required, min_n, clusters = idents, genes = genes)
    anyon <- any_on(b, genes = genes, clusters = idents)
    n_pass <- rowSums(pass)
    keep <- n_pass == k & rowSums(anyon) == n_pass
    hits <- which(pass[keep, , drop = FALSE], arr.ind = TRUE)
    if (nrow(hits) == 0) {
        return(data.frame(
            cluster = character(0), gene = character(0), offtargets = integer(0)
        ))
    }
    kept_genes <- genes[keep]
    out <- data.frame(
        cluster = idents[hits[, "col"]],
        gene = kept_genes[hits[, "row"]]
    )
    out$offtargets <- unname(offtarget_counts(b, out$gene, exclude = idents))
    out <- out[order(match(out$cluster, idents), out$offtargets, out$gene), ,
        drop = FALSE]
    rownames(out) <- NULL
    out
}

top_per_cluster <- function(tbl, n) {
    # First n rows per cluster of a distinct_within() table
    keep <- stats::ave(seq_len(nrow(tbl)), tbl$cluster, FUN = seq_along) <= n
    tbl[keep, , drop = FALSE]
}

top_per_cluster_capped <- function(tbl, n, max_genes) {
    # Like top_per_cluster(), but lowers n until at most `max_genes` distinct
    # genes are kept. Never goes below one gene per cluster, so every cluster
    # with a hit stays represented even if that exceeds `max_genes`.
    # Returns list(shown, n) with the n actually used.
    n <- max(1L, as.integer(n))
    repeat {
        shown <- top_per_cluster(tbl, n)
        if (n == 1L || length(unique(shown$gene)) <= max_genes) break
        n <- n - 1L
    }
    list(shown = shown, n = n)
}

empty_pairs <- function() {
    data.frame(gene1 = character(0), gene2 = character(0),
        offtargets = integer(0), separation = numeric(0), target_fold = numeric(0),
        worst_leak = numeric(0), ashman_D_lower = numeric(0),
        single_offtargets = integer(0))
}

pair_separation <- function(b, ident, cand, idx, chunk = 5000) {
    # Separation of each pair idx (rows: indices into cand) in cluster ident,
    # in log2 fold over ambient RNA (b$fold): target strength (the weaker
    # gene's fold, lowest over the splits where both genes are ON in the
    # target) minus the worst leak (over every other cluster and split where
    # the pair is not ON or ambiguous together, i.e. at least one gene is
    # OFF: the larger of the pair's smaller folds); and the weaker gene's
    # Ashman's D lower bound (lowest over splits). NA where the data is
    # missing.
    n <- nrow(idx)
    out <- data.frame(separation = rep(NA_real_, n), target_fold = NA_real_,
        worst_leak = NA_real_, ashman_D_lower = NA_real_)
    if (n == 0) return(out)
    if (!is.null(b$ashman_D_lower)) {
        d <- b$ashman_D_lower[cand, intersect(b$splits, colnames(b$ashman_D_lower)), drop = FALSE]
        d[is.na(d)] <- Inf
        d_gene <- apply(d, 1, min)
        d_gene[!is.finite(d_gene)] <- NA
        out$ashman_D_lower <- pmin(d_gene[idx[, 1]], d_gene[idx[, 2]])
    }
    if (is.null(b$fold)) return(out)
    n_cl <- length(b$clusters)
    n_col <- n_cl * length(b$splits)
    flat <- function(a) matrix(a[cand, b$clusters, b$splits, drop = FALSE], length(cand), n_col)
    fold <- flat(b$fold)
    on <- flat(b$on)
    nonoff <- on | flat(b$ambiguous)
    ti <- match(ident, b$clusters)
    target <- (seq_along(b$splits) - 1L) * n_cl + ti
    leak_col <- as.vector(b$present) & !seq_len(n_col) %in% target
    row_max <- function(m) {
        m[is.na(m)] <- -Inf
        v <- m[cbind(seq_len(nrow(m)), max.col(m, ties.method = "first"))]
        v[!is.finite(v)] <- NA
        v
    }
    for (start in seq(1, n, by = chunk)) {
        k <- start:min(n, start + chunk - 1L)
        g1 <- idx[k, 1]
        g2 <- idx[k, 2]
        m <- pmin(fold[g1, , drop = FALSE], fold[g2, , drop = FALSE])
        lk <- m
        lk[(nonoff[g1, , drop = FALSE] & nonoff[g2, , drop = FALSE]) |
            rep(!leak_col, each = length(k))] <- NA
        tg <- m[, target, drop = FALSE]
        tg[!(on[g1, target, drop = FALSE] & on[g2, target, drop = FALSE])] <- NA
        out$target_fold[k] <- -row_max(-tg)
        out$worst_leak[k] <- row_max(lk)
    }
    out$separation <- out$target_fold - out$worst_leak
    out
}

#' Rank gene pairs for a cluster
#'
#' Ranks pairs of genes that both pass the split criteria in `ident` by the
#' number of other clusters in which both genes are ON or ambiguous (in any
#' split): the off-targets of an intersectional (e.g. split-GAL4)
#' combination.
#'
#' Ties are broken by the binarizer 1.1 ranking data when `b` has it (from
#' [ds_binarize()] of a dataset binarized with [call_mixture()] tables):
#' first by **separation**, the pair's target strength minus its worst leak,
#' in log2 fold over ambient RNA (log2((UMIs + 1) / (expected ambient UMIs +
#' 1))): the target strength is the weaker gene's fold in `ident` (lowest
#' over the splits where both genes are ON there), and the worst leak is,
#' over the other clusters and splits where at least one of the two genes is
#' OFF, the largest of the pair's smaller folds (an intersectional driver
#' needs both genes). Then by the weaker gene's lower bound of Ashman's D
#' (10th percentile over the posterior), and last by the single-gene
#' off-targets. A leak of about 1 or less is a single stray UMI.
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param ident One cluster ID.
#' @param required Character vector of splits in which a gene must be ON
#'   (may be empty).
#' @param min_n Minimum number of splits in which a gene must be ON.
#' @param first Optional gene: only pairs containing it are returned.
#' @param limit Maximum number of pairs returned.
#' @param genes Genes to consider (default: every gene of `b`).
#'
#' @return A data.frame with columns `gene1`, `gene2`, `offtargets`,
#'   `separation`, `target_fold`, `worst_leak`, `ashman_D_lower` (`NA`
#'   without the ranking data) and `single_offtargets` (the off-targets of
#'   each gene on its own, summed), best first.
#' @seealso [find_pairs()] to rank pairs for every cluster of a dataset.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' pair_offtargets(b, ident = b$clusters[1], required = character(0),
#'     min_n = 1, limit = 5)
#' @export
pair_offtargets <- function(b, ident, required, min_n, first = NULL,
                            limit = 50, genes = b$genes) {
    cand <- target_genes(b, ident, required, min_n, genes = genes)
    if (!is.null(first)) {
        if (!first %in% cand) {
            return(empty_pairs())
        }
        cand <- c(first, setdiff(cand, first))
    }
    if (length(cand) < 2) {
        return(empty_pairs())
    }
    others <- setdiff(b$clusters, ident)
    a <- any_on(b, genes = cand, clusters = others) * 1L
    shared <- tcrossprod(a)

    # Keep only pairs that can make the top `limit` before sorting: there can be
    # millions of pairs, so drop the others before indexing them. Ties break
    # on how specific each gene is on its own (sum of single-gene
    # off-targets), then on gene position (alphabetical).
    if (!is.null(first)) {
        idx <- cbind(1L, seq(2, length(cand)))
        idx <- idx[can_rank(shared[idx], limit), , drop = FALSE]
    } else {
        pair <- upper.tri(shared)
        pair[pair] <- can_rank(shared[pair], limit)
        idx <- which(pair, arr.ind = TRUE)
    }
    vals <- as.integer(shared[idx])
    single <- diag(shared)
    single_sum <- as.integer(single[idx[, 1]] + single[idx[, 2]])
    sep <- pair_separation(b, ident, cand, idx)
    # missing separation or D sort last among pairs with as many off-targets
    ord <- head(order(vals, -sep$separation, -sep$ashman_D_lower, single_sum,
        idx[, 1], idx[, 2], method = "radix", na.last = TRUE), limit)
    round2 <- function(x) round(x, 2)
    data.frame(
        gene1 = cand[idx[ord, 1]],
        gene2 = cand[idx[ord, 2]],
        offtargets = vals[ord],
        separation = round2(sep$separation[ord]),
        target_fold = round2(sep$target_fold[ord]),
        worst_leak = round2(sep$worst_leak[ord]),
        ashman_D_lower = round2(sep$ashman_D_lower[ord]),
        single_offtargets = single_sum[ord]
    )
}

can_rank <- function(vals, limit) {
    # Which values can be among the `limit` smallest (all values tied with
    # the limit-th smallest are kept)
    if (length(vals) <= limit) return(rep(TRUE, length(vals)))
    vals <= sort(vals, partial = limit)[limit]
}

#' Clusters labeled by a gene combination
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param genes Gene names.
#' @param exclude Cluster IDs to leave out (typically the target).
#'
#' @return Character vector of the clusters, other than `exclude`, in which
#'   every gene of `genes` is ON or ambiguous in some split.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' combo_offtarget_clusters(b, b$genes[1:2])
#' @export
combo_offtarget_clusters <- function(b, genes, exclude = character(0)) {
    genes <- intersect(genes, b$genes)
    if (length(genes) == 0) {
        return(character(0))
    }
    a <- any_on(b, genes = genes, clusters = setdiff(b$clusters, exclude))
    colnames(a)[colSums(a) == length(genes)]
}

#' Per-gene ON/OFF states in long form
#'
#' @param b A binarized dataset from [ds_binarize()] or [binarize()].
#' @param genes Gene names.
#' @param clusters Cluster IDs (default: every cluster of `b`).
#'
#' @return A data.frame with columns `gene`, `cluster`, `split` and `state`,
#'   one of `"Expressed"`, `"Not expressed"`, `"Ambiguous"`, `"Non-bimodal"`
#'   or `"Absent"` (the cluster has no cells in the split).
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, band = 0.4)
#' head(state_table(b, b$genes[1], clusters = b$clusters[1:2]))
#' @export
state_table <- function(b, genes, clusters = b$clusters) {
    genes <- intersect(genes, b$genes)
    grid <- expand.grid(
        gene = genes, cluster = clusters, split = b$splits,
        stringsAsFactors = FALSE
    )
    gi <- match(grid$gene, b$genes)
    ci <- match(grid$cluster, b$clusters)
    si <- match(grid$split, b$splits)
    on <- b$on[cbind(gi, ci, si)]
    unsure <- b$ambiguous[cbind(gi, ci, si)]
    present <- b$present[cbind(ci, si)]
    bimodal <- b$bimodal[cbind(gi, si)]
    grid$state <- ifelse(
        !present, "Absent",
        ifelse(!bimodal, "Non-bimodal",
            ifelse(unsure, "Ambiguous",
                ifelse(on, "Expressed", "Not expressed")))
    )
    grid
}
