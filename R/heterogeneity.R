# Pseudobulk samples for the mixture model and per-gene cluster
# heterogeneity. make_pseudobulk() pools each cluster's cells into samples
# as Seurat::AggregateExpression() does (summed counts, then LogNormalize):
# random parts (sampling without replacement, which assumes the cluster is
# homogeneous), k-means subclusters, or the values of a cell metadata
# column. mixture_heterogeneity() then asks, gene by gene, whether the
# mixture model's calls disagree between the samples of a cluster more than
# the calls of genes expressed at a similar level (in any cluster of the
# split) do. Genes that do are heterogeneous in that cluster:
# add_probability() stores the scores, and binarize() marks those genes
# ambiguous there.

# A reference bin of gene x cluster pairs at a similar expression level needs
# at least this many pairs; splits with few scored pairs get fewer, wider bins
HET_MIN_BIN_GENES <- 10
# Layers holding the scores in a clusters file: z-score, spread of the calls
# and the flag (1 heterogeneous, 0 not, NaN not scored)
HET_LAYERS <- c("het_z", "het_spread", "het_flag")
# Floor of the reference SD, so that a bin whose genes all agree perfectly
# does not turn a small spread into an infinite z-score
HET_SD_FLOOR <- 0.01

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

cluster_pcs <- function(x, n_features, n_pcs) {
    # Principal components of one cluster's cells: x is genes x cells
    # (log-normalized), reduced to its n_features most variable genes,
    # centered and scaled. Returns cells x PCs.
    mu <- Matrix::rowMeans(x)
    v <- Matrix::rowMeans(x^2) - mu^2
    top <- head(order(v, decreasing = TRUE), n_features)
    top <- top[v[top] > 1e-8]
    if (length(top) == 0) return(NULL)
    m <- t(as.matrix(x[top, , drop = FALSE]))
    rank <- min(n_pcs, ncol(m), nrow(m) - 1)
    stats::prcomp(m, center = TRUE, scale. = TRUE, rank. = rank)$x
}

kmeans_groups <- function(coords, k, min_cells) {
    # Subcluster of each cell by k-means on coords (cells x dims); k is
    # lowered until every subcluster has min_cells cells. All 1 when it
    # reaches 1.
    while (k >= 2) {
        km <- tryCatch(stats::kmeans(coords, centers = k, nstart = 10, iter.max = 50),
            error = function(e) NULL)
        if (!is.null(km) && min(tabulate(km$cluster, k)) >= min_cells) {
            # number the subclusters by size, largest first
            sizes <- tabulate(km$cluster, k)
            return(match(km$cluster, order(sizes, decreasing = TRUE)))
        }
        k <- k - 1
    }
    rep(1L, nrow(coords))
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
#' With `method = "kmeans"` each cluster is split into `n_samples`
#' subclusters by k-means, on `embedding` when given or else on the first
#' `n_pcs` principal components of the cluster's `n_features` most variable
#' genes (of its log-normalized cells). The samples then follow the
#' structure within the cluster, and [mixture_heterogeneity()] can tell which
#' genes the subclusters disagree on: if the cluster can be split, the genes
#' that split it are the ones whose calls diverge.
#'
#' With `method = "metadata"` each cluster is split by the values of
#' `group`, a cell metadata column (e.g. library, batch, genotype): each
#' value is one sample, and also its replicate, the same in every cluster
#' (in factor level order, else sorted).
#' To pool libraries into two replicates, as the mixture_modeling_smk
#' pipeline did, map them to two groups first.
#'
#' A cluster gets fewer random parts when it has fewer than `min_cells`
#' cells per part (one sample below `2 * min_cells` cells), and k-means
#' lowers k until every subcluster has `min_cells` cells. With metadata, a
#' value with fewer than `min_cells` cells in a cluster is left out of it
#' (those cells are in no sample).
#'
#' @param counts Raw counts of one split's cells, genes x cells (a matrix or
#'   sparse `Matrix`, gene row names), e.g. the `counts` layer of a Seurat
#'   object.
#' @param cluster Cluster ID of each cell (column of `counts`).
#' @param n_samples Samples per cluster for `"random"`, k for `"kmeans"`.
#' @param method `"random"`, `"kmeans"` or `"metadata"` (see Details).
#' @param group For `"metadata"`, the value of each cell that splits the
#'   clusters (e.g. `obj$orig.ident`).
#' @param embedding For `"kmeans"`, optional cells x dimensions matrix (e.g.
#'   PCA coordinates), rows in the order of the columns of `counts`.
#' @param n_features,n_pcs For `"kmeans"` without `embedding`: the number of
#'   most variable genes and of principal components used per cluster.
#' @param min_cells Fewest cells per sample (see Details).
#' @param seed Random seed for the shuffle and k-means; the caller's random
#'   state is restored afterwards. `NULL` uses the current random state.
#'
#' @return A list: `lognorm` (genes x samples), `cluster` and `replicate` of
#'   each sample (the `group` value, or the sample's number within its
#'   cluster, largest subcluster first), `n_cells` per sample, `sample` (the
#'   sample of each cell, `NA` for cells left out) and `method`. Pass `lognorm`, `cluster` and
#'   `replicate` to [fit_mixture()].
#' @seealso [fit_mixture()], [mixture_heterogeneity()].
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(40 * 60, 2), 40, 60,
#'     dimnames = list(paste0("g", 1:40), paste0("c", 1:60)))
#' cluster <- rep(c("1", "2"), each = 30)
#' pb <- make_pseudobulk(counts, cluster, n_samples = 3, method = "kmeans",
#'     n_features = 20, n_pcs = 5, min_cells = 5)
#' pb$cluster
#' pb$n_cells
#' # one sample per library and cluster
#' pb <- make_pseudobulk(counts, cluster, method = "metadata",
#'     group = rep(c("L1", "L2"), 30))
#' pb$replicate
#' @export
make_pseudobulk <- function(counts, cluster, n_samples = 2,
                            method = c("random", "kmeans", "metadata"), group = NULL,
                            embedding = NULL, n_features = 500, n_pcs = 10,
                            min_cells = 10, seed = 42) {
    method <- match.arg(method)
    if (is.data.frame(counts)) counts <- as.matrix(counts)
    if (!(is.matrix(counts) && is.numeric(counts)) && !methods::is(counts, "Matrix")) {
        fail("counts must be a numeric matrix or a sparse Matrix, genes x cells")
    }
    counts <- methods::as(counts, "CsparseMatrix")
    if (is.null(rownames(counts))) fail("counts needs genes as row names")
    if (any(counts@x < 0)) fail("counts must not be negative")
    if (any(counts@x != round(counts@x))) {
        fail("counts must be raw counts (whole numbers), not normalized values")
    }
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
    if (!is.null(embedding)) {
        if (method != "kmeans") fail("embedding is used by method = \"kmeans\" only")
        embedding <- as.matrix(embedding)
        if (nrow(embedding) != n) fail("embedding has ", nrow(embedding), " rows for ", n, " cells")
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
            if (method == "random") return(rep_len(seq_len(k), length(idx))[sample(length(idx))])
            coords <- if (is.null(embedding)) {
                cluster_pcs(lognorm_by_cell(counts[, idx, drop = FALSE]), n_features, n_pcs)
            } else {
                embedding[idx, , drop = FALSE]
            }
            if (is.null(coords)) return(rep(1L, length(idx)))
            kmeans_groups(coords, k, min_cells)
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

het_reference_z <- function(spread, level, n_bins, trim) {
    # z-score of each spread against others at a similar level: they are
    # binned by level into n_bins quantile bins, and each bin's spreads,
    # Winsorized at the trim and 1 - trim quantiles, give its mean and SD
    z <- rep(NA_real_, length(spread))
    ok <- which(!is.na(spread) & !is.na(level))
    n <- length(ok)
    if (n < 2) return(z)
    nb <- max(1, min(n_bins, n %/% HET_MIN_BIN_GENES))
    bin <- ceiling(rank(level[ok], ties.method = "first") * nb / n)
    for (b in seq_len(nb)) {
        i <- ok[bin == b]
        x <- spread[i]
        lim <- stats::quantile(x, c(trim, 1 - trim), names = FALSE)
        w <- pmin(pmax(x, lim[1]), lim[2])
        z[i] <- (x - mean(w)) / max(stats::sd(w), HET_SD_FLOOR, na.rm = TRUE)
    }
    z
}

#' Find genes whose calls diverge within a cluster
#'
#' Scores each bimodal gene in each cluster by how much the mixture model's
#' calls disagree between the cluster's samples, ideally k-means subclusters
#' from [make_pseudobulk()]. In a homogeneous cluster with Gaussian noise,
#' genes expressed at a similar level should have calls that spread
#' similarly; a gene whose calls spread much more than theirs suggests the
#' cluster is a mix of cells that do and cells that do not express it.
#'
#' The spread of a gene in a cluster is the SD of P(On) of the cluster's
#' samples (`log_p_on_sample` of the bimodal fit). The reference pools every
#' bimodal gene x cluster pair of the split: the pairs are binned by the
#' gene's level in the cluster (`log1p(mean(expm1(x)))` of its samples) into
#' `n_bins` bins of equal size (fewer when there are fewer than 10 pairs per
#' bin). Each bin's spreads are Winsorized at the `trim` and `1 - trim`
#' quantiles, so that the heterogeneous pairs themselves weigh little, and
#' give a mean and SD (at least 0.01), against which each pair's spread is a
#' z-score. A gene is heterogeneous in the cluster when `z >= z_cutoff`,
#' `spread >= min_spread`, and its calls truly diverge: at least one sample
#' has P(On) `>= 1 - call_margin` and another `<= call_margin`. The last
#' condition drops clusters whose samples merely disagree on how sure the
#' call is, which a small, noisy subcluster can cause.
#'
#' The reference is pooled across clusters rather than made per cluster: a
#' gene ON in part of a cluster has a level between its ON and OFF levels,
#' so a mixed cluster's divergent genes crowd into the same level bins,
#' where they would make up too much of a per-cluster reference to be
#' Winsorized away.
#'
#' Clusters with a single sample cannot be scored (`NA`), and nor can genes
#' that are not bimodal (the same rule as [call_mixture()]).
#'
#' @param fit Fits from [fit_mixture()] or [merge_mixture_fits()], of every
#'   gene of the split (the reference uses them all).
#' @param n_bins Number of expression-level bins.
#' @param trim Winsorizing fraction at each end of a bin.
#' @param z_cutoff,min_spread A gene is heterogeneous in a cluster when its
#'   z-score is at least `z_cutoff` and its spread at least `min_spread`...
#' @param call_margin ...and some sample has P(On) at least
#'   `1 - call_margin` and another at most `call_margin`.
#'
#' @return A `scmarco_heterogeneity` list of genes x clusters matrices `z`,
#'   `spread`, `level`, `p_min` and `p_max` (the lowest and highest P(On) of
#'   the cluster's samples) and `heterogeneous` (logical; `FALSE` where not
#'   scored), `n_samples` per cluster, and `params`. Pass it to
#'   [add_probability()] (`heterogeneity`), list the flagged genes with
#'   [cluster_heterogeneity()].
#' @seealso [make_pseudobulk()], [cluster_heterogeneity()].
#' @examples
#' \dontrun{
#' pb <- make_pseudobulk(counts, cluster, n_samples = 4,
#'     method = "kmeans")
#' fit <- fit_mixture(pb$lognorm, pb$cluster, pb$replicate)
#' het <- mixture_heterogeneity(fit)
#' cluster_heterogeneity(het)
#' add_probability("data/mydata.clusters.h5ad", list(all = call_mixture(fit)),
#'     heterogeneity = list(all = het))
#' }
#' @export
mixture_heterogeneity <- function(fit, n_bins = 20, trim = 0.2, z_cutoff = 3,
                                  min_spread = 0.25, call_margin = 0.1) {
    if (!inherits(fit, "scmarco_mixture")) fail("fit must come from fit_mixture()")
    if (trim < 0 || trim >= 0.5) fail("trim must be in [0, 0.5)")
    if (call_margin <= 0 || call_margin >= 0.5) fail("call_margin must be in (0, 0.5)")
    genes <- names(fit$fits)
    has_samples <- vapply(fit$fits, function(g) !is.null(g$bimodal$log_p_on_sample), logical(1))
    if (!all(has_samples)) {
        fail("The fits lack log_p_on_sample; refit them with this version of fit_mixture()")
    }
    clusters <- unique(fit$cluster)
    by_cluster <- split(seq_along(fit$cluster), factor(fit$cluster, levels = clusters))
    p_on <- do.call(rbind, lapply(fit$fits, function(g) exp(g$bimodal$log_p_on_sample)))
    lognorm <- fit$lognorm[genes, , drop = FALSE]
    bimodal <- mixture_bimodal(fit)

    empty <- matrix(NA_real_, length(genes), length(clusters), dimnames = list(genes, clusters))
    spread <- empty
    level <- empty
    z <- empty
    p_min <- empty
    p_max <- empty
    for (cl in clusters) {
        idx <- by_cluster[[cl]]
        level[, cl] <- log1p(rowMeans(expm1(lognorm[, idx, drop = FALSE])))
        if (length(idx) < 2) next
        p <- p_on[, idx, drop = FALSE]
        spread[, cl] <- sqrt(rowSums((p - rowMeans(p))^2) / (length(idx) - 1))
        p_min[, cl] <- apply(p, 1, min)
        p_max[, cl] <- apply(p, 1, max)
    }
    z[bimodal, ] <- het_reference_z(as.vector(spread[bimodal, , drop = FALSE]),
        as.vector(level[bimodal, , drop = FALSE]), n_bins, trim)
    params <- list(n_bins = n_bins, trim = trim, z_cutoff = z_cutoff,
        min_spread = min_spread, call_margin = call_margin)
    # divergent calls: some sample confidently ON and another confidently OFF
    decisive <- !is.na(p_min) & p_min <= call_margin & p_max >= 1 - call_margin
    heterogeneous <- !is.na(z) & z >= z_cutoff & spread >= min_spread & decisive
    n_het <- colSums(heterogeneous)
    message(sum(n_het > 0), " of ", length(clusters), " clusters have heterogeneous genes",
        if (any(n_het > 0)) paste0(" (most: ", names(which.max(n_het)), ", ", max(n_het), ")"))
    structure(list(z = z, spread = spread, level = level, p_min = p_min, p_max = p_max,
        heterogeneous = heterogeneous, n_samples = setNames(lengths(by_cluster), clusters), params = params),
        class = "scmarco_heterogeneity")
}

#' List the genes that make clusters look heterogeneous
#'
#' Flags clusters that may be a mix of cell types, with the genes that
#' suggest it: those whose mixture-model calls diverge between the cluster's
#' subclusters (see [mixture_heterogeneity()]). Such genes are ambiguous in
#' that cluster in every search.
#'
#' @param x A result of [mixture_heterogeneity()], or a dataset from
#'   [ds_load()] whose probabilities were added with `heterogeneity`.
#' @param clusters Clusters to report (default: all).
#'
#' @return A data.frame with one row per heterogeneous gene x cluster
#'   (x split for a dataset): `cluster`, (`split`,) `gene`, `z` and
#'   `spread`, most divergent first within each cluster. Zero rows when
#'   nothing is flagged; `table(res$cluster)` counts the genes per cluster.
#' @seealso [mixture_heterogeneity()], [add_probability()].
#' @examples
#' toy <- system.file("extdata", "toy.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("toy", list(label = "Toy", clusters = toy))
#' cluster_heterogeneity(ds)
#' @export
cluster_heterogeneity <- function(x, clusters = NULL) {
    if (inherits(x, "scmarco_heterogeneity")) {
        if (is.null(clusters)) clusters <- colnames(x$z)
        clusters <- intersect(clusters, colnames(x$z))
        hit <- which(x$heterogeneous[, clusters, drop = FALSE], arr.ind = TRUE)
        out <- data.frame(cluster = clusters[hit[, "col"]],
            gene = rownames(x$z)[hit[, "row"]],
            z = x$z[, clusters, drop = FALSE][hit],
            spread = x$spread[, clusters, drop = FALSE][hit],
            stringsAsFactors = FALSE)
        out <- out[order(match(out$cluster, clusters), -out$z, out$gene), , drop = FALSE]
    } else if (inherits(x, "scmarco_dataset")) {
        if (is.null(clusters)) clusters <- x$clusters
        out <- data.frame(cluster = character(0), split = character(0),
            gene = character(0), z = numeric(0), spread = numeric(0))
        if (!x$has_heterogeneity) return(out)
        rows <- which(x$obs_cluster %in% clusters)
        hit <- which(x$heterogeneous[rows, , drop = FALSE], arr.ind = TRUE)
        if (nrow(hit) > 0) {
            z <- h5ad_layer(x$clusters_h, "het_z", obs = rows)
            spread <- h5ad_layer(x$clusters_h, "het_spread", obs = rows)
            out <- data.frame(cluster = x$obs_cluster[rows][hit[, "row"]],
                split = x$obs_split[rows][hit[, "row"]],
                gene = x$all_genes[hit[, "col"]],
                z = z[hit], spread = spread[hit], stringsAsFactors = FALSE)
            out <- out[order(match(out$cluster, x$clusters), match(out$split, x$splits),
                -out$z, out$gene), , drop = FALSE]
        }
    } else {
        fail("x must come from mixture_heterogeneity() or ds_load()")
    }
    rownames(out) <- NULL
    out
}

assemble_heterogeneity <- function(het, obs_split, obs_cluster, splits, genes) {
    # Per-split mixture_heterogeneity() results -> obs x genes matrices z,
    # spread and flag (1 heterogeneous, 0 not; NA where not scored) and the
    # shared params. Splits left out are not scored.
    if (inherits(het, "scmarco_heterogeneity")) {
        fail("heterogeneity must be named by split, e.g. list(Male = het_male)")
    }
    if (is.character(het)) het <- as.list(het)
    if (!is.list(het) || is.null(names(het)) || any(names(het) == "")) {
        fail("heterogeneity must be named by split, e.g. list(Male = het_male)")
    }
    het <- lapply(het, function(x) if (is.character(x) && length(x) == 1) readRDS(x) else x)
    unknown <- setdiff(names(het), splits)
    if (length(unknown) > 0) {
        fail("heterogeneity: split(s) not in the clusters file: ", first_few(unknown))
    }
    if (anyDuplicated(names(het))) fail("heterogeneity: split(s) given twice")
    z <- matrix(NA_real_, length(obs_split), length(genes))
    spread <- z
    flag <- z
    params <- NULL
    for (s in names(het)) {
        h <- het[[s]]
        where <- paste0("heterogeneity, split '", s, "': ")
        if (!inherits(h, "scmarco_heterogeneity")) {
            fail(where, "expected a result of mixture_heterogeneity()")
        }
        cutoffs <- c("z_cutoff", "min_spread", "call_margin")
        if (is.null(params)) params <- h$params
        if (!identical(h$params[cutoffs], params[cutoffs])) {
            fail("heterogeneity: every split must use the same z_cutoff, min_spread and call_margin")
        }
        extra <- setdiff(rownames(h$z), genes)
        if (length(extra) > 0) {
            fail(where, length(extra), " gene(s) are not in the clusters file, e.g. ",
                first_few(extra))
        }
        rows <- which(obs_split == s)
        present <- obs_cluster[rows]
        unknown <- setdiff(colnames(h$z), present)
        if (length(unknown) > 0) {
            fail(where, "cluster(s) with no cells in this split: ", first_few(unknown))
        }
        gi <- match(rownames(h$z), genes)
        ri <- rows[match(colnames(h$z), present)]
        z[ri, gi] <- t(h$z)
        spread[ri, gi] <- t(h$spread)
        flag[ri, gi] <- ifelse(is.na(t(h$z)), NA, t(h$heterogeneous) * 1)
    }
    list(z = z, spread = spread, flag = flag, params = params)
}
