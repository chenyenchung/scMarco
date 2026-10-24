# Combination-finding core
#
# Pure functions (no Shiny, no I/O) shared by every tab. They operate on a
# "binarized dataset" produced by binarize():
#
#   on[gene, cluster, split]  logical  confidently ON: prob >= 1 - margin
#                                      (FALSE when the cluster is absent from
#                                      the split or the gene is non-bimodal in
#                                      that split)
#   ambiguous[gene, cluster, split]
#                             logical  margin < prob < 1 - margin: neither
#                                      confidently ON nor confidently OFF
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

binarize <- function(prob, obs_split, obs_cluster, splits, clusters,
                     margin, bimodal = NULL) {
  #' prob: n_obs x n_genes numeric matrix (one row per split x cluster that
  #'   exists), colnames = genes
  #' obs_split, obs_cluster: character vectors of length n_obs
  #' splits, clusters: the full, ordered level sets
  #' margin: confidence margin in (0, 0.5]; ON when prob >= 1 - margin, OFF
  #'   when prob <= margin, ambiguous in between. At 0.5 there is no ambiguous
  #'   band and ON is prob > 0.5.
  #' bimodal: optional n_genes x n_splits logical matrix
  stopifnot(
    nrow(prob) == length(obs_split),
    length(obs_split) == length(obs_cluster),
    margin > 0, margin <= 0.5
  )
  genes <- colnames(prob)
  on <- array(
    FALSE,
    dim = c(length(genes), length(clusters), length(splits)),
    dimnames = list(genes, clusters, splits)
  )
  ambiguous <- on
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
  off <- is.na(prob) | prob <= margin + eps
  hit <- !off & prob >= 1 - margin - eps
  unsure <- !off & !hit
  for (i in seq_len(nrow(prob))) {
    on[, ci[i], si[i]] <- hit[i, ]
    ambiguous[, ci[i], si[i]] <- unsure[i, ]
  }

  if (is.null(bimodal)) {
    bimodal <- matrix(
      TRUE, length(genes), length(splits),
      dimnames = list(genes, splits)
    )
  } else {
    bimodal <- bimodal[genes, splits, drop = FALSE]
    for (s in seq_along(splits)) {
      on[!bimodal[, s], , s] <- FALSE
      ambiguous[!bimodal[, s], , s] <- FALSE
    }
  }

  structure(
    list(
      on = on, ambiguous = ambiguous, present = present, bimodal = bimodal,
      genes = genes, clusters = clusters, splits = splits, margin = margin
    ),
    class = "scmarco_binarized"
  )
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
any_on <- function(b, genes = b$genes, clusters = b$clusters) {
  maybe <- b$on[genes, clusters, , drop = FALSE] |
    b$ambiguous[genes, clusters, , drop = FALSE]
  out <- rowSums(maybe, dims = 2) > 0
  dim(out) <- c(length(genes), length(clusters))
  dimnames(out) <- list(genes, clusters)
  out
}

target_genes <- function(b, ident, required, min_n, genes = b$genes) {
  #' Genes that pass the split criteria in the single cluster `ident`
  if (length(ident) != 1 || !ident %in% b$clusters) {
    stop("ident must be a single known cluster")
  }
  genes <- intersect(genes, b$genes)
  pass <- pass_matrix(b, required, min_n, clusters = ident, genes = genes)
  genes[pass[, 1]]
}

offtarget_counts <- function(b, genes, exclude) {
  #' Number of other clusters in which each gene is ON (or ambiguous) in any
  #' split
  others <- setdiff(b$clusters, exclude)
  if (length(genes) == 0) {
    return(setNames(integer(0), character(0)))
  }
  counts <- rowSums(any_on(b, genes = genes, clusters = others))
  setNames(as.integer(counts), genes)
}

from_cluster <- function(b, ident, required, min_n, limit = 10,
                         genes = b$genes) {
  #' "From cluster": candidate markers of `ident`, fewest off-targets first
  cand <- target_genes(b, ident, required, min_n, genes = genes)
  counts <- offtarget_counts(b, cand, exclude = ident)
  counts <- counts[order(counts, names(counts))]
  head(counts, limit)
}

clusters_expressing <- function(b, genes, required, min_n,
                                mode = c("all", "any")) {
  #' "From gene": clusters in which all (or any) of `genes` pass
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

distinct_within <- function(b, idents, required, min_n, k = 1,
                            genes = b$genes) {
  #' "Find distinct genes within": genes that pass in exactly k of `idents`
  #' and are not ON or ambiguous (in any split) in any other cluster of
  #' `idents`.
  #' Returns a data.frame(cluster, gene, offtargets), where offtargets counts
  #' clusters OUTSIDE `idents` that also express the gene, sorted best first
  #' within each cluster; zero rows when nothing is found.
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
  #' First n rows per cluster of a distinct_within() table
  keep <- stats::ave(seq_len(nrow(tbl)), tbl$cluster, FUN = seq_along) <= n
  tbl[keep, , drop = FALSE]
}

empty_pairs <- function() {
  data.frame(gene1 = character(0), gene2 = character(0),
             offtargets = integer(0), single_offtargets = integer(0))
}

pair_offtargets <- function(b, ident, required, min_n, first = NULL,
                            limit = 50, genes = b$genes) {
  #' Rank gene pairs that both pass in `ident` by the number of OTHER
  #' clusters in which both genes are ON or ambiguous (any split). If `first`
  #' is given, only pairs containing it are returned.
  #' Returns data.frame(gene1, gene2, offtargets, single_offtargets).
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

  if (!is.null(first)) {
    idx <- cbind(1L, seq(2, length(cand)))
  } else {
    idx <- which(upper.tri(shared), arr.ind = TRUE)
  }
  vals <- as.integer(shared[idx])
  # Keep only pairs that can make the top `limit` before sorting: there can be
  # millions of pairs. Ties break on how specific each gene is on its own
  # (sum of single-gene off-targets), then on gene position (alphabetical).
  if (length(vals) > limit) {
    worst_kept <- sort(vals, partial = limit)[limit]
    keep <- vals <= worst_kept
    idx <- idx[keep, , drop = FALSE]
    vals <- vals[keep]
  }
  single <- diag(shared)
  single_sum <- as.integer(single[idx[, 1]] + single[idx[, 2]])
  ord <- head(order(vals, single_sum, idx[, 1], idx[, 2], method = "radix"), limit)
  data.frame(
    gene1 = cand[idx[ord, 1]],
    gene2 = cand[idx[ord, 2]],
    offtargets = vals[ord],
    single_offtargets = single_sum[ord]
  )
}

combo_offtarget_clusters <- function(b, genes, exclude = character(0)) {
  #' Clusters (other than `exclude`) in which ALL `genes` are ON (or
  #' ambiguous) in some split
  genes <- intersect(genes, b$genes)
  if (length(genes) == 0) {
    return(character(0))
  }
  a <- any_on(b, genes = genes, clusters = setdiff(b$clusters, exclude))
  colnames(a)[colSums(a) == length(genes)]
}

state_table <- function(b, genes, clusters = b$clusters) {
  #' Long table of per-gene states used by the tile plots:
  #' "Expressed", "Not expressed", "Ambiguous", "Non-bimodal", "Absent"
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
