het_fit <- function(n_genes = 200, k = 4, divergent = c("g1", "g2")) {
  #' Fits of three clusters (k samples each, one with a single sample) whose
  #' sample calls agree up to noise, except `divergent` genes in cluster "A",
  #' ON in half its samples and OFF in the others. Gene "flat" is not bimodal.
  set.seed(5)
  cluster <- c(rep(c("A", "B"), each = k), "C")
  genes <- c(paste0("g", seq_len(n_genes)), "flat")
  n <- length(cluster)
  level <- setNames(runif(length(genes), 0, 3), genes)
  lognorm <- t(vapply(genes, function(g) level[[g]] + rnorm(n, 0, 0.05), numeric(n)))
  colnames(lognorm) <- paste0("s", seq_len(n))
  fits <- lapply(genes, function(g) {
    base <- if (level[[g]] > 1.5) 0.97 else 0.03
    p <- pmin(pmax(base + rnorm(n, 0, 0.02), 1e-4), 1 - 1e-4)
    if (g %in% divergent) p[cluster == "A"] <- rep(c(0.99, 0.01), length.out = k)
    bimodal <- g != "flat"
    list(bimodal = list(log_p_on_cluster = log(c(0.5, 0.5, 0.5)), log_p_on_sample = log(p),
                        on_level = 2, off_level = 0, sd = 0.1, pi = 0.5),
         unimodal = list(level = 1, sd = 1),
         elpd = list(diff = if (bimodal) 30 else 1, diff_se = 2))
  })
  names(fits) <- genes
  structure(list(fits = fits, lognorm = lognorm, cluster = cluster,
                 replicate = c(rep(as.character(seq_len(k)), 2), "1")),
            class = "scmarco_mixture")
}

test_that("random pseudobulk deals each cell into one sample", {
  set.seed(1)
  x <- matrix(rpois(30 * 50, 2), 30, 50, dimnames = list(paste0("g", 1:30), paste0("c", 1:50)))
  cluster <- rep(c("1", "2", "3"), c(24, 20, 6))
  set.seed(9)
  before <- .Random.seed
  pb <- make_pseudobulk(Matrix::Matrix(x, sparse = TRUE), cluster, n_samples = 3, min_cells = 7)
  expect_identical(.Random.seed, before)
  expect_equal(pb$cluster, c("1", "1", "1", "2", "2", "3"))
  expect_equal(pb$replicate, c("1", "2", "3", "1", "2", "1"))
  expect_equal(unname(pb$n_cells), c(8, 8, 8, 10, 10, 6))
  expect_equal(sum(pb$n_cells), ncol(x))
  # as AggregateExpression(return.seurat = TRUE): summed counts, LogNormalized
  summed <- rowSums(x[, pb$sample == "2_1"])
  expect_equal(unname(pb$lognorm[, "2_1"]), unname(log1p(summed / sum(summed) * 1e4)))
  expect_identical(make_pseudobulk(x, cluster, n_samples = 3, min_cells = 7)$sample, pb$sample)
  expect_error(make_pseudobulk(log1p(x / 3), cluster), "raw counts")
})

test_that("metadata mode splits clusters by a cell column", {
  set.seed(4)
  x <- matrix(rpois(10 * 40, 3), 10, 40, dimnames = list(paste0("g", 1:10), NULL))
  cluster <- rep(c("A", "B"), each = 20)
  group <- rep(c("L1", "L2"), 20)
  group[1:5] <- "L3"          # 5 cells of A: too few, left out
  pb <- suppressMessages(make_pseudobulk(x, cluster, method = "metadata", group = group,
                                         min_cells = 6))
  expect_equal(pb$cluster, c("A", "A", "B", "B"))
  expect_equal(pb$replicate, c("L1", "L2", "L1", "L2"))
  expect_equal(unname(pb$n_cells), c(7, 8, 10, 10))
  expect_true(all(is.na(pb$sample[1:5])))
  summed <- rowSums(x[, which(pb$sample == "B_L2")])
  expect_equal(unname(pb$lognorm[, "B_L2"]), unname(log1p(summed / sum(summed) * 1e4)))
  expect_message(make_pseudobulk(x, cluster, method = "metadata", group = group,
                                 min_cells = 6), "left out")
  expect_error(make_pseudobulk(x, cluster, method = "metadata"), "needs group")
  expect_error(make_pseudobulk(x, cluster, group = group), "metadata")
})

test_that("k-means pseudobulk follows the structure within a cluster", {
  set.seed(2)
  n <- 80
  sub <- rep(1:2, each = n / 2)
  x <- matrix(rpois(50 * n, 3), 50, n, dimnames = list(paste0("g", 1:50), paste0("c", 1:n)))
  x[1:10, sub == 2] <- x[1:10, sub == 2] + 30
  pb <- make_pseudobulk(x, rep("A", n), n_samples = 2, method = "kmeans",
                        n_features = 30, n_pcs = 5)
  expect_equal(pb$cluster, c("A", "A"))
  expect_equal(max(table(pb$sample, sub)[1, ]), n / 2)
  # with an embedding, and lowering k when subclusters are too small
  emb <- cbind(c(rnorm(n - 3), 50, 50, 50), rnorm(n))
  pb <- make_pseudobulk(x, rep("A", n), n_samples = 3, method = "kmeans",
                        embedding = emb, min_cells = 5)
  expect_true(all(pb$n_cells >= 5))
  expect_error(make_pseudobulk(x, rep("A", n), embedding = emb), "kmeans")
})

test_that("partition_folds leaves a label with one sample alone", {
  set.seed(3)
  labels <- c(rep(1:2, each = 10), 3)
  folds <- partition_folds(21, 5, labels)
  expect_setequal(unique(folds), 1:5)
})

test_that("mixture_heterogeneity flags genes whose calls diverge in a cluster", {
  fit <- het_fit()
  het <- suppressMessages(mixture_heterogeneity(fit))
  expect_s3_class(het, "scmarco_heterogeneity")
  expect_equal(colnames(het$z), c("A", "B", "C"))
  expect_true(all(het$heterogeneous[c("g1", "g2"), "A"]))
  expect_equal(sum(het$heterogeneous), 2)
  # single-sample clusters and non-bimodal genes are not scored
  expect_true(all(is.na(het$z[, "C"])))
  expect_true(all(is.na(het$z["flat", ])))
  expect_equal(unname(het$n_samples), c(4, 4, 1))
  tab <- cluster_heterogeneity(het)
  expect_equal(tab$cluster, c("A", "A"))
  expect_setequal(tab$gene, c("g1", "g2"))
  expect_equal(nrow(cluster_heterogeneity(het, clusters = "B")), 0)
  # a stricter spread cutoff drops them
  expect_equal(sum(suppressMessages(mixture_heterogeneity(fit, min_spread = 0.9))$heterogeneous), 0)
  expect_equal(unname(het$p_min["g1", "A"]), 0.01)
  expect_equal(unname(het$p_max["g1", "A"]), 0.99)

  # samples that disagree only on how sure the call is are not heterogeneous
  unsure <- fit
  p <- exp(unsure$fits$g1$bimodal$log_p_on_sample)
  p[unsure$cluster == "A"] <- c(0.6, 0.01, 0.6, 0.01)
  unsure$fits$g1$bimodal$log_p_on_sample <- log(p)
  het <- suppressMessages(mixture_heterogeneity(unsure))
  expect_gte(het$spread["g1", "A"], 0.25)
  expect_false(het$heterogeneous["g1", "A"])
  expect_true(het$heterogeneous["g2", "A"])
  expect_true(suppressMessages(mixture_heterogeneity(unsure, call_margin = 0.45))$heterogeneous["g1", "A"])

  old <- fit
  old$fits$g3$bimodal$log_p_on_sample <- NULL
  expect_error(mixture_heterogeneity(old), "log_p_on_sample")
})

test_that("binarize makes heterogeneous genes ambiguous", {
  prob <- matrix(c(0.99, 0.01, 0.99, 0.99, 0.01, 0.99), 2,
                 dimnames = list(NULL, c("gA", "gB", "gC")))
  het <- matrix(c(TRUE, TRUE, FALSE, NA, TRUE, FALSE), 2)
  bimodal <- matrix(c(TRUE, TRUE, FALSE), 3, 1, dimnames = list(c("gA", "gB", "gC"), "s"))
  b <- binarize(prob, c("s", "s"), c("c1", "c2"), splits = "s", clusters = c("c1", "c2"),
                margin = 0.1, bimodal = bimodal, heterogeneous = het)
  expect_equal(unname(b$ambiguous[, , "s"]), rbind(c(TRUE, TRUE), c(FALSE, FALSE), c(FALSE, FALSE)))
  expect_equal(unname(b$on[, , "s"]), rbind(c(FALSE, FALSE), c(TRUE, TRUE), c(FALSE, FALSE)))
  # ON in c2 but ambiguous in c1: no longer a clean marker of c2
  expect_equal(unname(any_on(b, "gA")[1, ]), c(TRUE, TRUE))
  expect_equal(unname(b$heterogeneous[, , "s"]), rbind(c(TRUE, TRUE), c(FALSE, FALSE), c(FALSE, FALSE)))
  expect_equal(dim(restrict_clusters(b, "c2")$heterogeneous), c(3, 1, 1))
  plain <- binarize(prob, c("s", "s"), c("c1", "c2"), splits = "s", clusters = c("c1", "c2"),
                    margin = 0.1)
  expect_false(any(plain$heterogeneous))
})

test_that("heterogeneity reaches the dataset and the searches", {
  toy <- system.file("extdata", "toy.clusters.h5ad", package = "scMarco")
  path <- tempfile(fileext = ".clusters.h5ad")
  file.copy(toy, path)
  ds <- ds_load(paste0("toy", runif(1)), list(label = "Toy", clusters = path))
  present <- function(s) ds$clusters[ds$clusters %in% ds$obs_cluster[ds$obs_split == s]]
  tables <- lapply(c(Male = "Male", Female = "Female"), function(s) {
    tab <- as.data.frame(matrix(0.01, length(ds$all_genes), length(present(s)),
      dimnames = list(ds$all_genes, present(s))), check.names = FALSE)
    tab[1, ] <- 0.99
    tab$exprType <- "bimodal"
    tab
  })
  target <- present("Male")[1]
  gene <- ds$all_genes[1]
  shape <- matrix(NA_real_, length(ds$all_genes), length(present("Male")),
                  dimnames = list(ds$all_genes, present("Male")))
  z <- shape
  spread <- shape
  z[, ] <- 0
  spread[, ] <- 0.01
  z[gene, target] <- 8
  spread[gene, target] <- 0.6
  het <- structure(list(z = z, spread = spread,
    heterogeneous = !is.na(z) & z >= 3 & spread >= 0.25,
    params = list(n_bins = 20, trim = 0.2, z_cutoff = 3, min_spread = 0.25, call_margin = 0.1)),
    class = "scmarco_heterogeneity")

  suppressMessages(suppressWarnings(add_probability(path, tables, heterogeneity = list(Male = het))))
  ds <- ds_load(paste0("toy", runif(1)), list(label = "Toy", clusters = path))
  expect_true(ds$has_heterogeneity)
  b <- ds_binarize(ds, margin = 0.1)
  expect_false(b$on[gene, target, "Male"])
  expect_true(b$ambiguous[gene, target, "Male"])
  expect_true(b$on[gene, target, "Female"] || !b$present[target, "Female"])
  expect_equal(sum(b$heterogeneous), 1)
  tab <- cluster_heterogeneity(ds)
  expect_equal(tab[, c("cluster", "split", "gene")],
               data.frame(cluster = target, split = "Male", gene = gene))
  expect_equal(tab$z, 8)
  note <- heterogeneity_note(b, ds, target)
  expect_match(as.character(note), gene)
  expect_null(heterogeneity_note(b, ds, setdiff(ds$clusters, target)[1]))

  # binarizing again without heterogeneity drops the earlier scores
  suppressMessages(suppressWarnings(add_probability(path, tables)))
  ds <- ds_load(paste0("toy", runif(1)), list(label = "Toy", clusters = path))
  expect_false(ds$has_heterogeneity)
  expect_false(any(ds_binarize(ds, margin = 0.1)$heterogeneous))
  expect_equal(nrow(cluster_heterogeneity(ds)), 0)

  bad <- het
  bad$params$z_cutoff <- 2
  quietly <- function(x) suppressMessages(suppressWarnings(x))
  expect_error(quietly(add_probability(path, tables,
    heterogeneity = list(Male = het, Female = bad))), "same z_cutoff")
  expect_error(quietly(add_probability(path, tables, heterogeneity = list(Nope = het))), "Nope")
  expect_error(quietly(add_probability(path, tables, heterogeneity = het)), "named by split")
})
