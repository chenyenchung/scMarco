fake_gene_fit <- function(log_p_on, diff, diff_se, on_level = 1.23456789) {
  #' One gene's fit as fit_mixture() returns it, without running Stan
  list(
    bimodal = list(log_p_on_cluster = log_p_on, on_level = on_level,
                   off_level = 0.0123456, sd = 0.00312, pi = 0.987654),
    unimodal = list(level = 0.456, sd = 0.789),
    elpd = list(diff = diff, diff_se = diff_se)
  )
}

fake_mixture <- function() {
  #' Three clusters x two replicates; geneA and geneB bimodal, geneC not
  set.seed(1)
  cluster <- rep(c("1", "2", "10"), 2)
  lognorm <- rbind(
    geneA = c(2, 0.1, 0.1, 2.1, 0.2, 0.1),
    geneB = c(0.1, 1.8, 0.1, 0.2, 1.9, 0.1),
    geneC = c(0.3, 0.4, 0.2, 0.3, 0.4, 0.3)
  )
  colnames(lognorm) <- paste0("s", 1:6)
  fits <- list(
    geneA = fake_gene_fit(log(c(0.999, 0.004, 0.003)), diff = 30, diff_se = 5),
    geneB = fake_gene_fit(log(c(0.002, 0.995, 0.001)), diff = 25, diff_se = 4),
    geneC = fake_gene_fit(log(c(0.6, 0.7, 0.5)), diff = 3, diff_se = 2)
  )
  structure(list(fits = fits, lognorm = lognorm, cluster = cluster,
                 replicate = rep(c("1", "2"), each = 3)),
            class = "scmarco_mixture")
}

test_that("partition_folds fills every fold and splits every replicate", {
  set.seed(2)
  labels <- rep(1:2, each = 23)
  folds <- partition_folds(46, 10, labels)
  expect_setequal(unique(folds), 1:10)
  expect_equal(sum(folds == 10), 4 + 46 %% 10)
  for (l in 1:2) expect_gt(length(unique(folds[labels == l])), 1)
  set.seed(2)
  expect_identical(partition_folds(46, 10, labels), folds)
})

test_that("call_mixture calls bimodal genes and keeps exact fit parameters", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("matrixStats")
  tab <- suppressMessages(call_mixture(fake_mixture()))
  expect_equal(rownames(tab), c("geneA", "geneB", "geneC"))
  expect_equal(colnames(tab)[1:3], c("1", "2", "10"))
  expect_equal(tab$exprType, c("bimodal", "bimodal", "unimodal"))
  expect_equal(unname(unlist(tab["geneA", 1:3])), c(1, 0, 0))
  # a gene without bimodality gets one probability in every cluster
  expect_length(unique(unlist(tab["geneC", 1:3])), 1)
  expect_identical(tab$bimodal.on_level, rep(1.23456789, 3))
  expect_identical(tab$bimodal.sd, rep(0.00312, 3))
  expect_identical(tab$unimodal.level, rep(0.46, 3))
  expect_true(all(tab$exprType %in% "unimodal" | tab$diff.elpd > 2 * tab$diff.elpd.se))
})

test_that("call_mixture needs a bimodal gene", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("matrixStats")
  fit <- fake_mixture()
  fit$fits <- fit$fits["geneC"]
  expect_error(suppressMessages(call_mixture(fit)), "No bimodal gene")
})

test_that("merge_mixture_fits joins chunks of genes", {
  fit <- fake_mixture()
  a <- fit
  a$fits <- fit$fits[1:2]
  a$lognorm <- fit$lognorm[1:2, ]
  b <- fit
  b$fits <- fit$fits[3]
  b$lognorm <- fit$lognorm[3, , drop = FALSE]
  expect_identical(merge_mixture_fits(a, b), fit)
  expect_identical(merge_mixture_fits(list(a, b)), fit)
  expect_error(merge_mixture_fits(a, a), "fit twice")
  b$cluster <- rev(b$cluster)
  expect_error(merge_mixture_fits(a, b), "share their samples")
})

test_that("exact fit parameters reach the dataset unchanged", {
  toy <- system.file("extdata", "toy.clusters.h5ad", package = "scMarco")
  path <- tempfile(fileext = ".clusters.h5ad")
  file.copy(toy, path)
  ds <- ds_load(paste0("toy", runif(1)), list(label = "Toy", clusters = path))
  on_level <- c(Male = pi / 3, Female = exp(1) / 7)
  tables <- lapply(c(Male = "Male", Female = "Female"), function(s) {
    present <- ds$clusters[ds$clusters %in% ds$obs_cluster[ds$obs_split == s]]
    tab <- as.data.frame(matrix(0.01, length(ds$all_genes), length(present),
      dimnames = list(ds$all_genes, present)), check.names = FALSE)
    tab$bimodal.on_level <- on_level[[s]]
    tab$bimodal.off_level <- 1e-4
    tab$bimodal.sd <- 0.00123
    tab$bimodal.pi <- 0.9987
    tab$exprType <- "bimodal"
    tab
  })
  suppressMessages(suppressWarnings(add_probability(path, tables)))
  ds <- ds_load(paste0("toy", runif(1)), list(label = "Toy", clusters = path))
  expect_identical(unname(ds$mixture$on_level[1, c("Male", "Female")]),
                   unname(on_level))
  expect_identical(unname(ds$mixture$sd[1, "Male"]), 0.00123)
})

test_that("the P(On) curve floors only rounded-away parameters", {
  fit <- list(on_level = 1, off_level = 0, sd = 0.003, pi = 0.999)
  expect_identical(mixture_floored(fit), fit)
  expect_equal(mixture_p_on(0.5, fit), 0.999)
  rounded <- list(on_level = 1, off_level = 0, sd = 0, pi = 1)
  expect_equal(mixture_floored(rounded)[c("sd", "pi")],
               list(sd = MIXTURE_FLOOR, pi = 1 - MIXTURE_FLOOR))
  expect_true(is.finite(mixture_p_on(0.5, rounded)))
})

test_that("fit_mixture fits genes with Stan", {
  skip_on_cran()
  skip_on_bioc()
  skip_if(Sys.getenv("SCMARCO_TEST_STAN") == "",
          "set SCMARCO_TEST_STAN=true to compile and run the Stan models")
  skip_if_not_installed("rstan")
  set.seed(3)
  n <- 12
  cluster <- rep(as.character(1:6), 2)
  lognorm <- rbind(
    geneA = ifelse(cluster %in% c("1", "2"), 2, 0.05) + rnorm(n, 0, 0.05),
    geneB = rnorm(n, 0.3, 0.1)
  )
  colnames(lognorm) <- paste0("s", seq_len(n))
  fit <- suppressWarnings(suppressMessages(fit_mixture(
    lognorm, cluster, replicate = rep(c("1", "2"), each = 6), n_iter = 200,
    n_folds = 3
  )))
  expect_s3_class(fit, "scmarco_mixture")
  expect_named(fit$fits, c("geneA", "geneB"))
  expect_length(fit$fits$geneA$bimodal$log_p_on_cluster, 6)
  expect_length(fit$fits$geneA$elpd$bimodal_i, n)
  expect_gt(fit$fits$geneA$elpd$diff, 0)

  # seeded per gene: chunking does not matter, the caller's RNG is untouched
  set.seed(99)
  before <- .Random.seed
  alone <- suppressWarnings(suppressMessages(fit_mixture(
    lognorm, cluster, replicate = rep(c("1", "2"), each = 6), genes = "geneB",
    n_iter = 200, n_folds = 3
  )))
  expect_identical(.Random.seed, before)
  expect_identical(alone$fits$geneB, fit$fits$geneB)
})
