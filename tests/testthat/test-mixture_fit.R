fake_mixture <- function() {
  #' Three clusters x two replicates; geneA and geneB binarized, geneC's
  #' bimodal fit is no better than its unimodal one
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
    geneC = fake_gene_fit(log(c(0.6, 0.7, 0.5)), diff = 1, diff_se = 2)
  )
  new_mixture(fits, lognorm, cluster, replicate = rep(c("1", "2"), each = 3),
              settings = list(n_iter = 500, n_folds = 10, seed = 42, cores = 1))
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

test_that("call_mixture calls binarized genes and keeps exact fit parameters", {
  fit <- fake_mixture()
  # posterior means that differ from the medians, to tell the summaries apart
  fit$fits$geneA$bimodal$p_on_mean <- c(0.91, 0.05, 0.04)
  tab <- suppressMessages(call_mixture(fit))
  expect_equal(rownames(tab), c("geneA", "geneB", "geneC"))
  expect_equal(colnames(tab)[1:3], c("1", "2", "10"))
  expect_equal(tab$exprType, c("bimodal", "bimodal", "unimodal"))
  expect_equal(tab$not.binarized, c(NA, NA, "elpd"))
  expect_equal(unname(unlist(tab["geneA", 1:3])), c(0.91, 0.05, 0.04))
  # genes not binarized get no calls
  expect_true(all(is.na(unlist(tab["geneC", 1:3]))))
  # z, Benjamini-Hochberg q, Ashman's D
  expect_equal(tab$elpd.z, c(6, 6.25, 0.5))
  p <- pnorm(c(6, 6.25, 0.5), lower.tail = FALSE)
  expect_equal(tab$elpd.q, signif(p.adjust(p, "BH"), 3))
  expect_equal(tab$ashman.D, rep(round((1.23456789 - 0.0123456) / 0.00312, 3), 3))
  expect_equal(tab$ashman.D.lower, tab$ashman.D)
  expect_identical(tab$bimodal.on_level, rep(1.23456789, 3))
  expect_identical(tab$bimodal.sd, rep(0.00312, 3))
  expect_identical(tab$unimodal.level, rep(0.46, 3))
  expect_identical(attr(tab, "flags"),
                   matrix(0L, 3, 3, dimnames = list(rownames(tab), c("1", "2", "10"))))
  expect_true(all(is.na(attr(tab, "fold_ambient"))))
  b <- attr(tab, "binarization")
  expect_identical(b$summary, "mean")
  expect_identical(b$bh_q, 0.10)
  expect_false(b$ambient)
  # the median summary is binarizer 1.0's
  med <- suppressMessages(call_mixture(fit, summary = "median"))
  expect_equal(unname(unlist(med["geneA", 1:3])), c(1, 0, 0))
  # the cuts in use are printed
  msgs <- testthat::capture_messages(call_mixture(fit))
  expect_match(msgs, "2 of 3 genes are binarized", all = FALSE)
  expect_match(msgs, "false discovery rate of 0.1 .*set bh_q", all = FALSE)
  expect_match(msgs, "ambient check are skipped", all = FALSE)
})

test_that("call_mixture leaves out genes with a small Ashman's D or chains that disagree", {
  fit <- fake_mixture()
  # geneA: ON and OFF 1.5 SDs apart
  fit$fits$geneA <- fake_gene_fit(log(c(0.999, 0.004, 0.003)), diff = 30, diff_se = 5,
                                  on_level = 1, off_level = 0.25, sd = 0.5)
  # geneB: chains disagree on 1 of 3 clusters
  fit$fits$geneB$bimodal$chain_spread <- c(0, 0.7, 0)
  expect_warning(tab <- suppressMessages(call_mixture(fit)), "No gene is binarized")
  expect_equal(tab$ashman.D[1], 1.5)
  expect_equal(tab$not.binarized, c("ashman_D", "chains", "elpd"))
  expect_equal(tab$spread.share[2], round(1 / 3, 3))
  # looser cuts keep both; the cluster whose chains disagree is flagged
  tab <- suppressMessages(call_mixture(fit, min_D = 1, reject_spread_frac = 0.5))
  expect_equal(tab$exprType, c("bimodal", "bimodal", "unimodal"))
  expect_identical(unname(attr(tab, "flags")["geneB", ]), c(0L, CALL_FLAG_CHAINS, 0L))
  # Ashman's D lower bound: 10th percentile over draws
  d <- fit$fits$geneA$bimodal$draws
  d[, , 1] <- seq(0.6, 1.5, length.out = 10)
  fit$fits$geneA$bimodal$draws <- d
  D <- (seq(0.6, 1.5, length.out = 10) - 0.25) / 0.5
  tab <- suppressMessages(call_mixture(fit, min_D = 0))
  expect_equal(tab$ashman.D[1], round(median(D), 3))
  expect_equal(tab$ashman.D.lower[1], round(unname(quantile(D, 0.1)), 3))
})

test_that("call_mixture flags clusters below the detection floor or above ambient RNA", {
  fit <- fake_mixture()
  # 4 cells per cluster in one library; 10 more genes make up the library's
  # UMIs (the ambient profile)
  cl <- rep(c("1", "2", "10"), each = 4)
  counts <- matrix(5L, 13, 12, dimnames = list(c("geneA", "geneB", "geneC",
    paste0("bg", 1:10)), paste0("c", 1:12)))
  counts["geneA", ] <- c(40L, 40L, 40L, 40L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)
  counts["geneB", ] <- c(0L, 0L, 0L, 0L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L)
  amb <- ambient_summary(counts, cl, rho = 0.01)
  tab <- suppressMessages(call_mixture(fit, amb))
  flags <- attr(tab, "flags")
  # geneA is undetected in clusters 2 and 10, and above ambient in cluster 1,
  # where it is ON (binarize() applies the ambient flag to OFF calls only)
  expect_identical(unname(flags["geneA", ]),
                   c(CALL_FLAG_AMBIENT, CALL_FLAG_DETECTION, CALL_FLAG_DETECTION))
  # geneB's 120 UMIs in cluster 10 are far above 8 x its ambient share there
  expect_identical(unname(flags["geneB", ]),
                   c(CALL_FLAG_DETECTION, CALL_FLAG_AMBIENT, CALL_FLAG_AMBIENT))
  lam <- ambient_expected(amb, rownames(tab), c("1", "2", "10"))
  expect_equal(attr(tab, "fold_ambient"),
               round(log2((amb$umi[rownames(tab), c("1", "2", "10")] + 1) / (lam + 1)), 3))
  expect_true(attr(tab, "binarization")$ambient)
  expect_identical(attr(tab, "binarization")$rho, list(`1` = 0.01))
  # a permissive enough offset lets the UMIs pass as ambient
  tab <- suppressMessages(call_mixture(fit, amb, ambient_log2fc = 12))
  expect_identical(unname(attr(tab, "flags")["geneA", ]), c(0L, CALL_FLAG_DETECTION, CALL_FLAG_DETECTION))
  expect_identical(unname(attr(tab, "flags")["geneB", ]), c(CALL_FLAG_DETECTION, 0L, 0L))
  # ambient must cover the fit
  expect_error(suppressMessages(call_mixture(fit, ambient_summary(counts[-1, ], cl, rho = 0.01))),
               "Gene\\(s\\) not in ambient: geneA")
  expect_error(suppressMessages(call_mixture(fit, ambient_summary(counts[, 1:8], cl[1:8], rho = 0.01))),
               "Cluster\\(s\\) not in ambient: 10")
  expect_error(call_mixture(fit, list()), "must come from ambient_summary")
})

test_that("call_mixture checks its parameters", {
  fit <- fake_mixture()
  expect_error(call_mixture(fit, bh_q = 0), "bh_q must be one number in \\(0, 1\\]")
  expect_error(call_mixture(fit, bh_q = 1.5), "bh_q must be one number in \\(0, 1\\]")
  expect_error(call_mixture(fit, spread_cutoff = 2), "spread_cutoff must be")
  expect_error(call_mixture(fit, min_D = -1), "min_D must be one number of at least 0")
  expect_error(call_mixture(fit, summary = "mode"), "should be one of")
})

test_that("call_mixture calls binarizer 1.0 fits from their medians, with a warning", {
  fit <- fake_mixture()
  fit$fits <- lapply(fit$fits, function(g) {
    g$bimodal[c("draws", "p_on_mean", "chain_spread")] <- NULL
    g
  })
  fit$binarizer_version <- "1.0"
  expect_warning(msgs <- testthat::capture_messages(tab <- call_mixture(fit)),
                 "refitting with fit_mixture\\(\\) is recommended")
  expect_match(msgs, "posterior median and there is no chain-spread check", all = FALSE)
  expect_equal(unname(unlist(tab["geneA", 1:3])), c(1, 0, 0))
  expect_equal(tab$exprType, c("bimodal", "bimodal", "unimodal"))
  expect_equal(tab$ashman.D[1], round((1.23456789 - 0.0123456) / 0.00312, 3))
  expect_true(all(is.na(tab$ashman.D.lower)))
  expect_identical(attr(tab, "binarization")$summary, "median")
  # a 1.1 fit must keep its draws
  fit$binarizer_version <- SCMARCO_BINARIZER_VERSION
  expect_error(suppressMessages(call_mixture(fit)), "keeps no draws")
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
  suppressMessages(suppressWarnings(add_probability(path, tables, external = TRUE)))
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
  withr::local_options(scMarco.stan_cache = withr::local_tempdir())
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
  # binarizer 1.1: 8 chains on the bimodal fit, draws and summaries kept
  expect_identical(dim(fit$fits$geneA$bimodal$draws), c(100L, 8L, 4L))
  expect_length(fit$fits$geneA$bimodal$p_on_mean, 6)
  expect_length(fit$fits$geneA$bimodal$chain_spread, 6)
  expect_identical(fit$settings$prior, "range")
  expect_identical(fit$settings$n_chains_bimodal, 8L)

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

test_that("fits that did not converge are reported", {
  diag <- function(rhat, n_eff) cbind(Rhat = rhat, n_eff = n_eff)
  fit <- fake_mixture()
  # no diagnostics (fits made without Stan): convergence unknown, not reported
  expect_true(all(is.na(fit$convergence$converged)))
  for (g in names(fit$fits)) {
    fit$fits[[g]]$bimodal$diagnostics <- diag(c(1.001, 1.002), c(800, 900))
    fit$fits[[g]]$unimodal$diagnostics <- diag(c(1.001, 1.003), c(700, 750))
  }
  fit$fits$geneA$bimodal$diagnostics <- diag(c(1.62, 1.01), c(9.5, 500))
  fit$fits$geneC$unimodal$diagnostics <- diag(c(NaN, 1), c(800, 800))
  fit <- new_mixture(fit$fits, fit$lognorm, fit$cluster, fit$replicate, fit$settings)
  conv <- fit$convergence
  expect_equal(conv$gene, c("geneA", "geneB", "geneC"))
  expect_equal(conv$bimodal_rhat[1], 1.62)
  expect_equal(conv$bimodal_n_eff[1], 9.5)
  expect_equal(conv$converged, c(FALSE, TRUE, FALSE))
  expect_output(print(fit), "2 of 3 gene\\(s\\) did not converge.*geneA, geneC")
  expect_warning(suppressMessages(call_mixture(fit)),
                 "2 of 3 gene\\(s\\) did not converge.*Their calls are kept")
  # fits saved before convergence was recorded get it from their diagnostics
  fit$convergence <- NULL
  expect_equal(not_converged(fit), c("geneA", "geneC"))
})

test_that("compiled Stan models are cached on disk", {
  dir <- withr::local_tempdir()
  stan <- withr::local_tempfile(fileext = ".stan", lines = "model { }")
  compiled <- 0
  compile <- function() {
    compiled <<- compiled + 1
    list(model = "m")
  }
  usable <- function(m) identical(m, list(model = "m"))
  m <- suppressMessages(cached_stan_model("m", stan, compile, usable, dir = dir))
  expect_equal(compiled, 1)
  expect_true(file.exists(stan_cache_file("m", stan, dir)))
  # a later session loads it
  expect_message(m2 <- cached_stan_model("m", stan, compile, usable, dir = dir),
                 "Using the compiled m Stan model")
  expect_equal(compiled, 1)
  expect_identical(m2, m)
  # another Stan code gets another cache file
  writeLines("model { } ", stan)
  suppressMessages(cached_stan_model("m", stan, compile, usable, dir = dir))
  expect_equal(compiled, 2)
  # an unusable or unreadable cache is compiled again and replaced
  writeLines("not an rds", stan_cache_file("m", stan, dir))
  suppressMessages(cached_stan_model("m", stan, compile, usable, dir = dir))
  expect_equal(compiled, 3)
  expect_identical(readRDS(stan_cache_file("m", stan, dir)), m)
  # "" turns the cache off
  suppressMessages(cached_stan_model("m", stan, compile, usable, dir = ""))
  expect_equal(compiled, 4)
  expect_length(list.files(dir, "\\.tmp$"), 0)
})

test_that("the binarizer 1.1 prior is centred on the midpoint with range / 3", {
  p <- mixture_prior(c(0, 0.5, 3))
  expect_equal(p$prior_on_mean, 1.5)
  expect_equal(p$prior_off_mean, 1.5)
  expect_equal(p$prior_on_sd, 1)
  expect_equal(p$prior_off_sd, 1)
  # constant values: the SD floor
  expect_equal(mixture_prior(rep(0.2, 4))$prior_on_sd, MIXTURE_PRIOR_SD_FLOOR)
})

test_that("cluster summaries come from the draws: posterior mean and chain spread", {
  # 2 chains x 3 draws of (mu_on, mu_off, sd1, pi_on)
  d <- array(c(2, 2, 2, 2, 2, 2,  0, 0, 0, 0, 0, 0,  0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
               0.5, 0.5, 0.5, 0.5, 0.5, 0.5), dim = c(3, 2, 4))
  x <- c(2, 2.1, 0, 0.1)            # cluster 1 ON, cluster 2 OFF
  s <- cluster_summaries(d, x, c(1, 1, 2, 2), 2)
  expect_gt(s$p_on_mean[1], 0.99)
  expect_lt(s$p_on_mean[2], 0.01)
  expect_equal(unname(s$chain_spread), c(0, 0), tolerance = 1e-12)
  # identical draws in both chains -> no spread; recomputation matches the
  # closed form for one draw
  lp <- log_p_on_draws(matrix(d, ncol = 4)[1, , drop = FALSE], x, c(1, 1, 2, 2), 2)
  on <- sum(dnorm(x[1:2], 2, 0.3, log = TRUE)); off <- sum(dnorm(x[1:2], 0, 0.3, log = TRUE))
  expect_equal(unname(lp[1, 1]), log(0.5) + on - log(0.5 * exp(on) + 0.5 * exp(off)))
})

test_that("fit_mixture checks n_chains_bimodal", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("matrixStats")
  lognorm <- matrix(1, 1, 4, dimnames = list("g", paste0("s", 1:4)))
  expect_error(fit_mixture(lognorm, cluster = c("1", "1", "2", "2"), n_chains_bimodal = 1),
               "whole number of at least 2")
  expect_error(fit_mixture(lognorm, cluster = c("1", "1", "2", "2"), n_chains_bimodal = 2.5),
               "whole number of at least 2")
})

test_that("call_mixture with bh_q = 1 has no false-discovery filter", {
  fit <- fake_mixture()   # geneC: elpd difference 1 +- 2, fails the default q < 0.1
  strict <- suppressMessages(call_mixture(fit))
  off <- suppressMessages(call_mixture(fit, bh_q = 1))
  diff <- vapply(fit$fits, function(g) g$elpd$diff, numeric(1))
  # every gene with a positive elpd difference passes the evidence check
  expect_identical(is.na(off$not.binarized) | off$not.binarized != "elpd",
                   unname(diff > 0))
  expect_false(identical(off$not.binarized, strict$not.binarized))
  expect_identical(attr(off, "binarization")$bh_q, 1)
  msgs <- testthat::capture_messages(call_mixture(fit, bh_q = 1))
  expect_match(msgs, "no false discovery filter, bh_q = 1", all = FALSE)
})
