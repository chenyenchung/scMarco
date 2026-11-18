qc_fixture <- function() {
  # three clusters x two replicates; geneA ON in cluster 1, geneB in 2, geneC
  # not binarized; one library with background genes for the ambient profile
  cluster <- rep(c("1", "2", "10"), 2)
  lognorm <- rbind(geneA = c(2, 0.1, 0.1, 2.1, 0.2, 0.1),
                   geneB = c(0.1, 1.8, 0.1, 0.2, 1.9, 0.1),
                   geneC = c(0.3, 0.4, 0.2, 0.3, 0.4, 0.3))
  colnames(lognorm) <- paste0("s", 1:6)
  fits <- list(
    geneA = fake_gene_fit(log(c(0.999, 0.004, 0.003)), diff = 30, diff_se = 5),
    geneB = fake_gene_fit(log(c(0.002, 0.995, 0.001)), diff = 25, diff_se = 4,
                          chain_spread = c(0, 0, 0.8)),
    geneC = fake_gene_fit(log(c(0.6, 0.7, 0.5)), diff = 1, diff_se = 2))
  fit <- new_mixture(fits, lognorm, cluster, replicate = rep(c("1", "2"), each = 3),
                     settings = list(n_iter = 500, n_folds = 10, seed = 42, cores = 1,
                                     prior = "range", n_chains_bimodal = 8L, n_chains = 4L))
  cl <- rep(c("1", "2", "10"), each = 4)
  counts <- matrix(5L, 13, 12, dimnames = list(c("geneA", "geneB", "geneC",
    paste0("bg", 1:10)), paste0("c", 1:12)))
  counts["geneA", ] <- rep(c(40L, 0L, 0L), each = 4)
  counts["geneB", ] <- rep(c(0L, 30L, 30L), each = 4)
  amb <- ambient_summary(counts, cl, rho = 0.01)
  list(fit = fit, amb = amb,
       tab = suppressMessages(call_mixture(fit, amb, reject_spread_frac = 0.5)))
}

test_that("qc_report counts genes and calls as binarize() makes them", {
  x <- qc_fixture()
  qc <- qc_report(x$fit, x$tab, x$amb, dataset = "toy", split = "s")
  expect_identical(qc$qc_report_version, QC_REPORT_VERSION)
  expect_identical(qc$genes$total, 3L)
  expect_identical(qc$genes$binarized, 2L)
  expect_identical(qc$genes$not_binarized_evidence, 1L)
  # geneB's chains disagree on cluster 10 (spread 0.8), where it is above ambient too
  expect_identical(qc$convergence$genes_with_cluster_over_cutoff, 1L)
  expect_identical(qc$ambient$calls_changed_chains, 1L)
  # 2 binarized genes x 3 clusters: geneA ON in 1, geneB ON in 2, geneB in 10 ambiguous
  expect_identical(qc$calls$ON, 2L)
  expect_identical(qc$calls$ambiguous, 1L)
  expect_identical(qc$calls$OFF, 3L)
  expect_identical(qc$calls$no_call, 3L)
  expect_identical(qc$dataset$cells, 12L)
  expect_named(qc$dataset$umis_per_library, "1")
  expect_identical(qc$dataset$cells_per_library, list(`1` = 12L))
  expect_identical(qc$settings$n_chains_bimodal, 8L)
  expect_identical(qc$settings$band, 0.3)
  expect_identical(qc$ambient$rho, list(`1` = 0.01))
  expect_length(qc$warnings, 0)
  expect_true(length(qc$limitations) >= 4)
  # the band changes the counts
  wide <- qc_report(x$fit, x$tab, x$amb, band = 0.49)
  expect_identical(wide$calls$ON, 2L)
})

test_that("qc_report writes YAML and a self-contained HTML page without status", {
  x <- qc_fixture()
  base <- file.path(withr::local_tempdir(), "qc", "toy_s")
  expect_message(qc <- qc_report(x$fit, x$tab, x$amb, file = base, dataset = "toy", split = "s"),
                 "Wrote")
  y <- yaml::read_yaml(paste0(base, ".yaml"))
  expect_identical(y$qc_report_version, QC_REPORT_VERSION)
  expect_identical(y$genes$binarized, 2L)
  expect_named(y, c("qc_report_version", "generated", "dataset", "warnings", "settings", "genes",
                    "convergence", "ambient", "calls", "limitations", "provenance"),
               ignore.order = TRUE)
  html <- paste(readLines(paste0(base, ".html")), collapse = "\n")
  text <- gsub("data:image/png;base64,[A-Za-z0-9+/=]+", "", html)   # plots are random text
  for (h in c("Dataset", "Settings in use", "Genes", "Convergence", "Ambient RNA and detection",
              "Calls", "Limitations", "Provenance")) {
    expect_match(html, paste0("<h2>", h, "</h2>"), fixed = TRUE)
  }
  expect_match(html, "What is this", fixed = TRUE)
  expect_match(html, "data:image/png;base64,", fixed = TRUE)
  expect_false(grepl("status|work in progress|WIP", text, ignore.case = TRUE))
  expect_false(any(grepl("status", names(unlist(y)), ignore.case = TRUE)))
})

test_that("qc_report works without an ambient summary and checks its input", {
  x <- qc_fixture()
  tab <- suppressMessages(call_mixture(x$fit, reject_spread_frac = 0.5))
  qc <- qc_report(x$fit, tab)
  expect_null(qc$dataset$cells)
  expect_identical(qc$ambient$off_changed_ambient, 0L)
  expect_match(unlist(qc$warnings), "No ambient summary", all = FALSE)
  expect_error(qc_report(x$fit, tab[1:2, ]), "every gene of fit")
  bare <- tab
  attr(bare, "flags") <- NULL
  expect_error(qc_report(x$fit, bare), "must come from call_mixture")
  expect_error(qc_report(x$fit, tab, band = 0.5), "band must be")
  expect_error(qc_report(x$fit, tab, ambient = list()), "must come from ambient_summary")
})
