# The contract between fitting, the dataset files and the app: every fit
# object and binarized dataset records its schema and binarizer versions.

toy_copy <- function() {
  path <- tempfile(fileext = ".clusters.h5ad")
  file.copy(system.file("extdata", "toy.clusters.h5ad", package = "scMarco"), path)
  path
}

toy_tables <- function(path, version = NULL) {
  # One table per split of the toy dataset, tagged as call_mixture() does
  ch <- h5ad_open(path)
  tabs <- lapply(c(Male = "Male", Female = "Female"), function(s) {
    present <- unique(as.character(ch$obs$cluster[ch$obs$split == s]))
    tab <- as.data.frame(matrix(0.5, length(ch$genes), length(present),
      dimnames = list(ch$genes, present)), check.names = FALSE)
    tab$bimodal.on_level <- 2
    tab$bimodal.off_level <- 0
    tab$bimodal.sd <- 0.5
    tab$bimodal.pi <- 0.5
    tab$exprType <- "bimodal"
    tab
  })
  if (!is.null(version)) {
    # splits not named in `version` stay untagged, as tables not from call_mixture()
    for (s in intersect(names(tabs), names(version))) {
      attr(tabs[[s]], "scmarco_binarizer") <- version[[s]]
    }
  }
  tabs
}

set_version <- function(path, v) {
  obj <- "uns/scmarco/binarization/version"
  if (h5ad_exists(path, obj)) rhdf5::h5delete(path, obj)
  if (!is.null(v)) rhdf5::h5write(v, path, obj)
}

test_that("fit objects record the schema and binarizer versions", {
  fit <- new_mixture(list(), matrix(0, 0, 2), c("1", "2"), c("1", "1"),
                     settings = list(n_iter = 500, n_folds = 10, seed = 42, cores = 1))
  expect_identical(fit$schema_version, SCMARCO_FORMAT_VERSION)
  expect_identical(fit$binarizer_version, SCMARCO_BINARIZER_VERSION)
  expect_output(print(fit), paste0("binarizer ", SCMARCO_BINARIZER_VERSION))
  expect_silent(check_mixture(fit))

  old <- fit
  old$schema_version <- NULL
  expect_error(check_mixture(old), "before fits recorded their versions")
  old <- fit
  old$schema_version <- "9.0"
  expect_error(check_mixture(old), "update scMarco")
  # an older binarizer is read with a warning recommending a refit
  old <- fit
  old$binarizer_version <- "0.9"
  expect_warning(check_mixture(old), "fit by binarizer 0.9.*refitting with fit_mixture\\(\\) is recommended")
  # a newer or invalid one is refused
  new <- fit
  new$binarizer_version <- "9.9"
  expect_error(check_mixture(new), "newer than this scMarco's")
  bad <- fit
  bad$binarizer_version <- "x"
  expect_error(check_mixture(bad), "no valid binarizer version")
  # chunks from different binarizers are never merged
  expect_error(merge_mixture_fits(fit, old), "different binarizers cannot be merged")
})

test_that("merged fits must share their settings", {
  a <- new_mixture(list(g1 = 1), matrix(0, 1, 2, dimnames = list("g1", c("s1", "s2"))),
                   c("1", "2"), c("1", "1"),
                   settings = list(n_iter = 500, n_folds = 10, seed = 42, cores = 4))
  b <- a
  b$fits <- list(g2 = 1)
  rownames(b$lognorm) <- "g2"
  b$settings$cores <- 1
  expect_identical(merge_mixture_fits(a, b)$settings, a$settings)
  b$settings$n_iter <- 1000
  expect_error(merge_mixture_fits(a, b), "share their settings")
})

test_that("call_mixture tags its table with the binarizer version", {
  skip_if_not_installed("mclust")
  skip_if_not_installed("matrixStats")
  fit <- new_mixture(
    list(geneA = list(bimodal = list(log_p_on_cluster = log(c(0.99, 0.01)), on_level = 2,
                                     off_level = 0, sd = 0.1, pi = 0.5),
                      unimodal = list(level = 1, sd = 1), elpd = list(diff = 9, diff_se = 1))),
    rbind(geneA = c(s1 = 2, s2 = 0.1, s3 = 2.1, s4 = 0.2)), rep(c("1", "2"), 2),
    rep(c("1", "2"), each = 2),
    settings = list(n_iter = 500, n_folds = 10, seed = 42, cores = 1))
  tab <- suppressMessages(call_mixture(fit))
  expect_identical(attr(tab, "scmarco_binarizer"), SCMARCO_BINARIZER_VERSION)
  path <- tempfile(fileext = ".rds")
  saveRDS(list(Male = tab), path)
  expect_identical(table_binarizer(readRDS(path)$Male), SCMARCO_BINARIZER_VERSION)
})

test_that("add_probability records the tables' binarizer, not the package's", {
  path <- toy_copy()
  v <- SCMARCO_BINARIZER_VERSION
  suppressMessages(add_probability(path, toy_tables(path, c(Male = v, Female = v))))
  expect_identical(h5ad_open(path)$uns$scmarco$binarization$version, v)

  # tables without a binarizer are refused, and the dataset is left as it was
  expect_error(suppressMessages(add_probability(path, toy_tables(path))),
               "do not record a scMarco binarizer.*external = TRUE")
  expect_error(suppressMessages(add_probability(path, toy_tables(path, c(Male = v)))),
               "split\\(s\\) Female do not record")
  expect_identical(h5ad_open(path)$uns$scmarco$binarization$version, v)
  # ... unless asked for: every split is then external
  suppressMessages(add_probability(path, toy_tables(path), external = TRUE))
  expect_identical(h5ad_open(path)$uns$scmarco$binarization$version, BINARIZER_EXTERNAL)
  expect_warning(suppressMessages(
    add_probability(path, toy_tables(path, c(Male = v)), external = TRUE)),
    "although the tables of Male record binarizer")
  expect_identical(h5ad_open(path)$uns$scmarco$binarization$version, BINARIZER_EXTERNAL)

  expect_error(add_probability(path, toy_tables(path, c(Male = v, Female = "0.1"))),
               "different binarizers")
  expect_error(add_probability(path, toy_tables(path, c(Male = "99.0", Female = "99.0"))),
               "newer than this scMarco")
  expect_warning(suppressMessages(
    add_probability(path, toy_tables(path, c(Male = "0.1", Female = "0.1")))),
    "Consider refitting")
})

test_that("a prob layer cannot be written without a binarization record", {
  path <- tempfile(fileext = ".clusters.h5ad")
  args <- list(path, "F", "1", "F", "1", "a", layers = list(prob = matrix(0.5, 1, 1)))
  expect_error(do.call(write_clusters_h5ad, args), "binarization must be a list")
  expect_error(do.call(write_clusters_h5ad,
                       c(args, list(binarization = test_binarization("v1")))),
               "version number")
})

test_that("the validator and the app enforce the binarizer version", {
  path <- toy_copy()
  expect_length(validate_h5ad(path)$errors, 0)
  cfg <- list(toy = list(label = "Toy", clusters = path))

  set_version(path, NULL)
  expect_match(validate_h5ad(path)$errors, "lacks version", all = FALSE)
  expect_error(ds_check_config(cfg), "does not record its binarizer version")
  expect_error(ds_load(paste0("toy", runif(1)), cfg$toy), "binarizer version")

  set_version(path, "99.0")
  expect_match(validate_h5ad(path)$errors, "newer than this scMarco", all = FALSE)
  expect_error(ds_check_config(cfg), "update scMarco")

  set_version(path, "0.1")
  expect_match(validate_h5ad(path)$warnings, "consider refitting", all = FALSE)
  expect_warning(ds_check_config(cfg), "Consider refitting")

  set_version(path, SCMARCO_BINARIZER_VERSION)
  expect_silent(ds_check_config(cfg))
  ds <- ds_load(paste0("toy", runif(1)), cfg$toy)
  expect_identical(ds$binarization$version, SCMARCO_BINARIZER_VERSION)
})

test_that("migrated legacy datasets are recorded as external", {
  skip_if_not_installed("RSQLite")
  db <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  long <- expand.grid(gene = c("geneA", "geneB"), cluster = c("1", "2"),
                      stage = c("P15", "Adult"), type = c("prob", "lognorm"),
                      stringsAsFactors = FALSE)
  long$value <- 0.5
  DBI::dbWriteTable(con, "toy", long)
  DBI::dbDisconnect(con)
  out <- tempfile(fileext = ".clusters.h5ad")
  suppressMessages(sqlite_to_h5ad(db, "toy", c("P15", "Adult"), out = out))
  b <- h5ad_open(out)$uns$scmarco$binarization
  expect_identical(b$version, BINARIZER_EXTERNAL)
  expect_length(validate_h5ad(out)$errors, 0)
})
