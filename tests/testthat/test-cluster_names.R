root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))
suppressMessages(source(file.path(root, "script", "lib_write_h5ad.R")))

run_script <- function(script, args) {
  # some calls are expected to fail; their exit status is checked below
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(root, "script", script), args),
    stdout = TRUE, stderr = TRUE,
    env = "RENV_CONFIG_NAMESPACES_CHECK=FALSE"
  ))
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0, log = out)
}

make_legacy_pair <- function(dir) {
  # Two splits; labels embed IDs the legacy way, and one cluster is unnamed
  clusters <- c("12", "Mi4 (154)", "Tm27(Y) (9)")
  obs_split <- c("F", "F", "F", "M", "M")
  obs_cluster <- c("12", "Mi4 (154)", "Tm27(Y) (9)", "12", "Tm27(Y) (9)")
  prob <- matrix(seq(0.1, 1, length.out = 10), 5, 2)
  cpath <- file.path(dir, "toy.clusters.h5ad")
  write_clusters_h5ad(cpath, obs_split, obs_cluster, c("F", "M"), clusters,
                      c("a", "b"), layers = list(prob = prob))
  xpath <- file.path(dir, "toy.cells.h5ad")
  write_cells_h5ad(
    xpath, Matrix::Matrix(matrix(1, 4, 2), sparse = TRUE),
    barcodes = paste0("bc", 1:4), cell_split = c("F", "F", "M", "M"),
    cell_cluster = c("Mi4 (154)", "12", "12", "Tm27(Y) (9)"),
    splits = c("F", "M"), clusters = clusters, genes = c("a", "b")
  )
  list(clusters = cpath, cells = xpath, prob = prob)
}

test_that("legacy labels split into ID and name", {
  parts <- split_cluster_labels(c(
    "Mi4 (154)", "108", "Tm27(Y) (45)", "G/LQ1 (in pupae) (186)", "Apop"
  ))
  expect_equal(parts$id, c("154", "108", "45", "186", "Apop"))
  expect_equal(parts$name, c("Mi4", NA, "Tm27(Y)", "G/LQ1 (in pupae)", NA))
  expect_equal(order_cluster_ids(c("NB", "108", "11", "2", "Apop")),
               c("2", "11", "108", "Apop", "NB"))
})

test_that("labels show name (id), or the bare ID when unnamed", {
  labels <- format_cluster_labels(c("154", "12"), c("Mi4", NA))
  expect_equal(unname(labels), c("Mi4 (154)", "12"))
  ds <- list(cluster_labels = labels)
  expect_equal(cluster_label(ds, c("12", "154")), c("12", "Mi4 (154)"))
  expect_equal(cluster_choices(ds, "154"), c("Mi4 (154)" = "154"))
})

test_that("a dataset without names uses its IDs as labels", {
  fixture <- file.path(root, "tests", "fixtures", "example.clusters.h5ad")
  ds <- ds_load("example", list(label = "Example", clusters = fixture))
  expect_true(all(is.na(ds$cluster_names)))
  expect_equal(unname(ds$cluster_labels), ds$clusters)
})

test_that("--split-labels migrates both files and keeps the data", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))

  res <- run_script("set_cluster_names.R",
                    c("--clusters", p$clusters, "--cells", p$cells, "--split-labels", "TRUE"))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))

  ds <- ds_load("toy", list(label = "Toy", clusters = p$clusters, cells = p$cells))
  expect_equal(ds$clusters, c("9", "12", "154"))
  expect_equal(unname(ds$cluster_labels), c("Tm27(Y) (9)", "12", "Mi4 (154)"))
  expect_equal(rownames(ds$clusters_h$obs), c("F|12", "F|154", "F|9", "M|12", "M|9"))
  expect_equal(unname(ds$prob[, "a"]), p$prob[, 1])
  expect_equal(as.character(ds_cells_obs(ds)$cluster), c("154", "12", "12", "9"))

  res <- run_script("validate_h5ad.R", c("--clusters", p$clusters, "--cells", p$cells))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))
  expect_false(any(grepl("look like", res$log)))

  # a second run refuses: the names are already split out
  res <- run_script("set_cluster_names.R",
                    c("--clusters", p$clusters, "--split-labels", "TRUE"))
  expect_false(res$ok)
})

test_that("--names applies a table and --export reads it back", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))
  run_script("set_cluster_names.R", c("--clusters", p$clusters, "--split-labels", "TRUE"))

  tsv <- file.path(dir, "names.tsv")
  writeLines(c("cluster_id\tname", "12\tT4a", "154\t"), tsv)
  res <- run_script("set_cluster_names.R", c("--clusters", p$clusters, "--names", tsv))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))
  ds <- ds_load("toy", list(label = "Toy", clusters = p$clusters))
  # 12 named, 154 cleared, 9 untouched
  expect_equal(unname(ds$cluster_labels), c("Tm27(Y) (9)", "T4a (12)", "154"))

  out <- file.path(dir, "export.tsv")
  res <- run_script("set_cluster_names.R", c("--clusters", p$clusters, "--export", out))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))
  back <- read.delim(out, colClasses = "character", na.strings = character(0))
  expect_equal(back$cluster_id, c("9", "12", "154"))
  expect_equal(back$name, c("Tm27(Y)", "T4a", ""))

  writeLines(c("cluster_id\tname", "999\tFoo"), tsv)
  res <- run_script("set_cluster_names.R", c("--clusters", p$clusters, "--names", tsv))
  expect_false(res$ok)
  expect_true(any(grepl("999", res$log)))
})

test_that("names that differ between splits of a cluster are an error", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))
  rewrite_obs(p$clusters, function(obs) {
    obs$cluster_name <- factor(c("A", NA, NA, "B", NA))
    obs
  })
  expect_error(ds_load("toy2", list(label = "Toy", clusters = p$clusters)),
               "differs between splits")
})

test_that("clusters with no rows in any split keep only their ID", {
  dir <- withr::local_tempdir()
  cpath <- file.path(dir, "ghost.clusters.h5ad")
  suppressMessages(suppressWarnings(write_clusters_h5ad(
    cpath, c("F", "F"), c("Mi4 (154)", "12"), "F", c("12", "Ghost (77)", "Mi4 (154)"),
    c("a", "b"), layers = list(prob = matrix(0.5, 2, 2))
  )))
  res <- run_script("set_cluster_names.R", c("--clusters", cpath, "--split-labels", "TRUE"))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))
  expect_true(any(grepl("only their ID is kept: Ghost \\(77\\)", res$log)))
  ds <- ds_load("ghost", list(label = "Ghost", clusters = cpath))
  expect_equal(unname(ds$cluster_labels), c("12", "77", "Mi4 (154)"))

  tsv <- file.path(dir, "names.tsv")
  writeLines(c("cluster_id\tname", "77\tGhost"), tsv)
  res <- run_script("set_cluster_names.R", c("--clusters", cpath, "--names", tsv))
  expect_false(res$ok)
  expect_true(any(grepl("no rows in any split", res$log)))
})
