root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages(source(file.path(root, "script", "load.R"), chdir = TRUE))

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
  expect_equal(order_cluster_ids(c("NB", "108", "11", "2", "Apop", "1+2", "1")),
               c("1", "1+2", "2", "11", "108", "Apop", "NB"))
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

test_that("migrate_cluster_labels migrates both files and keeps the data", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))

  expect_message(migrate_cluster_labels(p$clusters, p$cells), "Split 2 labels")

  ds <- ds_load("toy", list(label = "Toy", clusters = p$clusters, cells = p$cells))
  expect_equal(ds$clusters, c("9", "12", "154"))
  expect_equal(unname(ds$cluster_labels), c("Tm27(Y) (9)", "12", "Mi4 (154)"))
  expect_equal(rownames(ds$clusters_h$obs), c("F|12", "F|154", "F|9", "M|12", "M|9"))
  expect_equal(unname(ds$prob[, "a"]), p$prob[, 1])
  expect_equal(as.character(ds_cells_obs(ds)$cluster), c("154", "12", "12", "9"))

  val <- validate_h5ad(p$clusters, p$cells)
  expect_length(val$errors, 0)
  expect_false(any(grepl("look like", val$warnings)))

  # a second run refuses: the names are already split out
  expect_error(migrate_cluster_labels(p$clusters), "already has cluster names")
})

test_that("set_cluster_names applies a table and get_cluster_names reads it back", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))
  suppressMessages(migrate_cluster_labels(p$clusters))

  tsv <- file.path(dir, "names.tsv")
  writeLines(c("cluster_id\tname", "12\tT4a", "154\t"), tsv)
  res <- suppressMessages(set_cluster_names(p$clusters, tsv))
  expect_equal(res$name, c("Tm27(Y)", "T4a", NA))
  ds <- ds_load("toy", list(label = "Toy", clusters = p$clusters))
  # 12 named, 154 cleared, 9 untouched
  expect_equal(unname(ds$cluster_labels), c("Tm27(Y) (9)", "T4a (12)", "154"))

  out <- file.path(dir, "export.tsv")
  names_df <- suppressMessages(get_cluster_names(p$clusters, out))
  expect_equal(names_df$name, c("Tm27(Y)", "T4a", NA))
  back <- read.delim(out, colClasses = "character", na.strings = character(0))
  expect_equal(back$cluster_id, c("9", "12", "154"))
  expect_equal(back$name, c("Tm27(Y)", "T4a", ""))

  # a named vector or data.frame works as well as a file
  suppressMessages(set_cluster_names(p$clusters, c(`154` = "Mi4", `12` = NA)))
  expect_equal(get_cluster_names(p$clusters)$name, c("Tm27(Y)", NA, "Mi4"))
  suppressMessages(set_cluster_names(p$clusters, back))
  expect_equal(get_cluster_names(p$clusters)$name, c("Tm27(Y)", "T4a", NA))
  expect_message(set_cluster_names(p$clusters, back), "No changes")

  writeLines(c("cluster_id\tname", "999\tFoo"), tsv)
  expect_error(set_cluster_names(p$clusters, tsv), "999")
  expect_error(set_cluster_names(p$clusters, data.frame(id = "9")), "cluster_id and name")
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
  expect_warning(suppressMessages(migrate_cluster_labels(cpath)),
                 "only their ID is kept: Ghost \\(77\\)")
  ds <- ds_load("ghost", list(label = "Ghost", clusters = cpath))
  expect_equal(unname(ds$cluster_labels), c("12", "77", "Mi4 (154)"))

  tsv <- file.path(dir, "names.tsv")
  writeLines(c("cluster_id\tname", "77\tGhost"), tsv)
  expect_error(set_cluster_names(cpath, tsv), "no rows in any split")
})

test_that("set_cluster_ids gives new IDs and names to both files and keeps the data", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_legacy_pair(dir)))
  tsv <- file.path(dir, "ids.tsv")
  writeLines(c("cluster_id\tnew_id\tname", "12\t5+6\tSPG", "Mi4 (154)\t154\tMi4",
               "Tm27(Y) (9)\t9\t"), tsv)
  expect_message(set_cluster_ids(p$clusters, tsv, cells = p$cells), "Gave 3 clusters new IDs")

  ds <- ds_load("toy_ids", list(label = "Toy", clusters = p$clusters, cells = p$cells))
  expect_equal(ds$clusters, c("5+6", "9", "154"))
  expect_equal(unname(ds$cluster_labels), c("SPG (5+6)", "9", "Mi4 (154)"))
  expect_equal(rownames(ds$clusters_h$obs), c("F|5+6", "F|154", "F|9", "M|5+6", "M|9"))
  expect_equal(unname(ds$prob[, "a"]), p$prob[, 1])
  expect_equal(as.character(ds_cells_obs(ds)$cluster), c("154", "5+6", "5+6", "9"))
  val <- validate_h5ad(p$clusters, p$cells)
  expect_length(val$errors, 0)

  # every cluster must be listed, and new IDs must be unique
  expect_error(set_cluster_ids(p$clusters, data.frame(cluster_id = "9", new_id = "1")),
               "Missing")
  expect_error(set_cluster_ids(p$clusters, data.frame(cluster_id = c("5+6", "9", "154"),
                                                      new_id = c("1", "1", "2"))),
               "Duplicate")
})
