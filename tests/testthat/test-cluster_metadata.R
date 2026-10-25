root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages(source(file.path(root, "script", "load.R"), chdir = TRUE))
source(file.path(root, "src", "ui_helpers.R"))

make_meta_toy <- function(dir) {
  # Four clusters in two splits; cluster 4 exists only in F
  cpath <- file.path(dir, "meta.clusters.h5ad")
  suppressMessages(suppressWarnings(write_clusters_h5ad(
    cpath, c("F", "F", "F", "F", "M", "M", "M"), c("1", "2", "3", "4", "1", "2", "3"),
    c("F", "M"), c("1", "2", "3", "4"), c("a", "b"),
    layers = list(prob = matrix(0.5, 7, 2))
  )))
  cpath
}

meta_tbl <- function() {
  data.frame(cluster_id = c("1", "2", "3", "4"),
             class = c("neuron", "neuron", "neuron", "glia"),
             subclass = c("Mi", "Tm", "Tm", ""))
}

test_that("set_cluster_metadata writes levels that get_cluster_metadata reads back", {
  dir <- withr::local_tempdir()
  cpath <- make_meta_toy(dir)
  expect_null(suppressMessages(get_cluster_metadata(cpath)))
  expect_message(set_cluster_metadata(cpath, meta_tbl()), "class > subclass")

  tsv <- file.path(dir, "meta.tsv")
  out <- suppressMessages(get_cluster_metadata(cpath, tsv))
  expect_equal(names(out), c("cluster_id", "class", "subclass"))
  expect_equal(out$subclass, c("Mi", "Tm", "Tm", NA))
  expect_equal(read.delim(tsv, colClasses = "character")$class, meta_tbl()$class)

  ds <- ds_load("meta", list(label = "Meta", clusters = cpath))
  expect_null(ds$cluster_meta)
  ds <- ds_load("meta_sub", list(label = "Meta", clusters = cpath,
                                 subset_by = list("subclass", "class")))
  expect_equal(names(ds$cluster_meta), c("subclass", "class"))
  expect_equal(rownames(ds$cluster_meta), c("1", "2", "3", "4"))
  expect_true(ds_check_config(list(meta = list(clusters = cpath, subset_by = "class"))))
  expect_error(ds_check_config(list(meta = list(clusters = cpath, subset_by = "nope"))),
               "subset_by")
  expect_length(validate_h5ad(cpath)$errors, 0)

  # A TSV with fewer clusters replaces the whole table; levels keep TSV order
  writeLines(c("cluster_id\tsubclass\tclass", "2\tTm\tneuron"), tsv)
  suppressMessages(set_cluster_metadata(cpath, tsv))
  out <- suppressMessages(get_cluster_metadata(cpath))
  expect_equal(names(out), c("cluster_id", "subclass", "class"))
  expect_equal(out$class, c(NA, "neuron", NA, NA))

  suppressMessages(set_cluster_metadata(cpath, NULL))
  expect_null(suppressMessages(get_cluster_metadata(cpath)))
})

test_that("bad metadata tables are rejected", {
  dir <- withr::local_tempdir()
  cpath <- make_meta_toy(dir)
  expect_error(set_cluster_metadata(cpath, data.frame(cluster_id = "9", class = "x")),
               "not in")
  expect_error(set_cluster_metadata(cpath, data.frame(cluster_id = c("1", "1"), class = "x")),
               "more than once")
  expect_error(set_cluster_metadata(cpath, data.frame(cluster_id = "1")),
               "at least one column")
  expect_error(set_cluster_metadata(cpath, data.frame(cluster_id = "1", `a/b` = "x",
                                                      check.names = FALSE)),
               "Column names")
})

test_that("metadata survives new cluster IDs", {
  dir <- withr::local_tempdir()
  cpath <- make_meta_toy(dir)
  suppressMessages(set_cluster_metadata(cpath, meta_tbl()))
  suppressMessages(set_cluster_ids(cpath, data.frame(cluster_id = c("1", "2", "3", "4"),
                                                     new_id = c("10", "20", "30", "40"))))
  out <- suppressMessages(get_cluster_metadata(cpath))
  expect_equal(out$cluster_id, c("10", "20", "30", "40"))
  expect_equal(out$subclass, c("Mi", "Tm", "Tm", NA))
})

test_that("filters cascade and define the search universe", {
  dir <- withr::local_tempdir()
  cpath <- make_meta_toy(dir)
  suppressMessages(set_cluster_metadata(cpath, meta_tbl()))
  ds <- ds_load("meta_cascade", list(label = "Meta", clusters = cpath,
                                     subset_by = c("class", "subclass")))

  none <- list(class = character(0), subclass = character(0))
  expect_equal(meta_universe(ds, none), ds$clusters)
  expect_equal(meta_level_choices(ds, none, "class"), c("glia", "neuron"))
  expect_equal(meta_level_choices(ds, none, "subclass"), c("Mi", "Tm"))

  glia <- list(class = "glia", subclass = character(0))
  expect_equal(meta_universe(ds, glia), "4")
  expect_equal(meta_level_choices(ds, glia, "subclass"), character(0))

  tm <- list(class = "neuron", subclass = "Tm")
  expect_equal(meta_universe(ds, tm), c("2", "3"))
  expect_equal(describe_scope(ds, tm),
               "Searching within class = neuron > subclass = Tm (2 of 4 clusters).")
  expect_null(describe_scope(ds, none))

  b <- restrict_clusters(ds_binarize(ds, 0.1), meta_universe(ds, tm))
  expect_equal(b$clusters, c("2", "3"))
})

test_that("metadata that differs between splits of a cluster is an error", {
  dir <- withr::local_tempdir()
  cpath <- make_meta_toy(dir)
  rewrite_obs(cpath, function(obs) {
    obs$meta_class <- factor(c("a", "a", "a", "a", "b", "a", "a"))
    obs
  })
  expect_error(ds_load("meta_bad", list(label = "Meta", clusters = cpath, subset_by = "class")),
               "differs between splits")
  expect_match(validate_h5ad(cpath)$errors, "meta_class", all = FALSE)
})
