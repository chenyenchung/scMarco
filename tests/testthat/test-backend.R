root <- file.path(testthat::test_path(), "..", "..")
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))
suppressMessages(source(file.path(root, "script", "lib_write_h5ad.R")))

fixture <- file.path(root, "tests", "fixtures", "example.clusters.h5ad")

test_that("the SQLite-derived fixture loads as a dataset", {
  ds <- ds_load("example", list(label = "Example", clusters = fixture))
  expect_equal(ds$splits, c("P15", "P30", "P40", "P50", "P70", "Adult"))
  expect_true(ds$ordered_splits)
  expect_equal(ds$split_label, "stage")
  expect_true(ds$has_mean)
  expect_false(ds$has_pct)
  expect_null(ds$cells_h)
  expect_equal(ncol(ds$prob), 300)
  expect_equal(nrow(ds$prob), length(ds$obs_split))
  # no bimodality information -> every gene is analysed
  expect_identical(ds$analysis_genes, ds$all_genes)

  b <- ds_binarize(ds, 0.5)
  expect_s3_class(b, "scmarco_binarized")
  expect_identical(ds_binarize(ds, 0.5), b)
})

test_that("config overrides win over uns", {
  ds <- ds_load("example2", list(
    label = "Example", clusters = fixture, split_label = "age",
    split_order = c("Adult", "P70", "P50", "P40", "P30", "P15")
  ))
  expect_equal(ds$split_label, "age")
  expect_equal(ds$splits[1], "Adult")
})

make_pair <- function(dir) {
  genes <- c("a", "b", "c")
  splits <- c("F", "M")
  clusters <- c("x", "y")
  obs_split <- c("F", "F", "M")
  obs_cluster <- c("x", "y", "x")
  prob <- rbind(c(0.9, 0.1, NA), c(0.2, 0.8, NA), c(0.7, 0.6, NA))
  bim <- matrix(c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE), 3, 2,
                dimnames = list(genes, splits))
  cpath <- file.path(dir, "toy.clusters.h5ad")
  write_clusters_h5ad(
    cpath, obs_split, obs_cluster, splits, clusters, genes,
    layers = list(prob = prob, pct_pos = prob / 2),
    is_bimodal = c(TRUE, TRUE, FALSE), bimodal_by_split = bim,
    n_cells = c(2, 1, 2), split_label = "sex"
  )
  x <- Matrix::sparseMatrix(
    i = c(1, 2, 4, 5, 3), j = c(1, 1, 1, 2, 3), x = c(1.5, 2, 0.5, 3, 1),
    dims = c(5, 3)
  )
  red <- matrix(c(1:5, 5:1), 5, 2)
  red_sub <- matrix(c(NaN, 1, 2, NaN, NaN, NaN, 3, 4, NaN, NaN), 5, 2)
  xpath <- file.path(dir, "toy.cells.h5ad")
  write_cells_h5ad(
    xpath, x, barcodes = paste0("bc", 1:5),
    cell_split = c("F", "F", "F", "M", "M"),
    cell_cluster = c("x", "x", "y", "x", "x"),
    splits = splits, clusters = clusters, genes = genes,
    reductions = list(umap_full = red, sub = red_sub)
  )
  list(clusters = cpath, cells = xpath, x = x)
}

test_that("a written pair round-trips through the reader", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_pair(dir)))
  ds <- ds_load("toy", list(label = "Toy", clusters = p$clusters, cells = p$cells))

  expect_equal(ds$splits, c("F", "M"))
  expect_equal(ds$split_label, "sex")
  expect_equal(ds$analysis_genes, c("a", "b"))
  expect_true(ds$has_pct)
  expect_equal(unname(ds$bimodal["c", ]), c(FALSE, FALSE))
  expect_true(is.na(ds$prob[1, "c"]))
  expect_equal(unname(ds$prob[, "a"]), c(0.9, 0.2, 0.7))

  long <- ds_cluster_layer_long(ds, "pct_pos", c("a", "b"), clusters = "x")
  expect_equal(nrow(long), 4)
  expect_equal(long$value[long$gene == "a" & long$split == "M"], 0.35)

  # per-gene CSC reads match the matrix we wrote
  for (j in 1:3) {
    expect_equal(ds_cell_expression(ds, c("a", "b", "c")[j]), p$x[, j])
  }
  expect_setequal(ds$reductions, c("X_umap_full", "X_sub"))
  sub <- ds_reduction(ds, "X_sub")
  expect_equal(dim(sub), c(5, 2))
  expect_equal(sum(!is.na(sub[, 1])), 2)
  cells_obs <- ds_cells_obs(ds)
  expect_equal(as.character(cells_obs$cluster), c("x", "x", "y", "x", "x"))

  b <- ds_binarize(ds, 0.5)
  expect_false(b$present["y", "M"])
  expect_false(any(b$on["c", , ]))
})
