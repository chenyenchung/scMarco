fixture <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")

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
  expect_false(any(grepl("unimodal", names(gene_choices(ds, ds$all_genes)))))

  b <- ds_binarize(ds, 0)
  expect_s3_class(b, "scmarco_binarized")
  expect_identical(ds_binarize(ds, 0), b)
})

test_that("config overrides win over uns", {
  ds <- ds_load("example2", list(
    label = "Example", clusters = fixture, split_label = "age",
    split_order = c("Adult", "P70", "P50", "P40", "P30", "P15")
  ))
  expect_equal(ds$split_label, "age")
  expect_equal(ds$splits[1], "Adult")
})

test_that("the trend view defaults to line for ordered splits and follows the config", {
  ds <- ds_load("example", list(label = "Example", clusters = fixture))
  expect_equal(ds$trend, "line")
  ds <- ds_load("example_unordered", list(label = "Example", clusters = fixture,
                                          ordered_splits = FALSE))
  expect_null(ds$trend)
  ds <- ds_load("example_both", list(label = "Example", clusters = fixture,
                                     ordered_splits = FALSE, trend = "both"))
  expect_equal(ds$trend, "both")
})

test_that("split labels get a plural that reads right", {
  expect_equal(plural_label("stage"), "stages")
  expect_equal(plural_label("sex"), "sexes")
  expect_equal(plural_label("species"), "species")
  expect_equal(plural_label("batch"), "batches")
  ds <- ds_load("example_plural", list(label = "Example", clusters = fixture))
  expect_equal(ds$split_plural, "stages")
  ds <- ds_load("example_plural2", list(label = "Example", clusters = fixture,
                                        split_label = "genotype",
                                        split_label_plural = "genotypes tested"))
  expect_equal(ds$split_plural, "genotypes tested")
})

test_that("the config check rejects trend views the files cannot back", {
  expect_silent(ds_check_config(list(a = list(clusters = fixture, trend = "line"))))
  expect_error(ds_check_config(list(a = list(clusters = fixture, trend = "lines"))),
               "trend must be one of")
  # The fixture has no cells file
  expect_error(ds_check_config(list(a = list(clusters = fixture, trend = "box"))),
               "needs a cells file")
})

test_that("the config check rejects a split_order that is not the file's splits", {
  splits <- c("P15", "P30", "P40", "P50", "P70", "Adult")
  expect_silent(ds_check_config(list(a = list(clusters = fixture, split_order = rev(splits)))))
  for (bad in list(splits[-1], c(splits, "Extra"), c(splits[1], splits))) {
    expect_error(ds_check_config(list(a = list(clusters = fixture, split_order = bad))),
                 "split_order must list each split")
  }
})

test_that("the config check rejects a file in another format_version", {
  for (v in c("1.0", "3.0")) {
    path <- tempfile(fileext = ".clusters.h5ad")
    file.copy(fixture, path)
    rhdf5::h5delete(path, "uns/scmarco/format_version")
    rhdf5::h5write(v, path, "uns/scmarco/format_version")
    expect_error(ds_check_config(list(a = list(clusters = path))),
                 if (v == "1.0") "rebuild the dataset" else "update scMarco")
  }
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
    n_cells = c(2, 1, 2), split_label = "sex", binarization = test_binarization()
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
  # unimodal genes stay selectable outside marker discovery, labeled as such
  expect_equal(gene_choices(ds, ds$all_genes),
               c(a = "a", b = "b", `c (unimodal)` = "c"))
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

  b <- ds_binarize(ds, 0)
  expect_false(b$present["y", "M"])
  expect_false(any(b$on["c", , ]))
})

test_that("unordered splits with a cells file get the box view by default", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_pair(dir)))
  ds <- ds_load("toy_view", list(label = "Toy", clusters = p$clusters, cells = p$cells))
  expect_false(ds$ordered_splits)
  expect_equal(ds$trend, "box")
  ds <- ds_load("toy_view2", list(label = "Toy", clusters = p$clusters, cells = p$cells,
                                  trend = "line"))
  expect_equal(ds$trend, "line")
})

test_that("cached() keeps the most recent entries and evaluates only on a miss", {
  env <- new.env(parent = emptyenv())
  n <- 0
  get <- function(key) cached(env, key, size = 2, {n <<- n + 1; toupper(key)})
  expect_equal(get("a"), "A")
  expect_equal(get("a"), "A")
  expect_equal(n, 1)
  get("b"); get("c")
  expect_setequal(ls(env), c("b", "c"))
  expect_equal(get("b"), "B")
  expect_equal(n, 3)
  # a failed computation stores nothing
  expect_error(cached(env, "d", 2, stop("no")), "no")
  expect_setequal(ls(env), c("b", "c"))
})

test_that("datasets keep two bands and their recent genes' cell expression", {
  dir <- withr::local_tempdir()
  suppressMessages(suppressWarnings(p <- make_pair(dir)))
  ds <- ds_load("toy-cache", list(label = "Toy", clusters = p$clusters, cells = p$cells))
  b1 <- ds_binarize(ds, 0.4)
  ds_binarize(ds, 0.3); ds_binarize(ds, 0.2)
  expect_setequal(ls(ds$binarized), format(c(0.2, 0.3), nsmall = 3))
  expect_identical(ds_binarize(ds, 0.4), b1)

  reads <- 0
  local_mocked_bindings(h5ad_x_gene = function(h, gene) {
    reads <<- reads + 1
    p$x[, match(gene, c("a", "b", "c"))]
  })
  for (i in 1:3) expect_equal(ds_cell_expression(ds, "b"), p$x[, 2])
  expect_equal(reads, 1)
  for (g in rep(c("a", "c"), 5)) ds_cell_expression(ds, g)
  expect_equal(reads, 3)
  expect_lte(length(ls(ds$cell_expression)), CELL_EXPRESSION_CACHE)
})
