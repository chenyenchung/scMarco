# Converters from the two legacy formats, on tiny synthetic inputs

write_h5v1 <- function(path) {
  # A scMarco-Ding-lab v1.0 db.h5 with one dataset "toy": 3 genes, clusters
  # a and b, splits Male and Female (b absent in Female), geneC non-bimodal
  # (-1) everywhere, and a UMAP per split
  rhdf5::h5createFile(path)
  rhdf5::h5createGroup(path, "toy")
  for (g in c("probability", "pct_pos", "barcodes", "clusters", "lognorm", "reduction")) {
    rhdf5::h5createGroup(path, paste0("toy/", g))
  }
  w <- function(obj, x) rhdf5::h5write(x, path, paste0("toy/", obj))
  w("gene_names", c("geneA", "geneB", "geneC"))
  w("cluster_names", c("a", "b"))
  w("splits", c("Male", "Female"))
  w("bimodal_genes", c(0L, 1L))
  cells <- list(Male = c("a", "a", "b", "b"), Female = c("a", "a"))
  prob <- list(Male = cbind(a = c(0.9, 0.1, -1), b = c(0.2, 0.95, -1)),
               Female = cbind(a = c(0.8, 0.05, -1), b = c(NA, NA, NA)))
  for (s in names(cells)) {
    w(paste0("probability/", s), prob[[s]])
    w(paste0("pct_pos/", s), pmax(prob[[s]], 0) / 2)
    n <- length(cells[[s]])
    w(paste0("barcodes/", s), paste0(s, "_cell", seq_len(n)))
    w(paste0("clusters/", s), cells[[s]])
    # CSR genes x cells: geneA in the a cells, geneB in the b cells
    x <- Matrix::sparseMatrix(i = seq_len(n), j = ifelse(cells[[s]] == "a", 1, 2),
                              x = seq_len(n) / 2, dims = c(n, 3))
    rhdf5::h5createGroup(path, paste0("toy/lognorm/", s))
    w(paste0("lognorm/", s, "/indices"), x@i)
    w(paste0("lognorm/", s, "/indptr"), x@p)
    w(paste0("lognorm/", s, "/data"), x@x)
    w(paste0("lognorm/", s, "/shape"), c(3, n))
    rhdf5::h5createGroup(path, paste0("toy/reduction/", s))
    rhdf5::h5createGroup(path, paste0("toy/reduction/", s, "/umap"))
    w(paste0("reduction/", s, "/umap/coords"), matrix(seq_len(2 * n), n))
  }
  path
}

test_that("a v1.0 db.h5 dataset converts and loads", {
  h5 <- write_h5v1(tempfile(fileext = ".h5"))
  out <- suppressWarnings(suppressMessages(
    h5v1_to_h5ad(h5, "toy", out_dir = tempfile(), split_label = "sex")
  ))
  res <- validate_h5ad(out$clusters, out$cells)
  expect_length(res$errors, 0)

  ds <- ds_load("v1", list(label = "v1", clusters = out$clusters, cells = out$cells))
  expect_setequal(ds$splits, c("Male", "Female"))
  expect_equal(ds$split_label, "sex")
  expect_equal(ds$analysis_genes, c("geneA", "geneB"))
  b <- ds_binarize(ds, 0.4)
  expect_false(b$present["b", "Female"])
  expect_true(b$on["geneA", "a", "Male"])
  expect_true(b$on["geneB", "b", "Male"])
  expect_equal(from_cluster(b, "b", character(0), 1), c(geneB = 0L))
  expect_equal(ds$reductions, "X_umap")
  expect_equal(sum(ds_cell_expression(ds, "geneA") > 0), 4)
})

test_that("a legacy SQLite table converts and loads", {
  skip_if_not_installed("RSQLite")
  db <- tempfile(fileext = ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  long <- expand.grid(gene = c("geneA", "geneB", "geneC"), cluster = c("1", "2"),
                      stage = c("P15", "Adult"), type = c("prob", "lognorm"),
                      stringsAsFactors = FALSE)
  long$value <- ifelse(long$type == "prob",
                       ifelse(long$gene == "geneA" & long$cluster == "1", 0.99, 0.01),
                       2)
  DBI::dbWriteTable(con, "toy", long)
  DBI::dbDisconnect(con)

  expect_error(sqlite_to_h5ad(db, "toy", stages = "P15", out_dir = tempdir()),
               "does not match")
  path <- suppressWarnings(suppressMessages(
    sqlite_to_h5ad(db, "toy", stages = c("P15", "Adult"), out_dir = tempfile())
  ))
  expect_length(validate_h5ad(path)$errors, 0)
  ds <- ds_load("sqlite", list(label = "sqlite", clusters = path))
  expect_equal(ds$splits, c("P15", "Adult"))
  expect_true(ds$ordered_splits)
  b <- ds_binarize(ds, 0.4)
  expect_equal(from_cluster(b, "1", c("P15", "Adult"), 2), c(geneA = 0L))
})
