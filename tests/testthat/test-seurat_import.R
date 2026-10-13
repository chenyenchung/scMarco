root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))

run_script <- function(script, args) {
  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(root, "script", script), args),
    stdout = TRUE, stderr = TRUE,
    env = "RENV_CONFIG_NAMESPACES_CHECK=FALSE"
  )
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0, log = out)
}

test_that("a Seurat object + mixture-model table imports and validates", {
  skip_if_not_installed("SeuratObject")
  suppressPackageStartupMessages(library(SeuratObject))
  dir <- withr::local_tempdir()

  set.seed(1)
  genes <- c("geneA", "geneB", "geneC", "geneD")
  n <- 40
  counts <- matrix(rpois(4 * n, 2), 4, n, dimnames = list(genes, paste0("bc", 1:n)))
  counts["geneA", 1:10] <- 20 # high in Male@x
  obj <- CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  ident <- rep(c("Male@x", "Male@y", "Female@x", "Female@y"), each = 10)
  ident[31:40] <- "Female@x" # Female@y is absent
  Idents(obj) <- factor(ident)
  emb <- matrix(rnorm(2 * n), n, 2, dimnames = list(Cells(obj), c("umap_1", "umap_2")))
  obj[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "umap_", assay = "RNA")
  saveRDS(obj, file.path(dir, "obj.rds"))

  prob <- data.frame(
    `Male@x` = c(0.9, 0.1, 0.2, 0.8), `Male@y` = c(0.1, 0.9, 0.2, 0.8),
    `Female@x` = c(0.9, 0.1, 0.2, 0.8),
    check.names = FALSE, row.names = genes
  )
  prob$exprType <- c("bimodal", "bimodal", "bimodal", "unimodal")
  for (qc in c("diff.elpd", "diff.elpd.se", "bimodal.on_level", "bimodal.off_level",
               "bimodal.sd", "bimodal.pi", "unimodal.level", "unimodal.sd")) {
    prob[[qc]] <- 0
  }
  saveRDS(prob, file.path(dir, "prob.rds"))

  res <- run_script("seurat_to_h5ad.R", c(
    "--seurat", file.path(dir, "obj.rds"), "--prob", file.path(dir, "prob.rds"),
    "--dataset", "toy", "--out-dir", dir, "--split-label", "sex",
    "--split-order", "Male,Female"
  ))
  expect_true(res$ok, info = paste(res$log, collapse = "\n"))

  clusters <- file.path(dir, "toy.clusters.h5ad")
  cells <- file.path(dir, "toy.cells.h5ad")
  val <- run_script("validate_h5ad.R", c("--clusters", clusters, "--cells", cells))
  expect_true(val$ok, info = paste(val$log, collapse = "\n"))

  ds <- ds_load("toy_seurat", list(label = "toy", clusters = clusters, cells = cells))
  expect_equal(ds$splits, c("Male", "Female"))
  expect_equal(ds$analysis_genes, c("geneA", "geneB", "geneC"))
  expect_equal(sort(rownames(ds$clusters_h$obs)),
               sort(c("Male|x", "Male|y", "Female|x")))
  b <- ds_binarize(ds, 0.5)
  expect_equal(from_cluster(b, "x", "Male", 1), c(geneA = 0L))
  expect_false(b$present["y", "Female"])
  expect_false(any(b$on["geneD", , ])) # unimodal -> never ON
  expect_equal(ds$reductions, "X_umap_full")
  # pct_pos matches the raw counts
  pct <- ds_cluster_layer_long(ds, "pct_pos", "geneA", clusters = "x")
  expected <- mean(counts["geneA", ident == "Male@x"] > 0)
  expect_equal(pct$value[pct$split == "Male"], expected)
})
