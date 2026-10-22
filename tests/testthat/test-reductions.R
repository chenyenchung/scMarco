root <- file.path(testthat::test_path(), "..", "..")
source(file.path(root, "src", "plot_core.R"))

test_that("per-split reductions are grouped by reduction name", {
  keys <- c("X_umap_full", "X_ANm", "X_P15_tsne", "X_Adult_tsne", "X_P15", "X_A_B_tsne")
  groups <- split_reductions(keys, c("P15", "Adult", "A", "A_B"))
  expect_equal(groups$tsne, c(P15 = "X_P15_tsne", Adult = "X_Adult_tsne", A_B = "X_A_B_tsne"))
  expect_equal(groups$umap, c(P15 = "X_P15"))
})

test_that("reduction choices fold per-split keys into one entry", {
  ds <- list(reductions = c("X_umap_full", "X_P15_tsne", "X_P30_tsne"),
             splits = c("P15", "P30"), split_label = "stage")
  ch <- reduction_choices(ds)
  expect_equal(unname(ch), c("X_umap_full", "split:tsne"))
  expect_equal(names(ch), c("All cells", "Each stage separately"))
  expect_equal(resolve_reduction(ds, "split:tsne"), c(P15 = "X_P15_tsne", P30 = "X_P30_tsne"))
  expect_equal(resolve_reduction(ds, "X_umap_full"), "X_umap_full")
  # Subset reductions that are not splits stay shared
  ds$reductions <- c("X_umap_full", "X_ANm")
  expect_equal(unname(reduction_choices(ds)), c("X_umap_full", "X_ANm"))
})
