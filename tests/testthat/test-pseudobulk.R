test_that("random pseudobulk deals each cell into one sample", {
  set.seed(1)
  x <- matrix(rpois(30 * 50, 2), 30, 50, dimnames = list(paste0("g", 1:30), paste0("c", 1:50)))
  cluster <- rep(c("1", "2", "3"), c(24, 20, 6))
  set.seed(9)
  before <- .Random.seed
  pb <- make_pseudobulk(Matrix::Matrix(x, sparse = TRUE), cluster, n_samples = 3, min_cells = 7)
  expect_identical(.Random.seed, before)
  expect_equal(pb$cluster, c("1", "1", "1", "2", "2", "3"))
  expect_equal(pb$replicate, c("1", "2", "3", "1", "2", "1"))
  expect_equal(unname(pb$n_cells), c(8, 8, 8, 10, 10, 6))
  expect_equal(sum(pb$n_cells), ncol(x))
  # as AggregateExpression(return.seurat = TRUE): summed counts, LogNormalized
  summed <- rowSums(x[, pb$sample == "2_1"])
  expect_equal(unname(pb$lognorm[, "2_1"]), unname(log1p(summed / sum(summed) * 1e4)))
  expect_identical(make_pseudobulk(x, cluster, n_samples = 3, min_cells = 7)$sample, pb$sample)
  expect_error(make_pseudobulk(log1p(x / 3), cluster), "raw counts")
})

test_that("metadata mode splits clusters by a cell column", {
  set.seed(4)
  x <- matrix(rpois(10 * 40, 3), 10, 40, dimnames = list(paste0("g", 1:10), NULL))
  cluster <- rep(c("A", "B"), each = 20)
  group <- rep(c("L1", "L2"), 20)
  group[1:5] <- "L3"          # 5 cells of A: too few, left out
  pb <- suppressMessages(make_pseudobulk(x, cluster, method = "metadata", group = group,
                                         min_cells = 6))
  expect_equal(pb$cluster, c("A", "A", "B", "B"))
  expect_equal(pb$replicate, c("L1", "L2", "L1", "L2"))
  expect_equal(unname(pb$n_cells), c(7, 8, 10, 10))
  expect_true(all(is.na(pb$sample[1:5])))
  summed <- rowSums(x[, which(pb$sample == "B_L2")])
  expect_equal(unname(pb$lognorm[, "B_L2"]), unname(log1p(summed / sum(summed) * 1e4)))
  expect_message(make_pseudobulk(x, cluster, method = "metadata", group = group,
                                 min_cells = 6), "left out")
  # by default every group is a sample, and small ones are named
  expect_message(pb <- make_pseudobulk(x, cluster, method = "metadata", group = group),
                 "3 sample\\(s\\) have fewer than 10 cells.*A_L1 \\(7\\), A_L2 \\(8\\), A_L3 \\(5\\)")
  expect_equal(pb$replicate, c("L1", "L2", "L3", "L1", "L2"))
  expect_false(anyNA(pb$sample))
  # random parts keep a minimum of 10 cells by default
  pb <- make_pseudobulk(x, rep(c("A", "B"), c(15, 25)), n_samples = 2)
  expect_equal(unname(pb$n_cells), c(15, 13, 12))
  expect_error(make_pseudobulk(x, cluster, method = "metadata"), "needs group")
  expect_error(make_pseudobulk(x, cluster, group = group), "metadata")
})

test_that("partition_folds leaves a label with one sample alone", {
  set.seed(3)
  labels <- c(rep(1:2, each = 10), 3)
  folds <- partition_folds(21, 5, labels)
  expect_setequal(unique(folds), 1:5)
})
