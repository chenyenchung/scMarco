test_that("find_pairs ranks pairs per cluster, the same on any number of cores", {
  ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
  serial <- find_pairs(ex, top = 3)
  expect_named(serial, c("cluster", "cluster_name", "gene1", "gene2", "offtargets",
    "single_offtargets"))
  expect_true(all(table(serial$cluster) <= 3))
  skip_on_os("windows")
  expect_identical(find_pairs(ex, top = 3, cores = 2), serial)
})

test_that("find_pairs writes the table it returns", {
  ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
  out <- tempfile(fileext = ".csv")
  expect_message(pairs <- find_pairs(ex, out = out, top = 2), "Wrote")
  back <- read.csv(out, colClasses = "character")
  expect_equal(nrow(back), nrow(pairs))
  expect_equal(back$gene1, pairs$gene1)
})
