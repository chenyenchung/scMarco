source(file.path(testthat::test_path(), "..", "..", "src", "combo_core.R"))

# Toy dataset: 5 genes, clusters A/B/C, splits s1/s2/s3; C is absent from s3.
# Values are probabilities; rows are split|cluster observations.
toy <- function(cut_off = 0.5, bimodal = NULL) {
  obs <- data.frame(
    split = c("s1", "s1", "s1", "s2", "s2", "s2", "s3", "s3"),
    cluster = c("A", "B", "C", "A", "B", "C", "A", "B")
  )
  #            g1   g2   g3   g4   g5
  prob <- rbind(
    c(0.9, 0.9, 0.1, 0.9, 0.5), # s1 A
    c(0.1, 0.9, 0.1, 0.1, 0.1), # s1 B
    c(0.1, 0.1, 0.9, 0.1, 0.1), # s1 C
    c(0.9, 0.9, 0.1, 0.1, 0.1), # s2 A
    c(0.1, 0.1, 0.1, 0.1, 0.1), # s2 B
    c(0.1, 0.1, 0.9, 0.9, 0.1), # s2 C
    c(0.9, 0.1, 0.1, 0.1, 0.1), # s3 A
    c(0.1, 0.1, 0.1, 0.1, 0.1)  # s3 B
  )
  colnames(prob) <- paste0("g", 1:5)
  binarize(prob, obs$split, obs$cluster,
           splits = c("s1", "s2", "s3"), clusters = c("A", "B", "C"),
           cut_off = cut_off, bimodal = bimodal)
}

test_that("binarize marks absent clusters and uses a strict threshold", {
  b <- toy()
  expect_false(b$present["C", "s3"])
  expect_true(all(b$present[c("A", "B"), ]))
  # g5 in s1 A is exactly 0.5: strict '>' means OFF (legacy OnPlot used >=)
  expect_false(b$on["g5", "A", "s1"])
  expect_true(b$on["g1", "A", "s3"])
})

test_that("target_genes applies required splits AND min_n", {
  b <- toy()
  expect_setequal(target_genes(b, "A", character(0), 1), c("g1", "g2", "g4"))
  expect_setequal(target_genes(b, "A", character(0), 2), c("g1", "g2"))
  expect_equal(target_genes(b, "A", character(0), 3), "g1")
  expect_setequal(target_genes(b, "A", c("s1", "s2"), 1), c("g1", "g2"))
  # One required split and min_n = 1 must not demand an extra split
  # (legacy scMarco From-gene HAVING bug)
  expect_setequal(target_genes(b, "A", "s1", 1), c("g1", "g2", "g4"))
  expect_error(target_genes(b, c("A", "B"), character(0), 1))
  expect_error(target_genes(b, "A", "nope", 1))
})

test_that("absent clusters cannot pass a split they are missing from", {
  b <- toy()
  expect_length(target_genes(b, "C", "s3", 1), 0)
  expect_setequal(target_genes(b, "C", character(0), 2), "g3")
})

test_that("from_cluster ranks by off-targets and never pads with NA", {
  b <- toy()
  res <- from_cluster(b, "A", character(0), 1, limit = 10)
  expect_equal(res, c(g1 = 0L, g2 = 1L, g4 = 1L))
  expect_false(anyNA(res))
  expect_length(from_cluster(b, "A", character(0), 1, limit = 2), 2)
  expect_length(from_cluster(b, "B", c("s1", "s2"), 1), 0)
})

test_that("clusters_expressing supports all/any across several genes", {
  b <- toy()
  expect_setequal(clusters_expressing(b, "g2", character(0), 1), c("A", "B"))
  expect_equal(clusters_expressing(b, c("g1", "g2"), character(0), 1), "A")
  expect_setequal(
    clusters_expressing(b, c("g1", "g3"), character(0), 1, mode = "any"),
    c("A", "C")
  )
  expect_length(clusters_expressing(b, c("g1", "g3"), character(0), 1), 0)
  expect_length(clusters_expressing(b, "unknown_gene", character(0), 1), 0)
})

test_that("distinct_within returns genes unique to one selected cluster", {
  b <- toy()
  res <- distinct_within(b, c("A", "B"), character(0), 1)
  expect_equal(res$gene, c("g1", "g4"))
  expect_equal(res$cluster, c("A", "A"))
  # ranked by clusters outside the selection that also express the gene
  expect_equal(res$offtargets, c(0L, 1L))
  expect_equal(top_per_cluster(res, 1)$gene, "g1")
  # g2 is ON in both A and B, so it is not distinct
  expect_false("g2" %in% res$gene)
  # nothing found -> zero rows, not a fake "Not found" row
  res_none <- distinct_within(b, c("A", "B"), c("s1", "s2", "s3"), 3,
                              genes = "g2")
  expect_equal(nrow(res_none), 0)
  expect_error(distinct_within(b, "A", character(0), 1))
})

test_that("pair_offtargets counts clusters where both genes are ON", {
  b <- toy()
  res <- pair_offtargets(b, "A", character(0), 1)
  key <- paste(res$gene1, res$gene2)
  expect_equal(res$offtargets[key == "g1 g2"], 0L)
  expect_equal(res$offtargets[key == "g1 g4"], 0L)
  expect_equal(res$offtargets[key == "g2 g4"], 0L)
  res_first <- pair_offtargets(b, "A", character(0), 1, first = "g2")
  expect_true(all(res_first$gene1 == "g2"))
  expect_equal(nrow(res_first), 2)
  expect_equal(nrow(pair_offtargets(b, "A", character(0), 1, first = "g3")), 0)
})

test_that("combo_offtarget_clusters checks a combination dataset-wide", {
  b <- toy()
  expect_equal(combo_offtarget_clusters(b, "g2", exclude = "A"), "B")
  expect_length(combo_offtarget_clusters(b, c("g1", "g2"), exclude = "A"), 0)
})

test_that("non-bimodal genes are never ON and are reported as such", {
  bim <- matrix(TRUE, 5, 3, dimnames = list(paste0("g", 1:5), paste0("s", 1:3)))
  bim["g1", "s3"] <- FALSE
  b <- toy(bimodal = bim)
  expect_false(b$on["g1", "A", "s3"])
  expect_equal(target_genes(b, "A", character(0), 3), character(0))
  st <- state_table(b, "g1")
  expect_equal(st$state[st$cluster == "A" & st$split == "s3"], "Non-bimodal")
  expect_equal(st$state[st$cluster == "C" & st$split == "s3"], "Absent")
  expect_equal(st$state[st$cluster == "A" & st$split == "s1"], "Expressed")
})
