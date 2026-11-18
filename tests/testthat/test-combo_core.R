# Toy dataset: 5 genes, clusters A/B/C, splits s1/s2/s3; C is absent from s3.
# Values are probabilities; rows are split|cluster observations.
toy <- function(band = 0, bimodal = NULL) {
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
           band = band, bimodal = bimodal)
}

test_that("binarize marks absent clusters and uses a strict threshold", {
  b <- toy()
  expect_false(b$present["C", "s3"])
  expect_true(all(b$present[c("A", "B"), ]))
  # g5 in s1 A is exactly 0.5: strict '>' means OFF (legacy OnPlot used >=)
  expect_false(b$on["g5", "A", "s1"])
  expect_true(b$on["g1", "A", "s3"])
})

test_that("binarize makes the calls flagged by the binarizer 1.1 ambiguous", {
  # one split, six clusters; columns: ON, OFF and intermediate calls, and no call
  prob <- matrix(c(0.95, 0.95, 0.95, 0.02, 0.02, 0.02, 0.5, NA), ncol = 1,
                 dimnames = list(NULL, "g"))
  flags <- matrix(c(CALL_FLAG_CHAINS, CALL_FLAG_DETECTION, CALL_FLAG_AMBIENT,
                    CALL_FLAG_CHAINS, CALL_FLAG_DETECTION, CALL_FLAG_AMBIENT, 0L,
                    CALL_FLAG_AMBIENT), ncol = 1)
  cl <- paste0("c", 1:8)
  b <- binarize(prob, rep("s", 8), cl, "s", cl, band = 0.3, flags = flags)
  # ON: chains disagree or below the detection floor -> ambiguous; above ambient stays ON
  # OFF: chains disagree or above ambient -> ambiguous; below the floor stays OFF
  expect_identical(unname(b$on["g", , "s"]),
                   c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(unname(b$ambiguous["g", , "s"]),
                   c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_true(b$flags_applied)
  # the flag that decided is kept for display
  expect_identical(unname(b$flagged["g", , "s"]),
                   c(CALL_FLAG_CHAINS, CALL_FLAG_DETECTION, 0L, CALL_FLAG_CHAINS, 0L,
                     CALL_FLAG_AMBIENT, 0L, 0L))
  st <- state_table(b, "g")
  expect_identical(st$reason, c("chains disagree", "below detection floor", "", "chains disagree",
                                "", "above ambient", "", ""))
  p <- onoff_plot(list(splits = "s", cluster_labels = setNames(cl, cl), has_pct = FALSE),
                  b, "g", cl)
  expect_match(p$labels$caption, "c: the fit's chains disagree.*d: too few cells.*a: more UMIs")
  # without flags only the band decides
  b0 <- binarize(prob, rep("s", 8), cl, "s", cl, band = 0.3)
  expect_identical(unname(b0$on["g", , "s"]), rep(c(TRUE, FALSE), c(3, 5)))
  expect_false(b0$flags_applied)
  expect_null(b0$flagged)
  expect_true(all(state_table(b0, "g")$reason == ""))
  expect_error(binarize(prob, rep("s", 8), cl, "s", cl, band = 0.3, flags = flags[-1, , drop = FALSE]))
})

test_that("pair_offtargets breaks ties by separation, then Ashman's D", {
  # g1-g3 are ON in the target T and OFF in A and B
  cl <- c("T", "A", "B")
  prob <- matrix(rep(c(0.95, 0.02, 0.02), 3), 3, dimnames = list(NULL, paste0("g", 1:3)))
  # log2 fold over ambient
  fold <- rbind(c(5, 4, 6),     # T
                c(2, 0.5, 1),   # A
                c(0, 3, 3))     # B
  b <- binarize(prob, rep("s", 3), cl, "s", cl, band = 0.3, fold = fold)
  res <- pair_offtargets(b, "T", character(0), 1)
  # g1+g3: target 5, worst leak max(min(2, 1), min(0, 3)) = 1 -> 4
  # g1+g2: target 4, leak max(0.5, 0) = 0.5 -> 3.5; g2+g3: 4 - 3 = 1
  expect_identical(paste(res$gene1, res$gene2), c("g1 g3", "g1 g2", "g2 g3"))
  expect_equal(res$separation, c(4, 3.5, 1))
  expect_equal(res$target_fold, c(5, 4, 4))
  expect_equal(res$worst_leak, c(1, 0.5, 3))
  expect_true(all(res$offtargets == 0))
  # an off-target cluster (both genes ON) is not a leak
  prob2 <- prob
  prob2[3, 1:2] <- 0.95
  res2 <- pair_offtargets(binarize(prob2, rep("s", 3), cl, "s", cl, band = 0.3, fold = fold),
                          "T", character(0), 1)
  expect_equal(res2$worst_leak[res2$gene1 == "g1" & res2$gene2 == "g2"], 0.5)
  expect_equal(res2$offtargets[res2$gene1 == "g1" & res2$gene2 == "g2"], 1)
  # without folds: the weaker gene's Ashman's D lower bound, then position
  D <- matrix(c(3, 1, 2), 3, dimnames = list(paste0("g", 1:3), "s"))
  b <- binarize(prob, rep("s", 3), cl, "s", cl, band = 0.3, ashman_D_lower = D)
  res <- pair_offtargets(b, "T", character(0), 1)
  expect_identical(paste(res$gene1, res$gene2), c("g1 g3", "g1 g2", "g2 g3"))
  expect_equal(res$ashman_D_lower, c(2, 1, 1))
  expect_true(all(is.na(res$separation)))
  # without either: as before (single-gene off-targets, then position)
  res <- pair_offtargets(binarize(prob, rep("s", 3), cl, "s", cl, band = 0.3),
                         "T", character(0), 1)
  expect_identical(paste(res$gene1, res$gene2), c("g1 g2", "g1 g3", "g2 g3"))
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

test_that("top_per_cluster_capped trims but keeps every cluster", {
  tbl <- data.frame(
    cluster = rep(c("A", "B", "C"), each = 4),
    gene = paste0("g", 1:12),
    offtargets = 0L
  )
  res <- top_per_cluster_capped(tbl, 4, 12)
  expect_equal(res$n, 4L)
  expect_equal(nrow(res$shown), 12)
  res <- top_per_cluster_capped(tbl, 4, 7)
  expect_equal(res$n, 2L)
  expect_equal(res$shown$gene, c("g1", "g2", "g5", "g6", "g9", "g10"))
  # one per cluster even when that exceeds the cap
  res <- top_per_cluster_capped(tbl, 4, 1)
  expect_equal(res$n, 1L)
  expect_setequal(res$shown$cluster, c("A", "B", "C"))
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

test_that("a non-bimodal gene makes co-expression undetermined, like NA", {
  bim <- matrix(TRUE, 5, 3, dimnames = list(paste0("g", 1:5), paste0("s", 1:3)))
  bim["g2", "s3"] <- FALSE
  b <- toy(bimodal = bim)
  co <- coexp_data(NULL, b, c("g1", "g2"))
  co <- co[co$gene == "Co-expressed", ]
  state <- function(cl, s) co$state[co$cluster == cl & co$split == s]
  expect_equal(state("A", "s1"), "Expressed")
  expect_equal(state("B", "s1"), "Not expressed")
  # gray wins over both ON (g1 in A) and OFF (g1 in B)
  expect_equal(state("A", "s3"), "Non-bimodal")
  expect_equal(state("B", "s3"), "Non-bimodal")
})

test_that("calls inside the ambiguity band are ambiguous and never help", {
  b <- toy(band = 0.3)
  # 0.9 is ON, 0.1 is OFF, 0.5 (g5 in s1 A) is ambiguous
  expect_true(b$on["g1", "A", "s1"])
  expect_false(b$on["g5", "A", "s1"])
  expect_true(b$ambiguous["g5", "A", "s1"])
  expect_false(any(b$ambiguous[c("g1", "g2", "g3", "g4"), , ]))
  # ambiguous is not ON in the target...
  expect_false("g5" %in% target_genes(b, "A", "s1", 1))
  # ...but counts as possibly ON off-target
  expect_equal(offtarget_counts(b, "g5", exclude = "B"), c(g5 = 1L))
  expect_equal(state_table(b, "g5", clusters = "A")$state,
               c("Ambiguous", "Not expressed", "Not expressed"))
  # the band edges are not ambiguous: band 0.4 makes 0.9 ON and 0.1 OFF
  b <- toy(band = 0.4)
  expect_true(b$on["g1", "A", "s1"])
  expect_false(b$ambiguous["g1", "B", "s1"])
})

test_that("an ambiguous gene makes co-expression undetermined", {
  b <- toy(band = 0.3)
  co <- coexp_data(NULL, b, c("g1", "g5"))
  co <- co[co$gene == "Co-expressed", ]
  state <- function(cl, s) co$state[co$cluster == cl & co$split == s]
  expect_equal(state("A", "s1"), "Ambiguous")
  expect_equal(state("A", "s2"), "Not expressed")
})

test_that("restrict_clusters hides the clusters left out from every search", {
  b <- restrict_clusters(toy(), c("C", "A"))
  expect_equal(b$clusters, c("A", "C"))
  expect_equal(dimnames(b$on)[[2]], c("A", "C"))
  expect_equal(rownames(b$present), c("A", "C"))
  # g2 is also ON in B, which is out of scope now
  expect_equal(from_cluster(toy(), "A", character(0), 1)[["g2"]], 1L)
  expect_equal(from_cluster(b, "A", character(0), 1)[["g2"]], 0L)
  expect_false("B" %in% clusters_expressing(b, "g2", character(0), 1))
  expect_error(target_genes(b, "B", character(0), 1))
  expect_identical(restrict_clusters(toy(), c("A", "B", "C")), toy())
})

test_that("any_on reads the matrix binarize() stores, also after restricting", {
  computed <- function(b, genes = b$genes, clusters = b$clusters) {
    b$anyon <- NULL
    any_on(b, genes = genes, clusters = clusters)
  }
  bimodal <- matrix(TRUE, 5, 3, dimnames = list(paste0("g", 1:5), c("s1", "s2", "s3")))
  bimodal["g5", "s1"] <- FALSE
  for (b in list(toy(band = 0.4), toy(band = 0.4, bimodal = bimodal))) {
    expect_identical(b$anyon, computed(b))
    expect_identical(any_on(b, "g5", c("C", "A")), computed(b, "g5", c("C", "A")))
    sub <- restrict_clusters(b, c("A", "C"))
    expect_identical(sub$anyon, computed(sub))
  }
  # g5 is ambiguous (0.5 at band 0.4) only in s1 A
  expect_identical(toy(band = 0.4)$anyon["g5", ], c(A = TRUE, B = FALSE, C = FALSE))
})

test_that("can_rank keeps every value tied with the limit-th smallest", {
  expect_identical(can_rank(c(3, 1, 2, 2, 5), 2), c(FALSE, TRUE, TRUE, TRUE, FALSE))
  expect_identical(can_rank(c(3, 1), 5), c(TRUE, TRUE))
  b <- toy()
  expect_identical(pair_offtargets(b, "A", character(0), 1, limit = 1),
    head(pair_offtargets(b, "A", character(0), 1, limit = 50), 1))
})

test_that("the deprecated margin gives the same calls as band = 0.5 - margin", {
  obs <- data.frame(split = c("s1", "s1"), cluster = c("A", "B"))
  prob <- matrix(c(0.95, 0.5, 0.05, 0.85), 2, dimnames = list(NULL, c("g1", "g2")))
  b1 <- binarize(prob, obs$split, obs$cluster, "s1", c("A", "B"), band = 0.4)
  expect_warning(b2 <- binarize(prob, obs$split, obs$cluster, "s1", c("A", "B"), margin = 0.1),
                 "`margin` is deprecated")
  expect_identical(b2, b1)
  expect_error(suppressWarnings(binarize(prob, obs$split, obs$cluster, "s1", c("A", "B"),
                                         band = 0.4, margin = 0.1)), "not both")
})
