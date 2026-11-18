toy_counts <- function() {
  set.seed(2)
  # 6 genes x 12 cells; clusters A (cells 1-6) and B (7-12); libraries L1/L2
  # alternate; g1 is expressed in A only, g6 is never detected
  m <- matrix(rpois(6 * 12, 3), 6, 12, dimnames = list(paste0("g", 1:6), paste0("c", 1:12)))
  m["g1", 7:12] <- 0
  m["g6", ] <- 0
  m
}

test_that("ambient_summary counts UMIs, detection and library totals per cluster", {
  m <- toy_counts()
  cl <- rep(c("A", "B"), each = 6)
  lib <- rep(c("L1", "L2"), 6)
  a <- ambient_summary(m, cl, lib, rho = 0.1)
  expect_s3_class(a, "scmarco_ambient")
  expect_equal(a$umi[, "A"], rowSums(m[, 1:6]))
  expect_equal(a$umi[, "B"], rowSums(m[, 7:12]))
  expect_equal(a$detect["g1", ], c(A = mean(m["g1", 1:6] > 0), B = 0))
  expect_equal(unname(a$detect["g6", ]), c(0, 0))
  expect_equal(a$n_cells, c(A = 6, B = 6))
  expect_equal(a$cluster_library_umi["A", "L1"], sum(m[, c(1, 3, 5)]))
  expect_equal(sum(a$cluster_library_umi), sum(m))
  expect_equal(a$cluster_library_cells, matrix(3, 2, 2, dimnames = list(c("A", "B"), c("L1", "L2"))))
  expect_equal(a$library_share[, "L2"], rowSums(m[, seq(2, 12, 2)]) / sum(m[, seq(2, 12, 2)]))
  expect_equal(a$rho, c(L1 = 0.1, L2 = 0.1))
  expect_identical(a$rho_source, "given")
  # expected ambient UMIs: sum over libraries of rho x share x cluster UMIs
  lam <- ambient_expected(a)
  expect_equal(lam["g2", "B"],
               sum(a$rho * a$library_share["g2", ] * a$cluster_library_umi["B", ]))
  expect_equal(unname(lam["g6", ]), c(0, 0))
  # genes keeps only some genes, with totals still from every gene
  b <- ambient_summary(m, cl, lib, rho = 0.1, genes = c("g2", "g3"))
  expect_identical(rownames(b$umi), c("g2", "g3"))
  expect_identical(b$cluster_library_umi, a$cluster_library_umi)
})

test_that("ambient_summary takes rho per library and checks its input", {
  m <- toy_counts()
  cl <- rep(c("A", "B"), each = 6)
  lib <- rep(c("L1", "L2"), 6)
  a <- ambient_summary(m, cl, lib, rho = c(L2 = 0.2, L1 = 0.05))
  expect_equal(a$rho, c(L1 = 0.05, L2 = 0.2))
  expect_error(ambient_summary(m, cl, lib, rho = c(L1 = 0.05)), "no value for library L2")
  expect_error(ambient_summary(m, cl, lib, rho = 1.5), "fractions in \\[0, 1\\)")
  expect_error(ambient_summary(m, cl[-1], lib, rho = 0.1), "cluster has 11 entries")
  expect_error(ambient_summary(m / 2, cl, lib, rho = 0.1), "raw counts")
  expect_error(ambient_summary(m, cl, lib, rho = 0.1, genes = "nope"), "not in counts")
  # one library when none is given
  expect_identical(ambient_summary(m, cl, rho = 0.1)$libraries, "1")
})

test_that("ambient_summary estimates each library's contamination with DecontX", {
  skip_if_not_installed("celda")
  # decontX computes a UMAP within each library, so it needs more cells than
  # its default 15 neighbours
  set.seed(3)
  m <- matrix(rpois(40 * 120, 1), 40, 120, dimnames = list(paste0("g", 1:40), paste0("c", 1:120)))
  m[1:10, 1:60] <- m[1:10, 1:60] + rpois(10 * 60, 8)        # genes ON in cluster A
  m[11:20, 61:120] <- m[11:20, 61:120] + rpois(10 * 60, 8)  # genes ON in cluster B
  a <- suppressWarnings(ambient_summary(m, rep(c("A", "B"), each = 60), rep(c("L1", "L2"), 60)))
  expect_named(a$rho, c("L1", "L2"))
  expect_true(all(a$rho >= 0 & a$rho < 1))
  expect_match(a$rho_source, "^decontX \\(celda ")
})

test_that("libraries with few cells per cluster are named in a warning", {
  # L1: 40 cells over 2 clusters; L2: 12 cells over 4 clusters
  cells <- matrix(c(20, 20, 0, 0, 3, 3, 3, 3), 4, 2,
                  dimnames = list(c("A", "B", "C", "D"), c("L1", "L2")))
  note <- small_library_note(cells)
  expect_match(note, "^1 of 2 libraries have fewer than 10 cells per cluster")
  expect_match(note, "L2: 12 cells over 4 clusters")
  expect_false(grepl("L1:", note))
  expect_match(note, "give the suspension as library")
  expect_null(small_library_note(cells[, "L1", drop = FALSE]))
})
