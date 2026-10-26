root <- file.path(testthat::test_path(), "..", "..")
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))
source(file.path(root, "src", "plot_core.R"))
source(file.path(root, "src", "ui_helpers.R"))

fixture <- file.path(root, "tests", "fixtures", "example.clusters.h5ad")

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")

test_that("cluster averages are lines across ordered splits, points otherwise", {
  ds <- ds_load("ordered", list(label = "Example", clusters = fixture))
  gene <- ds$all_genes[1]
  p <- trend_plot(ds, gene, interactive = FALSE)
  expect_true("GeomLine" %in% layer_geoms(p))
  expect_equal(p$labels$x, "Stage")

  ds <- ds_load("unordered", list(label = "Example", clusters = fixture,
                                  ordered_splits = FALSE))
  p <- trend_plot(ds, gene, interactive = FALSE)
  expect_true("GeomPoint" %in% layer_geoms(p))
  expect_false("GeomLine" %in% layer_geoms(p))
})

test_that("split colors follow whether splits are ordered", {
  scale_of <- function(splits, ordered) {
    split_color_scale(list(splits = splits, ordered_splits = ordered))
  }
  few <- c("Mel", "Yak", "San", "Tei")
  expect_true(inherits(scale_of(few, FALSE), "ScaleDiscrete"))
  # Set1 has 9 colors: more splits must not fall back to NA
  many <- paste0("c", 1:12)
  pal <- scale_of(many, FALSE)$palette(12)
  expect_false(anyNA(pal))
  expect_length(unique(pal), 12)
  expect_false(anyNA(scale_of(many, TRUE)$palette(12)))
})

test_that("the criteria text uses the singular or plural label", {
  expect_equal(describe_splits(character(0), 1, "sex", "sexes"),
               "Criteria: ON in at least 1 sex.")
  expect_equal(describe_splits("Mel", 2, "species", "species"),
               "Criteria: ON in Mel and ON in at least 2 species.")
})
