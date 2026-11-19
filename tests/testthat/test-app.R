toy_dir <- system.file("extdata", package = "scMarco")

toy_config <- function() {
  list(
    title = "Toy",
    datasets = list(
      toy = list(label = "Toy", clusters = "toy.clusters.h5ad",
                 cells = "toy.cells.h5ad", subset_by = list("class"))
    )
  )
}

test_that("relative config paths resolve against root", {
  cfg <- read_scmarco_config(toy_config(), root = toy_dir)
  expect_equal(cfg$datasets$toy$clusters, file.path(toy_dir, "toy.clusters.h5ad"))
  expect_true(file.exists(cfg$datasets$toy$cells))
  cfg <- toy_config()
  cfg$datasets$toy$clusters <- file.path(toy_dir, "toy.clusters.h5ad")
  expect_equal(read_scmarco_config(cfg, root = "elsewhere")$datasets$toy$clusters,
               file.path(toy_dir, "toy.clusters.h5ad"))
})

test_that("a YAML config is read", {
  path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(toy_config(), path)
  cfg <- read_scmarco_config(path, root = toy_dir)
  expect_equal(cfg$title, "Toy")
  expect_equal(names(cfg$datasets), "toy")
})

test_that("configs without datasets or with duplicate labels are refused", {
  expect_error(read_scmarco_config(list(title = "x")), "No datasets")
  cfg <- toy_config()
  cfg$datasets$again <- cfg$datasets$toy
  expect_error(read_scmarco_config(cfg, root = toy_dir), "unique label")
})

test_that("scmarco_app() builds an app", {
  app <- suppressMessages(scmarco_app(toy_config(), root = toy_dir))
  expect_s3_class(app, "shiny.appobj")
})

test_that("the server runs a From cluster search on the toy dataset", {
  config <- read_scmarco_config(toy_config(), root = toy_dir)
  # testServer() exposes the server function's own variables (ds, b, ...)
  server <- scmarco_server
  formals(server)$config <- config
  shiny::testServer(server, {
    session$setInputs(dataset = "toy", gene_group_to_use = "All", band_width = 80,
                      required_splits = character(0), min_n = 1)
    expect_equal(ds()$splits, c("Male", "Female"))
    expect_equal(length(universe()), 8)
    session$setInputs(meta_levels = "class", meta_1 = "glia")
    expect_equal(universe(), "8")
    expect_equal(b()$clusters, "8")
  })
})
