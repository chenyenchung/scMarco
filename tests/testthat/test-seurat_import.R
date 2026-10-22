root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))
source(file.path(root, "script", "lib_import.R"))

run_script <- function(script, args) {
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(file.path(root, "script", script), args),
    stdout = TRUE, stderr = TRUE,
    env = "RENV_CONFIG_NAMESPACES_CHECK=FALSE"
  ))
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0, log = out)
}
log_of <- function(res) paste(res$log, collapse = "\n")

GENES <- c("geneA", "geneB", "geneC", "geneD")
QC <- c("diff.elpd", "diff.elpd.se", "bimodal.on_level", "bimodal.off_level",
        "bimodal.sd", "bimodal.pi", "unimodal.level", "unimodal.sd")

prob_table <- function(cols, expr_type, genes = GENES) {
  #' One split's mixture-model output: cluster columns + QC columns
  df <- as.data.frame(cols, check.names = FALSE, row.names = genes)
  df$exprType <- expr_type
  for (qc in QC) df[[qc]] <- 0
  df
}

# Male has clusters 1, 2 and 10; Female only 1. geneB is bimodal in Male only.
male_prob <- function() prob_table(
  list(`1` = c(0.9, 0.1, 0.2, 0.8), `2` = c(0.1, 0.9, 0.2, 0.8),
       `10` = c(0.1, 0.1, 0.9, 0.8)),
  c("bimodal", "bimodal", "bimodal", "unimodal")
)
female_prob <- function() prob_table(
  list(`1` = c(0.9, 0.9, 0.2, 0.8)),
  c("bimodal", "unimodal", "bimodal", "unimodal")
)

toy_seurat <- function(dir, legacy_idents = FALSE) {
  suppressPackageStartupMessages(library(SeuratObject))
  set.seed(1)
  n <- 40
  counts <- matrix(rpois(4 * n, 2), 4, n, dimnames = list(GENES, paste0("bc", 1:n)))
  counts["geneA", 1:10] <- 20 # high in Male 1
  obj <- CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  split <- rep(c("Male", "Female"), each = 20)
  cluster <- c(rep("1", 10), rep("2", 5), rep("10", 5), rep("1", 20))
  if (legacy_idents) {
    Idents(obj) <- factor(paste(split, cluster, sep = "@"))
  } else {
    obj$scmarco_split <- split
    obj$scmarco_cluster <- cluster
  }
  emb <- matrix(rnorm(2 * n), n, 2, dimnames = list(Cells(obj), c("umap_1", "umap_2")))
  obj[["umap"]] <- CreateDimReducObject(embeddings = emb, key = "umap_", assay = "RNA")
  path <- file.path(dir, "obj.rds")
  saveRDS(obj, path)
  list(path = path, counts = counts, split = split, cluster = cluster)
}

expect_toy_binarized <- function(clusters, cells) {
  ds <- ds_load(paste0("toy", runif(1)), list(label = "toy", clusters = clusters, cells = cells))
  expect_equal(ds$splits, c("Male", "Female"))
  expect_equal(ds$clusters, c("1", "2", "10"))
  expect_equal(ds$analysis_genes, c("geneA", "geneB", "geneC"))
  expect_equal(unname(ds$bimodal["geneB", c("Male", "Female")]), c(TRUE, FALSE))
  b <- ds_binarize(ds, 0.5)
  expect_equal(from_cluster(b, "1", "Male", 1), c(geneA = 0L))
  expect_false(b$present["2", "Female"])
  expect_false(any(b$on["geneD", , ])) # unimodal everywhere -> never ON
  expect_false(b$on["geneB", "1", "Female"]) # 0.9, but not bimodal in Female
  ds
}

test_that("build from metadata columns, then binarize from per-split files", {
  skip_if_not_installed("SeuratObject")
  dir <- withr::local_tempdir()
  toy <- toy_seurat(dir)

  res <- run_script("seurat_to_h5ad.R", c(
    "--seurat", toy$path, "--dataset", "toy", "--out-dir", dir,
    "--split-label", "sex", "--split-order", "Male,Female"
  ))
  expect_true(res$ok, info = log_of(res))
  clusters <- file.path(dir, "toy.clusters.h5ad")
  cells <- file.path(dir, "toy.cells.h5ad")

  # built: valid as such, refused as binarized and by the app
  val <- run_script("validate_h5ad.R", c("--clusters", clusters, "--cells", cells,
                                         "--stage", "built"))
  expect_true(val$ok, info = log_of(val))
  val <- run_script("validate_h5ad.R", c("--clusters", clusters))
  expect_false(val$ok)
  expect_match(log_of(val), "add_probability_h5ad.R")
  expect_error(
    ds_check_config(list(toy = list(clusters = clusters),
                         gone = list(clusters = file.path(dir, "nope.h5ad")))),
    "toy: .*no probabilities.*\n.*\n- gone: clusters file not found"
  )

  # cluster names set before binarizing survive it
  names_tsv <- file.path(dir, "names.tsv")
  writeLines(c("cluster_id\tname", "1\tMi4"), names_tsv)
  res <- run_script("set_cluster_names.R", c("--clusters", clusters, "--names", names_tsv))
  expect_true(res$ok, info = log_of(res))

  saveRDS(male_prob(), file.path(dir, "male.rds"))
  saveRDS(female_prob(), file.path(dir, "female.rds"))
  res <- run_script("add_probability_h5ad.R", c(
    "--clusters", clusters,
    "--prob", sprintf("Male=%s,Female=%s", file.path(dir, "male.rds"),
                      file.path(dir, "female.rds"))
  ))
  expect_true(res$ok, info = log_of(res))
  expect_match(log_of(res), "Male 3, Female 2")

  val <- run_script("validate_h5ad.R", c("--clusters", clusters, "--cells", cells))
  expect_true(val$ok, info = log_of(val))
  expect_silent(ds_check_config(list(toy = list(clusters = clusters))))

  ds <- expect_toy_binarized(clusters, cells)
  expect_equal(ds$cluster_labels[["1"]], "Mi4 (1)")
  meta <- ds$clusters_h$uns$scmarco
  expect_equal(meta$binarization$method, "mixture model")
  expect_equal(meta$split_label, "sex")
  expect_equal(ds$reductions, "X_umap_full")
  pct <- ds_cluster_layer_long(ds, "pct_pos", "geneA", clusters = "1")
  in_cl <- toy$split == "Male" & toy$cluster == "1"
  expect_equal(pct$value[pct$split == "Male"], mean(toy$counts["geneA", in_cl] > 0))

  # binarizing again (named-list RDS) replaces the probabilities
  saveRDS(list(Male = male_prob(), Female = female_prob()), file.path(dir, "both.rds"))
  res <- run_script("add_probability_h5ad.R", c(
    "--clusters", clusters, "--prob", file.path(dir, "both.rds")
  ))
  expect_true(res$ok, info = log_of(res))
  expect_match(log_of(res), "replacing them")
  expect_toy_binarized(clusters, cells)
})

test_that("legacy identities and wide table still build in one step", {
  skip_if_not_installed("SeuratObject")
  dir <- withr::local_tempdir()
  toy <- toy_seurat(dir, legacy_idents = TRUE)
  # one exprType for every split, so geneB is bimodal in Female too
  wide <- data.frame(
    `Male@1` = c(0.9, 0.1, 0.2, 0.8), `Male@2` = c(0.1, 0.9, 0.2, 0.8),
    `Male@10` = c(0.1, 0.1, 0.9, 0.8), `Female@1` = c(0.9, 0.1, 0.2, 0.8),
    check.names = FALSE, row.names = GENES
  )
  wide$exprType <- c("bimodal", "bimodal", "bimodal", "unimodal")
  saveRDS(wide, file.path(dir, "prob.rds"))

  res <- run_script("seurat_to_h5ad.R", c(
    "--seurat", toy$path, "--prob", file.path(dir, "prob.rds"),
    "--dataset", "toy", "--out-dir", dir, "--split-order", "Male,Female"
  ))
  expect_true(res$ok, info = log_of(res))
  expect_match(log_of(res), "legacy '<split>@<cluster>' identities")
  clusters <- file.path(dir, "toy.clusters.h5ad")
  val <- run_script("validate_h5ad.R", c("--clusters", clusters))
  expect_true(val$ok, info = log_of(val))

  ds <- ds_load("toy_legacy", list(label = "toy", clusters = clusters))
  expect_true(all(ds$bimodal["geneB", ]))
  b <- ds_binarize(ds, 0.5)
  expect_equal(from_cluster(b, "1", "Male", 1), c(geneA = 0L))
})

# --- Input rules (no Rscript needed) -------------------------------------------

test_that("seurat_labels picks metadata columns, then legacy identities", {
  meta <- data.frame(scmarco_split = c("M", "F"), scmarco_cluster = c(1, 2))
  expect_equal(seurat_labels(meta, c("x", "y")),
               list(split = c("M", "F"), cluster = c("1", "2")))
  expect_equal(seurat_labels(data.frame(row.names = 1:2), c("M@1", "F@2")),
               list(split = c("M", "F"), cluster = c("1", "2")))
  expect_error(seurat_labels(data.frame(row.names = 1:2), c("M@1", "F")),
               "neither metadata columns")
  # a custom column that is missing is an error, never a silent fallback
  expect_error(seurat_labels(meta, c("M@1", "F@2"), cluster_col = "clust",
                             defaults = FALSE), "'clust' not found")
  expect_error(seurat_labels(meta[, "scmarco_cluster", drop = FALSE], c("a", "b")),
               "--split-col none")
  one <- seurat_labels(meta, c("a", "b"), split_col = "none", defaults = FALSE)
  expect_equal(one$split, c("all", "all"))
})

test_that("check_labels enforces the name rules", {
  expect_true(check_labels(c("M", "F"), c("ANm_12", "T4 (a).1")))
  expect_error(check_labels(c("M", NA), c("1", "2")), "empty or missing split")
  expect_error(check_labels(c("M", "F"), c("1", "")), "empty or missing cluster")
  expect_error(check_labels(c("M", "F"), c("1|2", "2")), "must not contain '\\|'")
  expect_error(check_labels(c("M ", "F"), c("1", "2")), "whitespace")
})

test_that("assemble_prob places per-split values and flags", {
  obs_split <- c("Male", "Male", "Male", "Female")
  obs_cluster <- c("1", "2", "10", "1")
  p <- assemble_prob(list(Female = female_prob(), Male = male_prob()),
                     obs_split, obs_cluster, c("Male", "Female"), GENES)
  expect_equal(p$prob[, 1], c(0.9, 0.1, 0.1, 0.9))       # geneA
  expect_equal(p$prob[4, ], c(0.9, NA, 0.2, NA))          # Female 1
  expect_equal(p$bimodal_by_split["geneB", ], c(Male = TRUE, Female = FALSE))
})

test_that("assemble_prob enforces the rules in the spec", {
  obs_split <- c("Male", "Male", "Male", "Female")
  obs_cluster <- c("1", "2", "10", "1")
  run <- function(tables, genes = GENES) {
    assemble_prob(tables, obs_split, obs_cluster, c("Male", "Female"), genes)
  }
  f <- female_prob()
  m <- male_prob()

  expect_error(run(list(Male = m, Female = f, Other = f)), "not in the clusters file: Other")
  expect_error(run(list(Male = m)), "No probability table for split\\(s\\): Female")
  expect_error(run(list(Male = m, Male = m, Female = f)), "given twice")
  bad <- m; bad$exprType <- NULL
  expect_error(run(list(Male = bad, Female = f)), "Male': no exprType")
  bad <- f; colnames(bad)[1] <- "2" # cluster 2 has no Female cells
  expect_error(run(list(Male = m, Female = bad)), "no cells in this split: 2")
  bad <- cbind(m[, "1", drop = FALSE], m)
  expect_error(run(list(Male = bad, Female = f)), "repeated: 1")
  expect_error(run(list(Male = m, Female = f), genes = GENES[1:3]),
               "not in the clusters file, e.g. geneD")
  bad <- f; bad[1, "1"] <- 1.5
  expect_error(run(list(Male = m, Female = bad)), "outside \\[0, 1\\]")
  bad <- f; bad[["1"]] <- as.character(bad[["1"]])
  expect_error(run(list(Male = m, Female = bad)), "must be numeric")
  bad <- f; rownames(bad) <- NULL
  expect_error(run(list(Male = m, Female = bad)), "row names")

  # the mixture model's raw character matrix is parsed; non-numbers still fail
  txt <- as.matrix(f)
  p <- run(list(Male = m, Female = txt))
  expect_equal(p$prob, run(list(Male = m, Female = f))$prob)
  txt[1, "1"] <- "high"
  expect_error(run(list(Male = m, Female = txt)), "must be numeric")

  # gaps are warnings: never ON / non-bimodal
  expect_warning(p <- run(list(Male = m[, colnames(m) != "10"], Female = f)),
                 "1 cluster\\(s\\) with cells have no column.*: 10")
  expect_true(is.na(p$prob[3, 1]))
  expect_warning(p <- run(list(Male = m, Female = f[1:3, ])), "1 gene\\(s\\) are not in")
  expect_false(p$bimodal_by_split["geneD", "Female"])

  # legacy hyphens match underscores
  p <- assemble_prob(list(S = prob_table(list(`a-1` = c(1, 0, 0, 0)), "bimodal")),
                     "S", "a_1", "S", GENES)
  expect_equal(p$prob[1, 1], 1)
})

test_that("legacy wide tables split into per-split tables", {
  wide <- data.frame(`M@1` = 1:2 / 4, `F@1` = 3:4 / 4, exprType = "bimodal",
                     check.names = FALSE, row.names = c("g1", "g2"))
  t <- legacy_wide_to_split_tables(wide)
  expect_equal(names(t), c("M", "F"))
  expect_equal(colnames(t$F), c("1", "exprType"))
  expect_equal(t$F[["1"]], c(0.75, 1))
  wide$oops <- 1
  expect_error(legacy_wide_to_split_tables(wide), "'<split>@<cluster>'.*oops")
})

test_that("build from one object per split with different genes and assays", {
  skip_if_not_installed("SeuratObject")
  dir <- withr::local_tempdir()
  toy <- toy_seurat(dir)
  obj <- readRDS(toy$path)
  # Female lacks geneD; both objects default to another assay than the counts
  for (s in c("Male", "Female")) {
    keep <- if (s == "Female") GENES[1:3] else GENES
    part <- CreateSeuratObject(
      Matrix::Matrix(toy$counts[keep, toy$split == s], sparse = TRUE)
    )
    part$cl <- toy$cluster[toy$split == s]
    part[["other"]] <- CreateAssayObject(counts = part[["RNA"]]$counts[1:2, ])
    DefaultAssay(part) <- "other"
    saveRDS(part, file.path(dir, paste0(s, ".rds")))
  }
  build_args <- c(
    "--seurat", sprintf("Male=%s,Female=%s", file.path(dir, "Male.rds"),
                        file.path(dir, "Female.rds")),
    "--cluster-col", "cl", "--assay", "RNA", "--dataset", "toy", "--out-dir", dir
  )
  # A missing reduction is an error, not a silently missing UMAP
  res <- run_script("seurat_to_h5ad.R", build_args)
  expect_false(res$ok)
  expect_match(log_of(res), "Male: no reduction named 'umap'. Reductions on the object: none", fixed = TRUE)
  res <- run_script("seurat_to_h5ad.R", c(build_args, "--reduction", "none"))
  expect_true(res$ok, info = log_of(res))
  clusters <- file.path(dir, "toy.clusters.h5ad")
  cells <- file.path(dir, "toy.cells.h5ad")
  val <- run_script("validate_h5ad.R", c("--clusters", clusters, "--cells", cells,
                                         "--stage", "built"))
  expect_true(val$ok, info = log_of(val))

  saveRDS(male_prob(), file.path(dir, "male.rds"))
  saveRDS(female_prob()[1:3, ], file.path(dir, "female.rds"))
  res <- run_script("add_probability_h5ad.R", c(
    "--clusters", clusters,
    "--prob", sprintf("Male=%s,Female=%s", file.path(dir, "male.rds"),
                      file.path(dir, "female.rds"))
  ))
  expect_true(res$ok, info = log_of(res))
  # splits keep the --seurat order
  ds <- expect_toy_binarized(clusters, cells)
  expect_equal(ds$all_genes, GENES)

  pct <- ds_cluster_layer_long(ds, "pct_pos", "geneA", clusters = "1")
  in_cl <- toy$split == "Male" & toy$cluster == "1"
  expect_equal(pct$value[pct$split == "Male"], mean(toy$counts["geneA", in_cl] > 0))
  pct <- ds_cluster_layer_long(ds, "pct_pos", "geneD", clusters = "1")
  expect_equal(pct$value[pct$split == "Female"], 0)
  # library size counts only the Female object's genes
  fem <- toy$counts[1:3, toy$split == "Female"]
  mean_a <- ds_cluster_layer_long(ds, "mean_lognorm", "geneA", clusters = "1")
  expect_equal(mean_a$value[mean_a$split == "Female"],
               mean(log1p(fem["geneA", ] / colSums(fem) * 1e4)), tolerance = 1e-6)

  expect_error(parse_seurat_arg("A=x.rds,A=y.rds"), "given twice")
  res <- run_script("seurat_to_h5ad.R", c(
    "--seurat", sprintf("Male=%s", file.path(dir, "Male.rds")),
    "--split-col", "x", "--dataset", "toy2", "--out-dir", dir
  ))
  expect_false(res$ok)
  expect_match(log_of(res), "--split-col cannot be used")
})
