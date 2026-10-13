# One-off parity check: new combination core vs. the two legacy forks' own
# implementations, on random (dataset, cluster(s), threshold, split selection)
# cases. Run from the project root after converting the data:
#
#   Rscript tests/parity/parity.R ../scMarco ../scMarco-Ding-lab
#
# Differences are expected only where a documented bug fix applies; any other
# mismatch is printed and makes the script exit non-zero.

args <- commandArgs(trailingOnly = TRUE)
legacy_sc <- normalizePath(if (length(args) >= 1) args[1] else "../scMarco")
legacy_ding <- normalizePath(if (length(args) >= 2) args[2] else "../scMarco-Ding-lab")
set.seed(20260928)

suppressMessages({
  source("src/combo_core.R")
  source("src/backend_h5ad.R")
  library(DBI)
  library(RSQLite)
})

load_legacy <- function(dir, files) {
  env <- new.env()
  old <- setwd(dir)
  on.exit(setwd(old))
  for (f in files) suppressMessages(sys.source(f, envir = env))
  env
}

mismatches <- list()
n_cases <- 0
report <- function(case, what, legacy, new) {
  n_cases <<- n_cases + 1
  if (!identical(legacy, new)) {
    mismatches[[length(mismatches) + 1]] <<- list(
      case = case, what = what,
      only_legacy = head(setdiff(legacy, new), 10),
      only_new = head(setdiff(new, legacy), 10)
    )
    cat("MISMATCH", what, case, "\n")
  } else {
    cat("ok      ", what, case, "\n")
  }
}

as_pairs <- function(counts) {
  # named off-target counts -> sorted "gene=count" strings (ties may reorder)
  if (length(counts) == 0) return(character(0))
  sort(paste0(names(counts), "=", counts))
}

# ---------------------------------------------------------------------------
# Ding-lab: HDF5 v1.0 backend vs. converted h5ad
# ---------------------------------------------------------------------------
ding <- load_legacy(legacy_ding, c("src/hdf5_backend.R", "src/utils.R"))
h5 <- ding$h5_connect(file.path(legacy_ding, "db.h5"))

for (ds_id in c("dmel_mixed", "four_species")) {
  ds <- ds_load(ds_id, list(
    label = ds_id,
    clusters = file.path("data", paste0(ds_id, ".clusters.h5ad"))
  ))
  legacy_genes <- ding$h5_get_bimodal_genes(h5, ds_id)
  stopifnot(setequal(legacy_genes, ds$analysis_genes))
  for (cut in c(0.3, 0.5, 0.8)) {
    b <- ds_binarize(ds, cut)
    for (i in 1:4) {
      splits <- sample(ds$splits, sample(seq_along(ds$splits), 1))
      ident <- sample(ds$clusters, 1)
      case <- sprintf("[%s cut=%.1f %s | %s]", ds_id, cut, ident,
                      paste(splits, collapse = "+"))

      legacy <- suppressMessages(ding$GetAllExpressedGenesForThis(
        h5_file = h5, dataset_name = ds_id,
        split = ds$splits %in% splits, splits = ds$splits,
        ident = ident, output_limit = 1e6, cut_off = cut,
        genes_list = legacy_genes
      ))
      new <- from_cluster(b, ident, splits, 1, limit = 1e6,
                          genes = ds$analysis_genes)
      report(case, "ding from_cluster", as_pairs(legacy), as_pairs(new))

      idents <- sample(ds$clusters, sample(2:5, 1))
      case <- sprintf("[%s cut=%.1f %s | %s]", ds_id, cut,
                      paste(idents, collapse = ","), paste(splits, collapse = "+"))
      legacy <- suppressMessages(ding$GetDistinctiveGenesMatrix(
        h5_file = h5, dataset_name = ds_id, selected_splits = splits,
        selected_clusters = idents, cut_off = cut, genes_list = legacy_genes,
        k = 1
      ))
      new <- distinct_within(b, idents, splits, 1, genes = ds$analysis_genes)
      report(case, "ding distinct   ",
             sort(paste(legacy$cluster, legacy$distinct_gene)),
             sort(paste(new$cluster, new$gene)))
    }
  }
}

# ---------------------------------------------------------------------------
# scMarco: SQLite long table vs. converted h5ad
# ---------------------------------------------------------------------------
sc <- load_legacy(legacy_sc, "src/utils.R")
db <- dbConnect(SQLite(), file.path(legacy_sc, "int/db.sqlite"), flags = SQLITE_RO)
sc_config <- yaml::read_yaml(file.path(legacy_sc, "config.yaml"))

for (tbl in c("ozel_2021_new", "kurmangaliyev_2020_new")) {
  ds <- ds_load(tbl, list(
    label = tbl, clusters = file.path("data", paste0(tbl, ".clusters.h5ad"))
  ))
  stages <- sc_config$databases[[tbl]]$stages
  stopifnot(identical(stages, ds$splits))
  for (cut in c(0.5, 0.8)) {
    message("Querying ", tbl, " at cut ", cut, " (full table scan)...")
    prob_mat <- dbGetQuery(db, paste(
      "SELECT * FROM", tbl, sc$wherecat(c(paste("value >", cut), "type = 'prob'"))
    ))
    b <- ds_binarize(ds, cut)
    for (i in 1:5) {
      n_sel <- sample(0:2, 1)
      sel <- sample(stages, n_sel)
      count <- sample(max(1, n_sel):min(length(stages), n_sel + 2), 1)
      ident <- sample(ds$clusters, 1)
      case <- sprintf("[%s cut=%.1f %s | req=%s n>=%d]", tbl, cut, ident,
                      paste(sel, collapse = "+"), count)

      legacy <- sc$GetAllExpressedGenesForThis(
        prob_mat = prob_mat, stage = stages %in% sel, count = count,
        stages = stages, ident = ident, output_limit = 1e6
      )
      legacy <- legacy[!is.na(legacy)] # legacy pads with NA (bug 11)
      new <- from_cluster(b, ident, sel, count, limit = 1e6)
      report(case, "sc from_cluster  ", as_pairs(legacy), as_pairs(new))

      idents <- sample(unique(prob_mat$cluster), sample(2:5, 1))
      case <- sprintf("[%s cut=%.1f %s | req=%s n>=%d]", tbl, cut,
                      paste(idents, collapse = ","), paste(sel, collapse = "+"), count)
      x <- sc$GetAllExpressedGenesDb(prob_mat, stage = stages %in% sel,
                                     count = count, stages = stages)
      legacy <- sc$FindDistinctiveGenes(x, idents = idents, max_presence = 1)
      legacy <- legacy[legacy$distinct_gene != "Not found", ] # bug 10
      new <- distinct_within(b, idents, sel, count)
      report(case, "sc distinct      ",
             sort(paste(legacy$cluster, legacy$distinct_gene)),
             sort(paste(new$cluster, new$gene)))
    }
  }
}
dbDisconnect(db)

cat(sprintf("\n%d comparisons, %d mismatches\n", n_cases, length(mismatches)))
if (length(mismatches) > 0) {
  str(mismatches)
  quit(status = 1)
}
