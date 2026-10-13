#!/usr/bin/env Rscript
# Convert one dataset of a scMarco-Ding-lab HDF5 "v1.0" database (db.h5) into
# the scMarco h5ad pair:
#   <out-dir>/<dataset>.clusters.h5ad  prob, pct_pos, mean_lognorm layers
#   <out-dir>/<dataset>.cells.h5ad     lognorm X (CSC), cluster labels, UMAPs
#
# Usage:
#   Rscript script/h5v1_to_h5ad.R --h5 db.h5 --dataset dmel_mixed \
#     --out-dir data --split-label sex [--ordered FALSE]

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
suppressPackageStartupMessages(library(rhdf5))

args <- parse_cli(list(
  h5 = NA, dataset = NA, `out-dir` = "data", `split-label` = "split",
  ordered = "FALSE"
))
ordered <- as.logical(args$ordered)
ds <- args$dataset
rd <- function(obj, ...) h5read(args$h5, paste0(ds, "/", obj), ...)

genes <- as.character(rd("gene_names"))
clusters <- as.character(rd("cluster_names"))
splits <- as.character(rd("splits"))
bimodal_idx <- as.integer(rd("bimodal_genes")) + 1L
message(sprintf("%s: %d genes, %d clusters, splits: %s",
                ds, length(genes), length(clusters), paste(splits, collapse = ", ")))

listing <- h5ls(args$h5)
has <- function(obj) any(paste0(listing$group, "/", listing$name) == paste0("/", ds, "/", obj))

obs_split <- character(0)
obs_cluster <- character(0)
prob_rows <- list()
pct_rows <- list()
mean_rows <- list()
n_cells <- integer(0)
bimodal_by_split <- matrix(
  FALSE, length(genes), length(splits), dimnames = list(genes, splits)
)

cell_x <- list()
cell_barcode <- character(0)
cell_split <- character(0)
cell_cluster <- character(0)
cell_offset <- c()

for (s in splits) {
  prob <- rd(paste0("probability/", s)) # genes x clusters
  dimnames(prob) <- list(genes, clusters)
  present <- colSums(!is.na(prob)) > 0
  nonbimodal <- rowSums(prob[, present, drop = FALSE] == -1, na.rm = TRUE) ==
    sum(present)
  bimodal_by_split[, s] <- !nonbimodal
  prob[prob == -1] <- NA

  pct <- if (has(paste0("pct_pos/", s))) {
    m <- rd(paste0("pct_pos/", s))
    dimnames(m) <- list(genes, clusters)
    m
  } else {
    NULL
  }

  # Cell level: CSR gene x cell == CSC cell x gene
  barcodes <- as.character(rd(paste0("barcodes/", s)))
  labels <- as.character(rd(paste0("clusters/", s)))
  shape <- as.numeric(rd(paste0("lognorm/", s, "/shape")))
  x <- sparseMatrix(
    i = as.integer(rd(paste0("lognorm/", s, "/indices"))) + 1L,
    p = as.integer(rd(paste0("lognorm/", s, "/indptr"))),
    x = as.numeric(rd(paste0("lognorm/", s, "/data"))),
    dims = c(shape[2], shape[1])
  )
  stopifnot(nrow(x) == length(barcodes), ncol(x) == length(genes))
  cell_offset[s] <- length(cell_barcode)
  cell_x[[s]] <- x
  cell_barcode <- c(cell_barcode, barcodes)
  cell_split <- c(cell_split, rep(s, length(barcodes)))
  cell_cluster <- c(cell_cluster, labels)

  for (cl in clusters[present]) {
    in_cl <- labels == cl
    obs_split <- c(obs_split, s)
    obs_cluster <- c(obs_cluster, cl)
    n_cells <- c(n_cells, sum(in_cl))
    prob_rows[[length(prob_rows) + 1]] <- prob[, cl]
    if (!is.null(pct)) pct_rows[[length(pct_rows) + 1]] <- pct[, cl]
    mean_rows[[length(mean_rows) + 1]] <- if (any(in_cl)) {
      Matrix::colMeans(x[in_cl, , drop = FALSE])
    } else {
      rep(NA_real_, length(genes))
    }
  }
  message(sprintf("  %s: %d cells, %d clusters present, %d non-bimodal genes",
                  s, length(barcodes), sum(present), sum(nonbimodal)))
}

# The legacy importer marks bimodality per gene (exprType) and zero-fills
# genes missing from the model; neither kind is bimodal in any split
is_bimodal <- seq_along(genes) %in% bimodal_idx
bimodal_by_split <- bimodal_by_split & is_bimodal

layers <- list(prob = do.call(rbind, prob_rows))
if (length(pct_rows) == length(prob_rows)) layers$pct_pos <- do.call(rbind, pct_rows)
layers$mean_lognorm <- do.call(rbind, mean_rows)

dir.create(args$`out-dir`, showWarnings = FALSE, recursive = TRUE)
clusters_path <- file.path(args$`out-dir`, paste0(ds, ".clusters.h5ad"))
write_clusters_h5ad(
  clusters_path, obs_split, obs_cluster, splits, sort(clusters), genes, layers,
  is_bimodal = is_bimodal,
  bimodal_by_split = bimodal_by_split, n_cells = n_cells,
  split_label = args$`split-label`, ordered_splits = ordered,
  extra_uns = list(
    source = paste0("h5v1:", basename(args$h5), ":", ds),
    mean_lognorm_method = "arithmetic mean of log-normalized values per cluster"
  )
)
message("Wrote ", clusters_path)

# Reductions: one obsm matrix per name, NaN for cells outside the subset
n_total <- length(cell_barcode)
reductions <- list()
reduction_meta <- list()
for (s in splits) {
  red_group <- paste0("/", ds, "/reduction/", s)
  names_s <- listing$name[listing$group == red_group & listing$otype == "H5I_GROUP"]
  for (r in names_s) {
    base <- paste0("reduction/", s, "/", r)
    coords <- rd(paste0(base, "/coords"))
    idx <- if (has(paste0(base, "/cell_indices"))) {
      as.integer(rd(paste0(base, "/cell_indices")))
    } else {
      seq_len(nrow(coords))
    }
    if (is.null(reductions[[r]])) {
      reductions[[r]] <- matrix(NaN, n_total, ncol(coords))
    }
    reductions[[r]][cell_offset[s] + idx, ] <- coords
    if (has(paste0(base, "/metadata"))) {
      reduction_meta[[r]] <- as.character(rd(paste0(base, "/metadata")))
    }
  }
}

cells_path <- file.path(args$`out-dir`, paste0(ds, ".cells.h5ad"))
write_cells_h5ad(
  cells_path, do.call(rbind, cell_x), cell_barcode, cell_split, cell_cluster,
  splits, sort(clusters), genes, reductions = reductions,
  reduction_meta = reduction_meta
)
message("Wrote ", cells_path, " with reductions: ",
        paste(names(reductions), collapse = ", "))

check_with_python(clusters_path)
check_with_python(cells_path)
