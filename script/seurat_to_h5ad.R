#!/usr/bin/env Rscript
# Build a new scMarco dataset from a Seurat object. Replaces the legacy
# seurat_to_h5.R + append_umaps_to_h5.R pipeline and writes the h5ad pair
# directly. The result is *built*: the app needs the mixture-model
# probabilities too, added here with --prob or later with
# script/add_probability_h5ad.R.
#
# Usage:
#   Rscript script/seurat_to_h5ad.R --seurat obj.rds --dataset dmel_mixed \
#     --out-dir data --split-label sex \
#     [--split-col scmarco_split] [--cluster-col scmarco_cluster] [--assay RNA] \
#     [--prob Male=male.rds,Female=female.rds] \
#     [--umaps dir_of_subset_objects] [--split-order Male,Female] [--ordered FALSE] \
#     [--reduction umap] [--exclude-clusters 0]
#
#   # one object per split: the split comes from the name, not --split-col
#   Rscript script/seurat_to_h5ad.R --seurat P15=P15.rds,P30=P30.rds ...
#
# Inputs (see docs/expected_input.md):
#   --seurat  Seurat object with raw counts in the default assay, normalized
#             here as Seurat's LogNormalize (log1p(counts / library size * 1e4)).
#             Each cell's split and cluster ID come from the metadata columns
#             --split-col / --cluster-col (--split-col none: one split, "all").
#             Without either default column, legacy "<split>@<cluster>"
#             identities are read instead. A "umap" reduction, if present,
#             becomes obsm X_umap_full.
#             Seurat v3 objects (saved before SeuratObject) are updated on read.
#             Or one object per split (Split1=a.rds,Split2=b.rds): the genes
#             are the union of the objects' genes (0 counts where a gene is
#             missing), and an object's "umap" reduction becomes obsm X_<split>.
#   --assay   assay holding the raw counts; defaults to the object's default
#             assay (e.g. pass RNA when the default is "integrated").
#   --reduction  reduction to store instead of "umap" (e.g. tsne); stored as
#             obsm X_<reduction>_full, or X_<split>_<reduction> per split.
#             An object without it is an error that lists its reductions;
#             --reduction none builds without one (the UMAP tab stays hidden).
#   --exclude-clusters  comma-separated cluster IDs whose cells are dropped
#             (e.g. 0 for cells discarded by the study).
#   --prob    optional; as for script/add_probability_h5ad.R.
#   --umaps   optional folder of Seurat objects (.rds) with a "umap" reduction
#             on a subset of the cells; the file name before the first "_"
#             names the reduction.
# Only SeuratObject is needed, not Seurat.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
source(file.path(get_script_dir(), "lib_import.R"))
suppressPackageStartupMessages(library(SeuratObject))

args <- parse_cli(list(
  seurat = NA, prob = "", dataset = NA, `out-dir` = "data",
  `split-col` = "scmarco_split", `cluster-col` = "scmarco_cluster", assay = "",
  `split-label` = "split", `split-order` = "", ordered = "FALSE", umaps = "",
  reduction = "umap", `exclude-clusters` = ""
))
red_name <- args$reduction
excluded <- if (nzchar(args$`exclude-clusters`)) strsplit(args$`exclude-clusters`, ",")[[1]] else character(0)

get_counts <- function(obj, assay = "") {
  if (!nzchar(assay)) assay <- DefaultAssay(obj)
  if (!assay %in% Assays(obj)) {
    fail("Assay '", assay, "' not found. Assays: ", paste(Assays(obj), collapse = ", "))
  }
  counts <- tryCatch(
    LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
  as(counts, "CsparseMatrix")
}

lognorm_by_cell <- function(counts) {
  #' genes x cells counts -> genes x cells LogNormalize values (genes sorted)
  counts <- counts[sort(rownames(counts)), , drop = FALSE]
  lib_size <- Matrix::colSums(counts)
  x <- counts %*% Matrix::Diagonal(x = 1e4 / pmax(lib_size, 1))
  x@x <- log1p(x@x)
  dimnames(x) <- dimnames(counts)
  as(x, "CsparseMatrix")
}

expand_genes <- function(x, genes) {
  #' genes x cells dgCMatrix with sorted genes -> the same on the sorted
  #' superset `genes`, with zeros for the added genes. Row indices map
  #' monotonically, so each column stays sorted.
  if (identical(rownames(x), genes)) return(x)
  idx <- match(rownames(x), genes)
  out <- x
  out@i <- as.integer(idx[x@i + 1] - 1)
  out@Dim <- c(length(genes), ncol(x))
  out@Dimnames <- list(genes, colnames(x))
  out
}

# --- Cells -------------------------------------------------------------------
inputs <- parse_seurat_arg(args$seurat)
# Read --prob before the slow steps so a bad path fails early
prob_tables <- if (nzchar(args$prob)) read_prob_arg(args$prob) else NULL
per_split <- !is.null(names(inputs))
if (per_split && args$`split-col` != "scmarco_split") {
  fail("--split-col cannot be used with one object per split: ",
       "the split comes from the names given in --seurat")
}

parts <- list()
for (k in seq_along(inputs)) {
  obj <- read_seurat(inputs[[k]])
  labels <- if (per_split) {
    lab <- seurat_labels(obj[[]], as.character(Idents(obj)), "none",
                         args$`cluster-col`, defaults = FALSE)
    lab$split <- rep(names(inputs)[k], length(lab$split))
    lab
  } else {
    seurat_labels(
      obj[[]], as.character(Idents(obj)), args$`split-col`, args$`cluster-col`,
      defaults = args$`split-col` == "scmarco_split" &&
        args$`cluster-col` == "scmarco_cluster"
    )
  }
  keep <- !labels$cluster %in% excluded
  if (!all(keep)) {
    message(sprintf("  dropping %d cells of excluded cluster(s)", sum(!keep)))
    labels <- lapply(labels, `[`, keep)
  }
  barcodes <- Cells(obj)[keep]
  red <- if (red_name != "none") {
    seurat_reduction(obj, red_name, barcodes,
                     if (per_split) paste0(names(inputs)[k], ": ") else "")
  }
  parts[[k]] <- list(
    x = lognorm_by_cell(get_counts(obj, args$assay)[, barcodes, drop = FALSE]),
    barcodes = barcodes, split = labels$split, cluster = labels$cluster, umap = red
  )
  if (per_split) {
    message(sprintf("  %s: %d cells, %d genes", names(inputs)[k],
                    length(barcodes), nrow(parts[[k]]$x)))
  }
  rm(obj)
  gc()
}

cell_split <- unlist(lapply(parts, `[[`, "split"))
cell_cluster <- unlist(lapply(parts, `[[`, "cluster"))
check_labels(cell_split, cell_cluster)
barcodes <- unlist(lapply(parts, `[[`, "barcodes"))
genes <- sort(unique(unlist(lapply(parts, function(p) rownames(p$x)))))
n_part <- vapply(parts, function(p) length(p$barcodes), 1L)

reductions <- list()
for (k in seq_along(parts)) {
  if (is.null(parts[[k]]$umap)) next
  if (!per_split) {
    reductions[[reduction_key(red_name)]] <- parts[[k]]$umap
    next
  }
  m <- matrix(NaN, length(barcodes), 2)
  m[sum(n_part[seq_len(k - 1)]) + seq_len(n_part[k]), ] <- parts[[k]]$umap
  reductions[[reduction_key(red_name, names(inputs)[k])]] <- m
}

x <- do.call(cbind, lapply(parts, function(p) expand_genes(p$x, genes)))
rm(parts)
gc()
x <- Matrix::t(x) # cells x genes
dimnames(x) <- list(barcodes, genes)

splits <- if (nzchar(args$`split-order`)) {
  strsplit(args$`split-order`, ",")[[1]]
} else if (per_split) {
  names(inputs)
} else {
  sort(unique(cell_split))
}
if (!setequal(splits, unique(cell_split))) {
  stop("--split-order must list exactly the splits: ", paste(sort(unique(cell_split)), collapse = ","))
}
clusters <- order_cluster_ids(unique(cell_cluster))
message(sprintf("%d cells, %d genes, %d splits, %d clusters",
                length(barcodes), length(genes), length(splits), length(clusters)))

# --- Cluster-level layers ----------------------------------------------------
# one row per split x cluster with cells, splits in order, then clusters
obs <- unique(data.frame(split = cell_split, cluster = cell_cluster))
obs <- obs[order(match(obs$split, splits), match(obs$cluster, clusters)), ]
obs_split <- obs$split
obs_cluster <- obs$cluster
cell_obs <- match(paste(cell_split, cell_cluster, sep = "|"),
                  paste(obs_split, obs_cluster, sep = "|"))
n_cells <- tabulate(cell_obs, nbins = nrow(obs))
# cluster sums as one product with an obs x cells indicator matrix; divide
# after summing so pct_pos stays within [0, 1] exactly
member <- Matrix::sparseMatrix(i = cell_obs, j = seq_along(cell_obs), x = 1,
                               dims = c(nrow(obs), length(cell_obs)))
x_pos <- x
x_pos@x <- as.numeric(x_pos@x > 0)
layers <- list(pct_pos = as.matrix(member %*% x_pos) / n_cells,
               mean_lognorm = as.matrix(member %*% x) / n_cells)
rm(x_pos)

dir.create(args$`out-dir`, showWarnings = FALSE, recursive = TRUE)
clusters_path <- file.path(args$`out-dir`, paste0(args$dataset, ".clusters.h5ad"))
write_clusters_h5ad(
  clusters_path, obs_split, obs_cluster, splits, clusters, genes, layers,
  n_cells = n_cells, split_label = args$`split-label`,
  ordered_splits = as.logical(args$ordered),
  extra_uns = list(
    source = paste0("seurat:", paste(basename(inputs), collapse = ",")),
    mean_lognorm_method = "arithmetic mean of log-normalized values per cluster"
  )
)
message("Wrote ", clusters_path)

# --- Reductions ----------------------------------------------------------------
if (nzchar(args$umaps)) {
  for (f in list.files(args$umaps, "\\.rds$", full.names = TRUE)) {
    name <- strsplit(sub("\\.rds$", "", basename(f)), "_")[[1]][1]
    sub_obj <- readRDS(f)
    if (!"umap" %in% Reductions(sub_obj)) {
      warning(basename(f), " has no umap reduction; skipped")
      next
    }
    emb <- Embeddings(sub_obj, "umap")
    m <- matrix(NaN, length(barcodes), 2)
    hit <- match(rownames(emb), barcodes)
    m[hit[!is.na(hit)], ] <- emb[!is.na(hit), 1:2]
    reductions[[name]] <- m
    message(sprintf("  reduction %s: %d cells", name, sum(!is.na(hit))))
  }
}

cells_path <- file.path(args$`out-dir`, paste0(args$dataset, ".cells.h5ad"))
write_cells_h5ad(
  cells_path, x, barcodes, cell_split, cell_cluster, splits, clusters, genes,
  reductions = reductions
)
message("Wrote ", cells_path)

if (!is.null(prob_tables)) {
  p <- assemble_prob(prob_tables, obs_split, obs_cluster, splits, genes)
  add_probability_layer(clusters_path, p$prob, p$bimodal_by_split, list(
    method = "mixture model", source = args$prob, date = format(Sys.Date())
  ))
  message(sprintf("Binarized %s: %d genes bimodal in at least one %s",
                  clusters_path, sum(rowSums(p$bimodal_by_split) > 0), args$`split-label`))
} else {
  message("Built without probabilities; add them before use with:\n  Rscript ",
          "script/add_probability_h5ad.R --clusters ", clusters_path,
          " --prob <split>=<file.rds>,...")
}
check_with_python(clusters_path)
check_with_python(cells_path)
