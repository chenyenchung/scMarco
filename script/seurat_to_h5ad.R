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
#     [--split-col scmarco_split] [--cluster-col scmarco_cluster] \
#     [--prob Male=male.rds,Female=female.rds] \
#     [--umaps dir_of_subset_objects] [--split-order Male,Female] [--ordered FALSE]
#
# Inputs (see docs/expected_input.md):
#   --seurat  Seurat object with raw counts in the default assay, normalized
#             here as Seurat's LogNormalize (log1p(counts / library size * 1e4)).
#             Each cell's split and cluster ID come from the metadata columns
#             --split-col / --cluster-col (--split-col none: one split, "all").
#             Without either default column, legacy "<split>@<cluster>"
#             identities are read instead. A "umap" reduction, if present,
#             becomes obsm X_umap_full.
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
  `split-col` = "scmarco_split", `cluster-col` = "scmarco_cluster",
  `split-label` = "split", `split-order` = "", ordered = "FALSE", umaps = ""
))

get_counts <- function(obj) {
  assay <- DefaultAssay(obj)
  counts <- tryCatch(
    LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
  as(counts, "CsparseMatrix")
}

# --- Cells -------------------------------------------------------------------
obj <- readRDS(args$seurat)
# Read --prob before the slow steps so a bad path fails early
prob_tables <- if (nzchar(args$prob)) read_prob_arg(args$prob) else NULL
labels <- seurat_labels(
  obj[[]], as.character(Idents(obj)), args$`split-col`, args$`cluster-col`,
  defaults = args$`split-col` == "scmarco_split" &&
    args$`cluster-col` == "scmarco_cluster"
)
check_labels(labels$split, labels$cluster)
cell_split <- labels$split
cell_cluster <- labels$cluster
barcodes <- Cells(obj)

counts <- get_counts(obj)
genes <- sort(rownames(counts))
counts <- counts[genes, , drop = FALSE]
lib_size <- Matrix::colSums(counts)
x <- counts %*% Matrix::Diagonal(x = 1e4 / pmax(lib_size, 1))
x@x <- log1p(x@x)
x <- Matrix::t(x) # cells x genes
dimnames(x) <- list(barcodes, genes)

splits <- if (nzchar(args$`split-order`)) {
  strsplit(args$`split-order`, ",")[[1]]
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
obs_split <- character(0)
obs_cluster <- character(0)
n_cells <- integer(0)
rows <- list(pct_pos = list(), mean_lognorm = list())
for (s in splits) {
  for (cl in clusters) {
    in_cl <- cell_split == s & cell_cluster == cl
    if (!any(in_cl)) next
    obs_split <- c(obs_split, s)
    obs_cluster <- c(obs_cluster, cl)
    n_cells <- c(n_cells, sum(in_cl))
    xs <- x[in_cl, , drop = FALSE]
    rows$pct_pos[[length(rows$pct_pos) + 1]] <- Matrix::colMeans(xs > 0)
    rows$mean_lognorm[[length(rows$mean_lognorm) + 1]] <- Matrix::colMeans(xs)
  }
}
layers <- lapply(rows, function(r) do.call(rbind, r))

dir.create(args$`out-dir`, showWarnings = FALSE, recursive = TRUE)
clusters_path <- file.path(args$`out-dir`, paste0(args$dataset, ".clusters.h5ad"))
write_clusters_h5ad(
  clusters_path, obs_split, obs_cluster, splits, clusters, genes, layers,
  n_cells = n_cells, split_label = args$`split-label`,
  ordered_splits = as.logical(args$ordered),
  extra_uns = list(
    source = paste0("seurat:", basename(args$seurat)),
    mean_lognorm_method = "arithmetic mean of log-normalized values per cluster"
  )
)
message("Wrote ", clusters_path)

# --- Reductions ----------------------------------------------------------------
reductions <- list()
if ("umap" %in% Reductions(obj)) {
  reductions$umap_full <- Embeddings(obj, "umap")[barcodes, 1:2, drop = FALSE]
}
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
