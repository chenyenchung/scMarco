#!/usr/bin/env Rscript
# Build a new scMarco dataset from a Seurat object and mixture-model output.
# Replaces the legacy seurat_to_h5.R + add_probability_to_h5.R +
# append_umaps_to_h5.R pipeline and writes the h5ad pair directly.
#
# Usage:
#   Rscript script/seurat_to_h5ad.R --seurat obj.rds --prob prob.rds \
#     --dataset dmel_mixed --out-dir data --split-label sex \
#     [--umaps dir_of_subset_objects] [--split-order Male,Female] [--ordered FALSE]
#
# Inputs (same conventions as the legacy scripts):
#   --seurat  Seurat object whose identities are "<split>@<cluster>"; raw counts
#             in the default assay. Normalized as Seurat's LogNormalize
#             (log1p(counts / library size * 1e4)). A "umap" reduction, if
#             present, becomes obsm X_umap_full.
#   --prob    data.frame, genes as row names, one column per "<split>@<cluster>"
#             plus the mixture-model QC columns (exprType marks bimodal genes).
#   --umaps   optional folder of Seurat objects (.rds) with a "umap" reduction
#             on a subset of the cells; the file name before the first "_"
#             names the reduction.
# Only SeuratObject is needed, not Seurat.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
suppressPackageStartupMessages(library(SeuratObject))

args <- parse_cli(list(
  seurat = NA, prob = NA, dataset = NA, `out-dir` = "data",
  `split-label` = "split", `split-order` = "", ordered = "FALSE", umaps = ""
))
QC_COLS <- c("diff.elpd", "diff.elpd.se", "bimodal.on_level", "bimodal.off_level",
             "bimodal.sd", "bimodal.pi", "unimodal.level", "unimodal.sd", "exprType")

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
ident <- as.character(Idents(obj))
parsed <- strsplit(ident, "@", fixed = TRUE)
if (any(lengths(parsed) != 2)) {
  stop("Identities must look like '<split>@<cluster>', e.g. ",
       paste(head(ident[lengths(parsed) != 2], 3), collapse = ", "))
}
cell_split <- vapply(parsed, `[`, "", 1)
cell_cluster <- vapply(parsed, `[`, "", 2)
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
clusters <- sort(unique(cell_cluster))
message(sprintf("%d cells, %d genes, %d splits, %d clusters",
                length(barcodes), length(genes), length(splits), length(clusters)))

# --- Mixture-model probabilities ---------------------------------------------
prob_df <- readRDS(args$prob)
prob_cols <- setdiff(colnames(prob_df), QC_COLS)
if (!"exprType" %in% colnames(prob_df)) stop("--prob lacks the exprType column")
col_parts <- strsplit(prob_cols, "@", fixed = TRUE)
if (any(lengths(col_parts) != 2)) stop("Probability columns must be '<split>@<cluster>'")
col_split <- vapply(col_parts, `[`, "", 1)
col_cluster <- vapply(col_parts, `[`, "", 2)
# the legacy pipeline wrote hyphens where the object has underscores
fix <- !col_cluster %in% clusters & gsub("-", "_", col_cluster) %in% clusters
col_cluster[fix] <- gsub("-", "_", col_cluster[fix])
unknown <- unique(paste(col_split, col_cluster, sep = "@")[
  !(col_split %in% splits & col_cluster %in% clusters)])
if (length(unknown) > 0) {
  stop("Probability columns without matching cells: ", paste(head(unknown, 5), collapse = ", "))
}
extra_genes <- setdiff(rownames(prob_df), genes)
if (length(extra_genes) > 0) {
  stop(length(extra_genes), " genes in --prob are not in the Seurat object, e.g. ",
       paste(head(extra_genes, 5), collapse = ", "))
}
is_bimodal <- genes %in% rownames(prob_df)[prob_df$exprType == "bimodal"]
message(sprintf("%d bimodal genes; %d genes absent from the model",
                sum(is_bimodal), length(setdiff(genes, rownames(prob_df)))))

# --- Cluster-level layers ----------------------------------------------------
obs_split <- character(0)
obs_cluster <- character(0)
n_cells <- integer(0)
rows <- list(prob = list(), pct_pos = list(), mean_lognorm = list())
gene_idx <- match(genes, rownames(prob_df))
for (s in splits) {
  for (cl in clusters) {
    in_cl <- cell_split == s & cell_cluster == cl
    if (!any(in_cl)) next
    obs_split <- c(obs_split, s)
    obs_cluster <- c(obs_cluster, cl)
    n_cells <- c(n_cells, sum(in_cl))
    col <- which(col_split == s & col_cluster == cl)
    p <- rep(NA_real_, length(genes))
    if (length(col) == 1) {
      p <- as.numeric(prob_df[[prob_cols[col]]])[gene_idx]
      p[!is_bimodal] <- NA
    } else {
      warning("No probability column for ", s, "@", cl, "; its genes will never be ON")
    }
    xs <- x[in_cl, , drop = FALSE]
    rows$prob[[length(rows$prob) + 1]] <- p
    rows$pct_pos[[length(rows$pct_pos) + 1]] <- Matrix::colMeans(xs > 0)
    rows$mean_lognorm[[length(rows$mean_lognorm) + 1]] <- Matrix::colMeans(xs)
  }
}
layers <- lapply(rows, function(r) do.call(rbind, r))
bimodal_by_split <- matrix(is_bimodal, length(genes), length(splits),
                           dimnames = list(genes, splits))

dir.create(args$`out-dir`, showWarnings = FALSE, recursive = TRUE)
clusters_path <- file.path(args$`out-dir`, paste0(args$dataset, ".clusters.h5ad"))
write_clusters_h5ad(
  clusters_path, obs_split, obs_cluster, splits, clusters, genes, layers,
  is_bimodal = is_bimodal, bimodal_by_split = bimodal_by_split,
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
check_with_python(clusters_path)
check_with_python(cells_path)
