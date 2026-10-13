# Shared writers for the scMarco h5ad schema (format_version 2.0).
# Sourced by the converters in this folder; see README.md for the schema.

suppressPackageStartupMessages({
  library(anndataR)
  library(Matrix)
})

SCMARCO_FORMAT_VERSION <- "2.0"

write_clusters_h5ad <- function(path, obs_split, obs_cluster, splits, clusters,
                                genes, layers, is_bimodal = NULL,
                                bimodal_by_split = NULL, n_cells = NULL,
                                split_label = "split", ordered_splits = FALSE,
                                extra_uns = list()) {
  #' layers: named list of n_obs x n_genes numeric matrices ("prob" required)
  #' bimodal_by_split: optional n_genes x n_splits logical matrix
  stopifnot("prob" %in% names(layers))
  obs <- data.frame(
    split = factor(obs_split, levels = splits, ordered = ordered_splits),
    cluster = factor(obs_cluster, levels = clusters),
    row.names = paste(obs_split, obs_cluster, sep = "|")
  )
  if (!is.null(n_cells)) obs$n_cells <- as.integer(n_cells)

  var <- data.frame(row.names = genes)
  if (!is.null(is_bimodal)) var$is_bimodal <- as.logical(is_bimodal)

  layers <- lapply(layers, function(m) {
    dimnames(m) <- list(rownames(obs), genes)
    m
  })

  varm <- list()
  uns <- list(
    format_version = SCMARCO_FORMAT_VERSION,
    split_label = split_label,
    split_order = splits,
    ordered_splits = ordered_splits
  )
  if (!is.null(bimodal_by_split)) {
    # anndataR writes logical matrices as flat vectors; store 0/1 integers
    m <- bimodal_by_split[genes, splits, drop = FALSE]
    storage.mode(m) <- "integer"
    varm$bimodal_by_split <- unname(m)
    uns$bimodal_splits <- splits
  }
  uns <- c(uns, extra_uns)

  ad <- AnnData(
    obs = obs, var = var, layers = layers,
    varm = if (length(varm) > 0) varm else NULL,
    uns = list(scmarco = uns)
  )
  write_h5ad(ad, path, mode = "w")
  rechunk_h5ad(
    path, n_obs = nrow(obs),
    dense = paste0("/layers/", names(layers)),
    gene_block = min(256, length(genes))
  )
  invisible(path)
}

write_cells_h5ad <- function(path, x, barcodes, cell_split, cell_cluster,
                             splits, clusters, genes, reductions = list(),
                             reduction_meta = list()) {
  #' x: n_cells x n_genes dgCMatrix (log-normalized). Written as csc_matrix,
  #'   which gives cheap per-gene reads.
  #' reductions: named list of n_cells x k matrices (NaN outside subsets)
  x <- as(x, "CsparseMatrix")
  obs <- data.frame(
    barcode = barcodes,
    split = factor(cell_split, levels = splits),
    cluster = factor(cell_cluster, levels = clusters),
    row.names = paste(cell_split, barcodes, sep = "|")
  )
  dimnames(x) <- list(rownames(obs), genes)
  obsm <- lapply(reductions, function(m) {
    m <- unname(m)
    storage.mode(m) <- "double"
    m
  })
  names(obsm) <- paste0("X_", names(reductions))

  ad <- AnnData(
    X = x, obs = obs, var = data.frame(row.names = genes),
    obsm = if (length(obsm) > 0) obsm else NULL,
    uns = list(scmarco = list(
      format_version = SCMARCO_FORMAT_VERSION,
      reduction_meta = reduction_meta
    ))
  )
  write_h5ad(ad, path, mode = "w")
  rechunk_h5ad(path, n_obs = nrow(obs), sparse_x = TRUE)
  invisible(path)
}

rechunk_h5ad <- function(path, n_obs, sparse_x = FALSE, dense = character(0),
                         gene_block = 256) {
  #' anndataR writes each array as a single chunk, so reading one gene would
  #' decompress everything. Re-chunk with h5repack (HDF5 command line tools):
  #' sparse X arrays in 64k blocks, dense obs x gene layers in gene blocks.
  if (Sys.which("h5repack") == "") {
    warning("h5repack not found; per-gene reads from ", path, " will be slow")
    return(invisible(FALSE))
  }
  layout <- character(0)
  if (sparse_x) {
    layout <- c(layout, "-l", "/X/data,/X/indices:CHUNK=65536")
  }
  for (obj in dense) {
    layout <- c(layout, "-l", sprintf("%s:CHUNK=%dx%d", obj, n_obs, gene_block))
  }
  if (length(layout) == 0) return(invisible(TRUE))
  tmp <- paste0(path, ".repack")
  status <- system2("h5repack", c("-f", "GZIP=4", layout, path, tmp))
  if (status != 0) {
    unlink(tmp)
    stop("h5repack failed on ", path)
  }
  file.rename(tmp, path)
  invisible(TRUE)
}

check_with_python <- function(path) {
  #' Confirm the file opens in Python anndata when a project venv exists
  py <- file.path(
    dirname(dirname(normalizePath(sys_script_path()))), ".venv", "bin", "python"
  )
  if (!file.exists(py)) {
    message("Skipping Python anndata check (no .venv)")
    return(invisible(NA))
  }
  out <- system2(py, c("-c", shQuote(sprintf(
    "import anndata; print(anndata.read_h5ad(%s))", shQuote(path, "sh")
  ))), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("Python anndata could not read ", path, ":\n", paste(out, collapse = "\n"))
  }
  message(paste(out, collapse = "\n"))
  invisible(TRUE)
}

sys_script_path <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg) == 0) return(file.path(getwd(), "script", "x.R"))
  sub("^--file=", "", arg[1])
}

parse_cli <- function(spec) {
  #' spec: named list of defaults (NA = required). Accepts --name value.
  args <- commandArgs(trailingOnly = TRUE)
  out <- spec
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(spec) || i == length(args)) {
      stop("Unknown or incomplete argument: ", args[i],
           "\nExpected: ", paste0("--", names(spec), collapse = " "))
    }
    out[[key]] <- args[i + 1]
    i <- i + 2
  }
  missing <- names(out)[vapply(out, function(v) identical(v, NA), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required argument(s): ", paste0("--", missing, collapse = " "))
  }
  out
}
