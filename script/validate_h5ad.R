#!/usr/bin/env Rscript
# Validate a scMarco dataset (format_version 2.0) before deploying it.
#
# Usage:
#   Rscript script/validate_h5ad.R --clusters data/x.clusters.h5ad \
#     [--cells data/x.cells.h5ad] [--stage built|binarized]
#
# --stage built accepts a dataset without probabilities (fresh from
# seurat_to_h5ad.R); the default, binarized, is what the app needs.
# Exits non-zero when any error is found; warnings are printed only.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(get_script_dir())
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(clusters = NA, cells = "", stage = "binarized"))
if (!args$stage %in% c("built", "binarized")) {
  stop("--stage must be built or binarized, not ", args$stage)
}
errors <- character(0)
warnings <- character(0)
err <- function(...) errors <<- c(errors, paste0(...))
warn <- function(...) warnings <<- c(warnings, paste0(...))

ch <- h5ad_open(args$clusters)
meta <- ch$uns$scmarco
if (is.null(meta)) {
  err("uns/scmarco is missing")
  meta <- list()
} else if (!identical(as.character(meta$format_version), SCMARCO_FORMAT_VERSION)) {
  err("format_version is ", meta$format_version, ", expected ", SCMARCO_FORMAT_VERSION)
}

for (col in c("split", "cluster")) {
  if (!is.factor(ch$obs[[col]])) err("obs/", col, " must be categorical")
}
pairs <- paste(ch$obs$split, ch$obs$cluster, sep = "|")
if (anyNA(ch$obs$split) || anyNA(ch$obs$cluster)) err("obs split/cluster contain missing values")
if (anyDuplicated(pairs)) err("obs has duplicated split x cluster rows")
if (!identical(rownames(ch$obs), pairs)) err("obs index is not \"<split>|<cluster>\"")
for (col in c("split", "cluster")) {
  lv <- levels(ch$obs[[col]])
  if (any(grepl("|", lv, fixed = TRUE))) err(col, " labels contain '|'")
  if (any(lv == "")) err(col, " labels include an empty string")
  if (any(lv != trimws(lv))) err(col, " labels start or end with whitespace")
}
names_ok <- tryCatch(
  read_cluster_names(ch$obs, levels(ch$obs$cluster)),
  error = function(e) {
    err("obs/cluster_name must be the same in every split of a cluster")
    NULL
  }
)
if (!is.null(names_ok) && all(is.na(names_ok))) {
  embedded <- split_cluster_labels(levels(ch$obs$cluster))
  if (any(!is.na(embedded$name))) {
    warn(sum(!is.na(embedded$name)), " cluster IDs look like \"name (id)\" labels, e.g. ",
         embedded$label[!is.na(embedded$name)][1], "; split them with ",
         "script/set_cluster_names.R --split-labels TRUE")
  }
}
empty <- setdiff(levels(ch$obs$cluster), as.character(ch$obs$cluster))
if (length(empty) > 0) {
  warn(length(empty), " cluster(s) have no rows in any split (always absent), e.g. ",
       paste(head(empty, 5), collapse = ", "))
}
if (anyDuplicated(ch$genes)) err("var index (genes) is not unique")
if (!is.null(meta$split_order) &&
    !setequal(as.character(meta$split_order), levels(ch$obs$split))) {
  err("uns split_order does not match the obs split categories")
}

binarized <- "prob" %in% ch$layers
if (!binarized && args$stage == "binarized") {
  err("layers/prob is missing; add it with script/add_probability_h5ad.R ",
      "(or validate with --stage built)")
}
if (binarized && !is.null(meta$binarization)) {
  missing <- setdiff(c("method", "source", "date"), names(meta$binarization))
  if (length(missing) > 0) {
    err("uns/scmarco/binarization lacks ", paste(missing, collapse = ", "))
  }
}
for (layer in intersect(c("prob", "pct_pos", "mean_lognorm"), ch$layers)) {
  m <- h5ad_layer(ch, layer)
  if (!identical(dim(m), c(nrow(ch$obs), length(ch$genes)))) {
    err("layers/", layer, " has shape ", paste(dim(m), collapse = "x"))
    next
  }
  if (layer %in% c("prob", "pct_pos")) {
    bad <- sum(!is.na(m) & (m < 0 | m > 1))
    if (bad > 0) err("layers/", layer, ": ", bad, " values outside [0, 1]")
  }
  if (all(is.na(m))) warn("layers/", layer, " is entirely missing")
}

if ("bimodal_by_split" %in% ch$varm) {
  bim <- h5ad_varm(ch, "bimodal_by_split")
  if (nrow(bim) != length(ch$genes) || ncol(bim) != length(meta$bimodal_splits)) {
    err("varm/bimodal_by_split shape does not match genes x uns bimodal_splits")
  }
  if (!setequal(as.character(meta$bimodal_splits), levels(ch$obs$split))) {
    err("uns bimodal_splits does not match the obs split categories")
  }
}
if (binarized && !"is_bimodal" %in% names(ch$var)) {
  warn("var/is_bimodal is missing: every gene will be used for marker finding")
}

if (nzchar(args$cells)) {
  xh <- h5ad_open(args$cells)
  if (!identical(xh$x_encoding, "csc_matrix")) {
    err("cells X is ", xh$x_encoding, "; per-gene reads need csc_matrix")
  }
  if (!identical(xh$genes, ch$genes)) {
    err("cells var (genes) differs from the clusters file")
  }
  if (!setequal(levels(xh$obs$split), levels(ch$obs$split))) {
    err("cells split categories differ from the clusters file")
  }
  cell_pairs <- unique(paste(xh$obs$split, xh$obs$cluster, sep = "|"))
  missing <- setdiff(cell_pairs, pairs)
  if (length(missing) > 0) {
    warn(length(missing), " split x cluster pairs have cells but no clusters row, e.g. ",
         paste(head(missing, 3), collapse = ", "))
  }
  for (key in xh$obsm) {
    m <- h5ad_obsm(xh, key)
    if (nrow(m) != nrow(xh$obs)) err("obsm/", key, " has ", nrow(m), " rows")
    if (all(is.na(m))) warn("obsm/", key, " has no coordinates")
  }
  if (length(xh$obsm) == 0) warn("cells file has no reductions: the UMAP tab stays hidden")
}

for (w in warnings) message("WARNING: ", w)
for (e in errors) message("ERROR:   ", e)
if (length(errors) > 0) {
  quit(status = 1)
}
message("OK (", if (binarized) "binarized" else "built", "): ", args$clusters,
        if (nzchar(args$cells)) paste(" +", args$cells))
