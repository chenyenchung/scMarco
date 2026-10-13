# h5ad backend
#
# A thin reader for the AnnData on-disk format (https://anndata.readthedocs.io
# /en/latest/fileformat-prose.html) built on rhdf5, so the app can read single
# genes without loading whole matrices and needs no Python at runtime.
#
# Each dataset is one required <ds>.clusters.h5ad (obs = split x cluster pairs
# that exist; layers prob / pct_pos / mean_lognorm) and one optional
# <ds>.cells.h5ad (X = log-normalized expression in CSC layout; obs = cells;
# obsm = reductions). See script/README.md for the full schema.

library(rhdf5)

# ---------------------------------------------------------------------------
# Low-level element readers
# ---------------------------------------------------------------------------

h5ad_attr <- function(path, obj, name) {
  attrs <- h5readAttributes(path, obj)
  attrs[[name]]
}

h5ad_exists <- function(path, obj) {
  fid <- H5Fopen(path, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  H5Lexists(fid, obj)
}

h5ad_list <- function(path, group) {
  if (!h5ad_exists(path, group)) {
    return(character(0))
  }
  fid <- H5Fopen(path, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))
  gid <- H5Gopen(fid, group)
  on.exit(H5Gclose(gid), add = TRUE)
  h5ls(gid, recursive = FALSE)$name
}

as_logical_enum <- function(x) {
  # h5py stores booleans as an enum that rhdf5 may return as a factor
  if (is.factor(x)) {
    return(as.logical(toupper(as.character(x))))
  }
  x
}

h5ad_read_column <- function(path, obj) {
  enc <- h5ad_attr(path, obj, "encoding-type")
  if (!is.null(enc) && enc == "categorical") {
    categories <- as.vector(h5read(path, paste0(obj, "/categories")))
    codes <- as.vector(h5read(path, paste0(obj, "/codes")))
    codes[codes < 0] <- NA
    ordered <- isTRUE(as_logical_enum(h5ad_attr(path, obj, "ordered")))
    return(factor(categories[codes + 1], levels = categories, ordered = ordered))
  }
  as.vector(as_logical_enum(h5read(path, obj)))
}

h5ad_read_df <- function(path, group) {
  index_name <- h5ad_attr(path, group, "_index")
  index <- as.character(h5read(path, paste0(group, "/", index_name)))
  col_order <- as.character(h5ad_attr(path, group, "column-order"))
  cols <- lapply(col_order, function(col) {
    h5ad_read_column(path, paste0(group, "/", col))
  })
  names(cols) <- col_order
  df <- if (length(cols) > 0) {
    as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  } else {
    data.frame(row.names = seq_along(index))
  }
  rownames(df) <- index
  df
}

h5ad_read_uns <- function(path, obj) {
  if (!h5ad_exists(path, obj)) {
    return(list())
  }
  enc <- h5ad_attr(path, obj, "encoding-type")
  if (!is.null(enc) && enc == "dict") {
    keys <- h5ad_list(path, obj)
    out <- lapply(keys, function(k) h5ad_read_uns(path, paste0(obj, "/", k)))
    names(out) <- keys
    return(out)
  }
  as.vector(as_logical_enum(h5read(path, obj)))
}

h5ad_read_dense <- function(path, obj, rows = NULL, cols = NULL) {
  # On disk: C order [n_rows, n_cols]; rhdf5 returns the transpose.
  mat <- h5read(path, obj, index = list(cols, rows))
  t(mat)
}

# ---------------------------------------------------------------------------
# File handles
# ---------------------------------------------------------------------------

h5ad_open <- function(path) {
  if (!file.exists(path)) {
    stop("h5ad file not found: ", path)
  }
  h <- list(
    path = path,
    mtime = file.info(path)$mtime,
    obs = h5ad_read_df(path, "obs"),
    var = h5ad_read_df(path, "var"),
    layers = h5ad_list(path, "layers"),
    obsm = h5ad_list(path, "obsm"),
    varm = h5ad_list(path, "varm"),
    has_x = h5ad_exists(path, "X"),
    uns = h5ad_read_uns(path, "uns")
  )
  h$genes <- rownames(h$var)
  if (h$has_x) {
    h$x_encoding <- h5ad_attr(path, "X", "encoding-type")
  }
  class(h) <- "h5ad_handle"
  h
}

h5ad_layer <- function(h, layer, genes = NULL, obs = NULL) {
  #' obs x genes matrix for a dense layer; subset by gene names / obs indices
  gi <- if (is.null(genes)) NULL else match(genes, h$genes)
  if (anyNA(gi)) {
    stop("Unknown gene(s): ", paste(genes[is.na(gi)], collapse = ", "))
  }
  mat <- h5ad_read_dense(h$path, paste0("layers/", layer), rows = obs, cols = gi)
  dimnames(mat) <- list(
    rownames(h$obs)[if (is.null(obs)) seq_len(nrow(h$obs)) else obs],
    if (is.null(genes)) h$genes else genes
  )
  mat
}

h5ad_varm <- function(h, key) {
  h5ad_read_dense(h$path, paste0("varm/", key))
}

h5ad_obsm <- function(h, key) {
  h5ad_read_dense(h$path, paste0("obsm/", key))
}

h5ad_x_gene <- function(h, gene) {
  #' Expression of one gene across all obs, from a CSC-encoded X
  if (!identical(h$x_encoding, "csc_matrix")) {
    stop("X must be stored as csc_matrix for per-gene reads")
  }
  j <- match(gene, h$genes)
  if (is.na(j)) {
    stop("Gene not found: ", gene)
  }
  indptr <- as.vector(h5read(h$path, "X/indptr", index = list(c(j, j + 1))))
  out <- numeric(nrow(h$obs))
  if (indptr[2] > indptr[1]) {
    span <- list(seq(indptr[1] + 1, indptr[2]))
    values <- as.vector(h5read(h$path, "X/data", index = span))
    rows <- as.vector(h5read(h$path, "X/indices", index = span))
    out[rows + 1] <- values
  }
  out
}

# ---------------------------------------------------------------------------
# Datasets: the object the app works with
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Cluster IDs and names
# ---------------------------------------------------------------------------
#
# obs/cluster holds a stable cluster ID, the key for every join. The optional
# obs/cluster_name holds an annotation (missing while unnamed). The app shows
# "name (id)", or the bare ID, but always works on IDs.

read_cluster_names <- function(obs, clusters, path = "") {
  #' Named character vector id -> name (NA when unnamed) from obs/cluster_name
  out <- setNames(rep(NA_character_, length(clusters)), clusters)
  if (!"cluster_name" %in% names(obs)) return(out)
  name <- trimws(as.character(obs$cluster_name))
  name[!is.na(name) & name == ""] <- NA
  id <- as.character(obs$cluster)
  per_id <- tapply(name, id, function(v) unique(v[!is.na(v)]), simplify = FALSE)
  conflict <- names(per_id)[lengths(per_id) > 1]
  if (length(conflict) > 0) {
    stop(path, ": obs/cluster_name differs between splits for cluster(s) ",
         paste(head(conflict, 5), collapse = ", "))
  }
  named <- names(per_id)[lengths(per_id) == 1]
  out[named] <- unlist(per_id[named])
  out
}

format_cluster_labels <- function(ids, names) {
  #' "name (id)" for named clusters, the bare ID otherwise
  setNames(ifelse(is.na(names), ids, paste0(names, " (", ids, ")")), ids)
}

cluster_label <- function(ds, ids) {
  #' Display labels of cluster IDs
  unname(ds$cluster_labels[as.character(ids)])
}

cluster_choices <- function(ds, ids = ds$clusters) {
  #' Choices for Shiny inputs: labels shown, IDs returned
  setNames(as.character(ids), cluster_label(ds, ids))
}

.ds_cache <- new.env(parent = emptyenv())

ds_load <- function(id, entry) {
  #' Load (and memoize per process) a dataset described by a config entry
  key <- paste(id, entry$clusters, file.info(entry$clusters)$mtime)
  if (!is.null(.ds_cache[[key]])) {
    return(.ds_cache[[key]])
  }

  ch <- h5ad_open(entry$clusters)
  if (!"prob" %in% ch$layers) {
    stop(entry$clusters, " has no 'prob' layer")
  }
  meta <- ch$uns$scmarco
  if (is.null(meta)) meta <- list()

  splits <- entry$split_order
  if (is.null(splits)) splits <- meta$split_order
  if (is.null(splits)) splits <- levels(ch$obs$split)
  clusters <- levels(ch$obs$cluster)
  if (is.null(clusters)) clusters <- sort(unique(ch$obs$cluster))
  clusters <- as.character(clusters)
  cluster_names <- read_cluster_names(ch$obs, clusters, entry$clusters)

  bimodal <- NULL
  if ("bimodal_by_split" %in% ch$varm) {
    bimodal <- h5ad_varm(ch, "bimodal_by_split") > 0
    dimnames(bimodal) <- list(ch$genes, as.character(meta$bimodal_splits))
  }
  analysis_genes <- ch$genes
  if ("is_bimodal" %in% names(ch$var)) {
    analysis_genes <- ch$genes[ch$var$is_bimodal %in% TRUE]
  }

  cells <- NULL
  if (!is.null(entry$cells)) {
    cells <- h5ad_open(entry$cells)
  }

  pick <- function(x, default) if (is.null(x)) default else x
  ds <- list(
    id = id,
    label = entry$label,
    clusters_h = ch,
    cells_h = cells,
    splits = as.character(splits),
    clusters = clusters,
    cluster_names = cluster_names,
    cluster_labels = format_cluster_labels(clusters, cluster_names),
    obs_split = as.character(ch$obs$split),
    obs_cluster = as.character(ch$obs$cluster),
    prob = h5ad_layer(ch, "prob"),
    bimodal = bimodal,
    all_genes = ch$genes,
    analysis_genes = analysis_genes,
    split_label = pick(entry$split_label, pick(meta$split_label, "split")),
    ordered_splits = isTRUE(pick(entry$ordered_splits, meta$ordered_splits)),
    has_pct = "pct_pos" %in% ch$layers,
    has_mean = "mean_lognorm" %in% ch$layers,
    reductions = if (is.null(cells)) character(0) else cells$obsm,
    binarized = new.env(parent = emptyenv())
  )
  class(ds) <- "scmarco_dataset"
  .ds_cache[[key]] <- ds
  ds
}

ds_binarize <- function(ds, cut_off) {
  #' Binarized view of the dataset at a threshold (memoized, a few per dataset)
  key <- format(cut_off, nsmall = 3)
  if (!is.null(ds$binarized[[key]])) {
    return(ds$binarized[[key]])
  }
  cached <- ls(ds$binarized)
  if (length(cached) >= 3) {
    rm(list = cached[1], envir = ds$binarized)
  }
  b <- binarize(
    ds$prob, ds$obs_split, ds$obs_cluster,
    splits = ds$splits, clusters = ds$clusters,
    cut_off = cut_off, bimodal = ds$bimodal
  )
  ds$binarized[[key]] <- b
  b
}

ds_cluster_layer_long <- function(ds, layer, genes, clusters = NULL) {
  #' Long data.frame(gene, cluster, split, value) of a cluster-level layer
  ch <- ds$clusters_h
  keep <- if (is.null(clusters)) {
    seq_len(nrow(ch$obs))
  } else {
    which(ds$obs_cluster %in% clusters)
  }
  genes <- intersect(genes, ds$all_genes)
  if (length(keep) == 0 || length(genes) == 0) {
    return(data.frame(
      gene = character(0), cluster = character(0),
      split = character(0), value = numeric(0)
    ))
  }
  mat <- h5ad_layer(ch, layer, genes = genes, obs = keep)
  data.frame(
    gene = rep(genes, each = length(keep)),
    cluster = rep(ds$obs_cluster[keep], times = length(genes)),
    split = rep(ds$obs_split[keep], times = length(genes)),
    value = as.vector(mat)
  )
}

ds_cells_obs <- function(ds) {
  if (is.null(ds$cells_h)) stop("Dataset has no cell-level file")
  ds$cells_h$obs
}

ds_cell_expression <- function(ds, gene) {
  if (is.null(ds$cells_h)) stop("Dataset has no cell-level file")
  h5ad_x_gene(ds$cells_h, gene)
}

ds_reduction <- function(ds, key) {
  if (is.null(ds$cells_h)) stop("Dataset has no cell-level file")
  coords <- h5ad_obsm(ds$cells_h, key)
  colnames(coords) <- paste0("dim", seq_len(ncol(coords)))
  coords
}
