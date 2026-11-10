# h5ad backend
#
# A thin reader for the AnnData on-disk format (https://anndata.readthedocs.io
# /en/latest/fileformat-prose.html) built on rhdf5, so the app can read single
# genes without loading whole matrices and needs no Python at runtime.
#
# Each dataset is one required <ds>.clusters.h5ad (obs = split x cluster pairs
# that exist; layers prob / pct_pos / mean_lognorm) and one optional
# <ds>.cells.h5ad (X = log-normalized expression in CSC layout; obs = cells;
# obsm = reductions). See the package README for the full schema.

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
    # obs x genes matrix for a dense layer; subset by gene names / obs indices
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
    # Expression of one gene across all obs, from a CSC-encoded X
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

read_cluster_column <- function(obs, col, clusters, path = "") {
    # Named character vector id -> value (NA when missing) of a per-cluster obs
    # column, which must be the same in every split of a cluster
    out <- setNames(rep(NA_character_, length(clusters)), clusters)
    if (!col %in% names(obs)) return(out)
    value <- trimws(as.character(obs[[col]]))
    value[!is.na(value) & value == ""] <- NA
    id <- as.character(obs$cluster)
    per_id <- tapply(value, id, function(v) unique(v[!is.na(v)]), simplify = FALSE)
    conflict <- names(per_id)[lengths(per_id) > 1]
    if (length(conflict) > 0) {
        stop(path, ": obs/", col, " differs between splits for cluster(s) ",
            paste(head(conflict, 5), collapse = ", "))
    }
    known <- intersect(names(per_id)[lengths(per_id) == 1], clusters)
    out[known] <- unlist(per_id[known])
    out
}

read_cluster_names <- function(obs, clusters, path = "") {
    # Named character vector id -> name (NA when unnamed) from obs/cluster_name
    read_cluster_column(obs, "cluster_name", clusters, path)
}

format_cluster_labels <- function(ids, names) {
    # "name (id)" for named clusters, the bare ID otherwise
    setNames(ifelse(is.na(names), ids, paste0(names, " (", ids, ")")), ids)
}

cluster_label <- function(ds, ids) {
    # Display labels of cluster IDs
    unname(ds$cluster_labels[as.character(ids)])
}

cluster_choices <- function(ds, ids = ds$clusters) {
    # Choices for Shiny inputs: labels shown, IDs returned
    setNames(as.character(ids), cluster_label(ds, ids))
}

# ---------------------------------------------------------------------------
# Cluster metadata
# ---------------------------------------------------------------------------
#
# Optional obs columns meta_<level>, one value per cluster (the same in every
# split; missing = unassigned). The config's subset_by lists the levels the app
# offers as filters, coarse to fine (e.g. class, subclass): a cascade that
# limits every search to a subset of clusters.

META_PREFIX <- "meta_"

read_cluster_meta <- function(obs, clusters, path = "", levels = NULL) {
    # data.frame, one row per cluster ID and one column per level (`levels`
    # order, default every level in obs order), or NULL when there is none
    cols <- names(obs)[startsWith(names(obs), META_PREFIX)]
    if (!is.null(levels)) cols <- paste0(META_PREFIX, unlist(levels))
    if (length(cols) == 0) return(NULL)
    out <- lapply(cols, function(col) unname(read_cluster_column(obs, col, clusters, path)))
    out <- as.data.frame(setNames(out, substring(cols, nchar(META_PREFIX) + 1)),
        optional = TRUE, stringsAsFactors = FALSE)
    rownames(out) <- clusters
    out
}

meta_match <- function(ds, filters, levels = names(ds$cluster_meta)) {
    # Logical per cluster: matches every non-empty filter among `levels`
    keep <- rep(TRUE, length(ds$clusters))
    for (level in intersect(levels, names(filters))) {
        chosen <- filters[[level]]
        if (length(chosen) == 0) next
        keep <- keep & ds$cluster_meta[[level]] %in% chosen
    }
    keep
}

meta_universe <- function(ds, filters) {
    # Cluster IDs inside the search scope: those matching every level with a
    # selection (unassigned never matches). All clusters without metadata or
    # filters.
    if (is.null(ds$cluster_meta)) return(ds$clusters)
    ds$clusters[meta_match(ds, filters)]
}

meta_level_choices <- function(ds, filters, level) {
    # Values of `level` among the clusters that pass the levels above it
    levels <- names(ds$cluster_meta)
    above <- levels[seq_len(match(level, levels) - 1)]
    values <- ds$cluster_meta[[level]][meta_match(ds, filters, above)]
    sort(unique(values[!is.na(values)]))
}

meta_active <- function(filters) {
    # The filters that restrict anything
    filters[lengths(filters) > 0]
}

.ds_cache <- new.env(parent = emptyenv())

# Expression views a config may ask for: cluster averages per split (line,
# needs mean_lognorm; drawn as lines only for ordered splits, as points
# otherwise), per-cell distributions (box, needs the cells file), or both in
# one panel
TREND_MODES <- c("line", "box", "both")

trend_mode <- function(entry, ordered_splits, has_mean, has_cells) {
    # The expression view the config asks for; without one, ordered splits
    # with cluster averages get the line plot and any other dataset with a
    # cells file the box plot. NULL: no expression view.
    if (!is.null(entry$trend)) return(entry$trend)
    if (ordered_splits && has_mean) return("line")
    if (has_cells) "box"
}

plural_label <- function(label) {
    # Default plural of a split label: "stage" -> "stages", "sex" -> "sexes";
    # a label already ending in "s" ("species") is kept as is
    if (grepl("s$", label)) return(label)
    if (grepl("(x|z|ch|sh)$", label)) return(paste0(label, "es"))
    paste0(label, "s")
}

format_problem <- function(path) {
    # NULL when the file has the current schema, else what is wrong and the fix
    obj <- "uns/scmarco/format_version"
    found <- if (h5ad_exists(path, obj)) as.character(h5read(path, obj)) else NA
    if (identical(found, SCMARCO_FORMAT_VERSION)) return(NULL)
    if (is.na(found)) {
        return(sprintf("%s has no format_version; it is not a scMarco dataset", path))
    }
    if (package_version(found) > package_version(SCMARCO_FORMAT_VERSION)) {
        return(sprintf("%s has format_version %s, newer than this scMarco reads (%s); update scMarco",
            path, found, SCMARCO_FORMAT_VERSION))
    }
    sprintf("%s has format_version %s, but this scMarco reads %s; rebuild the dataset",
        path, found, SCMARCO_FORMAT_VERSION)
}

binarizer_problem <- function(path) {
    # NULL when the dataset's binarizer version can be used, else what is
    # wrong and the fix. An older binarizer is usable: ds_check_config() warns.
    obj <- "uns/scmarco/binarization/version"
    v <- if (h5ad_exists(path, obj)) as.character(h5read(path, obj))
    switch(binarizer_status(v),
        missing = sprintf(paste0("%s does not record its binarizer version ",
            "(uns/scmarco/binarization/version); rerun add_probability(), or ",
            "rebuild the dataset"), path),
        invalid = sprintf("%s has an invalid binarizer version: %s", path, v),
        newer = sprintf("%s was binarized by binarizer %s, newer than this scMarco's (%s); update scMarco",
            path, v, SCMARCO_BINARIZER_VERSION),
        NULL)
}

#' Check the datasets of an app config
#'
#' Checks every configured dataset before any is loaded, so that a bad config
#' fails at once with one message listing each problem and its fix: missing
#' files, files in another `format_version`, datasets without probabilities
#' or without a usable binarizer version (one older than this scMarco's
#' binarizer only warns),
#' unknown `subset_by` levels, a `split_order` that does not match the file,
#' and `trend` views the data cannot support.
#'
#' @param datasets The `datasets` field of an app config: a named list of
#'   entries as described in [ds_load()].
#'
#' @return `TRUE`, invisibly; otherwise an error.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds_check_config(list(example = list(label = "Example", clusters = ex)))
#' @export
ds_check_config <- function(datasets) {
    problems <- character(0)
    for (id in names(datasets)) {
        path <- datasets[[id]]$clusters
        if (is.null(path) || !file.exists(path)) {
            problems <- c(problems, sprintf("- %s: clusters file not found: %s", id,
                if (is.null(path)) "(not set)" else path))
        } else if (!is.null(msg <- format_problem(path))) {
            problems <- c(problems, sprintf("- %s: %s", id, msg))
        } else if (!h5ad_exists(path, "layers/prob")) {
            problems <- c(problems, sprintf(
                "- %s: %s has no probabilities yet. Add them with\n    add_probability(\"%s\", c(<split> = \"<file.rds>\", ...))  # library(scMarco)",
                id, path, path
            ))
        } else if (!is.null(msg <- binarizer_problem(path))) {
            problems <- c(problems, sprintf("- %s: %s", id, msg))
        } else {
            v <- as.character(h5read(path, "uns/scmarco/binarization/version"))
            if (binarizer_status(v) == "older") {
                warning(sprintf("%s: binarized with binarizer %s; this scMarco's is %s. Consider refitting it.",
                    id, v, SCMARCO_BINARIZER_VERSION), call. = FALSE)
            }
            cells <- datasets[[id]]$cells
            if (!is.null(cells) && file.exists(cells) &&
                !is.null(msg <- format_problem(cells))) {
                problems <- c(problems, sprintf("- %s: %s", id, msg))
            }
            levels <- unlist(datasets[[id]]$subset_by)
            missing <- levels[!vapply(paste0("obs/", META_PREFIX, levels),
                function(obj) h5ad_exists(path, obj), logical(1))]
            if (length(missing) > 0) {
                problems <- c(problems, sprintf(
                    "- %s: subset_by lists level(s) %s has no metadata for: %s. Add them with set_cluster_metadata() or fix subset_by",
                    id, path, paste(missing, collapse = ", ")
                ))
            }
            if (anyDuplicated(levels)) {
                problems <- c(problems, sprintf("- %s: subset_by lists a level more than once", id))
            }
            order <- as.character(unlist(datasets[[id]]$split_order))
            if (length(order) > 0) {
                splits <- as.character(h5read(path, "obs/split/categories"))
                if (anyDuplicated(order) || !setequal(order, splits)) {
                    problems <- c(problems, sprintf(
                        "- %s: split_order must list each split of %s exactly once: %s",
                        id, path, paste(splits, collapse = ", ")
                    ))
                }
            }
            trend <- datasets[[id]]$trend
            if (!is.null(trend)) {
                if (!(length(trend) == 1 && trend %in% TREND_MODES)) {
                    problems <- c(problems, sprintf(
                        "- %s: trend must be one of %s", id, paste(TREND_MODES, collapse = ", ")
                    ))
                } else {
                    if (trend %in% c("line", "both") && !h5ad_exists(path, "layers/mean_lognorm")) {
                        problems <- c(problems, sprintf(
                            "- %s: trend: %s needs layers/mean_lognorm, which %s does not have; use trend: box",
                            id, trend, path
                        ))
                    }
                    cells <- datasets[[id]]$cells
                    if (trend %in% c("box", "both") && (is.null(cells) || !file.exists(cells))) {
                        problems <- c(problems, sprintf(
                            "- %s: trend: %s needs a cells file, but cells is %s",
                            id, trend, if (is.null(cells)) "not set" else paste("not found:", cells)
                        ))
                    }
                }
            }
        }
    }
    if (length(problems) > 0) {
        stop("Some configured datasets cannot be loaded:\n",
            paste(problems, collapse = "\n"), call. = FALSE)
    }
    invisible(TRUE)
}

#' Load an scMarco dataset
#'
#' Opens the h5ad files of a dataset and reads what every search needs (the
#' probability layer, bimodality, cluster names and metadata). Cell-level
#' data, cluster averages and reductions are read from disk on demand. The
#' result is cached per R process, keyed by `id` and the file's path and
#' modification time.
#'
#' @param id Dataset ID.
#' @param entry A list with the fields of a dataset in the app config:
#'   `label`, `clusters` (path of the `.clusters.h5ad` file, required),
#'   `cells` (path of the `.cells.h5ad` file, optional) and the optional
#'   overrides `split_order`, `split_label`, `split_label_plural`,
#'   `ordered_splits`, `trend` (`"line"`, `"box"` or `"both"`) and
#'   `subset_by` (cluster metadata levels offered as filters).
#'
#' @return An object of class `scmarco_dataset`: a list with, among others,
#'   `splits`, `clusters`, `cluster_names`, `cluster_labels`, `all_genes`,
#'   `analysis_genes` (genes used for marker finding), `split_label`,
#'   `reductions`, `binarization` (method, source, date and binarizer
#'   `version` of the probabilities) and the flags `has_pct`, `has_mean` and
#'   `has_mixture`.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, margin = 0.1)
#' ds$splits
#' head(ds$cluster_labels)
#' @export
ds_load <- function(id, entry) {
    key <- paste(id, entry$clusters, file.info(entry$clusters)$mtime)
    if (!is.null(.ds_cache[[key]])) {
        return(.ds_cache[[key]])
    }

    ch <- h5ad_open(entry$clusters)
    if (!"prob" %in% ch$layers) {
        stop(entry$clusters, " has no probabilities; add them with ",
            "scMarco::add_probability()", call. = FALSE)
    }
    meta <- ch$uns$scmarco
    if (is.null(meta)) meta <- list()
    for (msg in c(format_problem(entry$clusters), binarizer_problem(entry$clusters))) {
        stop(msg, call. = FALSE)
    }

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
    # Mixture-model fit per gene x split: component means, shared SD, ON weight
    mixture <- NULL
    mixture_keys <- paste0("mixture_", c("on_level", "off_level", "sd", "pi"))
    if (all(mixture_keys %in% ch$varm)) {
        mixture <- lapply(mixture_keys, function(key) {
            m <- h5ad_varm(ch, key)
            dimnames(m) <- list(ch$genes, as.character(meta$bimodal_splits))
            m
        })
        names(mixture) <- sub("^mixture_", "", mixture_keys)
    }
    # split -> the model it was fit in; splits sharing one were fit together
    mixture_models <- NULL
    if (!is.null(mixture)) {
        fit_splits <- as.character(meta$bimodal_splits)
        models <- if (is.null(meta$mixture_models)) fit_splits else
            as.character(meta$mixture_models)
        mixture_models <- setNames(models, fit_splits)
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
    ordered_splits <- isTRUE(pick(entry$ordered_splits, meta$ordered_splits))
    has_mean <- "mean_lognorm" %in% ch$layers
    split_label <- pick(entry$split_label, pick(meta$split_label, "split"))
    ds <- list(
        id = id,
        label = entry$label,
        clusters_h = ch,
        cells_h = cells,
        splits = as.character(splits),
        clusters = clusters,
        cluster_names = cluster_names,
        cluster_labels = format_cluster_labels(clusters, cluster_names),
        # Only the levels the config offers as filters; none without subset_by
        cluster_meta = if (length(entry$subset_by) > 0)
            read_cluster_meta(ch$obs, clusters, entry$clusters, entry$subset_by),
        obs_split = as.character(ch$obs$split),
        obs_cluster = as.character(ch$obs$cluster),
        prob = h5ad_layer(ch, "prob"),
        bimodal = bimodal,
        mixture = mixture,
        mixture_models = mixture_models,
        all_genes = ch$genes,
        analysis_genes = analysis_genes,
        unimodal_genes = setdiff(ch$genes, analysis_genes),
        split_label = split_label,
        split_plural = pick(entry$split_label_plural,
            pick(meta$split_label_plural, plural_label(split_label))),
        ordered_splits = ordered_splits,
        has_pct = "pct_pos" %in% ch$layers,
        has_mean = has_mean,
        trend = trend_mode(entry, ordered_splits, has_mean, !is.null(cells)),
        has_mixture = !is.null(mixture),
        binarization = meta$binarization,
        reductions = if (is.null(cells)) character(0) else cells$obsm,
        binarized = new.env(parent = emptyenv())
    )
    class(ds) <- "scmarco_dataset"
    .ds_cache[[key]] <- ds
    ds
}

#' Binarize a dataset at a confidence margin
#'
#' [binarize()] applied to a dataset. The two most recent margins are cached
#' with the dataset.
#'
#' @param ds A dataset from [ds_load()].
#' @param margin Confidence margin in (0, 0.5]: ON when the probability is at
#'   least `1 - margin`, OFF when it is at most `margin`.
#'
#' @return A binarized dataset (see [binarize()]).
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' ds <- ds_load("example", list(label = "Example", clusters = ex))
#' b <- ds_binarize(ds, margin = 0.1)
#' dim(b$on)
#' @export
ds_binarize <- function(ds, margin) {
    key <- format(margin, nsmall = 3)
    if (!is.null(ds$binarized[[key]])) {
        return(ds$binarized[[key]])
    }
    # Evict the oldest; .order is hidden from ls() and keeps insertion order
    order <- ds$binarized$.order
    if (length(order) >= 2) {
        rm(list = order[1], envir = ds$binarized)
        order <- order[-1]
    }
    b <- binarize(
        ds$prob, ds$obs_split, ds$obs_cluster,
        splits = ds$splits, clusters = ds$clusters,
        margin = margin, bimodal = ds$bimodal
    )
    ds$binarized[[key]] <- b
    ds$binarized$.order <- c(order, key)
    b
}

gene_choices <- function(ds, genes) {
    # Choices for Shiny inputs: "gene (unimodal)" shown for genes the mixture
    # model never called bimodal (kept out of marker discovery), the gene
    # returned
    label <- ifelse(genes %in% ds$unimodal_genes, paste0(genes, " (unimodal)"), genes)
    setNames(genes, label)
}

ds_cluster_layer_long <- function(ds, layer, genes, clusters = NULL) {
    # Long data.frame(gene, cluster, split, value) of a cluster-level layer
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

mixture_model_splits <- function(ds, model) {
    # Splits fit in one mixture model, in dataset order
    intersect(ds$splits, names(ds$mixture_models)[ds$mixture_models == model])
}

#' Mixture models of a dataset
#'
#' Splits are fit either on their own or together with other splits; each
#' fit is one mixture model.
#'
#' @param ds A dataset from [ds_load()].
#'
#' @return Named character vector: model names, named by the splits they
#'   cover (e.g. `"Mel + Yak"`). Empty when the dataset stores no fit
#'   parameters.
#' @examples
#' toy <- system.file("extdata", package = "scMarco")
#' ds <- ds_load("toy", list(
#'     label = "Toy", clusters = file.path(toy, "toy.clusters.h5ad"),
#'     cells = file.path(toy, "toy.cells.h5ad")
#' ))
#' mixture_model_choices(ds)
#' @export
mixture_model_choices <- function(ds) {
    if (is.null(ds$mixture_models)) return(character(0))
    models <- unique(ds$mixture_models[intersect(ds$splits, names(ds$mixture_models))])
    labels <- vapply(models, function(m) {
        paste(mixture_model_splits(ds, m), collapse = " + ")
    }, "")
    setNames(models, labels)
}

ds_mixture_fit <- function(ds, gene, model) {
    # Mixture-model fit of one gene in one model: list(on_level, off_level, sd,
    # pi, bimodal), or NULL when the dataset has no fit parameters. Splits fit
    # together share their parameters; bimodal is TRUE if it is in any of them.
    if (is.null(ds$mixture)) return(NULL)
    splits <- mixture_model_splits(ds, model)
    fit <- lapply(ds$mixture, function(m) m[gene, splits[1]])
    fit$bimodal <- any(ds$bimodal[gene, splits] %in% TRUE)
    fit
}

ds_cluster_average <- function(ds, gene, splits) {
    # data.frame(split, cluster, value): average expression of a gene in each
    # split x cluster pair of the given splits. With a cell-level file this is
    # log1p(mean(expm1(x))), the scale the mixture model was fit on (averaging
    # in linear space, as Seurat's AverageExpression); otherwise the mean of
    # log-normalized values. attr(, "scale") says which.
    if (!is.null(ds$cells_h)) {
        obs <- ds$cells_h$obs
        sel <- as.character(obs$split) %in% splits
        x <- ds_cell_expression(ds, gene)[sel]
        key <- data.frame(split = as.character(obs$split[sel]),
            cluster = as.character(obs$cluster[sel]))
        out <- aggregate(list(value = expm1(x)), key, mean)
        out$value <- log1p(out$value)
        attr(out, "scale") <- "log of mean"
        return(out)
    }
    if (!ds$has_mean) return(NULL)
    long <- ds_cluster_layer_long(ds, "mean_lognorm", gene)
    out <- long[long$split %in% splits, c("split", "cluster", "value")]
    attr(out, "scale") <- "mean of log"
    out
}
