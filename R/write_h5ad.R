# Shared writers for the scMarco h5ad schema (format_version 2.0).
# Used by build.R and migrate.R; see docs/expected_input.md for the contract.

# Versions of the contract between fitting, the dataset files and the app.
#
# scMarco schema version: the layout of the objects in the contract, i.e.
# the h5ad files (uns/scmarco/format_version, in both files) and the fit
# objects of fit_mixture() (schema_version). v1 was the legacy SQLite
# database, which migrate.R converts; v2 is still being consolidated, so
# required fields are added to it without a bump until it is released. Bump
# it when the layout changes after that; the app refuses any other version.
SCMARCO_FORMAT_VERSION <- "2.0"
# Version of the binarizer: fit_mixture() + call_mixture(). Bump it when a
# change alters the fits or the probabilities; datasets binarized with an
# older version then load with a warning to rerun the binarizer, and fits of
# a newer version are refused; fits of an older version are read and
# called with a warning, and fits of different versions are never merged. Recorded in fits, in
# call_mixture() tables and in uns/scmarco/binarization/version.
SCMARCO_BINARIZER_VERSION <- "1.1"
# binarization/version of probabilities that did not come from scMarco's
# binarizer: tables of the mixture_modeling_smk pipeline, hand-made tables,
# and migrated legacy databases
BINARIZER_EXTERNAL <- "external"

write_clusters_h5ad <- function(path, obs_split, obs_cluster, splits, clusters,
                                genes, layers, is_bimodal = NULL,
                                bimodal_by_split = NULL, n_cells = NULL,
                                split_label = "split", ordered_splits = FALSE,
                                binarization = NULL, extra_uns = list()) {
    # layers: named list of n_obs x n_genes numeric matrices. Without "prob"
    #   the file is in the *built* state; add_probability_layer() binarizes it.
    # bimodal_by_split: optional n_genes x n_splits logical matrix
    # binarization: list(method, source, date, version), required with a
    #   "prob" layer (see binarization_record())
    if ("prob" %in% names(layers)) {
        check_binarization_record(binarization)
    } else if (!is.null(binarization)) {
        stop("binarization given without a prob layer")
    }
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
    if (!is.null(binarization)) uns$binarization <- binarization
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
    write_h5ad_rechunked(ad, path)
    rechunk_h5ad(
        path, n_obs = nrow(obs),
        dense = paste0("/layers/", names(layers)),
        gene_block = min(256, length(genes))
    )
    invisible(path)
}

add_probability_layer <- function(path, prob, bimodal_by_split, binarization,
                                  mixture = NULL, mixture_models = NULL) {
    check_binarization_record(binarization)
    # Binarize a clusters file in place: write layers/prob, var/is_bimodal,
    # varm/bimodal_by_split, varm/mixture_* and uns/scmarco/binarization.
    # Replaces any earlier probabilities; obs, the other layers and uns are kept.
    # prob: n_obs x n_genes matrix in obs order
    # bimodal_by_split: n_genes x n_splits logical matrix, columns named by split
    # binarization: list(method, source, date, version), see
    #   binarization_record()
    # mixture: optional named list of n_genes x n_splits fit parameters
    #   (on_level, off_level, sd, pi), columns named by split
    # mixture_models: split -> the model it was fit in; splits sharing a model
    #   were fit together. Stored in uns/scmarco/mixture_models (in
    #   bimodal_splits order) alongside the parameters
    ad <- read_h5ad(path, as = "InMemoryAnnData")
    splits <- levels(ad$obs$split)
    genes <- rownames(ad$var)
    stopifnot(
        identical(dim(prob), c(nrow(ad$obs), length(genes))),
        setequal(colnames(bimodal_by_split), splits)
    )
    dimnames(prob) <- list(rownames(ad$obs), genes)
    layers <- ad$layers
    layers$prob <- prob
    ad$layers <- layers

    bim <- bimodal_by_split[genes, splits, drop = FALSE]
    var <- ad$var
    var$is_bimodal <- rowSums(bim) > 0
    ad$var <- var
    # anndataR writes logical matrices as flat vectors; store 0/1 integers
    storage.mode(bim) <- "integer"
    varm <- if (is.null(ad$varm)) list() else ad$varm
    varm$bimodal_by_split <- unname(bim)
    # drop parameters of an earlier run so they never outlive their probabilities
    varm <- varm[!startsWith(names(varm), "mixture_")]
    for (k in names(mixture)) {
        m <- mixture[[k]][genes, splits, drop = FALSE]
        storage.mode(m) <- "double"
        varm[[paste0("mixture_", k)]] <- unname(m)
    }
    ad$varm <- varm

    uns <- ad$uns
    uns$scmarco$bimodal_splits <- splits
    uns$scmarco$binarization <- binarization
    uns$scmarco$mixture_models <- NULL
    if (!is.null(mixture)) {
        if (is.null(mixture_models)) mixture_models <- setNames(splits, splits)
        uns$scmarco$mixture_models <- unname(as.character(mixture_models[splits]))
    }
    ad$uns <- uns

    tmp <- paste0(path, ".binarize")
    on.exit(unlink(tmp))
    write_h5ad_rechunked(ad, tmp)
    rechunk_h5ad(
        tmp, n_obs = nrow(ad$obs),
        dense = paste0("/layers/", names(layers)),
        gene_block = min(256, length(genes))
    )
    file.rename(tmp, path)
    invisible(path)
}

write_cells_h5ad <- function(path, x, barcodes, cell_split, cell_cluster,
                             splits, clusters, genes, reductions = list(),
                             reduction_meta = list()) {
    # x: n_cells x n_genes dgCMatrix (log-normalized). Written as csc_matrix,
    #   which gives cheap per-gene reads.
    # reductions: named list of n_cells x k matrices (NaN outside subsets)
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
    if (length(obsm) > 0) names(obsm) <- paste0("X_", names(reductions))

    ad <- AnnData(
        X = x, obs = obs, var = data.frame(row.names = genes),
        obsm = if (length(obsm) > 0) obsm else NULL,
        uns = list(scmarco = list(
            format_version = SCMARCO_FORMAT_VERSION,
            reduction_meta = reduction_meta
        ))
    )
    write_h5ad_rechunked(ad, path)
    rechunk_h5ad(path, n_obs = nrow(obs), sparse_x = TRUE)
    invisible(path)
}

binarization_record <- function(method, source, version) {
    # uns/scmarco/binarization of a binarized dataset
    list(method = method, source = source, date = format(Sys.Date()),
        version = version)
}

check_binarization_record <- function(b) {
    need <- c("method", "source", "date", "version")
    if (!is.list(b) || !all(need %in% names(b))) {
        stop("binarization must be a list with ", paste(need, collapse = ", "))
    }
    if (!valid_binarizer_version(b$version)) {
        stop("binarization version must be a version number or \"",
            BINARIZER_EXTERNAL, "\", not ", format(b$version))
    }
    invisible(TRUE)
}

valid_binarizer_version <- function(v) {
    length(v) == 1 && is.character(v) && !is.na(v) &&
        (v == BINARIZER_EXTERNAL || grepl("^[0-9]+(\\.[0-9]+)+$", v))
}

binarizer_status <- function(v) {
    # How a dataset's binarization/version relates to this scMarco's:
    # "missing", "invalid", "external", "current", "older" or "newer"
    if (is.null(v) || length(v) == 0) return("missing")
    v <- as.character(v)
    if (!valid_binarizer_version(v)) return("invalid")
    if (v == BINARIZER_EXTERNAL) return("external")
    cmp <- compareVersion(v, SCMARCO_BINARIZER_VERSION)
    c("older", "current", "newer")[cmp + 2]
}

rewrite_obs <- function(path, edit) {
    # Replace only obs of an h5ad file, in place: layers and X are untouched.
    # edit: function(obs) returning the new obs, same rows in the same order.
    ad <- suppressWarnings(read_h5ad(path, as = "HDF5AnnData", mode = "r+"))
    on.exit(ad$close())
    obs <- edit(ad$obs)
    stopifnot(nrow(obs) == ad$n_obs())
    ad$obs <- obs
    invisible(path)
}

split_cluster_labels <- function(labels) {
    # Split legacy "name (id)" labels into ID and name. A label without a
    # trailing "(id)" is an unnamed cluster whose ID is the label itself.
    m <- regmatches(labels, regexec("^(.*\\S)\\s+\\(([^()]+)\\)$", labels))
    parsed <- lengths(m) == 3
    id <- labels
    name <- rep(NA_character_, length(labels))
    id[parsed] <- vapply(m[parsed], `[`, "", 3)
    name[parsed] <- vapply(m[parsed], `[`, "", 2)
    data.frame(label = labels, id = id, name = name, stringsAsFactors = FALSE)
}

order_cluster_ids <- function(ids) {
    # IDs starting with a number first, by that number (so "1+2" follows "1"),
    # then the others alphabetically
    lead <- regmatches(ids, regexpr("^[0-9]+(\\.[0-9]+)?", ids))
    num <- rep(NA_real_, length(ids))
    num[grepl("^[0-9]", ids)] <- as.numeric(lead)
    ids[order(is.na(num), num, ids)]
}

write_h5ad_rechunked <- function(ad, path) {
    # anndataR's write_h5ad() for a file that rechunk_h5ad() re-chunks next.
    # rhdf5 notes, for every large array, that it is a single chunk; that is
    # what re-chunking fixes (rechunk_h5ad() warns when it cannot), so the
    # note is dropped. Other messages pass through.
    withCallingHandlers(
        write_h5ad(ad, path, mode = "w"),
        message = function(m) {
            if (startsWith(conditionMessage(m), "You created a large dataset with compression")) {
                invokeRestart("muffleMessage")
            }
        }
    )
}

rechunk_h5ad <- function(path, n_obs, sparse_x = FALSE, dense = character(0),
                         gene_block = 256) {
    # anndataR writes each array as a single chunk, so reading one gene would
    # decompress everything. Re-chunk with h5repack (HDF5 command line tools):
    # sparse X arrays in 64k blocks, dense obs x gene layers in gene blocks.
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
