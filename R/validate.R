# Validation of a scMarco dataset (format_version 2.0), used by
# add_probability() and by the tests.

#' Check an scMarco dataset
#'
#' Checks a clusters file, and optionally its cells file, against the scMarco
#' h5ad format (format_version 2.0): obs layout and labels, cluster names and
#' metadata, layer shapes and ranges, the binarization record and its
#' binarizer version, bimodality and mixture-model fields, and the agreement
#' of the cells file with the clusters file.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param cells Optional path of the `.cells.h5ad` file (`""` for none).
#' @param stage `"binarized"` (what the app needs) or `"built"`, which
#'   accepts a dataset without probabilities.
#'
#' @return A list with the character vectors `errors` and `warnings` and the
#'   flag `binarized`; print it with [report_validation()].
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' res <- validate_h5ad(ex)
#' report_validation(res, ex)
#' @export
validate_h5ad <- function(clusters, cells = "", stage = "binarized") {
    if (!stage %in% c("built", "binarized")) {
        stop("stage must be built or binarized, not ", stage)
    }
    errors <- character(0)
    warnings <- character(0)
    err <- function(...) errors <<- c(errors, paste0(...))
    warn <- function(...) warnings <<- c(warnings, paste0(...))

    ch <- h5ad_open(clusters)
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
                "migrate_cluster_labels()")
        }
    }
    for (col in names(ch$obs)[startsWith(names(ch$obs), META_PREFIX)]) {
        if (!is.factor(ch$obs[[col]])) err("obs/", col, " must be categorical")
        values <- tryCatch(
            read_cluster_column(ch$obs, col, levels(ch$obs$cluster)),
            error = function(e) {
                err("obs/", col, " must be the same in every split of a cluster")
                NULL
            }
        )
        lv <- levels(ch$obs[[col]])
        if (any(lv != trimws(lv))) err("obs/", col, " values start or end with whitespace")
        if (!is.null(values) && all(is.na(values))) warn("obs/", col, " is empty for every cluster")
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
    if (!binarized && stage == "binarized") {
        err("layers/prob is missing; add it with add_probability() ",
            "(or validate with stage = \"built\")")
    }
    if (binarized) {
        b <- meta$binarization
        missing <- setdiff(c("method", "source", "date", "version"), names(b))
        if (is.null(b)) {
            err("uns/scmarco/binarization is missing; rerun add_probability()")
        } else if (length(missing) > 0) {
            err("uns/scmarco/binarization lacks ", paste(missing, collapse = ", "))
        }
        switch(binarizer_status(b$version),
            invalid = err("uns/scmarco/binarization/version is not a version number or \"",
                BINARIZER_EXTERNAL, "\": ", b$version),
            newer = err("binarizer version ", b$version, " is newer than this scMarco's (",
                SCMARCO_BINARIZER_VERSION, "); update scMarco"),
            older = warn("binarized with binarizer ", b$version, "; this scMarco's is ",
                SCMARCO_BINARIZER_VERSION, ": consider refitting"),
            NULL)
    } else if (!is.null(meta$binarization)) {
        err("uns/scmarco/binarization is present but layers/prob is missing")
    }
    for (layer in intersect(c("prob", "pct_pos", "mean_lognorm", "call_flags", "fold_ambient"),
                            ch$layers)) {
        m <- h5ad_layer(ch, layer)
        if (!identical(dim(m), c(nrow(ch$obs), length(ch$genes)))) {
            err("layers/", layer, " has shape ", paste(dim(m), collapse = "x"))
            next
        }
        if (layer %in% c("prob", "pct_pos")) {
            bad <- sum(!is.na(m) & (m < 0 | m > 1))
            if (bad > 0) err("layers/", layer, ": ", bad, " values outside [0, 1]")
        }
        if (layer == "call_flags") {
            bad <- sum(is.na(m) | m != round(m) | m < 0 | m > 7)
            if (bad > 0) err("layers/call_flags: ", bad, " values that are not flags (0-7)")
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
    for (key in grep("^(mixture|call)_", ch$varm, value = TRUE)) {
        m <- h5ad_varm(ch, key)
        if (nrow(m) != length(ch$genes) || ncol(m) != length(meta$bimodal_splits)) {
            err("varm/", key, " shape does not match genes x uns bimodal_splits")
        }
    }
    if (!is.null(meta$mixture_models) &&
        length(meta$mixture_models) != length(meta$bimodal_splits)) {
        err("uns mixture_models does not match uns bimodal_splits")
    }
    # binarizer 1.1 data (call_mixture() tables)
    for (layer in intersect(c("call_flags", "fold_ambient"), ch$layers)) {
        if (!binarized) err("layers/", layer, " is present but layers/prob is missing")
    }
    if ("fold_ambient" %in% ch$layers && !"call_flags" %in% ch$layers) {
        err("layers/fold_ambient is present without layers/call_flags")
    }
    calls <- meta$binarization$calls
    if ("call_flags" %in% ch$layers && is.null(calls)) {
        err("layers/call_flags is present but uns/scmarco/binarization/calls is missing")
    }
    if (!is.null(calls)) {
        n_split <- length(meta$bimodal_splits)
        for (k in c(CALL_PARAMS, "rho_source")) {
            if (length(calls[[k]]) != n_split) {
                err("uns/scmarco/binarization/calls/", k, " does not have one value per split ",
                    "of uns bimodal_splits")
            }
        }
        n_rho <- length(calls$rho)
        if (length(calls$rho_split) != n_rho || length(calls$rho_library) != n_rho) {
            err("uns/scmarco/binarization/calls: rho, rho_split and rho_library differ in length")
        } else if (!all(as.character(calls$rho_split) %in% as.character(meta$bimodal_splits))) {
            err("uns/scmarco/binarization/calls/rho_split names splits not in uns bimodal_splits")
        }
    }
    if (binarized && binarizer_status(meta$binarization$version) == "current" &&
        !"call_flags" %in% ch$layers) {
        warn("layers/call_flags is missing: calls are made without the binarizer ",
            SCMARCO_BINARIZER_VERSION, " checks (chains, detection floor, ambient RNA); ",
            "rerun add_probability() with tables from call_mixture()")
    }
    if (binarized && !"is_bimodal" %in% names(ch$var)) {
        warn("var/is_bimodal is missing: every gene will be used for marker finding")
    }

    if (nzchar(cells)) {
        xh <- h5ad_open(cells)
        cv <- xh$uns$scmarco$format_version
        if (!identical(as.character(cv), SCMARCO_FORMAT_VERSION)) {
            err("cells format_version is ", if (is.null(cv)) "missing" else cv,
                ", expected ", SCMARCO_FORMAT_VERSION)
        }
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

    list(errors = errors, warnings = warnings, binarized = binarized)
}

#' Print a dataset check
#'
#' @param res A result of [validate_h5ad()].
#' @param clusters,cells The paths that were checked, for the message.
#'
#' @return `TRUE` when there are no errors, invisibly.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' report_validation(validate_h5ad(ex), ex)
#' @export
report_validation <- function(res, clusters, cells = "") {
    for (w in res$warnings) message("WARNING: ", w)
    for (e in res$errors) message("ERROR:   ", e)
    if (length(res$errors) > 0) return(invisible(FALSE))
    message("OK (", if (res$binarized) "binarized" else "built", "): ", clusters,
        if (nzchar(cells)) c(" + ", cells))
    invisible(TRUE)
}
