# Migrate the two supported legacy sources into the scMarco h5ad schema:
# a legacy scMarco SQLite table (sqlite_to_h5ad()) and a scMarco-Ding-lab
# HDF5 "v1.0" database (h5v1_to_h5ad()). Both write *binarized* datasets.

#' Convert a legacy scMarco SQLite table
#'
#' Converts one table of a legacy scMarco SQLite database, a long table with
#' columns `gene`, `cluster`, `stage`, `value` and `type` (`"prob"` or
#' `"lognorm"`), into `<out_dir>/<table>.clusters.h5ad` with the layers
#' `prob` and `mean_lognorm`. The table is streamed in chunks. Legacy labels
#' such as `"Mi4 (154)"` are written as they are; split them into ID and
#' name afterwards with [migrate_cluster_labels()]. Requires the DBI and
#' RSQLite packages.
#'
#' @param db Path of the SQLite database.
#' @param table Table name.
#' @param stages Every stage of the table, in display order.
#' @param out_dir Output directory.
#' @param split_label What a split is, shown in the app.
#' @param ordered The stages have a natural order.
#' @param n_genes Keep only the first `n_genes` genes (alphabetical), e.g. for
#'   test fixtures.
#' @param out Path to write instead of `<out_dir>/<table>.clusters.h5ad`.
#'
#' @return The path written, invisibly.
#' @examples
#' if (requireNamespace("RSQLite", quietly = TRUE)) {
#'     db <- tempfile(fileext = ".sqlite")
#'     con <- DBI::dbConnect(RSQLite::SQLite(), db)
#'     long <- expand.grid(gene = c("geneA", "geneB"), cluster = c("1", "2"),
#'         stage = c("P15", "Adult"), type = c("prob", "lognorm"),
#'         stringsAsFactors = FALSE)
#'     long$value <- runif(nrow(long))
#'     DBI::dbWriteTable(con, "toy", long)
#'     DBI::dbDisconnect(con)
#'     path <- sqlite_to_h5ad(db, "toy", stages = c("P15", "Adult"),
#'         out_dir = tempdir())
#'     validate_h5ad(path)$errors
#' }
#' @export
sqlite_to_h5ad <- function(db, table, stages, out_dir = "data",
                           split_label = "stage", ordered = TRUE,
                           n_genes = NULL, out = NULL) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    tbl <- DBI::dbQuoteIdentifier(con, table)

    message("Collecting genes and clusters...")
    genes <- sort(DBI::dbGetQuery(con, paste("SELECT DISTINCT gene FROM", tbl))$gene)
    if (!is.null(n_genes)) genes <- head(genes, n_genes)
    clusters <- sort(DBI::dbGetQuery(con, paste("SELECT DISTINCT cluster FROM", tbl))$cluster)
    db_stages <- DBI::dbGetQuery(con, paste("SELECT DISTINCT stage FROM", tbl))$stage
    if (!setequal(db_stages, stages)) {
        fail("stages (", paste(stages, collapse = ","), ") does not match the table (",
            paste(sort(db_stages), collapse = ","), ")")
    }
    message(sprintf("%d genes x %d clusters x %d stages",
        length(genes), length(clusters), length(stages)))

    dims <- c(length(genes), length(clusters), length(stages))
    prob <- array(NA_real_, dims)
    lognorm <- array(NA_real_, dims)

    res <- DBI::dbSendQuery(con, paste("SELECT gene, cluster, stage, value, type FROM", tbl))
    n_rows <- 0
    while (!DBI::dbHasCompleted(res)) {
        chunk <- DBI::dbFetch(res, n = 2e6)
        gi <- match(chunk$gene, genes)
        keep <- !is.na(gi)
        idx <- cbind(
            gi[keep],
            match(chunk$cluster[keep], clusters),
            match(chunk$stage[keep], stages)
        )
        is_prob <- chunk$type[keep] == "prob"
        prob[idx[is_prob, , drop = FALSE]] <- chunk$value[keep][is_prob]
        lognorm[idx[!is_prob, , drop = FALSE]] <- chunk$value[keep][!is_prob]
        n_rows <- n_rows + nrow(chunk)
        message(sprintf("  %.1fM rows", n_rows / 1e6))
    }
    DBI::dbClearResult(res)

    obs_split <- character(0)
    obs_cluster <- character(0)
    prob_rows <- list()
    mean_rows <- list()
    for (s in seq_along(stages)) {
        present <- colSums(!is.na(prob[, , s, drop = FALSE]), dims = 1) > 0
        present <- as.vector(present)
        for (ci in which(present)) {
            obs_split <- c(obs_split, stages[s])
            obs_cluster <- c(obs_cluster, clusters[ci])
            prob_rows[[length(prob_rows) + 1]] <- prob[, ci, s]
            mean_rows[[length(mean_rows) + 1]] <- lognorm[, ci, s]
        }
        message(sprintf("  %s: %d clusters present", stages[s], sum(present)))
    }

    layers <- list(prob = do.call(rbind, prob_rows))
    if (any(!is.na(lognorm))) layers$mean_lognorm <- do.call(rbind, mean_rows)
    rm(prob, lognorm)

    if (is.null(out)) {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        out <- file.path(out_dir, paste0(table, ".clusters.h5ad"))
    }
    write_clusters_h5ad(
        out, obs_split, obs_cluster, stages, clusters, genes, layers,
        split_label = split_label, ordered_splits = ordered,
        binarization = binarization_record("legacy scMarco mixture model",
            paste0("sqlite:", basename(db), ":", table), BINARIZER_EXTERNAL),
        extra_uns = list(
            source = paste0("sqlite:", basename(db), ":", table),
            mean_lognorm_method = "legacy scMarco cluster-average log-normalized value"
        )
    )
    message("Wrote ", out)
    invisible(out)
}

#' Convert a legacy scMarco-Ding-lab HDF5 dataset
#'
#' Converts one dataset of a scMarco-Ding-lab HDF5 "v1.0" database (`db.h5`)
#' into `<out_dir>/<dataset>.clusters.h5ad` (layers `prob`, `pct_pos` and
#' `mean_lognorm`) and `<out_dir>/<dataset>.cells.h5ad` (log-normalized
#' expression, cluster labels and UMAPs).
#'
#' @param h5 Path of the `db.h5` file.
#' @param dataset Name of the dataset (top-level group) in `h5`.
#' @param out_dir Output directory.
#' @param split_label What a split is, shown in the app.
#' @param ordered The splits have a natural order.
#'
#' @return The paths of the two files (`clusters`, `cells`), invisibly.
#' @examples
#' # A minimal db.h5: one split, two clusters, two genes
#' h5 <- tempfile(fileext = ".h5")
#' rhdf5::h5createFile(h5)
#' for (g in c("toy", "toy/probability", "toy/barcodes", "toy/clusters",
#'     "toy/lognorm", "toy/lognorm/adult")) {
#'     rhdf5::h5createGroup(h5, g)
#' }
#' w <- function(obj, x) rhdf5::h5write(x, h5, paste0("toy/", obj))
#' w("gene_names", c("geneA", "geneB"))
#' w("cluster_names", c("a", "b"))
#' w("splits", "adult")
#' w("bimodal_genes", c(0L, 1L)) # 0-based
#' w("probability/adult", cbind(a = c(0.9, 0.1), b = c(0.1, 0.9)))
#' w("barcodes/adult", paste0("cell", 1:4))
#' w("clusters/adult", c("a", "a", "b", "b"))
#' # log-normalized expression, CSR genes x cells
#' w("lognorm/adult/indptr", c(0L, 2L, 4L))
#' w("lognorm/adult/indices", c(0L, 1L, 2L, 3L))
#' w("lognorm/adult/data", c(1, 2, 1.5, 2.5))
#' w("lognorm/adult/shape", c(2, 4))
#' out <- h5v1_to_h5ad(h5, "toy", out_dir = tempdir())
#' validate_h5ad(out$clusters, out$cells)$errors
#' @export
h5v1_to_h5ad <- function(h5, dataset, out_dir = "data", split_label = "split",
                         ordered = FALSE) {
    ds <- dataset
    rd <- function(obj, ...) rhdf5::h5read(h5, paste0(ds, "/", obj), ...)

    genes <- as.character(rd("gene_names"))
    clusters <- as.character(rd("cluster_names"))
    splits <- as.character(rd("splits"))
    bimodal_idx <- as.integer(rd("bimodal_genes")) + 1L
    message(sprintf("%s: %d genes, %d clusters, splits: %s",
        ds, length(genes), length(clusters), paste(splits, collapse = ", ")))

    listing <- rhdf5::h5ls(h5)
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
        x <- Matrix::sparseMatrix(
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

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    clusters_path <- file.path(out_dir, paste0(ds, ".clusters.h5ad"))
    write_clusters_h5ad(
        clusters_path, obs_split, obs_cluster, splits, sort(clusters), genes, layers,
        is_bimodal = is_bimodal,
        bimodal_by_split = bimodal_by_split, n_cells = n_cells,
        split_label = split_label, ordered_splits = ordered,
        binarization = binarization_record("legacy scMarco-Ding-lab mixture model",
            paste0("h5v1:", basename(h5), ":", ds), BINARIZER_EXTERNAL),
        extra_uns = list(
            source = paste0("h5v1:", basename(h5), ":", ds),
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

    cells_path <- file.path(out_dir, paste0(ds, ".cells.h5ad"))
    write_cells_h5ad(
        cells_path, do.call(rbind, cell_x), cell_barcode, cell_split, cell_cluster,
        splits, sort(clusters), genes, reductions = reductions,
        reduction_meta = reduction_meta
    )
    message("Wrote ", cells_path, " with reductions: ",
        paste(names(reductions), collapse = ", "))
    invisible(list(clusters = clusters_path, cells = cells_path))
}
