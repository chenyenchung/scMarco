# Build a scMarco dataset and edit it in place: seurat_to_h5ad() writes a
# *built* dataset from Seurat objects, add_probability() binarizes it with the
# mixture-model output, and the cluster-name functions annotate or relabel
# its clusters. See docs/expected_input.md for the input formats.

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

#' Build an scMarco dataset from Seurat objects
#'
#' Writes `<out_dir>/<dataset>.clusters.h5ad` (one row per split x cluster,
#' with the layers `pct_pos` and `mean_lognorm` and the cell counts) and
#' `<out_dir>/<dataset>.cells.h5ad` (log-normalized expression per cell,
#' stored for fast per-gene reads, and reductions). The result is *built*:
#' before the app can use it, add the mixture-model probabilities with
#' `prob` or later with [add_probability()].
#'
#' Counts are normalized as Seurat's LogNormalize
#' (`log1p(counts / library size * 1e4)`). Seurat v3 objects saved before
#' SeuratObject existed are updated on read. Requires the SeuratObject
#' package; re-chunking the files for fast reads needs `h5repack` (HDF5
#' command-line tools) on the `PATH`.
#'
#' @param seurat One Seurat object (or `.rds` path) with raw counts, in which
#'   each cell's split and cluster ID come from the metadata columns
#'   `split_col` and `cluster_col`. Without either default column, legacy
#'   `"<split>@<cluster>"` identities are read instead. Or a named list or
#'   vector of objects or paths, one per split
#'   (`c(P15 = "a.rds", P30 = "b.rds")`); the genes are then the union of
#'   the objects' genes, with zero counts where a gene is missing.
#' @param dataset Dataset ID, used in the file names.
#' @param out_dir Output directory.
#' @param split_col,cluster_col Metadata columns holding each cell's split and
#'   cluster ID. `split_col = NULL` puts every cell in one split, `"all"`.
#' @param assay Assay holding the raw counts; `NULL` uses the default assay
#'   (pass e.g. `"RNA"` when the default is `"integrated"`).
#' @param split_label What a split is (`"stage"`, `"sex"`, ...), shown in the
#'   app.
#' @param split_order The splits in display order.
#' @param ordered The splits have a natural order (stages, doses); needs
#'   `split_order` when the splits come from one object.
#' @param prob Optional mixture-model output, as for [add_probability()].
#' @param umaps Optional folder of Seurat objects (`.rds`) with a `"umap"`
#'   reduction on a subset of the cells; the file name before the first `"_"`
#'   names the reduction.
#' @param reduction Reduction stored as `obsm/X_<reduction>_full`, or
#'   `X_<split>_<reduction>` per split (`X_<split>` for `"umap"`). `NULL`
#'   builds without one. A reduction named here must exist on every object;
#'   when the default `"umap"` is missing, the build warns and goes on
#'   without it.
#' @param exclude_clusters Cluster IDs whose cells are dropped.
#'
#' @return The paths of the two files (`clusters`, `cells`), invisibly.
#' @seealso [add_probability()], [validate_h5ad()].
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'     set.seed(1)
#'     counts <- matrix(rpois(20 * 60, 2), 20, 60,
#'         dimnames = list(paste0("gene", 1:20),
#'             paste0("cell", 1:60)))
#'     obj <- SeuratObject::CreateSeuratObject(counts)
#'     obj$scmarco_split <- rep(c("Male", "Female"), each = 30)
#'     obj$scmarco_cluster <- rep(c("1", "2", "3"), 20)
#'     out <- seurat_to_h5ad(obj, "demo", out_dir = tempdir(),
#'         split_label = "sex", reduction = NULL)
#'     validate_h5ad(out$clusters, stage = "built")$errors
#' }
#' @export
seurat_to_h5ad <- function(seurat, dataset, out_dir = "data",
                           split_col = "scmarco_split",
                           cluster_col = "scmarco_cluster", assay = NULL,
                           split_label = "split", split_order = NULL,
                           ordered = FALSE, prob = NULL, umaps = NULL,
                           reduction = "umap", exclude_clusters = character(0)) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        fail("seurat_to_h5ad() needs the SeuratObject package")
    }
    # a reduction the caller names must exist; the default may be missing
    reduction_required <- !missing(reduction)
    inputs <- seurat_inputs(seurat)
    # Read prob before the slow steps so a bad path fails early
    prob_tables <- if (!is.null(prob)) read_prob(prob) else NULL
    per_split <- !is.null(names(inputs))
    if (per_split && !identical(split_col, "scmarco_split")) {
        fail("split_col cannot be used with one object per split: ",
            "the split comes from the names of `seurat`")
    }

    # --- Cells -----------------------------------------------------------------
    parts <- list()
    for (k in seq_along(inputs)) {
        obj <- inputs[[k]]
        if (is.character(obj)) obj <- read_seurat(obj)
        ident <- as.character(SeuratObject::Idents(obj))
        labels <- if (per_split) {
            lab <- seurat_labels(obj[[]], ident, NULL, cluster_col, defaults = FALSE)
            lab$split <- rep(names(inputs)[k], length(lab$split))
            lab
        } else {
            seurat_labels(obj[[]], ident, split_col, cluster_col,
                defaults = identical(split_col, "scmarco_split") &&
                    identical(cluster_col, "scmarco_cluster"))
        }
        keep <- !labels$cluster %in% exclude_clusters
        if (!all(keep)) {
            message(sprintf("  dropping %d cells of excluded cluster(s)", sum(!keep)))
            labels <- lapply(labels, `[`, keep)
        }
        barcodes <- SeuratObject::Cells(obj)[keep]
        red <- if (!is.null(reduction)) {
            seurat_reduction(obj, reduction, barcodes,
                if (per_split) paste0(names(inputs)[k], ": ") else "",
                required = reduction_required)
        }
        parts[[k]] <- list(
            x = lognorm_by_cell(get_counts(obj, assay)[, barcodes, drop = FALSE]),
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
            reductions[[reduction_key(reduction)]] <- parts[[k]]$umap
            next
        }
        m <- matrix(NaN, length(barcodes), 2)
        m[sum(n_part[seq_len(k - 1)]) + seq_len(n_part[k]), ] <- parts[[k]]$umap
        reductions[[reduction_key(reduction, names(inputs)[k])]] <- m
    }

    x <- do.call(cbind, lapply(parts, function(p) expand_genes(p$x, genes)))
    rm(parts)
    gc()
    x <- Matrix::t(x) # cells x genes
    dimnames(x) <- list(barcodes, genes)

    splits <- if (!is.null(split_order)) {
        split_order
    } else if (per_split) {
        names(inputs)
    } else if (ordered) {
        fail("ordered = TRUE needs split_order, e.g. split_order = c(\"P15\", \"P30\", \"Adult\"); ",
            "sorting would put the splits in alphabetical order: ",
            paste(sort(unique(cell_split)), collapse = ","))
    } else {
        sort(unique(cell_split))
    }
    if (!setequal(splits, unique(cell_split))) {
        fail("split_order must list exactly the splits: ",
            paste(sort(unique(cell_split)), collapse = ","))
    }
    clusters <- order_cluster_ids(unique(cell_cluster))
    message(sprintf("%d cells, %d genes, %d splits, %d clusters",
        length(barcodes), length(genes), length(splits), length(clusters)))

    # --- Cluster-level layers --------------------------------------------------
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

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    clusters_path <- file.path(out_dir, paste0(dataset, ".clusters.h5ad"))
    write_clusters_h5ad(
        clusters_path, obs_split, obs_cluster, splits, clusters, genes, layers,
        n_cells = n_cells, split_label = split_label, ordered_splits = ordered,
        extra_uns = list(
            source = paste0("seurat:", paste(input_source(inputs), collapse = ",")),
            mean_lognorm_method = "arithmetic mean of log-normalized values per cluster"
        )
    )
    message("Wrote ", clusters_path)

    # --- Reductions ------------------------------------------------------------
    if (!is.null(umaps)) {
        for (f in list.files(umaps, "\\.rds$", full.names = TRUE)) {
            name <- strsplit(sub("\\.rds$", "", basename(f)), "_")[[1]][1]
            sub_obj <- readRDS(f)
            if (!"umap" %in% SeuratObject::Reductions(sub_obj)) {
                warning(basename(f), " has no umap reduction; skipped", call. = FALSE)
                next
            }
            emb <- SeuratObject::Embeddings(sub_obj, "umap")
            m <- matrix(NaN, length(barcodes), 2)
            hit <- match(rownames(emb), barcodes)
            m[hit[!is.na(hit)], ] <- emb[!is.na(hit), c(1, 2)]
            reductions[[name]] <- m
            message(sprintf("  reduction %s: %d cells", name, sum(!is.na(hit))))
        }
    }

    cells_path <- file.path(out_dir, paste0(dataset, ".cells.h5ad"))
    write_cells_h5ad(
        cells_path, x, barcodes, cell_split, cell_cluster, splits, clusters, genes,
        reductions = reductions
    )
    message("Wrote ", cells_path)

    if (!is.null(prob_tables)) {
        p <- assemble_prob(prob_tables, obs_split, obs_cluster, splits, genes)
        add_probability_layer(clusters_path, p$prob, p$bimodal_by_split, list(
            method = "mixture model", source = attr(prob_tables, "source"),
            date = format(Sys.Date())
        ), mixture = p$mixture, mixture_models = p$mixture_models)
        message(sprintf("Binarized %s: %d genes bimodal in at least one %s",
            clusters_path, sum(rowSums(p$bimodal_by_split) > 0), split_label))
    } else {
        message("Built without probabilities; add them before use with:\n  ",
            sprintf('add_probability("%s", c(<split> = "<file.rds>", ...))', clusters_path))
    }
    invisible(list(clusters = clusters_path, cells = cells_path))
}

seurat_inputs <- function(seurat) {
    # `seurat` -> list of objects or paths: unnamed of length 1 for one
    # object, named by split for one object per split
    if (is.character(seurat)) seurat <- as.list(seurat)
    if (!is.list(seurat) || inherits(seurat, "Seurat")) seurat <- list(seurat)
    per_split <- !is.null(names(seurat))
    if (!per_split && length(seurat) != 1) {
        fail("Several Seurat objects must be named by split, e.g. c(P15 = \"a.rds\", P30 = \"b.rds\")")
    }
    if (per_split) {
        if (any(names(seurat) == "")) fail("Every object of `seurat` needs a split name")
        if (anyDuplicated(names(seurat))) {
            fail("Split(s) given twice in `seurat`: ",
                first_few(unique(names(seurat)[duplicated(names(seurat))])))
        }
    }
    is_path <- vapply(seurat, is.character, logical(1))
    paths <- as.character(unlist(seurat[is_path]))
    missing <- paths[!file.exists(paths)]
    if (length(missing) > 0) fail("Seurat object(s) not found: ", first_few(missing))
    bad <- !is_path & !vapply(seurat, inherits, logical(1), "Seurat")
    if (any(bad)) fail("`seurat` must hold Seurat objects or .rds paths")
    seurat
}

input_source <- function(inputs) {
    # Provenance label of each input: the file name, or "<object>"
    vapply(inputs, function(x) if (is.character(x)) basename(x) else "<object>", "")
}

get_counts <- function(obj, assay = NULL) {
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(obj)
    if (!assay %in% SeuratObject::Assays(obj)) {
        fail("Assay '", assay, "' not found. Assays: ",
            paste(SeuratObject::Assays(obj), collapse = ", "))
    }
    a <- obj[[assay]]
    if (inherits(a, "Assay5")) {
        # Seurat v5 may keep the counts split by sample (counts.1, counts.2,
        # ...); LayerData() would return only the first of them
        layers <- grep("^counts(\\.|$)", SeuratObject::Layers(a), value = TRUE)
        if (length(layers) == 0) {
            fail("Assay '", assay, "' has no counts layer. Layers: ",
                paste(SeuratObject::Layers(a), collapse = ", "))
        }
        if (length(layers) > 1) {
            message(sprintf("  joining the %d counts layers of assay '%s' (%s)",
                length(layers), assay, first_few(layers)))
            a <- SeuratObject::JoinLayers(a, layers = "counts")
        }
        counts <- SeuratObject::LayerData(a, layer = "counts")
    } else {
        counts <- tryCatch(
            SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
            error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
        )
    }
    as(counts, "CsparseMatrix")
}

lognorm_by_cell <- function(counts) {
    # genes x cells counts -> genes x cells LogNormalize values (genes sorted)
    counts <- counts[sort(rownames(counts)), , drop = FALSE]
    lib_size <- Matrix::colSums(counts)
    x <- counts %*% Matrix::Diagonal(x = 1e4 / pmax(lib_size, 1))
    x@x <- log1p(x@x)
    dimnames(x) <- dimnames(counts)
    as(x, "CsparseMatrix")
}

expand_genes <- function(x, genes) {
    # genes x cells dgCMatrix with sorted genes -> the same on the sorted
    # superset `genes`, with zeros for the added genes. Row indices map
    # monotonically, so each column stays sorted.
    if (identical(rownames(x), genes)) return(x)
    idx <- match(rownames(x), genes)
    out <- x
    out@i <- as.integer(idx[x@i + 1] - 1)
    out@Dim <- c(length(genes), ncol(x))
    out@Dimnames <- list(genes, colnames(x))
    out
}

# ---------------------------------------------------------------------------
# Binarize
# ---------------------------------------------------------------------------

#' Add mixture-model probabilities to a dataset
#'
#' Binarizes a *built* dataset (see [seurat_to_h5ad()]) by writing the
#' mixture model's probabilities of expression, the bimodality of each gene
#' in each split and, when present, the fit parameters into its clusters
#' file, in place. Running it again replaces the earlier probabilities.
#'
#' Each split's table has genes as row names, one column per cluster ID with
#' the probability of expression, and QC columns, among them `exprType`
#' (`"bimodal"` marks the genes used in that split) and, for the mixture
#' model plot, `bimodal.on_level`, `bimodal.off_level`, `bimodal.sd` and
#' `bimodal.pi`.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param prob The mixture model's output, one table per split: a named
#'   vector of `.rds` paths (`c(Male = "male.rds", Female = "female.rds")`),
#'   a named list of data.frames, an `.rds` path holding such a list, or a
#'   legacy wide table with `"<split>@<cluster>"` columns (splits fit
#'   together, as one model).
#' @param validate Check the file afterwards with [validate_h5ad()] and stop
#'   on errors.
#'
#' @return The [validate_h5ad()] result (`NULL` without `validate`),
#'   invisibly.
#' @examples
#' toy <- system.file("extdata", "toy.clusters.h5ad", package = "scMarco")
#' path <- tempfile(fileext = ".clusters.h5ad")
#' file.copy(toy, path)
#' ds <- ds_load("toy", list(label = "Toy", clusters = path))
#' # one table per split: genes x cluster IDs, plus exprType
#' tables <- lapply(c(Male = "Male", Female = "Female"), function(s) {
#'     present <- ds$clusters[ds$clusters %in% ds$obs_cluster[ds$obs_split == s]]
#'     tab <- as.data.frame(matrix(0.01, length(ds$all_genes), length(present),
#'         dimnames = list(ds$all_genes, present)),
#'     check.names = FALSE)
#'     tab[1, ] <- 0.99
#'     tab$exprType <- "bimodal"
#'     tab
#' })
#' suppressWarnings(add_probability(path, tables))
#' @export
add_probability <- function(clusters, prob, validate = TRUE) {
    ch <- h5ad_open(clusters)
    if ("prob" %in% ch$layers) {
        message(clusters, " already has probabilities; replacing them")
    }
    splits <- levels(ch$obs$split)
    tables <- read_prob(prob)
    p <- assemble_prob(tables, as.character(ch$obs$split), as.character(ch$obs$cluster),
        splits, ch$genes)
    add_probability_layer(clusters, p$prob, p$bimodal_by_split, list(
        method = "mixture model", source = attr(tables, "source"), date = format(Sys.Date())
    ), mixture = p$mixture, mixture_models = p$mixture_models)
    groups <- split(names(p$mixture_models), p$mixture_models)
    message("Mixture models: ", paste(vapply(groups, paste, "", collapse = " + "),
        collapse = "; "))
    n_bim <- colSums(p$bimodal_by_split)
    message("Wrote probabilities to ", clusters, "; bimodal genes per split: ",
        paste(sprintf("%s %d", names(n_bim), n_bim), collapse = ", "))

    if (!validate) return(invisible(NULL))
    res <- validate_h5ad(clusters, stage = "binarized")
    if (!report_validation(res, clusters)) {
        fail("Validation failed after writing ", clusters)
    }
    invisible(res)
}

# ---------------------------------------------------------------------------
# Cluster names and IDs
# ---------------------------------------------------------------------------
# obs/cluster holds stable cluster IDs; obs/cluster_name holds the optional
# annotation, shown in the app as "name (id)". Only obs is rewritten, in
# place, so these take seconds even for large files.

#' Read or change cluster names
#'
#' Clusters have a stable ID (`obs/cluster`), used for every join, and an
#' optional name (`obs/cluster_name`), shown in the app as `"name (id)"`.
#' These functions rewrite only `obs`, in place, so they take seconds even
#' for large files.
#'
#' `set_cluster_names()` names the listed clusters (`NA` or `""` clears a
#' name) and keeps the names of the others.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param file Optional path: also write the table there as a TSV (empty =
#'   unnamed) that `set_cluster_names()` reads back.
#' @param names A data.frame or TSV path with columns `cluster_id` and
#'   `name`, or a named character vector `c(<id> = <name>)`.
#'
#' @return A data.frame with columns `cluster_id` and `name` (`NA` =
#'   unnamed); `set_cluster_names()` returns the names after the change,
#'   invisibly.
#' @seealso [set_cluster_ids()] to change the IDs.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' path <- tempfile(fileext = ".clusters.h5ad")
#' file.copy(ex, path)
#' head(get_cluster_names(path))
#' id <- get_cluster_names(path)$cluster_id[1]
#' set_cluster_names(path, setNames("my favorite cluster", id))
#' @export
get_cluster_names <- function(clusters, file = NULL) {
    ch <- h5ad_open(clusters)
    ids <- levels(ch$obs$cluster)
    current <- read_cluster_names(ch$obs, ids, clusters)
    out <- data.frame(cluster_id = ids, name = unname(current))
    if (!is.null(file)) {
        write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        message(sprintf("Wrote %d clusters (%d named) to %s",
            nrow(out), sum(!is.na(current)), file))
    }
    out
}

#' @rdname get_cluster_names
#' @export
set_cluster_names <- function(clusters, names) {
    ch <- h5ad_open(clusters)
    ids <- levels(ch$obs$cluster)
    current <- read_cluster_names(ch$obs, ids, clusters)
    # Names are stored per obs row, so a cluster absent from every split
    # cannot hold one
    no_rows <- setdiff(ids, as.character(ch$obs$cluster))

    tbl <- if (is.character(names) && !is.null(base::names(names))) {
        data.frame(cluster_id = base::names(names), name = unname(names))
    } else {
        read_id_table(names, c("cluster_id", "name"))
    }
    tbl$name[is.na(tbl$name)] <- ""
    unknown <- setdiff(tbl$cluster_id, ids)
    if (length(unknown) > 0) {
        fail(length(unknown), " cluster ID(s) not in ", clusters, ": ",
            paste(head(unknown, 10), collapse = ", "))
    }
    dup <- unique(tbl$cluster_id[duplicated(tbl$cluster_id)])
    if (length(dup) > 0) {
        fail("Cluster ID(s) listed more than once: ", paste(head(dup, 10), collapse = ", "))
    }
    homeless <- intersect(tbl$cluster_id[tbl$name != ""], no_rows)
    if (length(homeless) > 0) {
        fail("Cluster(s) with no rows in any split cannot be named: ",
            paste(head(homeless, 10), collapse = ", "))
    }
    updated <- current
    updated[tbl$cluster_id] <- ifelse(tbl$name == "", NA, tbl$name)
    result <- data.frame(cluster_id = ids, name = unname(updated))
    changed <- base::names(updated)[!mapply(identical, updated, current)]
    shared <- unique(updated[!is.na(updated) & duplicated(updated)])
    if (length(shared) > 0) {
        message("Note: name(s) shared by several clusters (kept apart by their IDs): ",
            paste(head(shared, 5), collapse = ", "))
    }
    if (length(changed) == 0) {
        message("No changes")
        return(invisible(result))
    }
    rewrite_obs(clusters, function(obs) names_obs(obs, updated))
    for (id in head(changed, 20)) {
        message(sprintf("  %s: %s -> %s", id,
            if (is.na(current[[id]])) "(unnamed)" else current[[id]],
            if (is.na(updated[[id]])) "(unnamed)" else updated[[id]]))
    }
    if (length(changed) > 20) message("  ... and ", length(changed) - 20, " more")
    message(sprintf("Updated %d cluster(s); %d of %d now named",
        length(changed), sum(!is.na(updated)), length(updated)))
    invisible(result)
}

#' Read or change cluster metadata
#'
#' Cluster metadata are per-cluster groupings (e.g. class and subclass),
#' stored as `obs/meta_<level>`, that let the app limit every search to part
#' of the dataset. The config's `subset_by` picks which levels the app
#' offers, coarse to fine.
#'
#' `set_cluster_metadata()` replaces all metadata. Unlisted clusters and
#' empty values are unassigned.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param file Optional path: also write the table there as a TSV (empty =
#'   unassigned) that `set_cluster_metadata()` reads back.
#' @param meta A data.frame or TSV path with a `cluster_id` column and one
#'   column per level, coarse to fine; `NULL` removes the metadata.
#'
#' @return A data.frame with `cluster_id` and one column per level (`NA` =
#'   unassigned), or `NULL` when there is none. `set_cluster_metadata()`
#'   returns it invisibly.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' path <- tempfile(fileext = ".clusters.h5ad")
#' file.copy(ex, path)
#' ids <- get_cluster_names(path)$cluster_id
#' set_cluster_metadata(path, data.frame(
#'     cluster_id = ids, half = rep(c("first", "second"), length.out = length(ids))
#' ))
#' head(get_cluster_metadata(path))
#' @export
get_cluster_metadata <- function(clusters, file = NULL) {
    ch <- h5ad_open(clusters)
    ids <- levels(ch$obs$cluster)
    meta <- read_cluster_meta(ch$obs, ids, clusters)
    if (is.null(meta)) {
        message("No cluster metadata in ", clusters)
        return(invisible(NULL))
    }
    out <- cbind(data.frame(cluster_id = ids), meta)
    rownames(out) <- NULL
    if (!is.null(file)) {
        write.table(out, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
        message(sprintf("Wrote %d clusters x %d levels to %s", nrow(out), ncol(meta), file))
    }
    out
}

#' @rdname get_cluster_metadata
#' @export
set_cluster_metadata <- function(clusters, meta) {
    ch <- h5ad_open(clusters)
    ids <- levels(ch$obs$cluster)
    if (is.null(meta)) {
        rewrite_obs(clusters, function(obs) meta_obs(obs, NULL))
        message("Removed the cluster metadata of ", clusters)
        return(invisible(NULL))
    }
    tbl <- read_id_table(meta, "cluster_id")
    levels <- setdiff(names(tbl), "cluster_id")
    if (length(levels) == 0) {
        fail("The table needs at least one column besides cluster_id ",
            "(to remove the metadata, pass meta = NULL)")
    }
    bad <- levels[levels == "" | levels != trimws(levels) | grepl("/", levels, fixed = TRUE)]
    if (length(bad) > 0 || anyDuplicated(levels)) {
        fail("Column names must be unique, non-empty, without '/' or surrounding spaces: ",
            paste(unique(c(bad, levels[duplicated(levels)])), collapse = ", "))
    }
    unknown <- setdiff(tbl$cluster_id, ids)
    if (length(unknown) > 0) {
        fail(length(unknown), " cluster ID(s) not in ", clusters, ": ",
            paste(head(unknown, 10), collapse = ", "))
    }
    dup <- unique(tbl$cluster_id[duplicated(tbl$cluster_id)])
    if (length(dup) > 0) {
        fail("Cluster ID(s) listed more than once: ", paste(head(dup, 10), collapse = ", "))
    }
    for (level in levels) tbl[[level]][!is.na(tbl[[level]]) & tbl[[level]] == ""] <- NA
    no_rows <- setdiff(ids, as.character(ch$obs$cluster))
    assigned <- rowSums(!is.na(tbl[levels])) > 0
    homeless <- intersect(tbl$cluster_id[assigned], no_rows)
    if (length(homeless) > 0) {
        fail("Cluster(s) with no rows in any split cannot hold metadata: ",
            paste(head(homeless, 10), collapse = ", "))
    }

    full <- data.frame(cluster_id = ids)
    for (level in levels) full[[level]] <- tbl[[level]][match(ids, tbl$cluster_id)]
    rewrite_obs(clusters, function(obs) meta_obs(obs, full))
    for (level in levels) {
        values <- full[[level]][!is.na(full[[level]])]
        message(sprintf("  %s: %d value(s) (%s), %d cluster(s) unassigned", level,
            length(unique(values)),
            paste(head(sort(unique(values)), 5), collapse = ", "),
            sum(is.na(full[[level]]))))
    }
    message(sprintf("Set %d metadata level(s) for %d clusters: %s", length(levels),
        length(ids), paste(levels, collapse = " > ")))
    invisible(full)
}

#' Give every cluster a new ID
#'
#' Renames cluster IDs (and sets names) in the clusters file and, if given,
#' the cells file. Rows and data are unchanged.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param ids A data.frame or TSV path with columns `cluster_id` (the current
#'   ID), `new_id` and optionally `name` (empty = unnamed). Every current
#'   cluster must be listed once.
#' @param cells Optional path of the dataset's `.cells.h5ad` file, so that its
#'   cluster IDs stay in sync.
#'
#' @return A data.frame with columns `label` (old ID), `id` and `name`,
#'   invisibly.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' path <- tempfile(fileext = ".clusters.h5ad")
#' file.copy(ex, path)
#' old <- get_cluster_names(path)$cluster_id
#' set_cluster_ids(path, data.frame(cluster_id = old, new_id = paste0("c", old)))
#' head(get_cluster_names(path))
#' @export
set_cluster_ids <- function(clusters, ids, cells = NULL) {
    ch <- h5ad_open(clusters)
    current_ids <- levels(ch$obs$cluster)
    tbl <- read_id_table(ids, c("cluster_id", "new_id"))
    if (!"name" %in% names(tbl)) tbl$name <- ""
    tbl$name[is.na(tbl$name)] <- ""
    unknown <- setdiff(tbl$cluster_id, current_ids)
    missing <- setdiff(current_ids, tbl$cluster_id)
    if (length(unknown) > 0 || length(missing) > 0 || anyDuplicated(tbl$cluster_id)) {
        fail("`ids` must list every cluster ID of ", clusters, " exactly once",
            if (length(unknown) > 0) paste0("\nUnknown: ", paste(head(unknown, 10), collapse = ", ")),
            if (length(missing) > 0) paste0("\nMissing: ", paste(head(missing, 10), collapse = ", ")))
    }
    if (any(tbl$new_id == "")) {
        fail("Empty new_id for: ", paste(tbl$cluster_id[tbl$new_id == ""], collapse = ", "))
    }
    tbl <- tbl[match(current_ids, tbl$cluster_id), ]
    parts <- data.frame(label = tbl$cluster_id, id = tbl$new_id,
        name = ifelse(tbl$name == "", NA, tbl$name), stringsAsFactors = FALSE)
    named <- relabel_clusters(clusters, cells, ch, parts)
    message(sprintf("Gave %d clusters new IDs; %d named, %d unnamed",
        nrow(parts), sum(named), sum(!named)))
    show_relabeled(parts, named)
}

#' Split legacy cluster labels into ID and name
#'
#' One-off migration of files whose cluster labels embed the ID:
#' `"Mi4 (154)"` becomes ID `154` with name `Mi4`; a label without a
#' trailing `"(id)"` is kept as an unnamed ID.
#'
#' @param clusters Path of the `.clusters.h5ad` file.
#' @param cells Optional path of the dataset's `.cells.h5ad` file, so that its
#'   cluster IDs stay in sync.
#'
#' @return A data.frame with columns `label`, `id` and `name`, invisibly.
#' @examples
#' ex <- system.file("extdata", "example.clusters.h5ad", package = "scMarco")
#' path <- tempfile(fileext = ".clusters.h5ad")
#' file.copy(ex, path)
#' # the example's clusters are already IDs, so nothing changes
#' migrate_cluster_labels(path)
#' @export
migrate_cluster_labels <- function(clusters, cells = NULL) {
    ch <- h5ad_open(clusters)
    labels <- levels(ch$obs$cluster)
    if (any(!is.na(read_cluster_names(ch$obs, labels, clusters)))) {
        fail(clusters, " already has cluster names; migrate_cluster_labels() is ",
            "for files whose cluster labels still embed them")
    }
    parts <- split_cluster_labels(labels)
    if (all(parts$id == parts$label)) {
        message("No label embeds an ID; nothing to change")
        return(invisible(parts))
    }
    named <- relabel_clusters(clusters, cells, ch, parts)
    message(sprintf("Split %d labels into ID + name; %d bare IDs kept as unnamed",
        sum(named), sum(!named)))
    show_relabeled(parts, named)
}

relabel_clusters <- function(clusters, cells, ch, parts) {
    # Rewrite both files' obs with new IDs and names.
    # parts: data.frame(label = current ID, id = new ID, name), one row per
    #   cluster. Returns which rows kept their name (logical).
    dup <- unique(parts$id[duplicated(parts$id)])
    if (length(dup) > 0) {
        fail("Duplicate new cluster IDs: ", paste(head(dup, 10), collapse = ", "),
            "\nOld IDs: ", paste(parts$label[parts$id %in% dup], collapse = ", "))
    }
    if (any(grepl("|", parts$id, fixed = TRUE))) fail("Cluster IDs must not contain '|'")
    no_rows <- setdiff(parts$label, as.character(ch$obs$cluster))
    lost <- parts$label[!is.na(parts$name) & parts$label %in% no_rows]
    if (length(lost) > 0) {
        warning(length(lost), " cluster(s) have no rows in any split, so only their ",
            "ID is kept: ", paste(lost, collapse = ", "), call. = FALSE, immediate. = TRUE)
    }
    id_of <- setNames(parts$id, parts$label)
    new_names <- setNames(parts$name, parts$id)
    id_levels <- order_cluster_ids(parts$id)

    if (!is.null(cells)) {
        xh <- h5ad_open(cells)
        stray <- setdiff(unique(as.character(xh$obs$cluster)), parts$label)
        if (length(stray) > 0) {
            fail("Cells file has clusters missing from the clusters file: ",
                paste(head(stray, 5), collapse = ", "))
        }
    }

    relabel <- function(obs) {
        obs$cluster <- factor(unname(id_of[as.character(obs$cluster)]), levels = id_levels)
        obs
    }
    rewrite_obs(clusters, function(obs) {
        obs <- names_obs(relabel(obs), new_names)
        rownames(obs) <- paste(obs$split, obs$cluster, sep = "|")
        obs
    })
    if (!is.null(cells)) rewrite_obs(cells, relabel)
    !is.na(parts$name) & !parts$label %in% lost
}

show_relabeled <- function(parts, named) {
    for (i in head(which(named), 5)) {
        message(sprintf("  %s -> ID %s, name %s", parts$label[i], parts$id[i], parts$name[i]))
    }
    invisible(parts)
}

names_obs <- function(obs, names) {
    # obs$cluster_name from an id -> name vector; dropped when nothing is named
    obs$cluster_name <- NULL
    if (any(!is.na(names))) {
        per_row <- unname(names[as.character(obs$cluster)])
        obs$cluster_name <- factor(per_row, levels = sort(unique(names[!is.na(names)])))
    }
    obs
}

meta_obs <- function(obs, meta) {
    # obs with its meta_<level> columns replaced by `meta` (data.frame with
    # cluster_id and one column per level; NULL drops them all)
    obs <- obs[!startsWith(names(obs), META_PREFIX)]
    if (is.null(meta)) return(obs)
    row <- match(as.character(obs$cluster), meta$cluster_id)
    for (level in setdiff(names(meta), "cluster_id")) {
        values <- meta[[level]]
        obs[[paste0(META_PREFIX, level)]] <- factor(values[row],
            levels = sort(unique(values[!is.na(values)])))
    }
    obs
}

read_id_table <- function(x, required) {
    # data.frame or TSV path -> data.frame of trimmed character columns
    if (is.character(x) && length(x) == 1) {
        tbl <- read.delim(x, colClasses = "character", na.strings = character(0),
            check.names = FALSE, quote = "")
        what <- x
    } else if (is.data.frame(x)) {
        tbl <- x
        tbl[] <- lapply(tbl, as.character)
        what <- "The table"
    } else {
        fail("Expected a data.frame or the path of a tab-separated table")
    }
    if (!all(required %in% names(tbl))) {
        fail(what, " needs the columns ", paste(required, collapse = " and "),
            if (what != "The table") " (tab-separated)")
    }
    tbl[] <- lapply(tbl, trimws)
    tbl
}
