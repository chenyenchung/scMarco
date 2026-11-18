# Input checks and readers for building a dataset from new data: split and
# cluster labels from a Seurat object, and mixture-model probability tables.
# Used by seurat_to_h5ad() and add_probability() in lib_build.R; see
# docs/expected_input.md for the formats and rules.

QC_COLS <- c("diff.elpd", "diff.elpd.se", "elpd.z", "elpd.q", "ashman.D", "ashman.D.lower",
    "spread.share", "bimodal.on_level", "bimodal.off_level", "bimodal.sd", "bimodal.pi",
    "unimodal.level", "unimodal.sd", "exprType", "not.binarized")

# Mixture-model fit parameters kept per gene x split (varm/mixture_<name>):
# the two Gaussian component means, their shared SD and the ON weight
MIXTURE_COLS <- c(on_level = "bimodal.on_level", off_level = "bimodal.off_level",
    sd = "bimodal.sd", pi = "bimodal.pi")

# Binarizer 1.1 per-gene values kept per gene x split (varm/call_<name>):
# evidence for two states (elpd z and its Benjamini-Hochberg q), Ashman's D
# (posterior median and 10th percentile) and the share of clusters whose
# chains disagree
CALL_STAT_COLS <- c(elpd_z = "elpd.z", elpd_q = "elpd.q", ashman_D = "ashman.D",
    ashman_D_lower = "ashman.D.lower", spread_share = "spread.share")
# Parameters of call_mixture() recorded per split in
# uns/scmarco/binarization/calls
CALL_PARAMS <- c("summary", "bh_q", "min_D", "spread_cutoff", "reject_spread_frac",
    "on_floor", "ambient_log2fc", "ambient_alpha", "ambient")

fail <- function(...) stop(..., call. = FALSE)

first_few <- function(x, n = 5) {
    paste0(paste(head(x, n), collapse = ", "), if (length(x) > n) ", ..." else "")
}

# ---------------------------------------------------------------------------
# Split and cluster labels
# ---------------------------------------------------------------------------

seurat_labels <- function(meta, ident, split_col = "scmarco_split",
                          cluster_col = "scmarco_cluster", defaults = TRUE) {
    # Split and cluster label of every cell.
    # meta: the object's meta.data; ident: as.character(Idents(obj)).
    # split_col NULL puts every cell in one split named "all".
    # defaults: TRUE when both column names are the defaults, so the legacy
    #   "<split>@<cluster>" identities may be used when neither column exists.
    has <- function(col) col %in% colnames(meta)
    if (!has(cluster_col) && !(!is.null(split_col) && has(split_col)) && defaults) {
        message("No '", cluster_col, "' / '", split_col, "' metadata columns; ",
            "reading legacy '<split>@<cluster>' identities")
        parsed <- strsplit(ident, "@", fixed = TRUE)
        bad <- lengths(parsed) != 2
        if (any(bad)) {
            fail("Found neither metadata columns '", split_col, "' and '", cluster_col,
                "' nor identities of the form '<split>@<cluster>' (e.g. ",
                first_few(unique(ident[bad]), 3), ").\nAdd the columns, e.g. ",
                "obj$", split_col, " <- obj$sex; obj$", cluster_col,
                " <- obj$seurat_clusters, or name yours with split_col / cluster_col.")
        }
        return(list(split = vapply(parsed, `[`, "", 1),
            cluster = vapply(parsed, `[`, "", 2)))
    }
    if (!has(cluster_col)) {
        fail("Metadata column '", cluster_col, "' not found. Columns: ",
            first_few(colnames(meta), 10), ". Name the cluster column with cluster_col.")
    }
    cluster <- meta[[cluster_col]]
    if (is.null(split_col)) {
        split <- rep("all", nrow(meta))
    } else if (!has(split_col)) {
        fail("Metadata column '", split_col, "' not found. Name the split column with ",
            "split_col, or pass split_col = NULL for a dataset with one split.")
    } else {
        split <- meta[[split_col]]
    }
    list(split = as.character(split), cluster = as.character(cluster))
}

retag_seurat_v3 <- function(x) {
    # Seurat v3 objects saved before SeuratObject existed name their classes'
    # package "Seurat"; point them to SeuratObject so it can load them without
    # Seurat. Walks S4 slots and plain lists only. UpdateSeuratObject() and
    # SeuratObject's internal UpdateClassPkg() both look the classes up first
    # and fail with "unable to load required package 'Seurat'" (SeuratObject 5.4.0).
    if (isS4(x)) {
        cl <- class(x)
        if (identical(attr(cl, "package"), "Seurat")) attr(cl, "package") <- "SeuratObject"
        slots <- attributes(x)
        slots$class <- NULL
        for (n in names(slots)) attr(x, n) <- retag_seurat_v3(slots[[n]])
        class(x) <- cl
    } else if (is.list(x) && identical(class(x), "list")) {
        for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- retag_seurat_v3(x[[i]])
    }
    x
}

read_seurat <- function(path) {
    # readRDS() a Seurat object; objects from before SeuratObject (Seurat v3)
    # are retagged and updated so SeuratObject alone can read them.
    obj <- readRDS(path)
    if (identical(attr(class(obj), "package"), "Seurat")) {
        message(basename(path), ": Seurat ", as.character(obj@version),
            " object; updating it for SeuratObject")
        obj <- suppressWarnings(suppressMessages(
            SeuratObject::UpdateSeuratObject(retag_seurat_v3(obj))
        ))
    }
    obj
}

seurat_reduction <- function(obj, red_name, cells, context = "", required = TRUE) {
    # First two dimensions of reduction `red_name` for `cells`. A missing
    # reduction lists the object's reductions: an error when `required` (the
    # caller named it), else a warning and NULL (no reduction stored).
    # context: prefix for the message, e.g. "P15: "
    if (!red_name %in% SeuratObject::Reductions(obj)) {
        found <- SeuratObject::Reductions(obj)
        found <- if (length(found) > 0) paste(found, collapse = ", ") else "none"
        if (!required) {
            warning(context, "no reduction named '", red_name, "' (reductions on the object: ",
                found, "), so none is stored and the UMAP tab stays hidden. Pass ",
                "reduction = \"<name>\" to store another one, or reduction = NULL to ",
                "build without one and silence this warning.", call. = FALSE)
            return(NULL)
        }
        fail(context, "no reduction named '", red_name, "'. Reductions on the object: ",
            found, ". Pass reduction = \"<name>\" to store one of them, or reduction = NULL to build without one.")
    }
    SeuratObject::Embeddings(obj, red_name)[cells, c(1, 2), drop = FALSE]
}

reduction_key <- function(red_name, split = NULL) {
    # obsm name (without "X_") of a reduction: "<red>_full" for one object,
    # "<split>" (umap) or "<split>_<red>" for one object per split
    if (is.null(split)) return(paste0(red_name, "_full"))
    if (red_name == "umap") split else paste0(split, "_", red_name)
}

check_labels <- function(split, cluster) {
    # Name rules for split and cluster labels (stop on the first broken rule)
    for (kind in c("split", "cluster")) {
        x <- if (kind == "split") split else cluster
        if (anyNA(x) || any(x == "")) {
            fail(sum(is.na(x) | x == ""), " cell(s) have an empty or missing ", kind,
                " label. Every cell needs one; subset the object to drop them.")
        }
        bad <- unique(x[grepl("|", x, fixed = TRUE)])
        if (length(bad) > 0) {
            fail(kind, " labels must not contain '|', which the files use as a separator: ",
                first_few(bad))
        }
        bad <- unique(x[x != trimws(x)])
        if (length(bad) > 0) {
            fail(kind, " labels must not start or end with whitespace: ",
                first_few(sprintf("'%s'", bad)))
        }
    }
    invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Mixture-model probabilities
# ---------------------------------------------------------------------------

read_prob <- function(prob) {
    # Mixture-model output -> named list split -> table (genes x clusters + QC
    # columns), with attr "models": split -> the mixture model it was fit in,
    # attr "binarizer": split -> the binarizer version its table records
    # (BINARIZER_EXTERNAL when it records none), and attr "source": where the
    # tables came from.
    # prob: a character vector of .rds paths or a list of tables / paths;
    #   each item is one of
    #   named, c(Male = "a.rds") or list(Male = df): one split's table,
    #     modeled on its own;
    #   unnamed, holding a named list of tables: one model per split;
    #   unnamed, holding a wide table with "<split>@<cluster>" columns: the
    #     splits in it were modeled together, as one model.
    # A single table (data.frame or matrix) is taken as a wide table.
    if (is.data.frame(prob) || is.matrix(prob)) prob <- list(prob)
    if (is.character(prob)) prob <- as.list(prob)
    if (!is.list(prob) || length(prob) == 0) {
        fail("prob must be .rds path(s) or a list of tables")
    }
    item_names <- if (is.null(names(prob))) rep("", length(prob)) else names(prob)
    is_path <- vapply(prob, function(x) is.character(x) && length(x) == 1, logical(1))
    paths <- as.character(unlist(prob[is_path]))
    missing <- paths[!file.exists(paths)]
    if (length(missing) > 0) fail("Probability file(s) not found: ", first_few(missing))

    tables <- list()
    models <- character(0)
    for (i in seq_along(prob)) {
        obj <- if (is_path[i]) readRDS(prob[[i]]) else prob[[i]]
        what <- if (is_path[i]) prob[[i]] else sprintf("prob item %d", i)
        if (item_names[i] != "") {
            add <- setNames(list(obj), item_names[i])
            model <- names(add)
        } else if (is.data.frame(obj) || is.matrix(obj)) {
            message("Reading ", what, " as a wide table with '<split>@<cluster>' ",
                "columns: its splits were modeled together")
            add <- legacy_wide_to_split_tables(obj)
            model <- rep(paste(names(add), collapse = "+"), length(add))
        } else if (is.list(obj) && !is.null(names(obj)) && all(names(obj) != "")) {
            add <- obj
            model <- names(add)
        } else {
            fail(what, " must hold a named list of per-split tables, e.g. ",
                "list(Male = male_df, Female = female_df), or a wide table with ",
                "'<split>@<cluster>' columns")
        }
        tables <- c(tables, add)
        models <- c(models, setNames(model, names(add)))
    }
    attr(tables, "models") <- models
    attr(tables, "binarizer") <- vapply(tables, table_binarizer, "")
    src <- ifelse(is_path, vapply(prob, function(x) if (is.character(x)) x[1] else "", ""),
        "<table>")
    attr(tables, "source") <- paste(ifelse(item_names == "", src, paste0(item_names, "=", src)),
        collapse = ",")
    tables
}

table_binarizer <- function(tab) {
    # Binarizer version a call_mixture() table records, else external
    v <- attr(tab, "scmarco_binarizer", exact = TRUE)
    if (is.null(v)) return(BINARIZER_EXTERNAL)
    if (!valid_binarizer_version(v)) fail("Table records an invalid binarizer version: ", format(v))
    v
}

prob_binarizer <- function(tables, external = FALSE) {
    # The one binarizer version to record for a dataset's tables. Tables that
    # record no scMarco binarizer (not from call_mixture()) are refused unless
    # `external` says so; then every split is recorded as external. Tables of
    # different binarizers would make probabilities that are not comparable
    # across splits.
    v <- attr(tables, "binarizer")
    unknown <- names(v)[v == BINARIZER_EXTERNAL]
    if (isTRUE(external)) {
        known <- v[v != BINARIZER_EXTERNAL]
        if (length(known) > 0) {
            warning("external = TRUE: recording every split as \"", BINARIZER_EXTERNAL,
                "\", although the tables of ", paste(names(known), collapse = ", "),
                " record binarizer ", paste(unique(known), collapse = ", "),
                ". Drop external = TRUE to keep their version.", call. = FALSE)
        }
        return(BINARIZER_EXTERNAL)
    }
    if (length(unknown) > 0) {
        fail("The table(s) of split(s) ", first_few(unknown), " do not record a scMarco ",
            "binarizer version, so where their probabilities came from is unknown. ",
            "Make them with call_mixture(), or, if they really come from elsewhere ",
            "(e.g. the mixture_modeling_smk pipeline), pass external = TRUE to record ",
            "every split as \"", BINARIZER_EXTERNAL, "\".")
    }
    if (length(unique(v)) > 1) {
        groups <- split(names(v), v)
        fail("The splits' tables come from different binarizers (",
            paste(sprintf("%s: %s", names(groups), vapply(groups, paste, "", collapse = ", ")),
                collapse = "; "),
            "). Make every table with the same binarizer, e.g. call_mixture() of this scMarco.")
    }
    status <- binarizer_status(unname(v[1]))
    if (status == "newer") {
        fail("The tables were made by binarizer ", v[1], ", newer than this scMarco's (",
            SCMARCO_BINARIZER_VERSION, "); update scMarco")
    }
    if (status == "older") {
        warning("The tables were made by binarizer ", v[1], "; this scMarco's is ",
            SCMARCO_BINARIZER_VERSION, ". Consider refitting.", call. = FALSE)
    }
    unname(v[1])
}

legacy_wide_to_split_tables <- function(df) {
    # Wide table (one model fit across splits) -> per-split tables. Its single
    # exprType (and QC columns) apply to every split.
    text_matrix <- is.matrix(df) && is.character(df)
    df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    qc <- intersect(colnames(df), QC_COLS)
    prob_cols <- colnames(df)[!colnames(df) %in% QC_COLS]
    parts <- strsplit(prob_cols, "@", fixed = TRUE)
    bad <- prob_cols[lengths(parts) != 2]
    if (length(bad) > 0) {
        fail("Legacy wide table columns must be '<split>@<cluster>' (or QC columns): ",
            first_few(bad))
    }
    col_split <- vapply(parts, `[`, "", 1)
    col_cluster <- vapply(parts, `[`, "", 2)
    tables <- lapply(unique(col_split), function(s) {
        out <- df[, c(prob_cols[col_split == s], qc), drop = FALSE]
        colnames(out) <- c(col_cluster[col_split == s], qc)
        # keep the raw model output's form, so numbers stored as text are parsed
        if (text_matrix) as.matrix(out) else out
    })
    names(tables) <- unique(col_split)
    tables
}

assemble_prob <- function(tables, obs_split, obs_cluster, splits, genes,
                          models = attr(tables, "models")) {
    # Per-split tables -> the prob layer (n_obs x n_genes, in obs order),
    # bimodal_by_split (n_genes x n_splits), mixture, a named list of
    # n_genes x n_splits fit parameters (MIXTURE_COLS; NULL when a table lacks
    # them, NA for genes missing from a split's table), and mixture_models,
    # split -> the model it was fit in (from `models`; each split its own model
    # by default), and calls, the binarizer 1.1 data of call_mixture() tables
    # (NULL for tables without it): flags and fold_ambient (n_obs x n_genes,
    # in obs order; flags 0 and fold NA where a table has no value), stats (a
    # named list of n_genes x n_splits CALL_STAT_COLS) and params (split ->
    # the table's "binarization" attribute). Applies the rules in
    # docs/expected_input.md: stops on errors, warns on gaps.
    given <- names(tables)
    if (anyDuplicated(given)) fail("Split(s) given twice: ", first_few(given[duplicated(given)]))
    unknown <- setdiff(given, splits)
    if (length(unknown) > 0) {
        fail("Split(s) not in the clusters file: ", first_few(unknown),
            ". Its splits are: ", paste(splits, collapse = ", "))
    }
    missing <- setdiff(splits, given)
    if (length(missing) > 0) {
        fail("No probability table for split(s): ", first_few(missing),
            ". Every split of the dataset needs one.")
    }

    prob <- matrix(NA_real_, length(obs_split), length(genes))
    bimodal <- matrix(FALSE, length(genes), length(splits), dimnames = list(genes, splits))
    mixture <- lapply(MIXTURE_COLS, function(col) {
        matrix(NA_real_, length(genes), length(splits), dimnames = list(genes, splits))
    })
    has_mixture <- TRUE
    # binarizer 1.1 data: stored only when every split's table has it
    with_calls <- vapply(tables[splits], function(t) !is.null(attr(t, "flags")), logical(1))
    if (any(with_calls) && !all(with_calls)) {
        warning("Only the tables of ", first_few(splits[with_calls]), " carry the call flags ",
            "of call_mixture(); no flags are stored, so the binarizer 1.1 checks are not ",
            "applied", call. = FALSE)
    }
    calls <- NULL
    if (all(with_calls)) {
        calls <- list(
            flags = matrix(0L, length(obs_split), length(genes)),
            fold = matrix(NA_real_, length(obs_split), length(genes)),
            stats = lapply(CALL_STAT_COLS, function(col) {
                matrix(NA_real_, length(genes), length(splits), dimnames = list(genes, splits))
            }),
            params = lapply(tables[splits], attr, "binarization"))
    }
    for (s in splits) {
        tab <- tables[[s]]
        flags_s <- attr(tab, "flags")
        fold_s <- attr(tab, "fold_ambient")
        where <- paste0("split '", s, "': ")
        if (!is.data.frame(tab) && !is.matrix(tab)) {
            fail(where, "expected a data.frame with genes as row names, got ", class(tab)[1])
        }
        text_matrix <- is.matrix(tab) && is.character(tab)
        tab <- as.data.frame(tab, check.names = FALSE, stringsAsFactors = FALSE)
        if (!"exprType" %in% colnames(tab)) {
            fail(where, "no exprType column, so bimodal genes cannot be told apart")
        }
        tab_genes <- rownames(tab)
        if (is.null(tab_genes) || all(tab_genes == as.character(seq_len(nrow(tab))))) {
            fail(where, "the table needs gene symbols as row names")
        }
        extra <- setdiff(tab_genes, genes)
        if (length(extra) > 0) {
            fail(where, length(extra), " gene(s) are not in the clusters file, e.g. ",
                first_few(extra))
        }

        is_cl <- !colnames(tab) %in% QC_COLS
        cl_cols <- colnames(tab)[is_cl]
        if (anyDuplicated(cl_cols)) {
            fail(where, "cluster column(s) repeated: ", first_few(unique(cl_cols[duplicated(cl_cols)])))
        }
        rows_s <- which(obs_split == s)
        present <- obs_cluster[rows_s]
        # the legacy pipeline wrote hyphens where the object has underscores
        target <- cl_cols
        fix <- !target %in% present & gsub("-", "_", target) %in% present
        target[fix] <- gsub("-", "_", target[fix])
        if (anyDuplicated(target)) {
            fail(where, "cluster column(s) repeated after matching '-' to '_': ",
                first_few(unique(target[duplicated(target)])))
        }
        unknown <- cl_cols[!target %in% present]
        if (length(unknown) > 0) {
            fail(where, "column(s) for cluster(s) with no cells in this split: ",
                first_few(unknown), ". Columns must be cluster IDs, as in the build step.")
        }
        values <- as.matrix(tab[, is_cl, drop = FALSE])
        if (is.logical(values) && all(is.na(values))) storage.mode(values) <- "double"
        if (text_matrix) {
            # the mixture model's raw output is a character matrix: exprType makes
            # every column text, so the numbers are parsed here
            num <- suppressWarnings(as.numeric(values))
            if (any(is.na(num) & !is.na(values))) {
                fail(where, "probability columns must be numeric")
            }
            values <- matrix(num, nrow(values), dimnames = dimnames(values))
        }
        if (length(values) > 0 && !is.numeric(values)) {
            fail(where, "probability columns must be numeric")
        }
        bad <- sum(!is.na(values) & (values < 0 | values > 1))
        if (bad > 0) fail(where, bad, " probability value(s) outside [0, 1]")

        no_col <- setdiff(present, target)
        if (length(no_col) > 0) {
            warning(where, length(no_col), " cluster(s) with cells have no column and are ",
                "never ON: ", first_few(no_col), call. = FALSE)
        }
        no_gene <- setdiff(genes, tab_genes)
        if (length(no_gene) > 0) {
            warning(where, length(no_gene), " gene(s) are not in the table and are ",
                "non-bimodal in this split", call. = FALSE)
        }

        gi <- match(tab_genes, genes)
        is_bim <- tab$exprType %in% "bimodal"
        bimodal[gi, s] <- is_bim
        if (all(MIXTURE_COLS %in% colnames(tab))) {
            for (k in names(MIXTURE_COLS)) {
                mixture[[k]][gi, s] <- suppressWarnings(as.numeric(tab[[MIXTURE_COLS[[k]]]]))
            }
        } else {
            has_mixture <- FALSE
        }
        values[!is_bim, ] <- NA
        for (j in seq_along(cl_cols)) {
            prob[rows_s[match(target[j], present)], gi] <- values[, j]
        }
        if (!is.null(calls)) {
            # by name: a table can be a selection of call_mixture()'s genes
            for (m in list(flags_s, fold_s)) {
                if (!is.matrix(m) || !all(tab_genes %in% rownames(m)) ||
                    !all(cl_cols %in% colnames(m))) {
                    fail(where, "the flags or fold_ambient attribute does not cover the ",
                        "table's genes and clusters")
                }
            }
            for (j in seq_along(cl_cols)) {
                row <- rows_s[match(target[j], present)]
                calls$flags[row, gi] <- as.integer(flags_s[tab_genes, cl_cols[j]])
                calls$fold[row, gi] <- fold_s[tab_genes, cl_cols[j]]
            }
            for (k in names(CALL_STAT_COLS)) {
                col <- CALL_STAT_COLS[[k]]
                if (col %in% colnames(tab)) calls$stats[[k]][gi, s] <- as.numeric(tab[[col]])
            }
        }
    }
    if (!has_mixture) {
        warning("Not every table has the columns ", paste(MIXTURE_COLS, collapse = ", "),
            "; the mixture-model plot will be unavailable", call. = FALSE)
        mixture <- NULL
    }
    mixture_models <- setNames(splits, splits)
    known <- intersect(names(models), splits)
    mixture_models[known] <- models[known]
    list(prob = prob, bimodal_by_split = bimodal, mixture = mixture,
        mixture_models = mixture_models, calls = calls)
}

calls_record <- function(params, splits) {
    # uns/scmarco/binarization/calls: call_mixture()'s parameters per split,
    # each a vector in `splits` order (split labels may hold characters that
    # HDF5 group names cannot), and the contamination per library in long
    # form (rho_split, rho_library, rho); rho_source is "" for splits called
    # without an ambient summary
    # params: split -> the "binarization" attribute of its table
    rec <- lapply(setNames(CALL_PARAMS, CALL_PARAMS), function(k) {
        v <- lapply(params[splits], `[[`, k)
        if (any(vapply(v, is.null, logical(1)))) fail("A call_mixture() table lacks the parameter ", k)
        unlist(v, use.names = FALSE)
    })
    rec$ambient <- as.integer(rec$ambient)     # stored as 0/1
    rec$rho_source <- vapply(params[splits], function(b) {
        if (is.null(b$rho_source)) "" else b$rho_source
    }, "", USE.NAMES = FALSE)
    rho <- lapply(splits, function(s) unlist(params[[s]]$rho))
    n <- lengths(rho)
    if (sum(n) > 0) {
        rec$rho_split <- rep(splits, n)
        rec$rho_library <- unlist(lapply(rho, names), use.names = FALSE)
        rec$rho <- unlist(rho, use.names = FALSE)
    }
    rec
}

write_prob <- function(path, p, source, version) {
    # Write assemble_prob()'s result into a clusters file (add_probability(),
    # seurat_to_h5ad()): probabilities, bimodality, fit parameters, the
    # binarization record and, from call_mixture() tables, the binarizer 1.1
    # data
    record <- binarization_record("mixture model", source, version)
    if (!is.null(p$calls)) record$calls <- calls_record(p$calls$params, colnames(p$bimodal_by_split))
    add_probability_layer(path, p$prob, p$bimodal_by_split, record,
        mixture = p$mixture, mixture_models = p$mixture_models, calls = p$calls)
}
