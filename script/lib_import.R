# Input checks and readers for building a dataset from new data: split and
# cluster labels from a Seurat object, and mixture-model probability tables.
# Sourced by seurat_to_h5ad.R and add_probability_h5ad.R; see
# docs/expected_input.md for the formats and rules.

QC_COLS <- c("diff.elpd", "diff.elpd.se", "bimodal.on_level", "bimodal.off_level",
             "bimodal.sd", "bimodal.pi", "unimodal.level", "unimodal.sd", "exprType")

fail <- function(...) stop(paste0(...), call. = FALSE)

first_few <- function(x, n = 5) {
  paste0(paste(head(x, n), collapse = ", "), if (length(x) > n) ", ..." else "")
}

# ---------------------------------------------------------------------------
# Split and cluster labels
# ---------------------------------------------------------------------------

seurat_labels <- function(meta, ident, split_col = "scmarco_split",
                          cluster_col = "scmarco_cluster", defaults = TRUE) {
  #' Split and cluster label of every cell.
  #' meta: the object's meta.data; ident: as.character(Idents(obj)).
  #' split_col "none" puts every cell in one split named "all".
  #' defaults: TRUE when both column names are the defaults, so the legacy
  #'   "<split>@<cluster>" identities may be used when neither column exists.
  has <- function(col) col %in% colnames(meta)
  if (!has(cluster_col) && !(split_col != "none" && has(split_col)) && defaults) {
    message("No '", cluster_col, "' / '", split_col, "' metadata columns; ",
            "reading legacy '<split>@<cluster>' identities")
    parsed <- strsplit(ident, "@", fixed = TRUE)
    bad <- lengths(parsed) != 2
    if (any(bad)) {
      fail("Found neither metadata columns '", split_col, "' and '", cluster_col,
           "' nor identities of the form '<split>@<cluster>' (e.g. ",
           first_few(unique(ident[bad]), 3), ").\nAdd the columns, e.g. ",
           "obj$", split_col, " <- obj$sex; obj$", cluster_col,
           " <- obj$seurat_clusters, or name yours with --split-col / --cluster-col.")
    }
    return(list(split = vapply(parsed, `[`, "", 1),
                cluster = vapply(parsed, `[`, "", 2)))
  }
  if (!has(cluster_col)) {
    fail("Metadata column '", cluster_col, "' not found. Columns: ",
         first_few(colnames(meta), 10), ". Name the cluster column with --cluster-col.")
  }
  cluster <- meta[[cluster_col]]
  if (split_col == "none") {
    split <- rep("all", nrow(meta))
  } else if (!has(split_col)) {
    fail("Metadata column '", split_col, "' not found. Name the split column with ",
         "--split-col, or pass --split-col none for a dataset with one split.")
  } else {
    split <- meta[[split_col]]
  }
  list(split = as.character(split), cluster = as.character(cluster))
}

retag_seurat_v3 <- function(x) {
  #' Seurat v3 objects saved before SeuratObject existed name their classes'
  #' package "Seurat"; point them to SeuratObject so it can load them without
  #' Seurat. Walks S4 slots and plain lists only. UpdateSeuratObject() and
  #' SeuratObject's internal UpdateClassPkg() both look the classes up first
  #' and fail with "unable to load required package 'Seurat'" (SeuratObject 5.4.0).
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
  #' readRDS() a Seurat object; objects from before SeuratObject (Seurat v3)
  #' are retagged and updated so SeuratObject alone can read them.
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

parse_seurat_arg <- function(arg) {
  #' --seurat value -> object path(s). "obj.rds": one object, unnamed.
  #' "P15=a.rds,P30=b.rds": one object per split, named by split.
  if (!grepl("=", arg, fixed = TRUE)) {
    if (!file.exists(arg)) fail("Seurat object not found: ", arg)
    return(arg)
  }
  pairs <- strsplit(strsplit(arg, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
  if (any(lengths(pairs) != 2)) {
    fail("--seurat must be one .rds or look like Split1=file1.rds,Split2=file2.rds; got ", arg)
  }
  files <- vapply(pairs, `[`, "", 2)
  names(files) <- vapply(pairs, `[`, "", 1)
  if (anyDuplicated(names(files))) {
    fail("Split(s) given twice in --seurat: ", first_few(unique(names(files)[duplicated(names(files))])))
  }
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) fail("Seurat object(s) not found: ", first_few(missing))
  files
}

check_labels <- function(split, cluster) {
  #' Name rules for split and cluster labels (stop on the first broken rule)
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

read_prob_arg <- function(arg) {
  #' --prob value -> named list split -> table (genes x clusters + QC columns).
  #' "Male=a.rds,Female=b.rds": one file per split. "x.rds": a named list of
  #' tables, or the legacy wide data.frame with "<split>@<cluster>" columns.
  if (grepl("=", arg, fixed = TRUE)) {
    pairs <- strsplit(strsplit(arg, ",", fixed = TRUE)[[1]], "=", fixed = TRUE)
    if (any(lengths(pairs) != 2)) {
      fail("--prob must look like Split1=file1.rds,Split2=file2.rds; got ", arg)
    }
    splits <- vapply(pairs, `[`, "", 1)
    files <- vapply(pairs, `[`, "", 2)
    missing <- files[!file.exists(files)]
    if (length(missing) > 0) fail("Probability file(s) not found: ", first_few(missing))
    tables <- lapply(files, readRDS)
    names(tables) <- splits
    return(tables)
  }
  if (!file.exists(arg)) fail("Probability file not found: ", arg)
  obj <- readRDS(arg)
  if (is.data.frame(obj)) {
    message("Reading a legacy wide table with '<split>@<cluster>' columns")
    return(legacy_wide_to_split_tables(obj))
  }
  if (!is.list(obj) || is.null(names(obj)) || any(names(obj) == "")) {
    fail(arg, " must hold a named list of per-split tables, e.g. ",
         "list(Male = male_df, Female = female_df), or a legacy wide data.frame")
  }
  obj
}

legacy_wide_to_split_tables <- function(df) {
  #' Legacy wide table -> per-split tables. Its single exprType (and QC
  #' columns) apply to every split.
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
    out
  })
  names(tables) <- unique(col_split)
  tables
}

assemble_prob <- function(tables, obs_split, obs_cluster, splits, genes) {
  #' Per-split tables -> the prob layer (n_obs x n_genes, in obs order) and
  #' bimodal_by_split (n_genes x n_splits). Applies the rules in
  #' docs/expected_input.md: stops on errors, warns on gaps.
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
  for (s in splits) {
    tab <- tables[[s]]
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
    values[!is_bim, ] <- NA
    for (j in seq_along(cl_cols)) {
      prob[rows_s[match(target[j], present)], gi] <- values[, j]
    }
  }
  list(prob = prob, bimodal_by_split = bimodal)
}
