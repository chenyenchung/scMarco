#!/usr/bin/env Rscript
# Read or change cluster names of a scMarco dataset.
#
# obs/cluster holds stable cluster IDs; obs/cluster_name holds the optional
# annotation, shown in the app as "name (id)". Only obs is rewritten, in place,
# so this takes seconds even for large files.
#
# Usage:
#   # write the current names to a table (cluster_id, name; empty = unnamed)
#   Rscript script/set_cluster_names.R --clusters data/x.clusters.h5ad --export names.tsv
#
#   # apply a table: listed clusters get the given name (empty clears it),
#   # unlisted clusters keep theirs
#   Rscript script/set_cluster_names.R --clusters data/x.clusters.h5ad --names names.tsv
#
#   # one-off migration of legacy labels: "Mi4 (154)" becomes ID 154 with name
#   # Mi4; a label without "(id)" is kept as an unnamed ID. Pass the cells file
#   # too so its cluster IDs stay in sync.
#   Rscript script/set_cluster_names.R --clusters data/x.clusters.h5ad \
#     [--cells data/x.cells.h5ad] --split-labels TRUE
#
#   # give every cluster a new ID (and name) from a table with the columns
#   # cluster_id (current ID), new_id and optional name (empty = unnamed);
#   # every current cluster must be listed. Rows and data are unchanged.
#   Rscript script/set_cluster_names.R --clusters data/x.clusters.h5ad \
#     [--cells data/x.cells.h5ad] --ids ids.tsv

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(get_script_dir())
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(
  clusters = NA, cells = "", names = "", export = "", ids = "", `split-labels` = "FALSE"
))
split_labels <- as.logical(args$`split-labels`)
modes <- c(export = nzchar(args$export), names = nzchar(args$names),
           split = isTRUE(split_labels), ids = nzchar(args$ids))
if (sum(modes) != 1) {
  stop("Give exactly one of --export, --names, --ids or --split-labels TRUE")
}
if (nzchar(args$cells) && !modes[["split"]] && !modes[["ids"]]) {
  stop("--cells is only used with --split-labels or --ids (cells files store IDs only)")
}

ch <- h5ad_open(args$clusters)
clusters <- levels(ch$obs$cluster)
current <- read_cluster_names(ch$obs, clusters, args$clusters)
# Names are stored per obs row, so a cluster absent from every split cannot
# hold one
no_rows <- setdiff(clusters, as.character(ch$obs$cluster))

names_obs <- function(obs, names) {
  #' obs$cluster_name from an id -> name vector; dropped when nothing is named
  obs$cluster_name <- NULL
  if (any(!is.na(names))) {
    per_row <- unname(names[as.character(obs$cluster)])
    obs$cluster_name <- factor(per_row, levels = sort(unique(names[!is.na(names)])))
  }
  obs
}

# --- Export ------------------------------------------------------------------
if (modes[["export"]]) {
  out <- data.frame(cluster_id = clusters, name = unname(current))
  write.table(out, args$export, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message(sprintf("Wrote %d clusters (%d named) to %s",
                  nrow(out), sum(!is.na(current)), args$export))
  quit(status = 0)
}

# --- Apply a names table -----------------------------------------------------
if (modes[["names"]]) {
  tbl <- read.delim(args$names, colClasses = "character", na.strings = character(0),
                    check.names = FALSE, quote = "")
  if (!all(c("cluster_id", "name") %in% names(tbl))) {
    stop(args$names, " needs the columns cluster_id and name (tab-separated)")
  }
  tbl$cluster_id <- trimws(tbl$cluster_id)
  tbl$name <- trimws(tbl$name)
  unknown <- setdiff(tbl$cluster_id, clusters)
  if (length(unknown) > 0) {
    stop(length(unknown), " cluster ID(s) not in ", args$clusters, ": ",
         paste(head(unknown, 10), collapse = ", "))
  }
  dup <- unique(tbl$cluster_id[duplicated(tbl$cluster_id)])
  if (length(dup) > 0) {
    stop("Cluster ID(s) listed more than once: ", paste(head(dup, 10), collapse = ", "))
  }
  homeless <- intersect(tbl$cluster_id[tbl$name != ""], no_rows)
  if (length(homeless) > 0) {
    stop("Cluster(s) with no rows in any split cannot be named: ",
         paste(head(homeless, 10), collapse = ", "))
  }
  updated <- current
  updated[tbl$cluster_id] <- ifelse(tbl$name == "", NA, tbl$name)
  changed <- names(updated)[!mapply(identical, updated, current)]
  shared <- unique(updated[!is.na(updated) & duplicated(updated)])
  if (length(shared) > 0) {
    message("Note: name(s) shared by several clusters (kept apart by their IDs): ",
            paste(head(shared, 5), collapse = ", "))
  }
  if (length(changed) == 0) {
    message("No changes")
    quit(status = 0)
  }
  rewrite_obs(args$clusters, function(obs) names_obs(obs, updated))
  for (id in head(changed, 20)) {
    message(sprintf("  %s: %s -> %s", id,
                    if (is.na(current[[id]])) "(unnamed)" else current[[id]],
                    if (is.na(updated[[id]])) "(unnamed)" else updated[[id]]))
  }
  if (length(changed) > 20) message("  ... and ", length(changed) - 20, " more")
  message(sprintf("Updated %d cluster(s); %d of %d now named",
                  length(changed), sum(!is.na(updated)), length(updated)))
  quit(status = 0)
}

# --- New IDs: split legacy "name (id)" labels, or apply an ID table --------
if (modes[["split"]]) {
  if (any(!is.na(current))) {
    stop(args$clusters, " already has cluster names; --split-labels is for ",
         "files whose cluster labels still embed them")
  }
  parts <- split_cluster_labels(clusters)
} else {
  tbl <- read.delim(args$ids, colClasses = "character", na.strings = character(0),
                    check.names = FALSE, quote = "")
  if (!all(c("cluster_id", "new_id") %in% names(tbl))) {
    stop(args$ids, " needs the columns cluster_id and new_id (tab-separated)")
  }
  if (!"name" %in% names(tbl)) tbl$name <- ""
  tbl[] <- lapply(tbl, trimws)
  unknown <- setdiff(tbl$cluster_id, clusters)
  missing <- setdiff(clusters, tbl$cluster_id)
  if (length(unknown) > 0 || length(missing) > 0 || anyDuplicated(tbl$cluster_id)) {
    stop(args$ids, " must list every cluster ID of ", args$clusters, " exactly once",
         if (length(unknown) > 0) paste0("\nUnknown: ", paste(head(unknown, 10), collapse = ", ")),
         if (length(missing) > 0) paste0("\nMissing: ", paste(head(missing, 10), collapse = ", ")))
  }
  if (any(tbl$new_id == "")) stop("Empty new_id for: ", paste(tbl$cluster_id[tbl$new_id == ""], collapse = ", "))
  tbl <- tbl[match(clusters, tbl$cluster_id), ]
  parts <- data.frame(label = tbl$cluster_id, id = tbl$new_id,
                      name = ifelse(tbl$name == "", NA, tbl$name), stringsAsFactors = FALSE)
}
dup <- unique(parts$id[duplicated(parts$id)])
if (length(dup) > 0) {
  stop("Duplicate new cluster IDs: ", paste(head(dup, 10), collapse = ", "),
       "\nOld IDs: ", paste(parts$label[parts$id %in% dup], collapse = ", "))
}
if (any(grepl("|", parts$id, fixed = TRUE))) stop("Cluster IDs must not contain '|'")
lost <- parts$label[!is.na(parts$name) & parts$label %in% no_rows]
if (length(lost) > 0) {
  warning(length(lost), " cluster(s) have no rows in any split, so only their ",
          "ID is kept: ", paste(lost, collapse = ", "), call. = FALSE, immediate. = TRUE)
}
id_of <- setNames(parts$id, parts$label)
new_names <- setNames(parts$name, parts$id)
id_levels <- order_cluster_ids(parts$id)

xh <- NULL
if (nzchar(args$cells)) {
  xh <- h5ad_open(args$cells)
  stray <- setdiff(unique(as.character(xh$obs$cluster)), clusters)
  if (length(stray) > 0) {
    stop("Cells file has clusters missing from the clusters file: ",
         paste(head(stray, 5), collapse = ", "))
  }
}
if (modes[["split"]] && all(parts$id == parts$label)) {
  message("No label embeds an ID; nothing to change")
  quit(status = 0)
}

relabel <- function(obs) {
  obs$cluster <- factor(unname(id_of[as.character(obs$cluster)]), levels = id_levels)
  obs
}
rewrite_obs(args$clusters, function(obs) {
  obs <- names_obs(relabel(obs), new_names)
  rownames(obs) <- paste(obs$split, obs$cluster, sep = "|")
  obs
})
if (!is.null(xh)) rewrite_obs(args$cells, relabel)

named <- !is.na(parts$name) & !parts$label %in% lost
if (modes[["split"]]) {
  message(sprintf("Split %d labels into ID + name; %d bare IDs kept as unnamed", sum(named), sum(!named)))
} else {
  message(sprintf("Gave %d clusters new IDs; %d named, %d unnamed", nrow(parts), sum(named), sum(!named)))
}
for (i in head(which(named), 5)) {
  message(sprintf("  %s -> ID %s, name %s", parts$label[i], parts$id[i], parts$name[i]))
}
check_with_python(args$clusters)
