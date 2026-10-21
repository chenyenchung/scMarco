#!/usr/bin/env Rscript
# Table for script/set_cluster_names.R --ids that gives the Ozel_2021 and
# Kurmangaliyev_2020 clusters the numeric IDs and annotations of the PNAS
# cluster table (pnas.2614314123.sd01.xlsx, sheet 1: Cluster, Annotation).
#
# Those datasets were built with the objects' Idents() (= MultiomeAnnotated)
# as cluster ID: annotation names, with bare numbers for unnamed clusters. The
# numbers are the MultiomeNN clusters, which the PNAS table numbers 1-259. A
# named ident is the NN cluster(s) whose cells it holds the majority of:
# usually one (Dm4 = 4), two for merged groups (G_PsCrtrg = 1+2), none for the
# P15 developmental groups (NB, GMC1, Im01, ...), which keep their label as ID.
# The grouping itself is unchanged, so the mixture-model probabilities stay valid.
#
# Usage:
#   Rscript script/misc/pnas_cluster_ids.R --pnas pnas.2614314123.sd01.xlsx \
#     --objects a.rds,b.rds,... --clusters data/x.clusters.h5ad --out ids.tsv
#
#   --objects  Seurat objects, or their meta.data saved as .rds (faster); the
#              NN mapping is pooled over all of them, so pass every object
#              the PNAS clusters were called on, not only those of one dataset.
#   Needs readxl.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(dirname(get_script_dir()))
source(file.path(root, "script", "lib_write_h5ad.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(pnas = NA, objects = NA, clusters = NA, out = NA))
is_num <- function(x) grepl("^[0-9]+$", x)

# --- PNAS annotation per cluster number --------------------------------------
tbl <- as.data.frame(readxl::read_excel(args$pnas, sheet = 1, col_types = "text",
                                       .name_repair = "minimal"))
ann <- trimws(tbl[[grep("^Annotation", names(tbl))]])
ann[is.na(ann)] <- ""
ann <- setNames(ann, trimws(tbl[[1]]))
stopifnot(all(is_num(names(ann))), !anyDuplicated(names(ann)))

# --- ident x MultiomeNN cell counts, pooled over objects ---------------------
counts <- NULL
for (f in strsplit(args$objects, ",")[[1]]) {
  md <- readRDS(f)
  if (!is.data.frame(md)) md <- md@meta.data
  for (col in c("MultiomeAnnotated", "MultiomeNN")) {
    if (!col %in% names(md)) stop(f, " has no ", col, " column")
  }
  tb <- table(ident = as.character(md[["MultiomeAnnotated"]]),
              nn = as.character(md[["MultiomeNN"]]))
  tb <- as.data.frame(tb, stringsAsFactors = FALSE)
  counts <- rbind(counts, tb[tb$Freq > 0, ])
  message(basename(f), ": ", nrow(md), " cells")
}
counts <- aggregate(Freq ~ ident + nn, counts, sum)
stray <- setdiff(counts$nn, names(ann))
if (length(stray) > 0) stop("MultiomeNN values missing from the PNAS table: ", paste(stray, collapse = ", "))

# each NN cluster belongs to the ident holding most of its cells
owner <- sapply(split(counts, counts$nn), function(d) d$ident[which.max(d$Freq)])
bad <- names(owner)[is_num(owner) & owner != names(owner)]
if (length(bad) > 0) stop("NN cluster(s) owned by another numeric ident: ", paste(bad, collapse = ", "))

# --- New ID and name per current cluster -------------------------------------
clusters <- levels(h5ad_open(args$clusters)$obs$cluster)
unknown <- setdiff(clusters, counts$ident)
if (length(unknown) > 0) stop("Clusters not found in the objects: ", paste(unknown, collapse = ", "))

rows <- lapply(clusters, function(id) {
  nns <- names(owner)[owner == id]
  nns <- nns[order(as.integer(nns))]
  if (length(nns) == 0) return(data.frame(cluster_id = id, new_id = id, name = ""))
  d <- counts[counts$ident == id, ]
  covered <- sum(d$Freq[d$nn %in% nns]) / sum(d$Freq)
  if (covered < 0.99) {
    stop(sprintf("%s: only %.1f%% of its cells are in NN cluster(s) %s",
                 id, 100 * covered, paste(nns, collapse = "+")))
  }
  nm <- unique(ann[nns][ann[nns] != ""])
  data.frame(cluster_id = id, new_id = paste(nns, collapse = "+"),
             name = paste(nm, collapse = "/"))
})
out <- do.call(rbind, rows)
write.table(out, args$out, sep = "\t", quote = FALSE, row.names = FALSE)

merged <- out[grepl("+", out$new_id, fixed = TRUE), ]
kept <- out[out$new_id == out$cluster_id & !is_num(out$cluster_id), ]
message(sprintf("Wrote %s: %d clusters, %d named; %d merged NN groups (%s); %d kept as is (%s)",
                args$out, nrow(out), sum(out$name != ""), nrow(merged),
                paste(merged$cluster_id, merged$new_id, sep = " = ", collapse = ", "),
                nrow(kept), paste(kept$cluster_id, collapse = ", ")))
