#!/usr/bin/env Rscript
# Rank marker gene pairs for every cluster of a dataset (replaces the legacy
# find_markers_h5.R brute-force search with one matrix product per cluster).
#
# Usage:
#   Rscript script/find_pairs.R --clusters data/x.clusters.h5ad \
#     --out results/x_pairs.csv [--cut 0.5] [--required Male,Female] \
#     [--min-n 1] [--top 20]
#
# Output columns: cluster (ID), cluster_name (empty when unnamed), gene1, gene2,
# offtargets (other clusters in which both genes are ON in any split),
# single_offtargets; rows sorted best first per cluster.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(get_script_dir())
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
source(file.path(root, "src", "combo_core.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(
  clusters = NA, out = NA, cut = "0.5", required = "", `min-n` = "1", top = "20"
))
ds <- ds_load("cli", list(label = "cli", clusters = args$clusters))
required <- if (nzchar(args$required)) strsplit(args$required, ",")[[1]] else character(0)
b <- ds_binarize(ds, as.numeric(args$cut))

results <- lapply(ds$clusters, function(cl) {
  res <- pair_offtargets(b, cl, required, as.integer(args$`min-n`),
                         limit = as.integer(args$top), genes = ds$analysis_genes)
  if (nrow(res) == 0) return(NULL)
  cbind(cluster = cl, cluster_name = unname(ds$cluster_names[cl]), res)
})
out <- do.call(rbind, results)
dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)
write.csv(out, args$out, row.names = FALSE)
message(sprintf("Wrote %d pairs for %d clusters to %s",
                if (is.null(out)) 0 else nrow(out),
                length(unique(out$cluster)), args$out))
