#!/usr/bin/env Rscript
# Names table for script/set_cluster_names.R --names from the cluster
# annotation of Simon et al. (2026, Nat. Neurosci., doi:10.1038/s41593-026-02256-6),
# Supplementary Table 41593_2026_2256_MOESM10_ESM.xlsx (sheet "Annotation":
# Cluster_number, Annotation, ...).
#
# The cluster numbers are the Ozel et al. (2021) cluster IDs used by the
# submission objects (FinalIdents) and the mixture-model tables, so every
# dataset cluster is looked up by ID. A cluster whose annotation is its own
# number stays unnamed; clusters missing from the table are an error.
#
# Usage:
#   Rscript script/misc/moesm10_cluster_names.R --xlsx 41593_2026_2256_MOESM10_ESM.xlsx \
#     --clusters data/x.clusters.h5ad --out names.tsv
#
#   Needs readxl.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(dirname(get_script_dir()))
source(file.path(root, "script", "lib_write_h5ad.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(xlsx = NA, clusters = NA, out = NA))

tbl <- as.data.frame(readxl::read_excel(args$xlsx, sheet = "Annotation", col_types = "text"))
stopifnot(all(c("Cluster_number", "Annotation") %in% names(tbl)))
id <- trimws(tbl$Cluster_number)
ann <- trimws(tbl$Annotation)
ann[is.na(ann)] <- ""
if (anyDuplicated(id)) stop("Cluster number(s) listed twice: ", paste(unique(id[duplicated(id)]), collapse = ", "))
ann <- setNames(ifelse(ann == id, "", ann), id)

clusters <- levels(h5ad_open(args$clusters)$obs$cluster)
missing <- setdiff(clusters, names(ann))
if (length(missing) > 0) stop("Clusters missing from ", args$xlsx, ": ", paste(missing, collapse = ", "))

out <- data.frame(cluster_id = clusters, name = unname(ann[clusters]))
write.table(out, args$out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %s: %d clusters, %d named; %d table rows not in the dataset (%s)",
                args$out, nrow(out), sum(out$name != ""),
                length(setdiff(names(ann), clusters)),
                paste(setdiff(names(ann), clusters), collapse = ", ")))
