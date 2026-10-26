# Names table for set_cluster_names() from the cluster
# annotation of Simon et al. (2026, Nat. Neurosci., doi:10.1038/s41593-026-02256-6),
# Supplementary Table 41593_2026_2256_MOESM10_ESM.xlsx (sheet "Annotation":
# Cluster_number, Annotation, ...).
#
# The cluster numbers are the Ozel et al. (2021) cluster IDs used by the
# submission objects (FinalIdents) and the mixture-model tables, so every
# dataset cluster is looked up by ID. A cluster whose annotation is its own
# number stays unnamed; clusters missing from the table are an error.
#
# Usage (needs readxl):
#   library(scMarco)
#   source("misc/moesm10_cluster_names.R")
#   names <- moesm10_cluster_names("41593_2026_2256_MOESM10_ESM.xlsx", "data/x.clusters.h5ad")
#   set_cluster_names("data/x.clusters.h5ad", names)

moesm10_cluster_names <- function(xlsx, clusters, out = NULL) {
  #' data.frame(cluster_id, name) for every cluster of `clusters` ("" =
  #' unnamed); with `out`, also written there as a TSV.
  tbl <- as.data.frame(readxl::read_excel(xlsx, sheet = "Annotation", col_types = "text"))
  stopifnot(all(c("Cluster_number", "Annotation") %in% names(tbl)))
  id <- trimws(tbl$Cluster_number)
  ann <- trimws(tbl$Annotation)
  ann[is.na(ann)] <- ""
  if (anyDuplicated(id)) stop("Cluster number(s) listed twice: ", paste(unique(id[duplicated(id)]), collapse = ", "))
  ann <- setNames(ifelse(ann == id, "", ann), id)

  ids <- levels(scMarco:::h5ad_open(clusters)$obs$cluster)
  missing <- setdiff(ids, names(ann))
  if (length(missing) > 0) stop("Clusters missing from ", xlsx, ": ", paste(missing, collapse = ", "))

  names_df <- data.frame(cluster_id = ids, name = unname(ann[ids]))
  if (!is.null(out)) write.table(names_df, out, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%d clusters, %d named; %d table rows not in the dataset (%s)",
                  nrow(names_df), sum(names_df$name != ""),
                  length(setdiff(names(ann), ids)),
                  paste(setdiff(names(ann), ids), collapse = ", ")))
  names_df
}
