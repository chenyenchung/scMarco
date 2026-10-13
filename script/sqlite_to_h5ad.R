#!/usr/bin/env Rscript
# Convert one table of a legacy scMarco SQLite database, a long table of
# (gene, cluster, stage, value, type in {prob, lognorm}), into
# <out-dir>/<table>.clusters.h5ad with layers prob and mean_lognorm.
#
# Usage:
#   Rscript script/sqlite_to_h5ad.R --db int/db.sqlite --table ozel_2021_new \
#     --stages P15,P30,P40,P50,P70,Adult --out-dir data \
#     [--split-label stage] [--ordered TRUE] [--genes N]
#
# --genes N keeps only the first N genes (alphabetical); use it for fixtures.
# The table is streamed in chunks, so the unindexed tables need only one scan.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})

args <- parse_cli(list(
  db = NA, table = NA, stages = NA, `out-dir` = "data",
  `split-label` = "stage", ordered = "TRUE", genes = "0", out = ""
))
stages <- strsplit(args$stages, ",")[[1]]
ordered <- as.logical(args$ordered)
max_genes <- as.integer(args$genes)

db <- dbConnect(SQLite(), args$db, flags = SQLITE_RO)
on.exit(dbDisconnect(db), add = TRUE)
tbl <- dbQuoteIdentifier(db, args$table)

message("Collecting genes and clusters...")
genes <- sort(dbGetQuery(db, paste("SELECT DISTINCT gene FROM", tbl))$gene)
if (max_genes > 0) genes <- head(genes, max_genes)
clusters <- sort(dbGetQuery(db, paste("SELECT DISTINCT cluster FROM", tbl))$cluster)
db_stages <- dbGetQuery(db, paste("SELECT DISTINCT stage FROM", tbl))$stage
if (!setequal(db_stages, stages)) {
  stop("--stages (", paste(stages, collapse = ","), ") does not match the table (",
       paste(sort(db_stages), collapse = ","), ")")
}
message(sprintf("%d genes x %d clusters x %d stages",
                length(genes), length(clusters), length(stages)))

dims <- c(length(genes), length(clusters), length(stages))
prob <- array(NA_real_, dims)
lognorm <- array(NA_real_, dims)

res <- dbSendQuery(db, paste("SELECT gene, cluster, stage, value, type FROM", tbl))
n_rows <- 0
while (!dbHasCompleted(res)) {
  chunk <- dbFetch(res, n = 2e6)
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
dbClearResult(res)

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

dir.create(args$`out-dir`, showWarnings = FALSE, recursive = TRUE)
out <- if (nzchar(args$out)) args$out else
  file.path(args$`out-dir`, paste0(args$table, ".clusters.h5ad"))
write_clusters_h5ad(
  out, obs_split, obs_cluster, stages, clusters, genes, layers,
  split_label = args$`split-label`, ordered_splits = ordered,
  extra_uns = list(
    source = paste0("sqlite:", basename(args$db), ":", args$table),
    mean_lognorm_method = "legacy scMarco cluster-average log-normalized value"
  )
)
message("Wrote ", out)
check_with_python(out)
