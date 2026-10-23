#!/usr/bin/env Rscript
# Binarize a built scMarco dataset: add the mixture-model probabilities to its
# clusters file, in place. Running it again replaces the earlier probabilities.
#
# Usage:
#   # one table per split
#   Rscript script/add_probability_h5ad.R --clusters data/x.clusters.h5ad \
#     --prob Male=male_prob.rds,Female=female_prob.rds
#
#   # one .rds holding list(Male = male_df, Female = female_df), or a wide
#   # table with "<split>@<cluster>" columns (splits modeled together)
#   Rscript script/add_probability_h5ad.R --clusters data/x.clusters.h5ad \
#     --prob prob.rds
#
#   # items mix freely: here Mel and Yak were modeled together, San alone
#   Rscript script/add_probability_h5ad.R --clusters data/x.clusters.h5ad \
#     --prob mel_yak_wide.rds,San=san_prob.rds
#
# Each table is the mixture model's output for one split: genes as row names,
# one column per cluster ID, plus QC columns; exprType marks the genes that
# are bimodal in that split. See docs/expected_input.md.

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
root <- dirname(get_script_dir())
source(file.path(get_script_dir(), "lib_write_h5ad.R"))
source(file.path(get_script_dir(), "lib_import.R"))
source(file.path(root, "src", "backend_h5ad.R"))

args <- parse_cli(list(clusters = NA, prob = NA))

ch <- h5ad_open(args$clusters)
if ("prob" %in% ch$layers) {
  message(args$clusters, " already has probabilities; replacing them")
}
splits <- levels(ch$obs$split)
tables <- read_prob_arg(args$prob)
p <- assemble_prob(tables, as.character(ch$obs$split), as.character(ch$obs$cluster),
                   splits, ch$genes)
add_probability_layer(args$clusters, p$prob, p$bimodal_by_split, list(
  method = "mixture model", source = args$prob, date = format(Sys.Date())
), mixture = p$mixture, mixture_models = p$mixture_models)
groups <- split(names(p$mixture_models), p$mixture_models)
message("Mixture models: ", paste(vapply(groups, paste, "", collapse = " + "),
                                  collapse = "; "))
n_bim <- colSums(p$bimodal_by_split)
message("Wrote probabilities to ", args$clusters, "; bimodal genes per split: ",
        paste(sprintf("%s %d", names(n_bim), n_bim), collapse = ", "))

out <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(file.path(get_script_dir(), "validate_h5ad.R"), "--clusters", args$clusters,
    "--stage", "binarized"),
  stdout = TRUE, stderr = TRUE
)
message(paste(out, collapse = "\n"))
if (!is.null(attr(out, "status")) && attr(out, "status") != 0) {
  stop("Validation failed after writing ", args$clusters, call. = FALSE)
}
check_with_python(args$clusters)
