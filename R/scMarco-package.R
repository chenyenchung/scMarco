#' scMarco: marker gene combinations from binarized single-cell expression
#'
#' scMarco helps choose marker gene combinations (for example split-GAL4
#' driver pairs) that label a cluster of a single-cell RNA-seq dataset and as
#' few other clusters as possible. Expression is binarized per cluster with a
#' mixture model, so a gene is ON, OFF or ambiguous in each cluster and in each
#' split of the dataset (developmental stage, sex, species, ...).
#'
#' The package has three parts:
#'
#' * an interactive Shiny app, [run_scmarco()], configured by a YAML file;
#' * the combination logic behind the app, usable from R on a dataset loaded
#'   with [ds_load()] and binarized with [ds_binarize()]: [from_cluster()],
#'   [clusters_expressing()], [distinct_within()], [pair_offtargets()],
#'   [find_pairs()] and the plots in [onoff_plot()];
#' * tools that build and edit datasets in the scMarco h5ad format:
#'   [seurat_to_h5ad()], [add_probability()], [validate_h5ad()],
#'   [set_cluster_names()], [set_cluster_metadata()] and the legacy
#'   converters [sqlite_to_h5ad()] and [h5v1_to_h5ad()].
#'
#' See `vignette("scMarco")` for a walk-through.
#'
#' @import shiny
#' @import ggplot2
#' @importFrom rhdf5 h5read h5readAttributes h5ls H5Fopen H5Fclose H5Gopen
#'   H5Gclose H5Lexists
#' @importFrom anndataR AnnData read_h5ad write_h5ad
#' @importFrom methods as
#' @importFrom stats setNames aggregate ave dnorm quantile plogis qlogis
#' @importFrom markdown mark_html
#' @importFrom utils head read.delim write.csv write.table
#' @keywords internal
"_PACKAGE"

# Columns referred to in ggplot2 aesthetics
utils::globalVariables(c(
    "cluster", "component", "dim1", "dim2", "label", "lower", "middle", "state",
    "upper", "value", "x", "y", "ymax", "ymin"
))
