# Add uns/scmarco/binarization (with its binarizer version) to binarized
# format_version 2.0 clusters files written before v2 required it. Only
# uns/scmarco is rewritten, in place.
#
# Every file binarized before the requirement got its tables from outside
# scMarco's binarizer (the mixture_modeling_smk pipeline, migrated legacy
# databases, hand-made fixtures), or from add_probability() stamping the
# package's version whatever the tables were, so the version recorded is
# "external". To record a real version, rerun add_probability() with tables
# from call_mixture() instead.
#
#   Rscript deploy/misc/stamp_binarizer_version.R data/*.clusters.h5ad

suppressMessages(library(scMarco))

stamp <- function(path) {
    ch <- scMarco:::h5ad_open(path)
    meta <- ch$uns$scmarco
    if (!identical(as.character(meta$format_version), scMarco:::SCMARCO_FORMAT_VERSION)) {
        stop(path, ": format_version ", meta$format_version, ", not ",
            scMarco:::SCMARCO_FORMAT_VERSION)
    }
    if (!"prob" %in% ch$layers) {
        message(path, ": not binarized, nothing to do")
        return(invisible(FALSE))
    }
    b <- meta$binarization
    if (!is.null(b$version) && !identical(as.character(b$version), "1.0")) {
        message(path, ": already records binarizer ", b$version)
        return(invisible(FALSE))
    }
    pick <- function(x, default) if (is.null(x)) default else as.character(x)
    record <- list(
        method = pick(b$method, "mixture model"),
        source = pick(b$source, pick(meta$source, "unknown")),
        date = pick(b$date, "unknown"),
        version = scMarco:::BINARIZER_EXTERNAL
    )
    ad <- suppressWarnings(anndataR::read_h5ad(path, as = "HDF5AnnData", mode = "r+"))
    on.exit(ad$close())
    uns <- ad$uns
    uns$scmarco$binarization <- record
    ad$uns <- uns
    message(path, ": binarizer version ", record$version, " (source ", record$source, ")")
    invisible(TRUE)
}

for (path in commandArgs(trailingOnly = TRUE)) stamp(path)
