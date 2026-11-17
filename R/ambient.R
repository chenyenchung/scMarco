# Ambient RNA for the binarizer 1.1 calling rules: per-cluster counts (UMIs,
# cells detecting each gene, UMIs per library) and each library's
# contamination, estimated with DecontX (celda) from the filtered counts and
# the cluster labels, so no empty droplets are needed. The expected ambient
# UMIs of a gene in a cluster follow from them:
#   lambda[g, c] = sum over libraries l of rho[l] * share[g, l] * umi[c, l]
# rho[l]: contamination of library l (share of its UMIs from ambient RNA);
# share[g, l]: gene g's share of library l's UMIs (the library's ambient
# profile); umi[c, l]: cluster c's UMIs from library l.

as_count_matrix <- function(counts) {
    # Raw counts, genes x cells, as a CsparseMatrix (checked)
    if (is.data.frame(counts)) counts <- as.matrix(counts)
    if (!(is.matrix(counts) && is.numeric(counts)) && !methods::is(counts, "Matrix")) {
        fail("counts must be a numeric matrix or a sparse Matrix, genes x cells")
    }
    counts <- methods::as(counts, "CsparseMatrix")
    if (is.null(rownames(counts))) fail("counts needs genes as row names")
    if (any(counts@x < 0)) fail("counts must not be negative")
    if (any(counts@x != round(counts@x))) {
        fail("counts must be raw counts (whole numbers), not normalized values")
    }
    counts
}

#' Summarize counts per cluster and estimate ambient RNA
#'
#' The per-cluster counts and the ambient-RNA estimate that [call_mixture()]
#' uses for its detection floor (ON needs enough cells detecting the gene) and
#' its ambient check (OFF needs the cluster's UMIs to be consistent with
#' ambient RNA). Each library's contamination is estimated with
#' [celda::decontX()] from the filtered counts and the cluster labels (one
#' batch per library), so no empty droplets are needed; or give it as `rho`,
#' e.g. from SoupX or CellBender.
#'
#' The expected ambient UMIs of gene g in cluster c are
#' \eqn{\sum_l \rho_l s_{gl} u_{cl}}{sum over libraries l of rho[l] x
#' share[g, l] x umi[c, l]}: library l's contamination, times gene g's share
#' of library l's UMIs (its ambient profile), times cluster c's UMIs from
#' library l.
#'
#' @param counts Raw counts of one split's cells, genes x cells (a matrix or
#'   sparse Matrix with gene row names). Use every gene, not only the ones
#'   fitted: the totals and DecontX need them all.
#' @param cluster Cluster ID of each cell (column of `counts`).
#' @param library 10x library (channel) of each cell; `NULL` puts every cell
#'   in one library.
#' @param rho Contamination per library, skipping DecontX: one number for
#'   every library, or a vector named by library. `NULL` estimates it.
#' @param genes Genes to keep the per-cluster counts of (default: all).
#'
#' @return An `scmarco_ambient` list: `umi` and `detect` (genes x clusters:
#'   UMIs, and share of the cluster's cells with at least one UMI), `n_cells`
#'   (per cluster), `cluster_library_umi` (clusters x libraries: UMIs of all
#'   genes), `library_share` (genes x libraries: the gene's share of the
#'   library's UMIs), `rho` (per library), `rho_source` (`"decontX"` with
#'   the celda version, or `"given"`), `genes`, `clusters` and `libraries`.
#' @seealso [call_mixture()], [make_pseudobulk()].
#' @references Yang S, Corbett SE, Koga Y, et al. (2020). Decontamination of
#'   ambient RNA in single-cell RNA-seq with DecontX. Genome Biology 21:57.
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(20 * 40, 2), 20, 40,
#'     dimnames = list(paste0("g", 1:20), paste0("c", 1:40)))
#' amb <- ambient_summary(counts, cluster = rep(c("A", "B"), each = 20),
#'     library = rep(c("L1", "L2"), 20), rho = 0.05)
#' amb$rho
#' dim(amb$umi)
#' @export
ambient_summary <- function(counts, cluster, library = NULL, rho = NULL, genes = NULL) {
    counts <- as_count_matrix(counts)
    n <- ncol(counts)
    if (length(cluster) != n) fail("cluster has ", length(cluster), " entries for ", n, " cells")
    if (anyNA(cluster)) fail("cluster has missing values")
    if (is.null(library)) library <- rep("1", n)
    if (length(library) != n) fail("library has ", length(library), " entries for ", n, " cells")
    if (anyNA(library)) fail("library has missing values")
    if (is.null(genes)) genes <- rownames(counts)
    unknown <- setdiff(genes, rownames(counts))
    if (length(unknown) > 0) fail("Gene(s) not in counts: ", first_few(unknown))
    cluster <- as.character(cluster)
    library <- as.character(library)
    clusters <- unique(cluster)
    libraries <- sort(unique(library))

    total <- Matrix::colSums(counts)                       # every gene
    ind <- Matrix::sparseMatrix(i = seq_len(n), j = match(cluster, clusters), x = 1,
        dims = c(n, length(clusters)))
    lib_ind <- Matrix::sparseMatrix(i = seq_len(n), j = match(library, libraries), x = 1,
        dims = c(n, length(libraries)))
    n_cells <- setNames(Matrix::colSums(ind), clusters)
    sub <- counts[genes, , drop = FALSE]
    umi <- as.matrix(sub %*% ind)
    detect <- sweep(as.matrix(((sub > 0) * 1) %*% ind), 2, n_cells, "/")
    dimnames(umi) <- dimnames(detect) <- list(genes, clusters)
    cluster_library_umi <- as.matrix(Matrix::crossprod(ind, lib_ind * total))
    dimnames(cluster_library_umi) <- list(clusters, libraries)
    lib_total <- Matrix::colSums(lib_ind * total)
    library_share <- sweep(as.matrix(sub %*% lib_ind), 2, lib_total, "/")
    dimnames(library_share) <- list(genes, libraries)

    if (is.null(rho)) {
        need_packages("celda", "ambient_summary()")
        keep <- Matrix::rowSums(counts) > 0
        res <- celda::decontX(counts[keep, , drop = FALSE], z = cluster, batch = library,
            verbose = FALSE)
        cont <- res$contamination
        rho <- vapply(libraries, function(l) {
            i <- library == l
            sum(cont[i] * total[i]) / sum(total[i])
        }, numeric(1))
        rho_source <- paste0("decontX (celda ", utils::packageVersion("celda"), ")")
    } else {
        if (!is.numeric(rho) || anyNA(rho) || any(rho < 0 | rho >= 1)) {
            fail("rho must be contamination fractions in [0, 1)")
        }
        if (length(rho) == 1 && is.null(names(rho))) {
            rho <- setNames(rep(rho, length(libraries)), libraries)
        }
        missing_lib <- setdiff(libraries, names(rho))
        if (length(missing_lib) > 0) fail("rho has no value for library ", first_few(missing_lib))
        rho <- rho[libraries]
        rho_source <- "given"
    }
    structure(list(umi = umi, detect = detect, n_cells = n_cells,
        cluster_library_umi = cluster_library_umi, library_share = library_share,
        rho = rho, rho_source = rho_source, genes = genes, clusters = clusters,
        libraries = libraries), class = "scmarco_ambient")
}

ambient_expected <- function(amb, genes = amb$genes, clusters = amb$clusters) {
    # Expected ambient UMIs, genes x clusters (see the header)
    lam <- amb$library_share[genes, , drop = FALSE] %*%
        (t(amb$cluster_library_umi[clusters, , drop = FALSE]) * amb$rho)
    dimnames(lam) <- list(genes, clusters)
    lam
}

#' @export
print.scmarco_ambient <- function(x, ...) {
    cat(sprintf("scMarco ambient summary: %d gene(s), %d cluster(s), %d library(ies)\n",
        length(x$genes), length(x$clusters), length(x$libraries)))
    cat("contamination per library (", x$rho_source, "): ",
        paste(sprintf("%s %.3f", names(x$rho), x$rho), collapse = ", "), "\n", sep = "")
    invisible(x)
}
