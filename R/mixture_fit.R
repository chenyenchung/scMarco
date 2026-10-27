# Mixture-model binarizer (version SCMARCO_BINARIZER_VERSION): fits each
# gene's expression across pseudobulk samples as one Gaussian or as an
# OFF/ON mixture of two with Stan, compares the fits by k-fold
# cross-validation, and calls each cluster's probability of expression.
# fit_mixture() fits genes, merge_mixture_fits() joins fits run in chunks,
# call_mixture() makes the per-split table that add_probability() reads.
#
# Attribution: the model, the Stan code (inst/stan/) and the calling logic
# are adapted from Fred P. Davis's opticlobe code
# (https://github.com/fredpdavis/opticlobe), published with Davis et al.
# (2020) eLife 9:e50901, via the mm_marker / mixture_modeling_smk pipeline
# (https://github.com/chenyenchung/mixture_modeling_smk), which this
# reimplements value for value. Unlike opticlobe, the held-out elpd averages
# each sample's likelihood over posterior draws
# (https://github.com/fredpdavis/opticlobe/pull/1).
#
# Copyright (c) 2018-2020 Fred P. Davis (original opticlobe code)
# Copyright (c) 2023-2026 Yen-Chung Chen (this adaptation)
# opticlobe carries no license; Fred P. Davis permitted this adaptation and
# its release under the MIT license (email, 2026-10-26). See inst/COPYRIGHTS.

# Stan settings of the pipeline
MIXTURE_CHAINS <- 4
MIXTURE_ADAPT_DELTA <- 0.98
# Bimodal genes' samples with P(On) above / below these define the global
# ON / OFF expression distributions that genes without bimodality are
# compared with
MIXTURE_SURE_ON <- 0.9
MIXTURE_SURE_OFF <- 0.1

stan_models <- new.env(parent = emptyenv())

need_packages <- function(pkgs, what) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
        fail(what, " needs the package(s) ", paste(missing, collapse = ", "),
            ": install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))")
    }
}

mixture_stan_model <- function(name) {
    # Compiled Stan model inst/stan/<name>.stan, compiled once per session
    if (is.null(stan_models[[name]])) {
        message("Compiling the ", name, " Stan model (once per session)")
        file <- system.file("stan", paste0(name, ".stan"), package = "scMarco")
        stan_models[[name]] <- rstan::stan_model(file, model_name = name,
            auto_write = FALSE)
    }
    stan_models[[name]]
}

run_stan <- function(name, data, n_iter, cores) {
    # Seeds come from R's RNG, so set.seed() makes a run reproducible. With
    # cores = 1 rstan takes more draws from it, so the seeds of later runs
    # differ from those of the pipeline, which ran its chains in parallel
    rstan::sampling(mixture_stan_model(name), data = data,
        chains = MIXTURE_CHAINS, iter = n_iter, cores = cores,
        control = list(adapt_delta = MIXTURE_ADAPT_DELTA), refresh = -1,
        verbose = FALSE)
}

partition_folds <- function(n, k, labels) {
    # Fold (1..k) of each of n samples for k-fold cross-validation: shuffled
    # into k folds of floor(n / k), the remainder joining fold k. A label
    # (replicate) found in one fold only swaps a sample with another label's
    # sample in another fold until every label spans two folds or more (a
    # label with one sample, e.g. a k-means subcluster number only one
    # cluster reaches, cannot and is left as it is).
    folds <- numeric(n)
    shuffled <- sample(seq_len(n), n)
    size <- floor(n / k)
    for (i in seq_len(k)) {
        start <- (i - 1) * size + 1
        folds[shuffled[start:(start + size - 1)]] <- i
    }
    if (n %% k > 0) folds[shuffled[(k * size + 1):n]] <- k

    repeat {
        balanced <- TRUE
        for (label in unique(labels)) {
            if (sum(labels == label) < 2) next
            parts <- unique(folds[labels == label])
            if (length(parts) == 1) {
                i <- sample(which(labels == label))[1]
                j <- sample(which(labels != label & folds != parts[1]))[1]
                folds[c(i, j)] <- folds[c(j, i)]
                balanced <- FALSE
            }
        }
        if (balanced) break
    }
    folds
}

log_mean_exp <- function(ll) {
    # elpd of each held-out sample: log of its mean likelihood over draws
    # ll: draws x samples log likelihoods
    vapply(seq_len(ncol(ll)), function(h) matrixStats::logSumExp(ll[, h]),
        numeric(1)) - log(nrow(ll))
}

fit_gene <- function(x, cluster_idx, replicate_idx, labels, n_folds, n_iter, cores) {
    # Bimodal and unimodal fits of one gene and their k-fold cross-validated
    # expected log predictive densities (elpd)
    # x: expression of each sample; cluster_idx, replicate_idx: 1-based
    #   indices of each sample's cluster and replicate; labels: replicate
    #   labels kept apart by partition_folds()
    n <- length(x)
    data <- list(n_samples = n, n_clusters = max(cluster_idx),
        n_replicates = max(replicate_idx), log_e = x, cluster = cluster_idx,
        replicate = replicate_idx)
    stats <- rstan::summary(run_stan("bimodal", data, n_iter, cores))$summary
    vec <- function(par, len) paste0(par, "[", seq_len(len), "]")
    bimodal <- list(
        log_p_on_cluster = unname(stats[vec("log_p_on_cluster", data$n_clusters), "50%"]),
        log_p_on_replicate = unname(stats[vec("log_p_on_replicate", data$n_replicates), "50%"]),
        log_p_on_sample = unname(stats[vec("log_p_on_sample", n), "50%"]),
        on_level = stats["mu_on", "50%"],
        off_level = stats["mu_off", "50%"],
        sd = stats["sd1", "50%"],
        pi = stats["pi_on", "50%"],
        diagnostics = stats[c("mu_on", "mu_off", "sd1", "pi_on"), c("Rhat", "n_eff")]
    )
    stats <- rstan::summary(run_stan("unimodal", data[c("n_samples", "log_e")],
        n_iter, cores))$summary
    unimodal <- list(
        level = stats["mu1", "50%"],
        sd = stats["sd1", "50%"],
        diagnostics = stats[c("mu1", "sd1"), c("Rhat", "n_eff")]
    )

    folds <- partition_folds(n, n_folds, labels)
    bimodal_i <- numeric(n)
    unimodal_i <- numeric(n)
    for (k in seq_len(n_folds)) {
        held <- which(folds == k)
        train <- which(folds != k)
        cv_data <- list(n_t = length(train), n_h = length(held),
            log_e_t = x[train], log_e_h = x[held])
        bim_fit <- run_stan("bimodal_cv", cv_data, n_iter, cores)
        uni_fit <- run_stan("unimodal_cv", cv_data, n_iter, cores)
        bimodal_i[held] <- log_mean_exp(as.matrix(bim_fit, pars = "log_lik_h"))
        unimodal_i[held] <- log_mean_exp(as.matrix(uni_fit, pars = "log_lik_h"))
    }
    elpd <- list(
        bimodal = sum(bimodal_i), unimodal = sum(unimodal_i),
        bimodal_se = sqrt(n * stats::var(bimodal_i)),
        unimodal_se = sqrt(n * stats::var(unimodal_i)),
        diff = sum(bimodal_i) - sum(unimodal_i),
        diff_se = sqrt(n * stats::var(bimodal_i - unimodal_i)),
        bimodal_i = bimodal_i, unimodal_i = unimodal_i
    )
    list(bimodal = bimodal, unimodal = unimodal, elpd = elpd)
}

#' Fit the mixture model of each gene
#'
#' Fits the log-normalized expression of each gene across pseudobulk samples
#' of one split twice with Stan: as a mixture of an OFF and an ON Gaussian
#' with a shared SD, and as a single Gaussian. Each fit is scored by k-fold
#' cross-validation. [call_mixture()] then calls bimodal genes and each
#' cluster's probability of expression. This is scMarco's binarizer, adapted
#' from Fred P. Davis's mixture model and code for Davis et al. (2020)
#' (<https://github.com/fredpdavis/opticlobe>), with the elpd fix of
#' <https://github.com/fredpdavis/opticlobe/pull/1>.
#'
#' Each gene takes 2 + 2 `n_folds` Stan runs, about a minute, so for a whole
#' dataset run chunks of genes in parallel jobs (`genes`) and join them with
#' merge_mixture_fits(). The Stan models are compiled once per R session.
#'
#' Fits are reproducible by default: each gene is fit after
#' `set.seed(seed + i)`, `i` being its row in `lognorm`, so a gene's fit does
#' not depend on how the genes are chunked as long as every chunk gets the
#' same `lognorm` and picks its genes with `genes`. The global random state is
#' restored afterwards. With `seed = NULL` the genes are fit one after another
#' from the current random state, as the mixture_modeling_smk pipeline did.
#' Stan's seeds come from R's generator, and rstan draws from it differently
#' when its chains run one after another, so `cores = 1` gives other values
#' than any `cores > 1`, which matches the pipeline.
#'
#' @param lognorm Log-normalized expression, genes x pseudobulk samples (a
#'   matrix or data.frame with gene row names). Each cluster needs samples in
#'   several replicates, e.g. its cells pooled into two random halves or
#'   into k-means subclusters with [make_pseudobulk()].
#' @param cluster Cluster ID of each sample (column of `lognorm`).
#' @param replicate Replicate of each sample, kept apart in cross-validation.
#'   `NULL` puts every sample in one replicate.
#' @param genes Genes to fit; `NULL` fits every row of `lognorm`.
#' @param n_iter Stan iterations per chain, half of them warm-up.
#' @param n_folds Number of cross-validation folds.
#' @param cores Cores for the 4 Stan chains (see Details for `cores = 1`).
#' @param seed Random seed (see Details); `NULL` uses the current random
#'   state.
#'
#' @return A `scmarco_mixture` list: `fits`, one fit per gene (`bimodal`
#'   and `unimodal` posterior medians, `elpd`), and the input `lognorm`,
#'   `cluster` and `replicate`.
#' @seealso [make_pseudobulk()] for the input, [call_mixture()],
#'   [mixture_heterogeneity()], [add_probability()].
#' @references Davis FP, Nern A, Picard S, et al. (2020). A genetic, genomic,
#'   and computational resource for exploring neural circuit function. eLife
#'   9:e50901. Code: <https://github.com/fredpdavis/opticlobe>.
#' @examples
#' \dontrun{
#' # one split: two pseudobulk replicates of each of 3 clusters
#' lognorm <- matrix(c(rnorm(6, 0.1, 0.05), rnorm(6, 2, 0.2)), 2, 6,
#'     byrow = TRUE, dimnames = list(c("geneA", "geneB"), paste0("s", 1:6)))
#' fit <- fit_mixture(lognorm, cluster = rep(c("1", "2", "3"), 2),
#'     replicate = rep(c("1", "2"), each = 3), n_folds = 2)
#' call_mixture(fit)
#' }
#' @export
fit_mixture <- function(lognorm, cluster, replicate = NULL, genes = NULL,
                        n_iter = 500, n_folds = 10,
                        cores = min(4, parallel::detectCores(), na.rm = TRUE),
                        seed = 42) {
    need_packages(c("rstan", "matrixStats"), "fit_mixture()")
    if (is.data.frame(lognorm)) {
        not_num <- colnames(lognorm)[!vapply(lognorm, is.numeric, logical(1))]
        if (length(not_num) > 0) fail("lognorm columns must be numeric: ", first_few(not_num))
        lognorm <- as.matrix(lognorm)
    }
    if (!is.numeric(lognorm)) fail("lognorm must be a numeric matrix")
    if (is.null(rownames(lognorm))) fail("lognorm needs genes as row names")
    if (length(cluster) != ncol(lognorm)) {
        fail("cluster has ", length(cluster), " entries for ", ncol(lognorm), " samples")
    }
    if (is.null(replicate)) replicate <- rep("1", ncol(lognorm))
    if (length(replicate) != ncol(lognorm)) {
        fail("replicate has ", length(replicate), " entries for ", ncol(lognorm), " samples")
    }
    if (is.null(genes)) genes <- rownames(lognorm)
    unknown <- setdiff(genes, rownames(lognorm))
    if (length(unknown) > 0) fail("Gene(s) not in lognorm: ", first_few(unknown))
    cluster <- as.character(cluster)
    replicate <- as.character(replicate)

    cluster_idx <- as.integer(factor(cluster, levels = unique(cluster)))
    replicate_idx <- as.integer(factor(replicate, levels = unique(replicate)))
    labels <- as.integer(as.factor(replicate))
    if (!is.null(seed)) {
        # leave the caller's random state as it was
        had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
        if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(if (had_seed) {
            assign(".Random.seed", old_seed, envir = globalenv())
        } else {
            rm(".Random.seed", envir = globalenv())
        }, add = TRUE)
    }
    fits <- list()
    for (i in seq_along(genes)) {
        message(sprintf("Fitting %s (%d/%d)", genes[i], i, length(genes)))
        if (!is.null(seed)) set.seed(seed + match(genes[i], rownames(lognorm)))
        fits[[genes[i]]] <- fit_gene(as.numeric(lognorm[genes[i], ]), cluster_idx,
            replicate_idx, labels, n_folds, n_iter, cores)
    }
    structure(list(fits = fits, lognorm = lognorm[genes, , drop = FALSE],
        cluster = cluster, replicate = replicate), class = "scmarco_mixture")
}

#' @rdname fit_mixture
#' @param ... `scmarco_mixture` fits of the same samples, or one list of
#'   them, e.g. chunks of genes fit in parallel.
#' @export
merge_mixture_fits <- function(...) {
    parts <- list(...)
    if (length(parts) == 1 && !inherits(parts[[1]], "scmarco_mixture")) parts <- parts[[1]]
    if (!all(vapply(parts, inherits, logical(1), "scmarco_mixture"))) {
        fail("merge_mixture_fits() takes fits from fit_mixture()")
    }
    first <- parts[[1]]
    for (p in parts[-1]) {
        if (!identical(colnames(p$lognorm), colnames(first$lognorm)) ||
            !identical(p$cluster, first$cluster) || !identical(p$replicate, first$replicate)) {
            fail("Fits to merge must share their samples, clusters and replicates")
        }
    }
    fits <- do.call(c, lapply(parts, `[[`, "fits"))
    if (anyDuplicated(names(fits))) {
        fail("Gene(s) fit twice: ", first_few(unique(names(fits)[duplicated(names(fits))])))
    }
    structure(list(fits = fits,
        lognorm = do.call(rbind, lapply(parts, `[[`, "lognorm")),
        cluster = first$cluster, replicate = first$replicate), class = "scmarco_mixture")
}

mixture_bimodal <- function(fit) {
    # Named logical per gene: the bimodal fit beats the unimodal one by more
    # than two standard errors of the elpd difference
    diff <- vapply(fit$fits, function(g) g$elpd$diff, numeric(1))
    diff_se <- vapply(fit$fits, function(g) g$elpd$diff_se, numeric(1))
    diff > 0 & diff > 2 * diff_se
}

#' Call expression from mixture-model fits
#'
#' Turns the fits of [fit_mixture()] into one split's probability table for
#' [add_probability()]. A gene is bimodal when its bimodal fit beats the
#' unimodal one by more than two standard errors of the elpd difference; its
#' probability of expression in a cluster is then the bimodal model's. For
#' the other genes, expression in all of the split's samples is compared
#' with the ON and OFF distributions pooled from the bimodal genes (samples
#' with P(On) > 0.9 and < 0.1), which gives one probability for every
#' cluster. Needs at least one bimodal gene.
#'
#' @param fit Fits from [fit_mixture()] or [merge_mixture_fits()], of
#'   every gene of the split (the pooled distributions use them all).
#'
#' @return A data.frame, genes x clusters probabilities of expression
#'   (rounded to 2 decimals), followed by `diff.elpd`, `diff.elpd.se`
#'   (rounded), the bimodal fit's `bimodal.on_level`, `bimodal.off_level`,
#'   `bimodal.sd` and `bimodal.pi` (posterior medians, unrounded so that the
#'   mixture-model plot is exact), `unimodal.level`, `unimodal.sd` (rounded)
#'   and `exprType` (`"bimodal"` or `"unimodal"`).
#' @seealso [fit_mixture()], [add_probability()].
#' @references Davis FP, Nern A, Picard S, et al. (2020). A genetic, genomic,
#'   and computational resource for exploring neural circuit function. eLife
#'   9:e50901. Code: <https://github.com/fredpdavis/opticlobe>.
#' @examples
#' \dontrun{
#' tables <- list(Male = call_mixture(male_fit), Female = call_mixture(female_fit))
#' add_probability("data/mydata.clusters.h5ad", tables)
#' }
#' @export
call_mixture <- function(fit) {
    need_packages(c("mclust", "matrixStats"), "call_mixture()")
    if (!inherits(fit, "scmarco_mixture")) fail("fit must come from fit_mixture()")
    genes <- names(fit$fits)
    lognorm <- fit$lognorm[genes, , drop = FALSE]
    clusters <- unique(fit$cluster)
    pick <- function(f) vapply(fit$fits, f, numeric(1))
    diff <- pick(function(g) g$elpd$diff)
    diff_se <- pick(function(g) g$elpd$diff_se)
    bimodal <- mixture_bimodal(fit)
    p_on <- do.call(rbind, lapply(fit$fits, function(g) exp(g$bimodal$log_p_on_cluster)))
    dimnames(p_on) <- list(genes, clusters)
    message(sum(bimodal), " of ", length(genes), " genes are bimodal")
    if (!any(bimodal)) fail("No bimodal gene to learn the ON and OFF distributions from")

    # Global ON and OFF distributions: samples of bimodal genes whose cluster
    # is surely ON or OFF, read gene by gene
    bim_lognorm <- t(lognorm[bimodal, , drop = FALSE])
    bim_p_on <- t(p_on[bimodal, fit$cluster, drop = FALSE])
    normal_fit <- function(x) {
        # what Mclust(x, G = 1) fits
        f <- mclust::mvnX(x)$parameters
        c(mean = unname(f$mean), sd = sqrt(f$variance$sigmasq))
    }
    on_e <- bim_lognorm[bim_p_on > MIXTURE_SURE_ON]
    off_e <- bim_lognorm[bim_p_on < MIXTURE_SURE_OFF]
    on <- normal_fit(on_e)
    off <- normal_fit(off_e)
    pi_on <- length(on_e) / (length(on_e) + length(off_e))
    message(sprintf("Global ON: mean %.3f, SD %.3f; OFF: mean %.3f, SD %.3f; ON weight %.3f",
        on[["mean"]], on[["sd"]], off[["mean"]], off[["sd"]], pi_on))

    # P(On) of genes without bimodality: all samples ON vs all samples OFF
    log_on <- numeric(length(genes))
    log_off <- numeric(length(genes))
    for (s in seq_len(ncol(lognorm))) {
        log_on <- log_on + (log(pi_on) +
            dnorm(lognorm[, s], mean = on[["mean"]], sd = on[["sd"]], log = TRUE))
        log_off <- log_off + (log1p(-1 * pi_on) +
            dnorm(lognorm[, s], mean = off[["mean"]], sd = off[["sd"]], log = TRUE))
    }
    log_p_unimodal <- log_on - vapply(seq_along(genes), function(i) {
        matrixStats::logSumExp(c(log_off[i], log_on[i]))
    }, numeric(1))
    p_on[!bimodal, ] <- exp(log_p_unimodal[!bimodal])

    round2 <- function(x) round(x, digits = 2)
    out <- data.frame(round2(p_on), check.names = FALSE)
    out$diff.elpd <- round2(diff)
    out$diff.elpd.se <- round2(diff_se)
    out$bimodal.on_level <- pick(function(g) g$bimodal$on_level)
    out$bimodal.off_level <- pick(function(g) g$bimodal$off_level)
    out$bimodal.sd <- pick(function(g) g$bimodal$sd)
    out$bimodal.pi <- pick(function(g) g$bimodal$pi)
    out$unimodal.level <- round2(pick(function(g) g$unimodal$level))
    out$unimodal.sd <- round2(pick(function(g) g$unimodal$sd))
    out$exprType <- ifelse(bimodal, "bimodal", "unimodal")
    out
}
