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
# Binarizer 1.1 changes the fit (specification in the scMarco_test design
# record): a per-gene prior on the component means instead of normal(7, 5),
# 8 chains for the full bimodal fit, and the parameter draws, posterior mean
# P(On) and chain spread kept per gene. Binarizer 1.0 is the original
# implementation (equal to mm_marker within seed noise).
#
# Copyright (c) 2018-2020 Fred P. Davis (original opticlobe code)
# Copyright (c) 2023-2026 Yen-Chung Chen (this adaptation)
# opticlobe carries no license; Fred P. Davis permitted this adaptation and
# its release under the MIT license (email, 2026-10-26). See inst/COPYRIGHTS.

# Stan settings of the pipeline
MIXTURE_CHAINS <- 4
MIXTURE_ADAPT_DELTA <- 0.98
# Binarizer 1.1: chains of the full bimodal fit (its chains' agreement on
# each cluster is the convergence check; the other 21 Stan runs per gene keep
# MIXTURE_CHAINS), and the per-gene prior on the component means,
# normal(midpoint, range / MIXTURE_PRIOR_SD_DIV) of the gene's values with an
# SD floor (replacing the original normal(7, 5), set for bulk log1p(TPM))
MIXTURE_CHAINS_BIMODAL <- 8
MIXTURE_PRIOR_SD_DIV <- 3
MIXTURE_PRIOR_SD_FLOOR <- 0.01
# Bimodal genes' samples with P(On) above / below these define the global
# ON / OFF expression distributions that genes without bimodality are
# compared with
MIXTURE_SURE_ON <- 0.9
MIXTURE_SURE_OFF <- 0.1

# A fit whose Stan run has a parameter with Rhat above MIXTURE_MAX_RHAT or an
# effective sample size below MIXTURE_MIN_N_EFF is reported as not converged
MIXTURE_MAX_RHAT <- 1.05
MIXTURE_MIN_N_EFF <- 100

stan_models <- new.env(parent = emptyenv())

need_packages <- function(pkgs, what) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) {
        fail(what, " needs the package(s) ", paste(missing, collapse = ", "),
            ": install.packages(c(", paste0('"', missing, '"', collapse = ", "), "))")
    }
}

mixture_stan_model <- function(name) {
    # Compiled Stan model inst/stan/<name>.stan: kept for the session, and
    # cached on disk so that later sessions (e.g. the other jobs of a chunked
    # fit) load it instead of compiling it again
    if (is.null(stan_models[[name]])) {
        file <- system.file("stan", paste0(name, ".stan"), package = "scMarco")
        stan_models[[name]] <- cached_stan_model(name, file,
            compile = function() {
                message("Compiling the ", name, " Stan model (needs ~12 GB of memory)")
                tryCatch(
                    rstan::stan_model(file, model_name = name, auto_write = FALSE),
                    error = function(e) fail("Compiling the ", name, " Stan model failed. ",
                        "The C++ compiler needs ~12 GB of memory; when it is killed for ",
                        "lack of memory the error is unclear, e.g. 'invalid connection'. ",
                        "Original error: ", conditionMessage(e)))
            })
    }
    stan_models[[name]]
}

stan_cache_dir <- function() {
    # Where compiled Stan models are cached: option scMarco.stan_cache, else
    # SCMARCO_STAN_CACHE, else the user cache directory. "" turns caching off.
    dir <- getOption("scMarco.stan_cache", Sys.getenv("SCMARCO_STAN_CACHE", NA))
    if (is.na(dir)) dir <- tools::R_user_dir("scMarco", which = "cache")
    dir
}

stan_cache_file <- function(name, file, dir = stan_cache_dir()) {
    # Cache file of model `name`: keyed by its Stan code and by what the
    # compiled code depends on, so an edited model or another rstan,
    # StanHeaders or R gets its own file
    v <- function(p) as.character(utils::packageVersion(p))
    key <- paste(unname(tools::md5sum(file)), v("rstan"), v("StanHeaders"),
        R.version$major, R.version$minor, R.version$platform)
    key_file <- tempfile()   # md5sum() hashes files only
    on.exit(unlink(key_file), add = TRUE)
    writeLines(key, key_file)
    tag <- substr(unname(tools::md5sum(key_file)), 1, 16)
    file.path(dir, sprintf("%s-%s.rds", name, tag))
}

stanmodel_usable <- function(m) {
    # A cached stanmodel can be used when it holds its compiled code for this
    # system (rstan's own check for auto_write)
    methods::is(m, "stanmodel") && isTRUE(m@dso@dso_saved) &&
        identical(m@dso@system, R.version$system)
}

cached_stan_model <- function(name, file, compile, usable = stanmodel_usable,
                              dir = stan_cache_dir()) {
    # Load model `name` from the disk cache, or compile() it and cache it. A
    # cache that cannot be read or written only costs a compilation.
    if (!nzchar(dir)) return(compile())
    path <- stan_cache_file(name, file, dir)
    if (file.exists(path)) {
        m <- tryCatch(readRDS(path), error = function(e) NULL)
        if (!is.null(m) && usable(m)) {
            message("Using the compiled ", name, " Stan model from ", path)
            return(m)
        }
    }
    m <- compile()
    # write to a temporary file and rename it, so parallel jobs that compile
    # at once never read a half-written cache
    tmp <- tempfile(paste0(name, "-"), tmpdir = dir, fileext = ".tmp")
    ok <- tryCatch({
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        saveRDS(m, tmp)
        file.rename(tmp, path)
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) {
        message("Cached the compiled ", name, " Stan model in ", path)
    } else {
        unlink(tmp)
        message("Could not cache the compiled ", name, " Stan model in ", dir,
            "; set options(scMarco.stan_cache = <writable folder>) to reuse it")
    }
    m
}

# The draws fit_gene() reads from each model. Stan saves only these: the
# others (e.g. the CV models' log_lik_t, bimodal's per-cluster sums) cost
# memory and summary() time, and dropping them does not change the draws.
MIXTURE_SAVED_PARS <- list(
    bimodal = c("mu_on", "mu_off", "sd1", "pi_on", "log_p_on_sample",
        "log_p_on_replicate", "log_p_on_cluster"),
    unimodal = c("mu1", "sd1"),
    bimodal_cv = "log_lik_h",
    unimodal_cv = "log_lik_h"
)

run_stan <- function(name, data, n_iter, cores, chains = MIXTURE_CHAINS) {
    # Seeds come from R's RNG, so set.seed() makes a run reproducible. With
    # cores = 1 rstan takes more draws from it, so the seeds of later runs
    # differ from those of the pipeline, which ran its chains in parallel
    rstan::sampling(mixture_stan_model(name), data = data,
        pars = MIXTURE_SAVED_PARS[[name]],
        chains = chains, iter = n_iter, cores = cores,
        control = list(adapt_delta = MIXTURE_ADAPT_DELTA), refresh = -1,
        verbose = FALSE)
}

partition_folds <- function(n, k, labels) {
    # Fold (1..k) of each of n samples for k-fold cross-validation: shuffled
    # into k folds of floor(n / k), the remainder joining fold k. A label
    # (replicate) found in one fold only swaps a sample with another label's
    # sample in another fold until every label spans two folds or more (a
    # label with one sample, e.g. a metadata group only one cluster has
    # enough cells in, cannot and is left as it is).
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

mixture_prior <- function(x) {
    # Binarizer 1.1 prior on both component means (data items of bimodal.stan
    # and bimodal_cv.stan): normal(midpoint, range / 3) of the values x, the
    # original normal(7, 5) on 0-14 expressed per gene (Richardson & Green
    # 1997, data-range priors for mixtures)
    lo <- min(x)
    hi <- max(x)
    sdv <- max((hi - lo) / MIXTURE_PRIOR_SD_DIV, MIXTURE_PRIOR_SD_FLOOR)
    list(prior_on_mean = (lo + hi) / 2, prior_on_sd = sdv,
        prior_off_mean = (lo + hi) / 2, prior_off_sd = sdv)
}

log_p_on_draws <- function(par, x, cl, n_cl) {
    # Per-draw log P(On) of each cluster from the parameter draws, as
    # bimodal.stan computes log_p_on_cluster
    # par: draws x 4 (mu_on, mu_off, sd1, pi_on); returns draws x n_cl
    on <- sapply(seq_len(nrow(par)), function(d) stats::dnorm(x, par[d, 1], par[d, 3], log = TRUE))
    off <- sapply(seq_len(nrow(par)), function(d) stats::dnorm(x, par[d, 2], par[d, 3], log = TRUE))
    con <- rowsum(on, cl, reorder = TRUE)[seq_len(n_cl), , drop = FALSE]
    coff <- rowsum(off, cl, reorder = TRUE)[seq_len(n_cl), , drop = FALSE]
    a <- sweep(con, 2, log(par[, 4]), "+")
    b <- sweep(coff, 2, log1p(-par[, 4]), "+")
    m <- pmax(a, b)
    t(a - m - log(exp(a - m) + exp(b - m)))
}

cluster_summaries <- function(draws, x, cl, n_cl) {
    # Posterior mean of each cluster's P(On), and the spread (max - min) of
    # the chains' means: chains that disagree on a cluster's state
    # draws: iter x chain x 4 array of mu_on, mu_off, sd1, pi_on
    par <- matrix(draws, ncol = 4)
    chain <- rep(seq_len(dim(draws)[2]), each = dim(draws)[1])
    p <- exp(log_p_on_draws(par, x, cl, n_cl))
    cm <- rowsum(p, chain) / as.vector(table(chain))
    list(p_on_mean = colMeans(p), chain_spread = apply(cm, 2, max) - apply(cm, 2, min))
}

log_mean_exp <- function(ll) {
    # elpd of each held-out sample: log of its mean likelihood over draws
    # ll: draws x samples log likelihoods
    vapply(seq_len(ncol(ll)), function(h) matrixStats::logSumExp(ll[, h]),
        numeric(1)) - log(nrow(ll))
}

fit_gene <- function(x, cluster_idx, replicate_idx, labels, n_folds, n_iter, cores,
                     n_chains_bimodal = MIXTURE_CHAINS_BIMODAL) {
    # Bimodal and unimodal fits of one gene and their k-fold cross-validated
    # expected log predictive densities (elpd)
    # x: expression of each sample; cluster_idx, replicate_idx: 1-based
    #   indices of each sample's cluster and replicate; labels: replicate
    #   labels kept apart by partition_folds()
    # Binarizer 1.1: the bimodal models take the per-gene prior
    # (mixture_prior(): the gene's values for the full fit, the training
    # fold's for each CV fit); the full bimodal fit runs n_chains_bimodal
    # chains and keeps its parameter draws, the posterior mean of each
    # cluster's P(On) and the spread of the chains' means
    n <- length(x)
    data <- list(n_samples = n, n_clusters = max(cluster_idx),
        n_replicates = max(replicate_idx), log_e = x, cluster = cluster_idx,
        replicate = replicate_idx)
    bim <- run_stan("bimodal", c(data, mixture_prior(x)), n_iter, cores,
        chains = n_chains_bimodal)
    stats <- rstan::summary(bim)$summary
    draws <- as.array(bim, pars = c("mu_on", "mu_off", "sd1", "pi_on"))
    vec <- function(par, len) paste0(par, "[", seq_len(len), "]")
    bimodal <- list(
        log_p_on_cluster = unname(stats[vec("log_p_on_cluster", data$n_clusters), "50%"]),
        log_p_on_replicate = unname(stats[vec("log_p_on_replicate", data$n_replicates), "50%"]),
        log_p_on_sample = unname(stats[vec("log_p_on_sample", n), "50%"]),
        on_level = stats["mu_on", "50%"],
        off_level = stats["mu_off", "50%"],
        sd = stats["sd1", "50%"],
        pi = stats["pi_on", "50%"],
        diagnostics = stats[c("mu_on", "mu_off", "sd1", "pi_on"), c("Rhat", "n_eff")],
        draws = draws
    )
    bimodal <- c(bimodal, cluster_summaries(draws, x, cluster_idx, data$n_clusters))
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
        bim_fit <- run_stan("bimodal_cv", c(cv_data, mixture_prior(x[train])), n_iter, cores)
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
#' merge_mixture_fits().
#'
#' The four Stan models are compiled the first time they are needed, which
#' takes a few minutes and ~12 GB of memory, and cached in
#' `tools::R_user_dir("scMarco", "cache")`, so later sessions and jobs load
#' them instead. Set `options(scMarco.stan_cache = <folder>)` (or the
#' environment variable `SCMARCO_STAN_CACHE`) to cache them elsewhere, e.g. a
#' folder that every node of a cluster sees, or `""` to compile in every
#' session. Before submitting many jobs, fit one gene once so the jobs find
#' the models compiled. A new rstan, StanHeaders or R version compiles them
#' again.
#'
#' A gene whose bimodal or unimodal fit has a parameter with Rhat above 1.05
#' or an effective sample size below 100 is reported as not converged, with
#' a warning, in `$convergence` and by [call_mixture()].
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
#'   several replicates, e.g. its cells pooled into two random halves with
#'   [make_pseudobulk()].
#' @param cluster Cluster ID of each sample (column of `lognorm`).
#' @param replicate Replicate of each sample, kept apart in cross-validation.
#'   `NULL` puts every sample in one replicate.
#' @param genes Genes to fit; `NULL` fits every row of `lognorm`.
#' @param n_iter Stan iterations per chain, half of them warm-up.
#' @param n_folds Number of cross-validation folds.
#' @param cores Cores for the Stan chains (see Details for `cores = 1`).
#' @param seed Random seed (see Details); `NULL` uses the current random
#'   state.
#' @param n_chains_bimodal Chains of the full bimodal fit (default 8). Their
#'   agreement on each cluster is the convergence check; with fewer than 4,
#'   chains in different modes are often missed (a warning says so). The
#'   unimodal and cross-validation fits run 4 chains.
#'
#' @return A `scmarco_mixture` list: `fits`, one fit per gene (`bimodal`
#'   and `unimodal` posterior medians, `elpd`; the bimodal fit also keeps
#'   its parameter `draws`, each cluster's posterior mean P(On) `p_on_mean`
#'   and the spread of its chains' means `chain_spread`), the input
#'   `lognorm`, `cluster` and `replicate`, the fit `settings` (`n_iter`,
#'   `n_folds`, `seed`, `cores`, `prior`, `n_chains_bimodal`, `n_chains`), `convergence` (one row per gene: the largest Rhat
#'   and smallest effective sample size of the bimodal and unimodal fits'
#'   parameters, and `converged`), and the versions that make up its contract:
#'   `schema_version` (the scMarco schema, as `format_version` of the
#'   dataset files) and `binarizer_version`
#'   (of the binarizer that fit it). [merge_mixture_fits()] and
#'   [call_mixture()] refuse a newer schema or binarizer; fits of an older
#'   binarizer (e.g. 1.0) are read and called with a warning recommending a
#'   refit, and fits of different binarizers are never merged.
#' @seealso [make_pseudobulk()] for the input, [call_mixture()],
#'   [add_probability()].
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
                        seed = 42, n_chains_bimodal = 8) {
    need_packages(c("rstan", "matrixStats"), "fit_mixture()")
    if (!is.numeric(n_chains_bimodal) || length(n_chains_bimodal) != 1 ||
        is.na(n_chains_bimodal) || n_chains_bimodal < 2 ||
        n_chains_bimodal != round(n_chains_bimodal)) {
        fail("n_chains_bimodal must be a whole number of at least 2")
    }
    n_chains_bimodal <- as.integer(n_chains_bimodal)
    if (n_chains_bimodal < 4) {
        warning("n_chains_bimodal = ", n_chains_bimodal, ": with fewer than 4 chains, ",
            "chains that settle in different modes are often missed, so fewer ",
            "unresolved clusters are flagged (default ", MIXTURE_CHAINS_BIMODAL, ")",
            call. = FALSE)
    }
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
            replicate_idx, labels, n_folds, n_iter, cores, n_chains_bimodal)
    }
    out <- new_mixture(fits, lognorm[genes, , drop = FALSE], cluster, replicate,
        settings = list(n_iter = n_iter, n_folds = n_folds,
            seed = if (is.null(seed)) NA_real_ else seed, cores = cores,
            prior = "range", n_chains_bimodal = n_chains_bimodal,
            n_chains = MIXTURE_CHAINS))
    bad <- not_converged(out)
    if (length(bad) > 0) warning(convergence_note(bad, length(genes)), call. = FALSE)
    out
}

new_mixture <- function(fits, lognorm, cluster, replicate, settings) {
    structure(list(
        schema_version = SCMARCO_FORMAT_VERSION,
        binarizer_version = SCMARCO_BINARIZER_VERSION,
        fits = fits, lognorm = lognorm, cluster = cluster,
        replicate = replicate, settings = settings,
        convergence = mixture_convergence(fits)
    ), class = "scmarco_mixture")
}

mixture_convergence <- function(fits) {
    # One row per gene: worst Rhat and smallest effective sample size of the
    # bimodal and unimodal fits' parameters, and whether both converged (NA
    # for fits without diagnostics)
    worst <- function(d, col, f) {
        if (is.null(d)) return(NA_real_)
        x <- d[, col]
        if (anyNA(x)) return(NaN)   # e.g. Rhat of a constant chain
        f(x)
    }
    stat <- function(model, col, f) vapply(fits, function(g) {
        part <- if (is.list(g)) g[[model]]
        worst(if (is.list(part)) part$diagnostics, col, f)
    }, numeric(1))
    out <- data.frame(
        gene = names(fits),
        bimodal_rhat = stat("bimodal", "Rhat", max),
        bimodal_n_eff = stat("bimodal", "n_eff", min),
        unimodal_rhat = stat("unimodal", "Rhat", max),
        unimodal_n_eff = stat("unimodal", "n_eff", min),
        row.names = NULL
    )
    ok <- function(rhat, n_eff) rhat <= MIXTURE_MAX_RHAT & n_eff >= MIXTURE_MIN_N_EFF
    out$converged <- ok(out$bimodal_rhat, out$bimodal_n_eff) &
        ok(out$unimodal_rhat, out$unimodal_n_eff)
    # diagnostics recorded but not computable (NaN): not converged
    nan <- is.nan(out$bimodal_rhat) | is.nan(out$bimodal_n_eff) |
        is.nan(out$unimodal_rhat) | is.nan(out$unimodal_n_eff)
    out$converged[nan] <- FALSE
    out
}

not_converged <- function(fit) {
    # Genes whose fit did not converge; fits saved before fit_mixture()
    # recorded convergence get it computed from their diagnostics
    conv <- fit$convergence
    if (is.null(conv)) conv <- mixture_convergence(fit$fits)
    conv$gene[conv$converged %in% FALSE]
}

convergence_note <- function(genes, n) {
    sprintf(paste0("%d of %d gene(s) did not converge (Rhat > %s or effective sample ",
        "size < %d in the bimodal or unimodal fit): %s. See `$convergence`; ",
        "refitting them with a larger n_iter may help."),
        length(genes), n, MIXTURE_MAX_RHAT, MIXTURE_MIN_N_EFF, first_few(genes))
}

check_mixture <- function(x, what = "fit") {
    # Stop unless x is a scmarco_mixture this scMarco can use: the current
    # object schema and binarizer. Fits are not upgraded: a fit of another
    # binarizer would give probabilities neither version makes.
    if (!inherits(x, "scmarco_mixture")) fail(what, " must come from fit_mixture()")
    schema <- x$schema_version
    if (is.null(schema)) {
        fail(what, " was made before fits recorded their versions; ",
            "refit it with fit_mixture()")
    }
    if (!identical(schema, SCMARCO_FORMAT_VERSION)) {
        fail(what, " has schema version ", format(schema), " but this scMarco reads ",
            SCMARCO_FORMAT_VERSION, "; ", if (compareVersion(as.character(schema),
                SCMARCO_FORMAT_VERSION) > 0) "update scMarco" else
                "refit it with fit_mixture()")
    }
    # Fits of an older binarizer (e.g. 1.0) are read and called with a
    # warning; newer or invalid versions are refused
    v <- x$binarizer_version
    status <- binarizer_status(v)
    if (status %in% c("missing", "invalid", "external")) {
        fail(what, " records no valid binarizer version (", format(v), "); refit it ",
            "with fit_mixture()")
    }
    if (status == "newer") {
        fail(what, " was fit by binarizer ", v, ", newer than this scMarco's (",
            SCMARCO_BINARIZER_VERSION, "); update scMarco")
    }
    if (status == "older") {
        warning(what, " was fit by binarizer ", v, "; this scMarco's is ",
            SCMARCO_BINARIZER_VERSION, ". It is read and called as fit, but ",
            "refitting with fit_mixture() is recommended.", call. = FALSE)
    }
    invisible(TRUE)
}

#' @export
print.scmarco_mixture <- function(x, ...) {
    cat(sprintf("scMarco mixture-model fit: %d gene(s), %d sample(s), %d cluster(s)\n",
        length(x$fits), length(x$cluster), length(unique(x$cluster))))
    cat(sprintf("binarizer %s, schema %s\n", format(x$binarizer_version),
        format(x$schema_version)))
    bad <- not_converged(x)
    if (length(bad) > 0) cat(convergence_note(bad, length(x$fits)), "\n")
    invisible(x)
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
    # chunks from different binarizers are never merged
    versions <- vapply(parts, function(p) format(p$binarizer_version), "")
    if (length(unique(versions)) > 1) {
        fail("Fits from different binarizers cannot be merged (",
            paste(unique(versions), collapse = ", "), "); refit every chunk with ",
            "the same scMarco")
    }
    for (i in seq_along(parts)) check_mixture(parts[[i]], sprintf("fit %d", i))
    first <- parts[[1]]
    for (p in parts[-1]) {
        if (!identical(colnames(p$lognorm), colnames(first$lognorm)) ||
            !identical(p$cluster, first$cluster) || !identical(p$replicate, first$replicate)) {
            fail("Fits to merge must share their samples, clusters and replicates")
        }
        # cores may differ between jobs; the rest may not
        keys <- c("n_iter", "n_folds", "seed", "prior", "n_chains_bimodal", "n_chains")
        if (!identical(p$settings[keys], first$settings[keys])) {
            fail("Fits to merge must share their settings (n_iter, n_folds, seed, ",
                "prior, n_chains_bimodal, n_chains)")
        }
    }
    fits <- do.call(c, lapply(parts, `[[`, "fits"))
    if (anyDuplicated(names(fits))) {
        fail("Gene(s) fit twice: ", first_few(unique(names(fits)[duplicated(names(fits))])))
    }
    new_mixture(fits, do.call(rbind, lapply(parts, `[[`, "lognorm")),
        first$cluster, first$replicate, settings = first$settings)
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
#' cluster. Needs at least one bimodal gene. Genes whose fit did not converge
#' (see [fit_mixture()]) are named in a warning; their calls are kept.
#'
#' @param fit Fits from [fit_mixture()] or [merge_mixture_fits()], of
#'   every gene of the split (the pooled distributions use them all).
#'
#' @return A data.frame, genes x clusters probabilities of expression
#'   (rounded to 2 decimals), followed by `diff.elpd`, `diff.elpd.se`
#'   (rounded), the bimodal fit's `bimodal.on_level`, `bimodal.off_level`,
#'   `bimodal.sd` and `bimodal.pi` (posterior medians, unrounded so that the
#'   mixture-model plot is exact), `unimodal.level`, `unimodal.sd` (rounded)
#'   and `exprType` (`"bimodal"` or `"unimodal"`). Its attribute
#'   `scmarco_binarizer` holds the binarizer version, which
#'   [add_probability()] records in the dataset; it survives [saveRDS()],
#'   but not subsetting, after which the table counts as external.
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
    check_mixture(fit)
    genes <- names(fit$fits)
    lognorm <- fit$lognorm[genes, , drop = FALSE]
    clusters <- unique(fit$cluster)
    pick <- function(f) vapply(fit$fits, f, numeric(1))
    diff <- pick(function(g) g$elpd$diff)
    diff_se <- pick(function(g) g$elpd$diff_se)
    bimodal <- diff > 0 & diff > 2 * diff_se
    p_on <- do.call(rbind, lapply(fit$fits, function(g) exp(g$bimodal$log_p_on_cluster)))
    dimnames(p_on) <- list(genes, clusters)
    message(sum(bimodal), " of ", length(genes), " genes are bimodal")
    bad <- not_converged(fit)
    if (length(bad) > 0) {
        warning(convergence_note(bad, length(genes)), " Their calls are kept.", call. = FALSE)
    }
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
    attr(out, "scmarco_binarizer") <- fit$binarizer_version
    out
}
