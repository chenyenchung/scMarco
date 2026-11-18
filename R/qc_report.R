# Binarization QC report, version QC_REPORT_VERSION: a structured YAML file
# (stable keys, for comparing runs) and a self-contained HTML page (base R
# graphics embedded as base64 PNG) for one split's call_mixture() table.
# The version changes only after release, when a key or its meaning changes.

QC_REPORT_VERSION <- "1"

# Row labels of the HTML page: what each value is, not its key
QC_LABELS <- list(
    dataset = c(name = "Dataset", split = "Split", cells = "Cells", clusters = "Clusters",
        pseudobulk_samples = "Pseudobulk samples (cluster x replicate)",
        single_sample_clusters = "Clusters with only one pseudobulk sample", genes = "Genes"),
    genes = c(total = "Genes tested",
        binarized = "Genes binarized (their ON/OFF calls are used downstream)",
        not_binarized_evidence = "Not binarized: evidence for two states too weak at the chosen false discovery rate",
        not_binarized_ashman_D = "Not binarized: ON and OFF levels too close (Ashman's D)",
        not_binarized_chains = "Not binarized: the model's chains disagree on too many clusters",
        bimodality_cut_z = "Evidence cut (z = elpd difference / SE) at the chosen false discovery rate",
        near_cut = "Genes within 0.5 of the cut (may switch between binarized and not between runs)"),
    convergence = c(
        global_not_converged = "Genes whose fit failed the standard convergence check (R-hat > 1.05 or effective sample size < 100)",
        genes_with_cluster_over_cutoff = "Genes whose chains disagree on at least one cluster"),
    ambient = c(off_changed_ambient = "OFF calls made ambiguous: more UMIs than ambient RNA explains",
        on_changed_detection = "ON calls made ambiguous: too few cells detect the gene",
        calls_changed_chains = "Calls made ambiguous: the chains disagree on the cluster",
        genes_off_changed = "Genes with at least one OFF call made ambiguous (ambient check)",
        genes_on_changed = "Genes with at least one ON call made ambiguous (detection check)"),
    calls = c(ON = "ON calls (gene x cluster)", OFF = "OFF calls", ambiguous = "Ambiguous calls",
        no_call = "Entries without a call (gene not binarized)"))

QC_SETTING_INFO <- c(
    binarizer = "Binarizer that produced the calls",
    prior = "Prior on the ON and OFF levels: centred on each gene's midpoint, spread = its range / 3",
    n_iter = "MCMC iterations per chain (half are warm-up)",
    n_folds = "Folds of the cross-validation that decides whether a gene has two states",
    seed = "Random seed (each gene uses seed + its row number)",
    n_chains_bimodal = "Chains for the two-state model; disagreement between chains is the convergence check",
    bh_q = "False discovery rate for binarizing a gene (Benjamini-Hochberg on the evidence for two states)",
    min_D = "Minimum Ashman's D (distance between the ON and OFF levels relative to their spread); genes below are not binarized",
    spread_cutoff = "Chain spread above which a cluster's call is ambiguous",
    reject_spread_frac = "Share of clusters above the spread cutoff at which a gene is not binarized",
    summary = "How each cluster's P(ON) is summarised from the posterior",
    band = "Half-width of the ambiguity band around P(ON) = 0.5: ON above 0.5 + band, OFF below 0.5 - band",
    on_floor = "Minimum share of a cluster's cells with at least one UMI for an ON call",
    ambient_log2fc = "OFF becomes ambiguous only above this log2 fold over the expected ambient RNA (deliberately maximally permissive)",
    ambient_alpha = "Level of the test that a cluster's UMIs exceed that offset over ambient RNA",
    rho_source = "Where the contamination (ambient RNA) estimates come from")

QC_LIMITATIONS <- c(
    "Small clusters: the tests and the separation ranking have less power, so their calls are less certain, as in any single-cell analysis.",
    "Broad genes: when a gene is expressed in most clusters, whether a cluster is ON or OFF can stay unresolved; such calls are ambiguous or the gene is not binarized.",
    "Heterogeneous clusters: a cluster that mixes cell types is not flagged; the ambiguity band covers part of it.",
    "mRNA is not driver expression: a weak but real expression that a GAL4 driver may amplify can be called OFF; the ambient check is deliberately permissive and the leak ranking of the pair search shows smaller excesses.")

# Share of genes near the cut that changed status between two seeds (Ozel
# et al. 2021 P40, binarizer 1.1, seeds 42 and 43), shown next to near_cut
QC_NEAR_CUT_REFERENCE <- "about 19% in a two-seed comparison (Ozel et al. 2021, P40)"

qc_num <- function(x) {
    # integer quantities stay integers (no "2.0"), names kept
    if (is.numeric(x) && length(x) > 0 && all(is.finite(x)) && all(abs(x - round(x)) < 1e-9)) {
        setNames(as.integer(round(x)), names(x))
    } else {
        x
    }
}

qc_quant <- function(x, p = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)) {
    # named quantiles; observed values (type 1) for integers, so they stay integers
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NULL)
    is_int <- all(abs(x - round(x)) < 1e-9)
    v <- unname(stats::quantile(x, p, type = if (is_int) 1 else 7))
    v <- if (is_int) as.integer(v) else signif(v, 4)
    as.list(setNames(v, paste0("q", sprintf("%02d", round(100 * p)))))
}

#' Report on the quality of a binarization
#'
#' Summarizes one split's binarization (binarizer 1.1) for people and for
#' comparing runs: the dataset, the settings in use, which genes were
#' binarized and why the others were not, convergence (chain agreement), the
#' detection and ambient-RNA checks, the calls, known limitations and
#' provenance. With `file`, writes a structured YAML file (stable keys,
#' `qc_report_version`) and a self-contained HTML page (plots embedded).
#'
#' The calls are counted as [binarize()] makes them: the ambiguity band
#' `band` and the flags of [call_mixture()].
#'
#' @param fit Fits of the split from [fit_mixture()] or
#'   [merge_mixture_fits()].
#' @param table The split's table from [call_mixture()] of `fit`, with its
#'   attributes (not subset).
#' @param ambient The split's [ambient_summary()], as given to
#'   [call_mixture()]; `NULL` leaves out the cell counts and the observed vs
#'   expected ambient UMIs.
#' @param file Optional path without extension: writes `<file>.yaml` and
#'   `<file>.html`.
#' @param band Ambiguity band at which the calls are counted (the app's
#'   default is 0.3).
#' @param dataset,split Names shown in the report.
#'
#' @return The report as a list (the content of the YAML file), invisibly
#'   when `file` is given.
#' @seealso [call_mixture()], [ambient_summary()].
#' @examples
#' \dontrun{
#' amb <- ambient_summary(counts, cluster = obj$scmarco_cluster, library = obj$orig.ident)
#' tab <- call_mixture(fit, amb)
#' qc_report(fit, tab, amb, file = "qc/mydata_Male", dataset = "mydata", split = "Male")
#' }
#' @export
qc_report <- function(fit, table, ambient = NULL, file = NULL, band = 0.3,
                      dataset = "", split = "") {
    check_mixture(fit)
    params <- attr(table, "binarization")
    flags <- attr(table, "flags")
    if (is.null(params) || is.null(flags)) {
        fail("table must come from call_mixture() of binarizer ", SCMARCO_BINARIZER_VERSION,
            ", with its attributes (not subset)")
    }
    if (!is.null(ambient) && !inherits(ambient, "scmarco_ambient")) {
        fail("ambient must come from ambient_summary()")
    }
    if (!is.numeric(band) || length(band) != 1 || is.na(band) || band < 0 || band >= 0.5) {
        fail("band must be one number in [0, 0.5)")
    }
    genes <- rownames(table)
    if (!identical(genes, names(fit$fits))) {
        fail("table must have every gene of fit, in its order (not a selection): ",
            "the genes binarized and the evidence cut depend on all of them")
    }
    clusters <- setdiff(colnames(table), QC_COLS)
    binarized <- table$exprType %in% "bimodal"
    reason <- table$not.binarized

    # calls as binarize() makes them
    prob <- t(as.matrix(table[, clusters, drop = FALSE]))
    b <- binarize(prob, rep("s", length(clusters)), clusters, "s", clusters, band = band,
        bimodal = matrix(binarized, dimnames = list(genes, "s")),
        flags = t(flags[genes, clusters, drop = FALSE]))
    one <- function(a) matrix(a, length(genes), length(clusters), dimnames = list(genes, clusters))
    on <- one(b$on)
    amb <- one(b$ambiguous)
    flagged <- one(b$flagged)
    off <- binarized & !on & !amb
    n_called <- sum(binarized) * length(clusters)

    # per gene and per cluster
    z <- table$elpd.z
    z_cut <- if (any(binarized | reason %in% c("ashman_D", "chains"))) {
        min(z[!(reason %in% "elpd")])
    } else NA_real_
    spread <- do.call(rbind, lapply(fit$fits, function(g) {
        s <- g$bimodal$chain_spread
        if (is.null(s)) rep(NA_real_, length(clusters)) else s
    }))
    max_spread <- apply(spread, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    conv <- fit$convergence
    if (is.null(conv)) conv <- mixture_convergence(fit$fits)
    n_rep <- tapply(fit$replicate, fit$cluster, function(r) length(unique(r)))[clusters]
    n_cells <- if (!is.null(ambient)) ambient$n_cells[clusters]
    st <- fit$settings
    chains <- st$n_chains_bimodal

    warnings <- c(
        if (any(n_rep == 1)) sprintf("%d cluster(s) have a single pseudobulk sample; their fits are the least certain.", sum(n_rep == 1)),
        if (!is.null(chains) && chains < MIXTURE_CHAINS_BIMODAL) {
            sprintf("The two-state model ran %d chains (default %d); chain disagreement may be missed.", chains, MIXTURE_CHAINS_BIMODAL)
        },
        if (binarizer_status(fit$binarizer_version) == "older") {
            sprintf("The fit is from binarizer %s; refitting with binarizer %s is recommended.", fit$binarizer_version, SCMARCO_BINARIZER_VERSION)
        },
        if (!isTRUE(params$ambient)) "No ambient summary was given to call_mixture(): the detection and ambient checks did not run.")

    rho <- params$rho
    settings <- list(binarizer = fit$binarizer_version, prior = st$prior %||% "",
        n_iter = st$n_iter, n_folds = st$n_folds, seed = st$seed, n_chains_bimodal = chains,
        bh_q = params$bh_q, min_D = params$min_D, spread_cutoff = params$spread_cutoff,
        reject_spread_frac = params$reject_spread_frac, summary = params$summary, band = band,
        on_floor = params$on_floor, ambient_log2fc = params$ambient_log2fc,
        ambient_alpha = params$ambient_alpha, rho_source = params$rho_source %||% "")
    settings <- settings[!vapply(settings, is.null, logical(1))]

    blocked_ratio <- NULL
    detect_on <- NULL
    if (!is.null(ambient)) {
        obs <- ambient$umi[genes, clusters, drop = FALSE]
        lambda <- ambient_expected(ambient, genes, clusters)
        hit <- flagged == CALL_FLAG_AMBIENT
        blocked_ratio <- (obs / pmax(lambda, 1e-12))[hit]
        detect_on <- ambient$detect[genes, clusters, drop = FALSE][on]
    }
    n_on_cl <- colSums(on)
    n_amb_cl <- colSums(amb)
    amb_frac <- n_amb_cl / max(sum(binarized), 1)
    high_amb <- clusters[amb_frac > stats::quantile(amb_frac, 0.75) + 3 * stats::IQR(amb_frac)]
    n_on_g <- rowSums(on)[binarized]

    qc <- list(
        qc_report_version = QC_REPORT_VERSION,
        generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        dataset = list(name = dataset, split = split,
            cells = if (!is.null(n_cells)) qc_num(sum(n_cells)),
            clusters = length(clusters), pseudobulk_samples = length(fit$cluster),
            replicates_per_cluster = qc_quant(as.numeric(n_rep), c(0, 0.5, 1)),
            single_sample_clusters = sum(n_rep == 1),
            cells_per_cluster = if (!is.null(n_cells)) qc_quant(as.numeric(n_cells), c(0, 0.5, 1)),
            genes = length(genes),
            umis_per_library = if (!is.null(ambient)) as.list(qc_num(colSums(ambient$cluster_library_umi)))),
        warnings = as.list(warnings),
        settings = lapply(settings, qc_num),
        genes = list(total = length(genes), binarized = sum(binarized),
            not_binarized_evidence = sum(reason %in% "elpd"),
            not_binarized_ashman_D = sum(reason %in% "ashman_D"),
            not_binarized_chains = sum(reason %in% "chains"),
            bimodality_cut_z = signif(z_cut, 3),
            near_cut = if (is.na(z_cut)) 0L else sum(abs(z - z_cut) < 0.5)),
        convergence = list(global_not_converged = sum(conv$converged %in% FALSE),
            genes_with_cluster_over_cutoff = sum(max_spread > params$spread_cutoff, na.rm = TRUE),
            max_chain_spread_per_gene = qc_quant(max_spread, c(0.5, 0.9, 0.99, 1)),
            not_binarized_chains_first = as.list(head(genes[reason %in% "chains"], 20))),
        ambient = list(rho = if (!is.null(rho)) lapply(rho, signif, 4),
            off_changed_ambient = sum(flagged == CALL_FLAG_AMBIENT),
            on_changed_detection = sum(flagged == CALL_FLAG_DETECTION),
            calls_changed_chains = sum(flagged == CALL_FLAG_CHAINS),
            genes_off_changed = sum(rowSums(flagged == CALL_FLAG_AMBIENT) > 0),
            genes_on_changed = sum(rowSums(flagged == CALL_FLAG_DETECTION) > 0),
            changed_off_observed_over_expected = qc_quant(blocked_ratio, c(0.1, 0.5, 0.9)),
            on_call_detection = qc_quant(detect_on, c(0.01, 0.1, 0.5))),
        calls = list(ON = sum(on), OFF = sum(off), ambiguous = sum(amb),
            no_call = length(genes) * length(clusters) - n_called,
            on_genes_per_cluster = qc_quant(n_on_cl), ambiguous_genes_per_cluster = qc_quant(n_amb_cl),
            clusters_high_ambiguity = as.list(high_amb),
            on_clusters_per_binarized_gene = if (length(n_on_g) > 0) as.list(table(cut(n_on_g,
                c(-1, 0, 1, 2, 5, 20, 100, Inf),
                labels = c("0", "1", "2", "3-5", "6-20", "21-100", ">100"))))),
        limitations = as.list(QC_LIMITATIONS),
        provenance = list(R = R.version.string,
            scMarco = as.character(utils::packageVersion("scMarco")),
            rstan = if (requireNamespace("rstan", quietly = TRUE)) as.character(utils::packageVersion("rstan")),
            contamination = params$rho_source %||% ""))
    drop_null <- function(x) {
        if (!is.list(x)) return(x)
        x <- x[!vapply(x, is.null, logical(1))]
        lapply(x, drop_null)
    }
    qc <- drop_null(qc)
    if (is.null(file)) return(qc)

    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    yaml::write_yaml(qc, paste0(file, ".yaml"))
    plots <- qc_plots(z, z_cut, max_spread, binarized, table, flagged, on, off, amb, n_cells,
        params)
    qc_write_html(qc, plots, paste0(file, ".html"))
    message("Wrote ", file, ".yaml and ", file, ".html")
    invisible(qc)
}

qc_png <- function(expr, width = 640, height = 400) {
    # base R plot -> data URI of a PNG
    f <- tempfile(fileext = ".png")
    on.exit(unlink(f))
    type <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")
    grDevices::png(f, width = width, height = height, res = 96, type = type)
    op <- graphics::par(mar = c(4.2, 4.2, 2.5, 1))
    tryCatch(force(expr), finally = {
        graphics::par(op)
        grDevices::dev.off()
    })
    paste0("data:image/png;base64,", base64enc::base64encode(f))
}

qc_plots <- function(z, z_cut, max_spread, binarized, table, flagged, on, off, amb, n_cells,
                     params) {
    out <- list()
    if (any(is.finite(z))) out$z <- qc_png({
        zz <- pmin(pmax(z[is.finite(z)], -5), 15)
        graphics::hist(zz, breaks = 80, col = "grey70", border = NA, main = "",
            xlab = "z = elpd difference / SE (clipped to [-5, 15])", ylab = "Genes")
        if (is.finite(z_cut)) graphics::abline(v = z_cut, col = "red", lwd = 2)
        graphics::abline(v = 0, lty = 2)
    })
    if (any(is.finite(max_spread))) out$spread <- qc_png({
        h <- graphics::hist(max_spread[is.finite(max_spread)], breaks = seq(0, 1, 0.025), plot = FALSE)
        plot(h$mids, pmax(h$counts, 0.8), log = "y", type = "h", lwd = 6, lend = 1, col = "grey50",
            xlab = "Largest chain spread over clusters", ylab = "Genes (log scale)")
        graphics::abline(v = params$spread_cutoff, col = "red", lwd = 2)
    })
    if (isTRUE(params$ambient) && any(binarized)) out$levels <- qc_png({
        x <- table$bimodal.off_level[binarized]
        y <- table$bimodal.on_level[binarized]
        share <- function(changed, calls) {
            n <- rowSums(calls)[binarized] + rowSums(changed)[binarized]
            ifelse(n > 0, rowSums(changed)[binarized] / pmax(n, 1), NA)
        }
        s_on <- share(flagged == CALL_FLAG_DETECTION, on)
        s_off <- share(flagged == CALL_FLAG_AMBIENT, off)
        pal <- grDevices::colorRampPalette(c("#d9d9d9", "#fdae61", "#d7191c"))(101)
        col_of <- function(v) ifelse(is.na(v), "#d9d9d9", pal[1 + round(100 * pmin(pmax(v, 0), 1))])
        graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
        for (k in 1:2) {
            v <- if (k == 1) s_on else s_off
            o <- order(v, na.last = FALSE)
            plot(x[o], y[o], pch = 16, cex = 0.45, col = col_of(v[o]),
                xlab = "OFF level (log1p CP10k)", ylab = "ON level (log1p CP10k)",
                main = if (k == 1) "Detection check: share of the gene's\nON calls made ambiguous" else
                    "Ambient check: share of the gene's\nOFF calls made ambiguous", cex.main = 0.85)
            graphics::abline(0, 1, lty = 2, col = "grey40")
        }
        graphics::legend("bottomright", c("0%", "50%", "100%"), col = pal[c(1, 51, 101)], pch = 16,
            bty = "n", cex = 0.8, title = "share changed")
    }, width = 900, height = 460)
    if (!is.null(n_cells) && any(binarized)) out$clusters <- qc_png({
        graphics::par(mfrow = c(1, 2))
        plot(n_cells, colSums(on), log = "x", pch = 16, cex = 0.6, xlab = "Cells in cluster", ylab = "ON genes")
        plot(n_cells, colSums(amb), log = "x", pch = 16, cex = 0.6, xlab = "Cells in cluster", ylab = "Ambiguous genes")
    }, width = 760, height = 360)
    out
}

htmlesc <- function(s) {
    s <- gsub("&", "&amp;", s, fixed = TRUE)
    s <- gsub("<", "&lt;", s, fixed = TRUE)
    gsub(">", "&gt;", s, fixed = TRUE)
}

qc_fmt <- function(v) {
    # one value: integers without decimals, others to 3 significant digits
    if (is.null(v)) return("")
    if (is.numeric(v) && length(v) == 1) {
        if (is.finite(v) && abs(v - round(v)) < 1e-9 && abs(v) < 1e15) {
            return(formatC(round(v), format = "d", big.mark = ","))
        }
        return(format(signif(v, 3)))
    }
    paste(unlist(v), collapse = ", ")
}

qc_table <- function(rows, labels = NULL, header = c("", "Value"), extra = NULL) {
    # rows: named list of values; labels: key -> label; extra: key -> third column
    if (length(rows) == 0) return("")
    keys <- names(rows)
    lab <- if (is.null(labels)) keys else ifelse(keys %in% names(labels), labels[keys], keys)
    h <- paste0("<tr>", paste0("<th>", header, "</th>", collapse = ""), "</tr>")
    r <- vapply(seq_along(keys), function(i) {
        paste0("<tr><td>", htmlesc(lab[i]), "</td><td>", htmlesc(qc_fmt(rows[[i]])), "</td>",
            if (!is.null(extra)) paste0("<td>", htmlesc(if (keys[i] %in% names(extra)) extra[[keys[i]]] else ""), "</td>"),
            "</tr>")
    }, "")
    paste0("<table>", h, paste(r, collapse = ""), "</table>")
}

qc_write_html <- function(qc, plots, file) {
    title <- "scMarco binarization QC"
    sec <- function(h, body) paste0("<section><h2>", h, "</h2>", body, "</section>")
    img <- function(k, cap) {
        if (is.null(plots[[k]])) return("")
        paste0("<figure><img src='", plots[[k]], "' alt='", htmlesc(cap), "'><figcaption>",
            htmlesc(cap), "</figcaption></figure>")
    }
    triple <- function(q) paste(vapply(q, qc_fmt, ""), collapse = " / ")
    warn <- if (length(qc$warnings) > 0) {
        paste0("<h3>Warnings</h3><ul class='warn'>",
            paste0("<li>", htmlesc(unlist(qc$warnings)), "</li>", collapse = ""), "</ul>")
    } else "<p>No warnings.</p>"
    ds <- qc$dataset
    ds_rows <- ds[intersect(names(QC_LABELS$dataset), names(ds))]
    ds_rows <- ds_rows[vapply(ds_rows, function(v) !identical(v, ""), logical(1))]
    if (!is.null(ds$cells_per_cluster)) {
        ds_rows[["Cells per cluster (min / median / max)"]] <- triple(ds$cells_per_cluster)
    }
    ds_rows[["Replicates per cluster (min / median / max)"]] <- triple(ds$replicates_per_cluster)
    s <- qc$settings
    g <- qc$genes
    a <- qc$ambient
    checks <- if (is.null(s$on_floor)) "" else paste0(
        "<p>Two checks run after the model. <strong>ON needs detection:</strong> a cluster is ",
        "called ON only if at least ", qc_fmt(100 * s$on_floor), "% of its cells have one or more ",
        "UMIs of the gene. <strong>OFF needs absence:</strong> every cluster carries some ambient ",
        "RNA, the free-floating RNA that ends up in all droplets. For each gene and cluster, the ",
        "expected ambient UMIs are the contamination of each library times the gene's share of that ",
        "library's UMIs times the cluster's UMIs from that library, summed over libraries. An OFF ",
        "call becomes ambiguous only if the cluster's UMIs are significantly more than ",
        qc_fmt(2^s$ambient_log2fc), " times that expectation (log2 fold > ", qc_fmt(s$ambient_log2fc),
        "; Poisson test at level ", qc_fmt(s$ambient_alpha), "). This threshold is deliberately ",
        "maximally permissive: it only catches clear expression, and smaller excesses are left to ",
        "the leak ranking of the pair search. Calls that fail either check become ambiguous, which ",
        "counts as a possible off-target downstream. Large clusters give these tests more power ",
        "than small ones.</p>")
    body <- paste0(
        "<header><h1>", title, "</h1><p>", htmlesc(paste(ds$name %||% "", ds$split %||% "")),
        " &middot; QC report version ", qc$qc_report_version, " &middot; ", qc$generated,
        "</p></header>",
        sec("Dataset", paste0(qc_table(ds_rows, QC_LABELS$dataset),
            if (!is.null(ds$umis_per_library)) paste0("<h3>UMIs per library</h3>",
                qc_table(ds$umis_per_library, header = c("Library", "UMIs"))),
            warn)),
        sec("Settings in use", qc_table(s, header = c("Setting", "Value", "What is this"),
            extra = QC_SETTING_INFO)),
        sec("Genes", paste0(
            qc_table(g[intersect(names(QC_LABELS$genes), names(g))], QC_LABELS$genes),
            "<p>Genes near the cut can switch between binarized and not binarized from one run to the next: ",
            QC_NEAR_CUT_REFERENCE, ".</p>",
            img("z", "Evidence for two states per gene (z = elpd difference / SE); red: the cut at the chosen false discovery rate"))),
        sec("Convergence", paste0(
            qc_table(qc$convergence[intersect(names(QC_LABELS$convergence), names(qc$convergence))],
                QC_LABELS$convergence),
            img("spread", "Largest chain spread per gene (log scale); red: the spread cutoff"),
            if (length(qc$convergence$not_binarized_chains_first) > 0) paste0(
                "<p>Genes not binarized because their chains disagree (first 20): ",
                htmlesc(paste(unlist(qc$convergence$not_binarized_chains_first), collapse = ", ")), "</p>"))),
        sec("Ambient RNA and detection", paste0(checks,
            if (!is.null(a$rho)) paste0("<h3>Contamination per library (share of UMIs from ambient RNA)</h3>",
                qc_table(a$rho, header = c("Library", "Contamination"))),
            qc_table(a[intersect(names(QC_LABELS$ambient), names(a))], QC_LABELS$ambient),
            img("levels", "Each binarized gene's fitted OFF and ON levels, coloured by the share of its ON calls (left) or OFF calls (right) made ambiguous"))),
        sec("Calls", paste0(
            qc_table(qc$calls[intersect(names(QC_LABELS$calls), names(qc$calls))], QC_LABELS$calls),
            if (!is.null(qc$calls$on_clusters_per_binarized_gene)) paste0("<h3>ON clusters per binarized gene</h3>",
                qc_table(qc$calls$on_clusters_per_binarized_gene, header = c("ON clusters", "Genes"))),
            img("clusters", "Per cluster: ON genes and ambiguous genes against cluster size"),
            if (length(qc$calls$clusters_high_ambiguity) > 0) paste0(
                "<p>Clusters with unusually many ambiguous calls: ",
                htmlesc(paste(unlist(qc$calls$clusters_high_ambiguity), collapse = ", ")), "</p>"))),
        sec("Limitations", paste0("<ul>", paste0("<li>", htmlesc(unlist(qc$limitations)), "</li>",
            collapse = ""), "</ul>")),
        sec("Provenance", qc_table(qc$provenance)))
    css <- paste0(":root{--bg:#fff;--fg:#1d1d1f;--muted:#666;--line:#ddd;--warn:#a33}",
        "@media (prefers-color-scheme: dark){:root{--bg:#141414;--fg:#eee;--muted:#aaa;--line:#333;--warn:#e88}}",
        "body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,sans-serif;max-width:960px;margin:0 auto;padding:16px}",
        "h1{font-size:1.5em}h2{border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:2em}",
        "header p{color:var(--muted)}table{border-collapse:collapse;margin:8px 0}",
        "th{text-align:left;color:var(--muted);font-weight:600;padding:3px 10px;border-bottom:2px solid var(--line)}",
        "td{border-bottom:1px solid var(--line);padding:3px 10px;vertical-align:top}",
        "td:nth-child(2){font-variant-numeric:tabular-nums;white-space:nowrap}",
        "td:nth-child(3){color:var(--muted)}img{max-width:100%;height:auto;background:#fff}figure{margin:12px 0}",
        "figcaption{color:var(--muted);font-size:.9em}.warn li{color:var(--warn)}")
    html <- paste0("<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'><title>", title,
        "</title><style>", css, "</style></head><body>", body, "</body></html>")
    writeLines(html, file)
    invisible(file)
}
