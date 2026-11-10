# Generates inst/extdata/toy.clusters.h5ad and toy.cells.h5ad: a small
# synthetic dataset with every optional part of the scMarco format (cells
# file, reduction, cluster names and metadata, mixture-model fit
# parameters), used by the examples and the vignette.
#
#   Rscript inst/scripts/make_toy_data.R   # from the package root, with
#                                          # scMarco installed
#
# Two splits (Male, Female), eight clusters (cluster 8 has no Female cells),
# 40 genes. Counts are Poisson: each cluster expresses its own set of marker
# genes, so the combination searches have something to find. The
# probabilities come from a two-component Gaussian mixture fit per gene and
# split on cluster averages, as the mixture-model pipeline would produce.

library(scMarco)
set.seed(20261026)
out_dir <- file.path("inst", "extdata")

genes <- sort(c(
    "acj6", "ap", "ase", "bsh", "brp", "cas", "ChAT", "D", "dac", "dsx",
    "elav", "erm", "ey", "fru", "Gad1", "grh", "hb", "hth", "Isl", "Kr",
    "Lim3", "nSyb", "Optix", "pdm2", "ple", "pros", "repo", "Rx", "side",
    "svp", "Tdc2", "tll", "toy", "Trh", "VGlut", "vvl", "dpr8", "beat-Ia",
    "DIP-beta", "CG4328"
))
clusters <- as.character(1:8)
cluster_names <- c("Mi1", "Tm3", "T4", "Pm2", "Lawf1", "C3", "Dm8", "glia")
splits <- c("Male", "Female")

# Marker genes of each cluster: on in those clusters, off elsewhere
on_in <- list(
    Kr = c("1", "2"), hb = "1", pdm2 = c("2", "3"), cas = "3", grh = c("3", "4"),
    acj6 = "4", ey = c("5", "6"), toy = "5", Lim3 = c("6", "7"), Isl = "7",
    vvl = c("1", "5"), bsh = c("1", "6"), ap = c("2", "7"), dac = c("3", "5"),
    Optix = "6", Rx = c("4", "7"), hth = c("1", "4", "6"), tll = "2",
    svp = c("5", "7"), D = c("3", "6"), erm = "4", ase = c("2", "5"),
    repo = "8", Gad1 = c("2", "4", "6"), VGlut = c("1", "7"),
    ChAT = c("3", "5"), Tdc2 = "7", ple = "6", Trh = "1", side = c("2", "3"),
    dpr8 = c("4", "5"), `beat-Ia` = c("6", "7"), `DIP-beta` = c("1", "3")
)
# Sex-specific genes: dsx and fru only in male neurons
male_only <- list(dsx = c("2", "5", "7"), fru = c("1", "2", "5"))
# Pan-neuronal genes, and genes that are never bimodal (no clear ON/OFF)
pan_neuronal <- c("elav", "nSyb", "brp", "pros")

cells_per <- 30
cell_split <- character(0)
cell_cluster <- character(0)
for (s in splits) {
    for (cl in clusters) {
        if (s == "Female" && cl == "8") next
        cell_split <- c(cell_split, rep(s, cells_per))
        cell_cluster <- c(cell_cluster, rep(cl, cells_per))
    }
}
n_cells <- length(cell_split)
barcodes <- sprintf("cell%04d", seq_len(n_cells))

rate <- matrix(0.05, n_cells, length(genes), dimnames = list(barcodes, genes))
for (g in names(on_in)) rate[cell_cluster %in% on_in[[g]], g] <- 4
for (g in names(male_only)) {
    rate[cell_cluster %in% male_only[[g]] & cell_split == "Male", g] <- 4
}
rate[cell_cluster != "8", pan_neuronal] <- 6
rate[, "CG4328"] <- 1 # flat: unimodal everywhere
counts <- matrix(rpois(length(rate), rate * runif(n_cells, 0.6, 1.4)),
                 n_cells, dimnames = dimnames(rate))
lib_size <- rowSums(counts) + 200 # unmeasured genes
lognorm <- log1p(counts / lib_size * 1e4)
x <- Matrix::Matrix(lognorm, sparse = TRUE)

# Cluster-level layers: one row per split x cluster with cells
obs <- unique(data.frame(split = cell_split, cluster = cell_cluster))
key <- paste(cell_split, cell_cluster, sep = "|")
obs_key <- paste(obs$split, obs$cluster, sep = "|")
pct_pos <- t(vapply(obs_key, function(k) colMeans(counts[key == k, ] > 0),
                    numeric(length(genes))))
mean_lognorm <- t(vapply(obs_key, function(k) colMeans(lognorm[key == k, ]),
                         numeric(length(genes))))
n_obs_cells <- as.integer(table(factor(key, levels = obs_key)))

clusters_path <- file.path(out_dir, "toy.clusters.h5ad")
cells_path <- file.path(out_dir, "toy.cells.h5ad")
scMarco:::write_clusters_h5ad(
    clusters_path, obs$split, obs$cluster, splits, clusters, genes,
    layers = list(pct_pos = unname(pct_pos), mean_lognorm = unname(mean_lognorm)),
    n_cells = n_obs_cells, split_label = "sex", ordered_splits = FALSE,
    extra_uns = list(source = "synthetic: inst/scripts/make_toy_data.R")
)

# A UMAP: one blob per cluster
centers <- cbind(cos(seq_along(clusters) * 2 * pi / 8),
                 sin(seq_along(clusters) * 2 * pi / 8)) * 6
ci <- match(cell_cluster, clusters)
umap <- centers[ci, ] + matrix(rnorm(2 * n_cells), n_cells)
scMarco:::write_cells_h5ad(
    cells_path, x, barcodes, cell_split, cell_cluster, splits, clusters, genes,
    reductions = list(umap_full = umap)
)

# Mixture-model output, one table per split: P(On) per gene and cluster from
# a two-Gaussian fit on the log of the mean expression of each cluster
mixture_table <- function(s) {
    rows <- which(obs$split == s)
    avg <- vapply(rows, function(r) {
        log1p(colMeans(expm1(lognorm[key == obs_key[r], , drop = FALSE])))
    }, numeric(length(genes)))
    colnames(avg) <- obs$cluster[rows]
    tab <- matrix(NA_real_, length(genes), length(rows),
                  dimnames = list(genes, obs$cluster[rows]))
    na <- rep(NA_real_, length(genes))
    qc <- data.frame(row.names = genes, exprType = rep("unimodal", length(genes)),
                     bimodal.on_level = na, bimodal.off_level = na,
                     bimodal.sd = na, bimodal.pi = na)
    for (g in genes) {
        v <- avg[g, ]
        on <- v > (max(v) + min(v)) / 2
        if (max(v) - min(v) < 1 || all(on) || !any(on)) next
        sd <- max(sqrt((sum((v[on] - mean(v[on]))^2) +
                        sum((v[!on] - mean(v[!on]))^2)) / length(v)), 0.1)
        fit <- list(on_level = mean(v[on]), off_level = mean(v[!on]), sd = sd,
                    pi = mean(on))
        d_on <- fit$pi * dnorm(v, fit$on_level, sd)
        d_off <- (1 - fit$pi) * dnorm(v, fit$off_level, sd)
        tab[g, ] <- d_on / (d_on + d_off)
        qc[g, ] <- list("bimodal", fit$on_level, fit$off_level, sd, fit$pi)
    }
    cbind(as.data.frame(tab, check.names = FALSE), qc)
}
# made here, not by call_mixture(): recorded as an external binarizer
add_probability(clusters_path, lapply(c(Male = "Male", Female = "Female"),
                                      mixture_table), external = TRUE)

set_cluster_names(clusters_path, setNames(cluster_names, clusters))
set_cluster_metadata(clusters_path, data.frame(
    cluster_id = clusters,
    class = c(rep("neuron", 7), "glia"),
    neuropil = c("medulla", "medulla", "lobula plate", "medulla", "lamina",
                 "medulla", "medulla", NA)
))
report_validation(validate_h5ad(clusters_path, cells_path), clusters_path,
                  cells_path)
