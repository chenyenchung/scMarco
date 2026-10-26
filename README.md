# scMarco

An R package and Shiny app for selecting marker gene combinations from
mixture-model binarized scRNA-seq expression. This codebase merges the two
earlier forks:

- `scMarco` (optic lobe): SQLite backend, ordered stages, expression-trend plot;
- `scMarco-Ding-lab` (dsx neurons): HDF5 backend, sex/species splits, UMAPs,
  % positive overlays, non-bimodal gene handling.

Every lab runs the same code with its own `config.yaml` and data folder. Tabs
appear only when a dataset has the data they need. The searches and plots
behind the app can also be used from R; see `vignette("scMarco")`.

## Installing

```r
# install.packages("BiocManager")
BiocManager::install("chenyenchung/scMarco")
```

For development, `renv::restore()` restores the pinned library and
`R CMD INSTALL .` installs the package from a checkout.

## Running

```r
library(scMarco)
run_scmarco("config.yaml")          # paths in the config relative to the working directory
run_scmarco("lab.yaml", root = "/srv/scmarco")
```

The lab instances live in `deploy/` (not part of the package): `app.R`,
the configs, their README files, gene lists and the (git-ignored) data.

```sh
cd deploy
SCMARCO_CONFIG=configs/optic_lobe.yaml Rscript -e 'shiny::runApp(".")'
SCMARCO_CONFIG=configs/dinglab.yaml    Rscript -e 'shiny::runApp(".")'
Rscript -e 'shiny::runApp(".")'        # config.yaml: every dataset
```

For deployment (e.g. shinyapps.io), copy the lab's config to `config.yaml`
and ship `app.R`, `readme/`, `gene_lists/` and its `data/*.h5ad` files;
`app.R` is `scMarco::scmarco_app()`, so the server needs the package.

## Config

```yaml
title: "scMarco"
theme: minty                 # any bootswatch theme
readme:                      # Markdown files shown in "Read Me", in order,
  - readme/optic_lobe.md     # followed by the packaged usage guide
usage_guide: true            # false: leave the usage guide out
datasets:
  ozel_2021:                 # id
    label: "Optic lobe, multiome annotation: ..."
    clusters: data/ozel_2021.clusters.h5ad       # required
    cells: data/ozel_2021.cells.h5ad             # optional
    default: true            # selected at startup
    split_order: [...]       # optional overrides of the file's uns/scmarco
    split_label: stage
    ordered_splits: true
    trend: line              # optional: line (needs mean_lognorm), box (needs
                             # cells) or both; default line for ordered splits
    subset_by: [class, subclass]  # optional: obs/meta_<level> columns users can
                                  # subset clusters by, coarse to fine
gene_list:
  flybase_tf: {path: gene_lists/flybase_tfs.txt, label: "FlyBase TFs"}
```

## Data format (scMarco h5ad, format_version 2.0)

Each dataset is one or two [AnnData](https://anndata.readthedocs.io) files.

**`<id>.clusters.h5ad`** (required): one observation per split × cluster
that exists.

| slot | content |
|---|---|
| `obs` | index `"<split>\|<cluster>"`; `split` (categorical, ordered when splits are ordered), `cluster` (categorical cluster ID), optional `cluster_name` (categorical; missing = unnamed), optional `meta_<level>` (categorical cluster metadata; the config's `subset_by` picks which levels limit the searches), optional `n_cells` |
| `var` | index = genes (sorted); optional `is_bimodal` (bool): genes used for marker finding |
| `varm/bimodal_by_split` | optional genes × splits 0/1, column order in `uns/scmarco/bimodal_splits` |
| `varm/mixture_on_level`, `mixture_off_level`, `mixture_sd`, `mixture_pi` | optional genes × splits mixture-model fit (ON and OFF component means, shared SD, ON weight), same column order; enables the Mixture model tab |
| `layers/prob` | mixture-model probability of expression (NaN = not modeled / non-bimodal) |
| `layers/pct_pos` | optional fraction of cells with > 0 UMIs |
| `layers/mean_lognorm` | optional cluster-average log-normalized expression (trend plot) |
| `uns/scmarco` | `format_version`, `split_label`, `split_order`, `ordered_splits`, `source`, `binarization` (`method`, `source`, `date`, `version`) once binarized |

A cluster absent from a split simply has no row, so missing values never mean
two things.

`cluster` is a **stable ID** (e.g. `154`): every join uses it, and it never
changes when a cluster is annotated. The annotation lives in `cluster_name`
(the same in every split of a cluster), and the app shows `name (id)`, or the
bare ID while a cluster is unnamed. Change names with
`set_cluster_names()` rather than editing the file: it checks the table
and rewrites only `obs`.

**`<id>.cells.h5ad`** (optional): enables the UMAP tab and the true
co-expression percentage.

| slot | content |
|---|---|
| `X` | log-normalized expression, cells × genes, **CSC** (fast per-gene reads) |
| `obs` | index `"<split>\|<barcode>"`; `barcode`, `split`, `cluster` (ID; names come from the clusters file) |
| `obsm/X_<name>` | reductions, cells × 2; NaN for cells outside a subset UMAP; `X_umap_full` is shown as "All cells" |

Files are written with anndataR and then re-chunked with `h5repack`, so
per-gene reads touch only small blocks. The app reads them with a small rhdf5
reader (`R/backend_h5ad.R`); Python is not needed at runtime.

## Building datasets

New datasets are built and edited with the package's functions:

| function | use |
|---|---|
| `seurat_to_h5ad()` | new data: Seurat object (split and cluster ID in metadata columns), or one object per split (`c(P15 = "a.rds", P30 = "b.rds")`), + optional subset UMAPs → a *built* dataset; `prob` also binarizes it |
| `add_probability()` | binarize a built dataset with the mixture-model output, one table per split; the app needs this step |
| `get_cluster_names()` / `set_cluster_names()` | read, export or apply a `cluster_id`/`name` table |
| `get_cluster_metadata()` / `set_cluster_metadata()` | read, export or replace a `cluster_id`/`<level>`/... table of cluster groupings (e.g. class, subclass) that lets users limit searches to part of the dataset |
| `set_cluster_ids()` | give every cluster a new ID (and name) from a `cluster_id`/`new_id`/`name` table |
| `migrate_cluster_labels()` | turn legacy `name (id)` labels into ID + name |
| `validate_h5ad(clusters, cells, stage)` | check a dataset before deploying (`stage = "built"` for one without probabilities yet); `add_probability()` runs it after writing |
| `sqlite_to_h5ad()` | migrate a legacy scMarco SQLite table |
| `h5v1_to_h5ad()` | migrate a legacy scMarco-Ding-lab `db.h5` dataset |
| `find_pairs()` | offline ranking of marker pairs for every cluster (table with `cluster` ID and `cluster_name`; `out` also writes a CSV) |

Each function has a help page (`?seurat_to_h5ad`). Inputs can be file paths
or objects already in memory. `deploy/misc/` holds one-off helpers for
specific datasets; `source()` one after `library(scMarco)` when needed.
Building needs the HDF5 command-line tools (`h5repack`) for re-chunking.

## Code layout

| path | role |
|---|---|
| `R/combo_core.R` | pure combination logic (From cluster, From gene, distinct genes, pairs) |
| `R/backend_h5ad.R` | h5ad reader and per-process dataset cache |
| `R/plot_core.R` | On/Off, co-expression, trend and UMAP plots |
| `R/module_*.R` | one Shiny module per tab |
| `R/app.R` | config, sidebar (dataset, genes, clusters, threshold, split criteria), tab gating |
| `R/build.R`, `R/import.R`, `R/write_h5ad.R`, `R/validate.R` | building, binarizing, annotating and checking datasets |
| `R/migrate.R`, `R/pairs.R` | legacy converters; offline pair ranking |
| `inst/extdata/` | example datasets (`inst/scripts/` regenerates them) |
| `vignettes/` | the package vignette |
| `deploy/` | lab instances, data and data-dependent checks (not in the package) |
| `docs/` | the user guide website (Quarto) |

## Tests

```sh
Rscript -e 'testthat::test_local()'                 # unit tests (package fixtures only)
R CMD build . && R CMD check scMarco_*.tar.gz       # full package check
Rscript -e 'BiocCheck::BiocCheck("scMarco_0.99.0.tar.gz")'
cd deploy                                           # needs the converted data:
SCMARCO_CONFIG=configs/all.yaml Rscript tests/integration/app_flow.R
Rscript tests/parity/parity.R ../../scMarco ../../scMarco-Ding-lab
```

## Behavior changes from the legacy forks

- Split criteria are global: "must be ON in these splits" and "…in at least N".
- Thresholds are strict (`prob > cut`) everywhere; legacy scMarco's On/Off plot used `>=`.
- From gene with several genes means *all* genes (or explicitly *any*), not a
  pooled count; legacy scMarco errored and legacy Ding-lab pooled genes across splits.
- The co-expression % is the share of cells expressing *all* genes (was the
  mean of per-gene percentages).
- Clusters absent from a split are drawn as a cross, not as OFF.
- Find distinct ranks hits by off-targets outside the selection, plots the top N
  per cluster and checks the combination against the whole dataset.
