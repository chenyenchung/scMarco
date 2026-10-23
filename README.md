# scMarco (harmonized)

A Shiny app for selecting marker gene combinations from mixture-model
binarized scRNA-seq expression. This codebase merges the two earlier forks:

- `scMarco` (optic lobe): SQLite backend, ordered stages, expression-trend plot;
- `scMarco-Ding-lab` (dsx neurons): HDF5 backend, sex/species splits, UMAPs,
  % positive overlays, non-bimodal gene handling.

Every lab runs the same code with its own `config.yaml` and data folder. Tabs
appear only when a dataset has the data they need.

## Running

```sh
# once: restore the R library (renv) and, optionally, Python for checks
Rscript -e 'renv::restore()'
uv venv .venv && uv pip install --python .venv/bin/python anndata h5py

# run an instance
SCMARCO_CONFIG=configs/optic_lobe.yaml Rscript -e 'shiny::runApp(".")'
SCMARCO_CONFIG=configs/dinglab.yaml    Rscript -e 'shiny::runApp(".")'
Rscript -e 'shiny::runApp(".")'        # config.yaml: every dataset
```

For deployment (e.g. shinyapps.io), copy the lab's config to `config.yaml`
and ship its `data/*.h5ad` files alongside the app.

## Config

```yaml
title: "scMarco"
theme: minty                 # any bootswatch theme
readme:                      # Markdown files shown in "Read Me", in order
  - readme/optic_lobe.md
  - readme/usage.md
datasets:
  ozel_2021:                 # id
    label: "Optic lobe, multiome annotation: ..."
    clusters: data/ozel_2021.clusters.h5ad       # required
    cells: data/ozel_2021.cells.h5ad             # optional
    default: true            # selected at startup
    split_order: [...]       # optional overrides of the file's uns/scmarco
    split_label: stage
    ordered_splits: true
gene_list:
  flybase_tf: {path: gene_lists/flybase_tfs.txt, label: "FlyBase TFs"}
```

## Data format (scMarco h5ad, format_version 2.0)

Each dataset is one or two [AnnData](https://anndata.readthedocs.io) files.

**`<id>.clusters.h5ad`** (required): one observation per split × cluster
that exists.

| slot | content |
|---|---|
| `obs` | index `"<split>\|<cluster>"`; `split` (categorical, ordered when splits are ordered), `cluster` (categorical cluster ID), optional `cluster_name` (categorical; missing = unnamed), optional `n_cells` |
| `var` | index = genes (sorted); optional `is_bimodal` (bool): genes used for marker finding |
| `varm/bimodal_by_split` | optional genes × splits 0/1, column order in `uns/scmarco/bimodal_splits` |
| `varm/mixture_on_level`, `mixture_off_level`, `mixture_sd`, `mixture_pi` | optional genes × splits mixture-model fit (ON and OFF component means, shared SD, ON weight), same column order; enables the Mixture model tab |
| `layers/prob` | mixture-model probability of expression (NaN = not modeled / non-bimodal) |
| `layers/pct_pos` | optional fraction of cells with > 0 UMIs |
| `layers/mean_lognorm` | optional cluster-average log-normalized expression (trend plot) |
| `uns/scmarco` | `format_version`, `split_label`, `split_order`, `ordered_splits`, `source` |

A cluster absent from a split simply has no row, so missing values never mean
two things.

`cluster` is a **stable ID** (e.g. `154`): every join uses it, and it never
changes when a cluster is annotated. The annotation lives in `cluster_name`
(the same in every split of a cluster), and the app shows `name (id)`, or the
bare ID while a cluster is unnamed. Change names with
`script/set_cluster_names.R` rather than editing the file: it checks the table
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
reader (`src/backend_h5ad.R`); Python is not needed at runtime.

## Building datasets (`script/`)

| script | use |
|---|---|
| `seurat_to_h5ad.R` | new data: Seurat object (split and cluster ID in metadata columns), or one object per split (`--seurat P15=a.rds,P30=b.rds`), + optional subset UMAPs → a *built* dataset; `--prob` also binarizes it |
| `add_probability_h5ad.R` | binarize a built dataset with the mixture-model output, one table per split; the app needs this step |
| `sqlite_to_h5ad.R` | migrate a legacy scMarco SQLite table |
| `h5v1_to_h5ad.R` | migrate a legacy scMarco-Ding-lab `db.h5` dataset |
| `set_cluster_names.R` | export or apply a `cluster_id`/`name` table; `--ids` gives every cluster a new ID (and name) from a `cluster_id`/`new_id`/`name` table; `--split-labels TRUE` turns legacy `name (id)` labels into ID + name |
| `validate_h5ad.R` | check a dataset before deploying (`--stage built` for one without probabilities yet) |
| `find_pairs.R` | offline ranking of marker pairs for every cluster (CSV with `cluster` ID and `cluster_name`) |

Each script prints its usage at the top of the file. The builders need the HDF5
command-line tools (`h5repack`) for re-chunking.

## Code layout

| file | role |
|---|---|
| `src/combo_core.R` | pure combination logic (From cluster, From gene, distinct genes, pairs); unit tested |
| `src/backend_h5ad.R` | h5ad reader and per-process dataset cache |
| `src/plot_core.R` | On/Off, co-expression, trend and UMAP plots |
| `src/module_*.R` | one Shiny module per tab |
| `app.R` | config, sidebar (dataset, genes, clusters, threshold, split criteria), tab gating |

## Tests

```sh
Rscript -e 'testthat::test_dir("tests/testthat")'              # unit tests (fixtures only)
SCMARCO_CONFIG=configs/all.yaml Rscript tests/integration/app_flow.R  # full workflow per dataset
Rscript tests/parity/parity.R ../scMarco ../scMarco-Ding-lab    # vs. legacy implementations
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
