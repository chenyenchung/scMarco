# To-do: finish the harmonized scMarco

The repo is an R package (`R/`, `man/`, `inst/`, `vignettes/`,
`tests/testthat/`). The lab instances, converted data, `misc/` helpers and the
data-dependent integration/parity scripts live in `deploy/`. The input formats
are specified in [`docs/expected_input.md`](docs/expected_input.md).

## Bioconductor submission

- [x] Package skeleton: DESCRIPTION (0.99.0, biocViews), NAMESPACE via
  roxygen, help pages with runnable examples for every export, NEWS.md,
  BiocStyle vignette, testthat suite, small example data with generating
  scripts in `inst/scripts/`.
- [x] `R CMD check`: Status OK. BiocCheck: no warnings.
- [ ] Register `chenyenjung@gmail.com` on https://support.bioconductor.org and
  add the `scmarco` watched tag (the one BiocCheck ERROR), and subscribe to
  bioc-devel.
- [x] **Binarizer license**: Fred P. Davis permitted the adaptation of his
  opticlobe code and its MIT licensing (email, 2026-10-26, archived). He is
  `ctb`, `cph` in `Authors@R`; notices in the file headers, `LICENSE`,
  `LICENSE.md` and `inst/COPYRIGHTS`.
- [ ] Optional: ORCID iD in `Authors@R`; `fnd` role for funding.
- [ ] Consider an interface from `SingleCellExperiment` (reviewers tend to ask
  for interoperability with Bioconductor classes): e.g. build a dataset from
  an SCE as `seurat_to_h5ad()` does from Seurat.
- [ ] Remaining BiocCheck NOTEs: long functions (Shiny servers), lines > 80
  characters (most in `R/build.R`, `R/import.R`, `R/plot_core.R`), `<<-` in
  `validate_h5ad()` (`R/validate.R`) and the app server (`R/app.R`).
- [ ] Update the Quarto user guide (`docs/src/config.qmd` still describes the
  legacy SQLite config) and publish the GitHub repo, then submit at
  https://github.com/Bioconductor/Contributions.

## Migration (two paths only)

Only the Ding-lab production HDF5 v1.0 (`db.h5`) and the legacy scMarco SQLite
are supported; the internal `dinglab.sqlite` and the pre-v1 HDF5 are out of
scope. Both converters are in `R/migrate.R`.

- [ ] **Legacy SQLite** (`sqlite_to_h5ad()`): take
  `legacy_config = "../scMarco/config.yaml"` so stages and labels come from the
  old config instead of a manual `stages`. Add a way to migrate every table
  listed there.
- [ ] **Ding-lab HDF5 v1.0** (`h5v1_to_h5ad()`): handle every dataset in
  `db.h5`, taking `split_label` from the legacy Ding config. Treat `lognorm` as
  optional, as the v1 spec does: with no `lognorm`, write no cells file and no
  `mean_lognorm` (`pct_pos` and `reduction` are already optional). Fix the
  per-split "non-bimodal" log line so it reports the value after the
  `is_bimodal` fix.
- [ ] `migrate_legacy()`: one function that takes `sqlite = <db>` or
  `h5 = <db.h5>`, plus `legacy_config = <yaml>`. It converts the data, copies
  the gene lists the old config references into `gene_lists/`, writes
  `configs/<name>.yaml`, and runs the validator.
- [ ] **Value round-trip check** (`deploy/tests/parity/roundtrip.R`), alongside
  the logic check `deploy/tests/parity/parity.R`, which already passes. Sample
  random cells and entries and compare legacy and new values for:
  - `prob` and `pct_pos`;
  - SQLite `lognorm` vs `mean_lognorm`;
  - v1 CSR vs `cells.h5ad` per-gene expression;
  - UMAP coordinates, including `cell_indices` subsets;
  - per-cell cluster labels.
- [ ] Converter unit test for a subset UMAP with `cell_indices`
  (`tests/testthat/test-migrate.R` covers the rest).
- [ ] Document known migration facts: the optic-lobe data has no bimodality
  information (every gene is used, as before), and `ozel_2021` is deprecated but
  kept.

## Wrap-up

- [ ] Link `docs/expected_input.md` from `README.md`.
- [ ] Re-run all checks:
  - `Rscript -e 'testthat::test_local()'`
  - from `deploy/`: `SCMARCO_CONFIG=configs/all.yaml Rscript tests/integration/app_flow.R`
  - from `deploy/`: `Rscript tests/parity/parity.R ../../scMarco ../../scMarco-Ding-lab`
  - from `deploy/`: `Rscript tests/parity/roundtrip.R ../../scMarco ../../scMarco-Ding-lab`
- [ ] Click through in a browser once for each lab config.
- [ ] Decide where the converted `deploy/data/*.h5ad` files live (they are
  git-ignored).

## Done when

- The unit tests pass.
- `migrate_legacy()` on `../scMarco` and `../scMarco-Ding-lab` reproduces the
  current `deploy/data/*.h5ad`, and both parity scripts report 0 mismatches.
- `tests/integration/app_flow.R` passes on every lab config.

## Done

- Expected-input spec (`docs/expected_input.md`), kept in sync with
  `R/import.R` and `R/build.R`.
- Two-step build: `seurat_to_h5ad()` writes a *built* dataset (optional
  `prob`, `split_col`/`cluster_col`, `<split>@<cluster>` fallback); `add_probability()` binarizes it from per-split tables, a
  named-list RDS or the legacy wide RDS, then validates with
  `stage = "binarized"`.
- `validate_h5ad(stage = "built" | "binarized")` with the Seurat-derived checks.
- The app checks every configured dataset at startup and lists each one
  missing `prob` with the `add_probability()` fix.
- Tests for the build/binarize/validate/load flow, bad inputs, and the
  converters on tiny synthetic SQLite and v1 `db.h5` inputs.

## Backlog (future)

- [ ] **Backend updater** (like `UpdateSeuratObject()`), needed with the first
  bump of `SCMARCO_FORMAT_VERSION`: rewrite an older file to the current
  schema in place. Then point the "rebuild the dataset" message of
  `ds_check_config()` (`format_problem()` in `R/backend_h5ad.R`) at it.

- [ ] **Binarizer input**: the mixture model is now in the package
  (`fit_mixture()`, `call_mixture()`), but its input, log-normalized
  pseudobulk replicates of each cluster, is still made outside it
  (mm_marker's `get_input.R`: `Seurat::AggregateExpression()` on each cluster's
  cells split into two halves by library). `make_pseudobulk()` now pools a
  counts matrix as `AggregateExpression()` does (random parts, metadata groups
  or k-means subclusters); still to add: taking a Seurat object directly, and a
  way to run gene chunks as parallel jobs.
- [ ] **Cluster heterogeneity**: try `mixture_heterogeneity()`'s cutoffs
  (`z_cutoff = 3`, `min_spread = 0.25`) on real data, e.g. clusters known to
  merge two cell types.
