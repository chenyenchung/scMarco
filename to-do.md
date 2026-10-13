# To-do: finish the harmonized scMarco

Status (2026-09-28): the app runs end to end on all five converted datasets.
The logic parity check against both legacy forks passes 88/88, and the unit and
integration tests are green. This list covers what is left for:

1. **Migration**, from the two supported legacy sources only: the Ding-lab
   production HDF5 v1.0 (`db.h5`) and the legacy scMarco SQLite. The internal
   `dinglab.sqlite` and the pre-v1 HDF5 are out of scope.
2. **Building backend files without the binarization data** from a Seurat
   object, with the mixture-model probabilities added in a separate step. A
   dataset without `layers/prob` makes the app **error out**.

The input formats are specified in [`docs/expected_input.md`](docs/expected_input.md).

## 0. Expected input spec

- [x] Write `docs/expected_input.md` (Seurat object, probability tables, file states).
- [ ] Keep it in sync with `check_seurat_input()` and `add_probability_h5ad.R` as they are implemented.

## 1. Shared plumbing

- [ ] `script/lib_write_h5ad.R`: make `prob` optional in `write_clusters_h5ad()`.
  Add `add_probability_layer(path, prob, is_bimodal, bimodal_by_split, meta)`,
  which reads with `anndataR::read_h5ad`, adds the layer, `var` and `varm`,
  rewrites the file, then reuses `rechunk_h5ad()`.
- [ ] `script/validate_h5ad.R`: add `--stage built|binarized` (default
  `binarized`). A missing `prob` is an error for `binarized` and fine for
  `built`. Add the Seurat-derived checks from the spec (name rules, obs
  index separator, `uns/scmarco/binarization`).
- [ ] `src/backend_h5ad.R` `ds_load()` and `app.R`: at startup, check **every**
  configured dataset, then stop with one message that lists each dataset
  missing `prob` and the fix (`Rscript script/add_probability_h5ad.R ...`).
  Today `ds_load()` stops on the first such dataset with a terse message.

## 2. Migration (two paths only)

- [ ] **Legacy SQLite** (`script/sqlite_to_h5ad.R`): take
  `--legacy-config ../scMarco/config.yaml` so stages and labels come from the
  old config instead of a manual `--stages`. Add `--all` to migrate every table
  listed there.
- [ ] **Ding-lab HDF5 v1.0** (`script/h5v1_to_h5ad.R`): add `--all` for every
  dataset in `db.h5`, taking `split_label` from the legacy Ding config. Treat
  `pct_pos`, `lognorm` and `reduction` as optional, as the v1 spec does: with no
  `lognorm`, write no cells file and no `mean_lognorm`. Fix the per-split
  "non-bimodal" log line so it reports the value after the `is_bimodal` fix.
- [ ] `script/migrate_legacy.R`: one wrapper that takes
  `--sqlite <db> --legacy-config <yaml>` or `--h5 <db.h5> --legacy-config <yaml>`.
  It converts the data, copies the gene lists the old config references into
  `gene_lists/`, writes `configs/<name>.yaml`, and runs the validator.
- [ ] **Value round-trip check** (`tests/parity/roundtrip.R`), alongside the logic
  check `tests/parity/parity.R`, which already passes. Sample random cells and
  entries and compare legacy and new values for:
  - `prob` and `pct_pos`;
  - SQLite `lognorm` vs `mean_lognorm`;
  - v1 CSR vs `cells.h5ad` per-gene expression;
  - UMAP coordinates, including `cell_indices` subsets;
  - per-cell cluster labels.
- [ ] Converter unit tests on tiny synthetic inputs built inside the tests: a
  3-gene SQLite, and a v1 `db.h5` written with rhdf5 that has absent clusters,
  `-1` non-bimodal values and a subset UMAP.
- [ ] Document known migration facts: the optic-lobe data has no bimodality
  information (every gene is used, as before), and `ozel_2021` is deprecated but
  kept.

## 3. Build from Seurat without binarization

- [ ] Split `script/seurat_to_h5ad.R` into a **build** step that makes `--prob`
  optional. The build step writes `clusters.h5ad` with `pct_pos`,
  `mean_lognorm` and `n_cells`, plus `cells.h5ad`, and leaves the file in the
  *built* state.
- [ ] Support `--split-col` / `--cluster-col` (defaults `scmarco_split` /
  `scmarco_cluster`), with the `<split>@<cluster>` identities as fallback. Put
  input checks in `check_seurat_input()` with actionable errors, following the spec.
- [ ] New `script/add_probability_h5ad.R`: reads the canonical long table or the
  legacy wide RDS and validates it against the clusters file. It writes
  `prob`, `is_bimodal`, `bimodal_by_split` and `uns/scmarco/binarization`, then
  validates with `--stage binarized`.
- [ ] Tests, extending `tests/testthat/test-seurat_import.R`:
  1. build, then validate (`built`);
  2. check that the app refuses the dataset;
  3. add probabilities, from both the long table and the legacy RDS;
  4. validate (`binarized`);
  5. check that `ds_load()` works and `from_cluster()` gives the expected result.

  Add one bad-input test per rule in the spec.

## 4. Wrap-up

- [ ] Update `README.md` with the two-step build, the migration commands and a link
  to `docs/expected_input.md`.
- [ ] Re-run all checks:
  - `Rscript -e 'testthat::test_dir("tests/testthat")'`
  - `SCMARCO_CONFIG=configs/all.yaml Rscript tests/integration/app_flow.R`
  - `Rscript tests/parity/parity.R ../scMarco ../scMarco-Ding-lab`
  - `Rscript tests/parity/roundtrip.R ../scMarco ../scMarco-Ding-lab`
- [ ] Click through in a browser once for each lab config.
- [ ] Make the first git commit, once the user OKs it, and decide where the
  converted `data/*.h5ad` files live (they are git-ignored).

## Done when

- The unit tests pass, including the new converter and two-step build tests.
- `script/migrate_legacy.R` on `../scMarco` and `../scMarco-Ding-lab`
  reproduces the current `data/*.h5ad`, and both parity scripts report 0
  mismatches.
- An app started with a *built*-only dataset fails at startup with the
  listing message. After `add_probability_h5ad.R`, it starts and
  `tests/integration/app_flow.R` passes.

## Backlog (future)

- [ ] **smk pipeline helper** for
  [mixture_modeling_smk](https://github.com/chenyenchung/mixture_modeling_smk):
  export the pipeline's inputs from a *built* dataset, run or submit the
  pipeline, and import its output through `add_probability_h5ad.R`. The long
  table in `docs/expected_input.md` is canonical so this helper needs no new
  import format. The export format still has to be confirmed against the pipeline.
