# scMarco 0.99.0

* First version as an R package, prepared for Bioconductor submission.
* The Shiny app runs with `run_scmarco()` or `scmarco_app()`; relative paths
  in the config are resolved against `root`. The usage guide now ships with
  the package and is appended to the configured README files
  (`usage_guide: false` turns it off).
* The combination searches, plots and dataset tools are exported and
  documented: `ds_load()`, `ds_binarize()`, `from_cluster()`,
  `pair_offtargets()`, `distinct_within()`, `onoff_plot()`, `umap_plot()`,
  `seurat_to_h5ad()`, `add_probability()`, `validate_h5ad()`, ...
* The mixture-model binarizer is part of the package: `fit_mixture()`,
  `merge_mixture_fits()` and `call_mixture()` port the mixture_modeling_smk
  pipeline and give the same values. They adapt Fred P. Davis's mixture
  model code (https://github.com/fredpdavis/opticlobe; Davis et al. 2020), released
  under the MIT license with his permission. Its fit parameters are kept unrounded,
  so the mixture-model plot draws the exact fitted curves.
* Fits, `call_mixture()` tables and datasets record the scMarco schema
  version (2.0) and the binarizer version (1.0), and each step checks them:
  fits of another version cannot be merged or called, `add_probability()`
  records the tables' binarizer, refuses tables from outside scMarco's
  binarizer unless `external = TRUE` (recorded as `"external"`) and refuses
  mixed ones, and the validator and the app
  require the binarizer version of a binarized dataset.
* `make_pseudobulk()` pools each cluster's raw counts into pseudobulk samples
  for `fit_mixture()` as `Seurat::AggregateExpression()` does (summed counts,
  then LogNormalize): random parts (sampling without replacement) or the
  values of a cell metadata column.
* `make_pseudobulk(method = "metadata")` has no minimum sample size by
  default: every library (or other group) in a cluster is a sample, and a
  message names the samples with fewer than 10 cells. Random parts keep a
  minimum of 10 cells. Give `min_cells` to set either.
* Fixed: `fit_mixture()` no longer loops forever when a replicate label has
  a single sample.
* The confidence margin is now the *ambiguity band* `band`, the half-width
  around P(On) = 0.5 (ON at P(On) >= 0.5 + band, OFF at <= 0.5 - band), in
  `binarize()`, `ds_binarize()`, `mixture_plot()`, `find_pairs()` and the
  app's slider. The app's default is 30% (ON at 0.8, OFF at 0.2, Fred P.
  Davis's thresholds) instead of a 10% margin. `margin`, and `cut` in
  `find_pairs()`, still work with a deprecation warning.
* Binarizer 1.1 calls (`call_mixture()`): genes are binarized at a false
  discovery rate of 10% (Benjamini-Hochberg on the elpd z, instead of two
  standard errors), with Ashman's D >= 2 and chains that agree; genes not
  binarized get no calls, and the dataset-wide ON/OFF distributions are gone.
  P(On) is the posterior mean. Flags mark calls whose chains disagree, ON
  calls below a detection floor and OFF calls above ambient RNA
  (`ambient_summary()`, DecontX per library: the cell suspension each 10x
  channel was loaded from; it warns about libraries with fewer than 10 cells
  per cluster, whose contamination DecontX tends to underestimate). Every cut is an argument and is
  printed. mclust is no longer needed. `add_probability()` stores the flags,
  the log2 fold over ambient RNA, the per-gene values and the parameters with
  the dataset, and `validate_h5ad()` checks them. `binarize()` (`flags`) and
  `ds_binarize()` make the flagged calls ambiguous on top of the band.
  `pair_offtargets()` and `find_pairs()` break ties among pairs with as many
  off-targets by *separation* (the pair's target strength minus its worst
  leak, in log2 fold over ambient RNA), then by the weaker gene's lower bound
  of Ashman's D; each is its own column (`separation`, `target_fold`,
  `worst_leak`, `ashman_D_lower`). In the app, a letter on a gray tile says which check
  made a call ambiguous, a worst leak of 1 or less shows as noise, the
  ranking columns appear only for datasets that have the data, and the
  sidebar says whether the dataset's call checks are applied.
* `qc_report()` reports on one split's binarization (QC report version 1): a
  YAML file with stable keys for comparing runs and a self-contained HTML
  page with the dataset, the settings in use (with what each is), the genes
  binarized and why the others were not, chain agreement, the detection and
  ambient checks, the calls, known limitations and provenance.
* `umap_plot()` labels only the clusters in `highlight` (the app's clusters
  of interest) unless `label_all = TRUE` (the app's "Label all clusters"
  checkbox); the legend, when shown, is the key to every cluster's color.
* `mixture_plot()` takes several models and stacks them top to bottom; in the
  app, checkboxes pick the stages, sexes or species to stack (all by default).
* New synthetic example dataset `toy` (inst/extdata) with a cells file,
  reduction, cluster names and metadata, and mixture-model fits.
* Fixed: a config without `gene_list` no longer fails to build the UI.
