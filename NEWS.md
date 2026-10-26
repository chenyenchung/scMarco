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
* New synthetic example dataset `toy` (inst/extdata) with a cells file,
  reduction, cluster names and metadata, and mixture-model fits.
* Fixed: a config without `gene_list` no longer fails to build the UI.
