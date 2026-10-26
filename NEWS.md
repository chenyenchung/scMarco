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
* New synthetic example dataset `toy` (inst/extdata) with a cells file,
  reduction, cluster names and metadata, and mixture-model fits.
* Fixed: a config without `gene_list` no longer fails to build the UI.
