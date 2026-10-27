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
* `make_pseudobulk()` pools each cluster's raw counts into pseudobulk samples
  for `fit_mixture()` as `Seurat::AggregateExpression()` does (summed counts,
  then LogNormalize): random parts (sampling without replacement), the values
  of a cell metadata column, or k-means subclusters.
* Cluster heterogeneity: `mixture_heterogeneity()` scores how much the calls
  of a cluster's k-means subclusters disagree for each gene, against genes
  expressed at a similar level (Winsorized reference per level bin, pooled
  across the split's clusters), and flags those whose subclusters are
  confidently ON and OFF. `add_probability(heterogeneity = )` stores them;
  genes heterogeneous in a cluster are ambiguous there in every search.
  `cluster_heterogeneity()` lists them, and the From cluster search warns
  when the chosen cluster may be heterogeneous.
* Fixed: `fit_mixture()` no longer loops forever when a replicate label has
  a single sample.
* New synthetic example dataset `toy` (inst/extdata) with a cells file,
  reduction, cluster names and metadata, and mixture-model fits.
* Fixed: a config without `gene_list` no longer fails to build the UI.
