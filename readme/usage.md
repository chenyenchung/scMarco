## How to use this tool

This tool streamlines selecting marker gene combinations with a Bayesian
mixture model that binarizes gene expression in scRNA-seq
([Davis et al. (2018)](https://elifesciences.org/articles/50901#s4); we maintain
[a fork](https://github.com/chenyenchung/mixture_modeling_smk) that patches
[a bug in elpd estimation](https://github.com/fredpdavis/opticlobe/pull/1)).
It was originally designed for [Chen et al. (2023)](https://doi.org/10.1101/2023.02.03.527019).

We recommend starting **from a cluster (From cluster)** or **from a gene that
marks several clusters (From gene)**. Either tab shows which other clusters a
single marker would also label, and hands those clusters to **Find distinct
genes within** so you can pick a second marker that separates them.

Once you have a combination, **Co-expression plot** shows how every cluster in
the dataset expresses it.

Binarization greatly reduces complexity, but expression is a spectrum and some
clusters are hard to call ON or OFF. Where the dataset provides them, also check
the log-normalized values: **Plot expression trend** for datasets with ordered
stages, and **UMAP** for datasets with cell-level data.

### Shared options

#### Which dataset to use?

See the dataset descriptions above.

#### Which genes to consider?

By default every gene the mixture model called bimodal is used for marker
finding (every gene, when the dataset has no bimodality information). You can
restrict the analysis to a gene list. Gene symbols change between genome
releases, so **check that your list uses the same naming scheme as the data**.

#### Cluster(s) of interest

The cluster to start from (**From cluster**, exactly one), the clusters to tell
apart (**Find distinct genes within**, two or more), or the clusters to
highlight (**Plot expression trend**). Annotated clusters are shown as
`name (id)`, e.g. `Mi4 (154)`, and unannotated ones by their ID; you can type
either to search.

#### Probability threshold

A gene is **ON** in a cluster when its probability of being expressed is
greater than this threshold.

#### Must be ON in these … / ON in at least this many …

These two controls apply to every tab. A gene counts as a marker of a cluster
when it is ON in **every checked** stage, sex or species, **and** ON in at least
the number chosen on the slider. With nothing checked and the slider at 1, a
gene expressed in any one of them qualifies. Off-target clusters are always
counted if the gene is ON in *any* of them.

### From cluster

Pick one cluster and click **Show genes**. Candidate markers are ranked by how
many other clusters also express them. Tick **Also rank gene pairs** to see the
pairs of candidates that share the fewest off-target clusters. Select a gene to
send its off-target clusters to **Find distinct genes within**.

### From gene

Pick one or more genes and click **Show clusters**. With several genes, choose
whether a cluster must express **all** of them (a combination) or **any** of
them. Select the clusters to distinguish and move on to **Find distinct genes
within**.

### Find distinct genes within

Among the clusters of interest, find genes that are ON in exactly one of them
(or the number you choose) and OFF in the others. The ON/OFF plot shows each gene
per cluster: black is ON, white is OFF, gray is a gene that is not bimodal there,
and a cross marks a cluster that does not exist in that stage, sex or species.
When you arrived here from **From cluster** or **From gene**, a table checks
each combination against the **whole dataset**, not just the selected clusters.

### Co-expression plot

Select genes and see, for every cluster that expresses any of them, where each
gene and the full combination are ON. With cell-level data, the percentage on
the **Co-expressed** panel is the share of cells that express **all** the genes.

### Plot expression trend

An interactive line plot of cluster-average log-normalized expression across
ordered stages. Hover over a line to see its cluster. Clusters of interest are
highlighted, and get distinct colors in the downloaded static version.

### UMAP

Plot any stored reduction colored by a gene's log-normalized expression, by
cluster, or by stage, sex or species. Cells are limited to the checked items of
**Must be ON in these …**; with none checked, every cell is shown.
