## How to use this tool

### Shared options

#### Which genes to consider?

By default, every gene that is detected in the dataset is used for analysis and
visualization. You can provide your own gene list to focus on them. Please note
that gene symbols could change over time, and you **need to examine if your
provided gene list utilize the same naming scheme as your data**.

#### Clusters of interest

The clusters to focus on (find distinct genes) or to highlight (plot expression
trend).

#### Probability threshold to be considered as expressed

The tool utilizes binomial probability of a gene being expressed in a cluster/
cell type. A cutoff is selected here to decide the probability threshold to
consider a gene **ON**.

### Find distinct genes

Among your clusters of interest, find genes that are only **ON** in one cluster
but not others. Use the slide bar to define **ON** as _expressed in at least N
stages_. Alternatively, use the check boxes to select the stages a gene must be
expressed to be considered as **ON**.

Click **Find genes** to visualize genes that distinguishes between these
clusters.

### Co-expression plot

Select genes that you are interested in, and visualize their co-expression
pattern through different stages.

### Plot expression trend

While binarization is powerful and simplistic, it has its own limitations.
Therefore, one would often want to also see how the log-normalized expression
changes of a gene over different stages and in contrast to other cell types.

This tool generates an interactive line plot of log-transformed cluster average
normalized expression. Cluster label will be shown when you point your cursor to
a point. Clusters of interest will be highlighted (and coded with distinct
colors when you download a static version of the plot).
