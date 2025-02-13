## How to use this tool?

This tool accompanies [Chen et al. (2023)](https://doi.org/10.1101/2023.02.03.527019)
and aims to streamline the art and science of selecting marker combinations with
a Bayesian mixture model to binarize gene expression in scRNA-seq ([Davis et al. (2018)](https://elifesciences.org/articles/50901#s4);
we maintain [a fork](https://github.com/chenyenchung/mixture_modeling_smk) that patched
[a bug in elpd estimation](https://github.com/fredpdavis/opticlobe/pull/1)).

We recommend users to start **from a cluster (From cluster)** or **from a gene
that marks several clusters (From gene)** of interest. Either tab will help you
identify other clusters that could be marked when you only use one marker, and
help you select another marker to specifically label it (**Find distinct genes
within**).

Once marker combinations are found, you might want to see how other clusters
in the same dataset express this combination (**Co-exprssion plot**).

We acknowledge that while binarization greatly reduces computational complexity,
it has its limitation: Gene expression is a spectrum, and there would be
clusters that are hard to categorize to either ON or OFF state. We thus 
encourage users to also examine log-normalized expression value and its
dynamic (**Plot expression trend**) and decide whether a satisfactory pair of
markers have been identified.

### Shared options

#### Which dataset to use?

This tool retrieves information stored in a [SQLite3](https://sqlite.org/index.html)
database, which is preprocessed by Bayesian mixture modeling from single-cell
transcriptomics data.

Multiple datasets could co-exist in the database, and this switch allows users
to explore what intrigues them the most.

#### Which genes to consider?

By default, every gene that is detected in the dataset is used for analysis and
visualization. You can provide your own gene list to focus on them. Please note
that gene symbols used in reference genomes could change over time, and you 
**need to examine if your provided gene list utilize the same naming scheme as
your data**.

#### Clusters of interest

The clusters to focus on (Find distinct genes within) or to highlight (Plot 
expression trend).

#### Probability threshold to be considered as expressed

Binomial probability of a gene being expressed in a cluster/cell type . A cutoff is selected here to decide the probability threshold to
consider a gene **ON**.

### From cluster

Starting from one cluster, define how consistent among different stages/
conditions you desire for a marker. After clicking **Show genes**, the 
consistently expressing genes would be demonstrated and ranked by how many
other clusters they are also detected in.

You can select a gene, and move forward to **Find distinct genes within**.


### From gene

Similar to **From clustser**, but starting from a gene. You define how 
consistent among different stages/conditions you desire for a marker. After
clicking **Show clusters**, all clusters that consistently express this gene
would be shown, and you can choose the clusters that you want to distinguish
and move forward to **Find distinct genes within**.
other clusters they are also detected in.

### Find distinct genes within

Among your clusters of interest, find genes that are only **ON** in one cluster
but not others. Use the slide bar to define **ON** as _expressed in at least N
stages_. Alternatively, use the check boxes to select the stages a gene must be
expressed to be considered as **ON**.

Click **Find genes** to visualize how well you can distinguish these clusters
with the addition of another gene. You could pick a gene from this list to
use together with what you chose in **From cluster** or **From gene**.

You might want to visualize this gene pair to see how they express in other
clusters in the dataset with **Co-expression plot**.

### Co-expression plot

Select genes that you are interested in, and visualize their co-expression
pattern through different stages.

### Plot expression trend

This tool generates an interactive line plot of log-transformed cluster average
normalized expression. Cluster label will be shown when you point your cursor to
a point. Clusters of interest will be highlighted (and coded with distinct
colors when you download a static version of the plot).
