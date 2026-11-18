# Expected input

This spec covers building a new scMarco dataset. There are two steps, both R
functions of the scMarco package (`library(scMarco)`):

1. **Build** (`seurat_to_h5ad()`): a Seurat object →
   `<id>.clusters.h5ad` (no probabilities yet) + `<id>.cells.h5ad`.
2. **Binarize** (`add_probability()`): mixture-model output →
   `layers/prob` and the bimodality flags in `<id>.clusters.h5ad`.

```r
seurat_to_h5ad("obj.rds", "dmel_mixed", split_label = "sex")
# tables from call_mixture(), saved with saveRDS()
add_probability("data/dmel_mixed.clusters.h5ad",
                c(Male = "male_prob.rds", Female = "female_prob.rds"))
```

Both take paths or objects already in memory.

The app **only accepts binarized datasets**. If any configured dataset has no
`layers/prob`, the app stops at startup, lists the offending datasets and names
the command that fixes them.

Legacy data (scMarco SQLite, Ding-lab HDF5 v1.0) is migrated with its own
scripts and produces binarized files directly; see the README.

---

## 1. Seurat object (build step)

A Seurat v4 or v5 object, or an `.rds` file readable by `readRDS()` holding
one. Only the SeuratObject package is needed; full Seurat isn't required.
Seurat v3 objects saved before SeuratObject existed (e.g. the Ozel et al. 2021
submission objects) are updated for SeuratObject on read.

When each split was saved as its own object (e.g. one per stage), pass them
all, named by split: `seurat = c(P15 = "P15.rds", P30 = "P30.rds", ...)` (or a
named list of objects). The split is the name, and `split_col` is not used. The genes are the union of the objects' genes; a
gene missing from an object counts as 0 in its cells. Splits keep the order
given unless `split_order` says otherwise. An object's `umap` reduction is
stored as `obsm/X_<split>`.

### Expression

| Requirement | Details |
|---|---|
| Assay | The default assay (usually `RNA`); another with `assay`, e.g. `assay = "RNA"` when the default is `integrated`. |
| Counts | **Raw counts**: non-negative integers, genes × cells. In v5, the `counts` layer (layers split by sample, `counts.1`, `counts.2`, ..., are joined by the builder); in v4, the `counts` slot. |
| Normalization | Done by the builder, as Seurat's LogNormalize: `log1p(counts / library size × 1e4)`. Any `data` layer already present is ignored. |
| Genes | Unique, non-empty names; the builder sorts them. They must use the same symbols as the mixture-model output. |
| Cells | Unique barcodes (Seurat enforces this). |

### Split and cluster labels

Each cell needs a **split** (stage, sex, species, …) and a **cluster ID**.
Use the ID the clustering gave (e.g. `154`), not an annotation: the ID is the
key every file and table joins on, so it must not change when the cluster is
named.

| Source | How |
|---|---|
| Preferred | Metadata columns `scmarco_split` and `scmarco_cluster`; other names with `split_col` / `cluster_col`. `split_col = NULL` puts every cell in one split named `all`. A column named with either argument must exist. |
| Fallback (legacy) | `Idents(obj)` of the form `<split>@<cluster>`, exactly one `@`. Used only when neither default column exists and neither option is given. |

`exclude_clusters = "0"` drops the cells of the listed cluster IDs before
anything else, e.g. cells the study discarded but kept in the object.

Name rules, checked with actionable errors (`check_labels()` in `R/import.R`;
`validate_h5ad()` in `R/validate.R` checks them again):

- no empty or missing value in either label;
- no `|`, because the files use `"<split>|<cluster>"` and `"<split>|<barcode>"`
  as observation names;
- no leading or trailing whitespace;
- spaces, parentheses, `-`, `_` and `.` are allowed (e.g. `151`, `ANm_12`);
- a split × cluster pair with no cells simply does not exist. It is shown as
  "absent", not OFF.

### Cluster names

Names are added after the build, and changed at any time, from a
tab-separated table:

```
cluster_id	name
151	T4b
154	Mi4
```

```r
set_cluster_names("data/<id>.clusters.h5ad", "names.tsv")  # or a data.frame, or c(`151` = "T4b")
get_cluster_names("data/<id>.clusters.h5ad", "names.tsv")  # current names (and the table)
```

Listed clusters get the given name; an empty name clears it; unlisted clusters
keep theirs. Unknown or repeated IDs are errors. The app shows `T4b (151)`
for a named cluster and `151` for an unnamed one.

To change the IDs themselves (e.g. a dataset built with annotation names as
IDs), give every cluster a new ID, and optionally a name, in one step. Pass the
cells file too so its IDs stay in sync; the data are not touched:

```
cluster_id	new_id	name
Dm4	4	Dm4
G_PsCrtrg	1+2	PCG
```

```r
set_cluster_ids("data/<id>.clusters.h5ad", "ids.tsv", cells = "data/<id>.cells.h5ad")
```

### Cluster metadata (optional)

A table of groupings lets app users limit the searches to part of the dataset,
e.g. "among Tm neurons only". Its first column is `cluster_id`; every other
column is one level (e.g. class, then subclass):

```
cluster_id	class	subclass
151	neuron	T4/T5
154	neuron	Mi
3	glia	astrocyte
```

```r
set_cluster_metadata("data/<id>.clusters.h5ad", "meta.tsv")   # or a data.frame
get_cluster_metadata("data/<id>.clusters.h5ad", "meta.tsv")   # current table (and the TSV)
set_cluster_metadata("data/<id>.clusters.h5ad", NULL)         # remove it
```

The table replaces any earlier one. Unlisted clusters and empty values are
unassigned, and an unassigned cluster is left out whenever its level is
filtered. Unknown or repeated IDs are errors, and column names must be unique
and must not contain `/`. Each level is stored as `obs/meta_<level>`; the
metadata is kept when `set_cluster_ids()` renames the clusters.

The app only offers the levels listed under the dataset's `subset_by` in the
config, in that order (coarse to fine); without it, there is no subsetting:

```yaml
datasets:
  my_data:
    clusters: data/<id>.clusters.h5ad
    subset_by: [class, subclass]
```

In the app, users pick which of these levels to subset by ("Subset clusters
by", none by default = all clusters); each picked level gets a "Search within"
filter. Choosing values at a level narrows the choices of the picked levels
below it. The cluster picker, From cluster, From gene, Find distinct and the
co-expression plot then only see the matching clusters. Expression trend, UMAP
and Mixture model still show every cluster.

### Split order

| Argument | Meaning |
|---|---|
| `split_order = c("Mel", "Yak", "San", "Tei")` | Display order; must list every split exactly once. Defaults to alphabetical. |
| `ordered = TRUE` | The splits are a sequence (stages). This enables the expression-trend plot. |
| `split_label = "stage"` | Word used in the UI ("Must be ON in these stages"). |

### Reductions (optional)

| Source | Stored as |
|---|---|
| A reduction named `umap` on the object (all cells) | `obsm/X_umap_full`, shown as "All cells" |
| `umaps = "<folder>"`: one `.rds` Seurat object per subset, each with a `umap` reduction | `obsm/X_<name>`, where `<name>` is the file name before the first `_` (e.g. `ANm_umap.rds` → `X_ANm`) |

`reduction = "tsne"` stores that reduction instead of `umap`, as
`obsm/X_tsne_full` (one object) or `obsm/X_<split>_tsne` (one object per
split); the plots label its axes tSNE. The build stops if an object lacks a
reduction named with `reduction`, and lists the reductions it has (e.g.
`umap.rna` or `wnn.umap`). When the default `umap` is missing, the build
warns with the same list and goes on without a reduction. Pass the right
name as `reduction`, or `reduction = NULL` to build without one.

Reductions named after a split (`X_<split>` or `X_<split>_<reduction>`) are
treated as computed separately for each split. The UMAP tab shows each group as
one entry, "Each <split label> separately". It replaces "Split by" with a
multi-select of the splits to show, and draws each split as its own panel with
its own axes.

Subset cells are matched to the main object **by barcode**. Cells outside a
subset get NaN coordinates, and unmatched barcodes are ignored. All reductions
should come from the same embedding when they are shown side by side, e.g. one
UMAP split by sex.

---

## 2. Mixture-model probabilities (binarize step)

One probability per gene × split × cluster: the posterior that the gene is
in the "ON" component of the mixture model, plus whether each gene is bimodal.

The binarizer ships with the package (Suggests: rstan, matrixStats).
Per split, `fit_mixture()` fits each gene's log-normalized expression across
pseudobulk samples (each cluster's cells pooled into a few replicates) and
`call_mixture()` turns the fits into the table described in 2a. Fit chunks of
genes in parallel jobs and join them with `merge_mixture_fits()`. Fits are
seeded per gene (`seed = 42`), so they do not depend on the chunking:

```r
fit <- fit_mixture(lognorm, cluster = meta$cluster, replicate = meta$replicate)
amb <- ambient_summary(counts, cluster = meta_cells$cluster, library = meta_cells$library)
add_probability("data/<id>.clusters.h5ad", list(Male = call_mixture(fit, amb), ...))
```

`call_mixture()` binarizes a gene when its bimodal fit beats the unimodal one
at a false discovery rate of 10% (Benjamini-Hochberg on the elpd z), its
Ashman's D is at least 2 and its chains agree on at least 90% of clusters;
other genes get no calls. Each cluster's P(On) is the posterior mean, and
flags mark calls to treat as ambiguous: chains that disagree, ON with fewer
than 5% of the cluster's cells detecting the gene, OFF with UMIs above 8x the
ambient RNA expected (DecontX contamination per library, from
`ambient_summary()` on the split's raw counts and cell labels). Every cut is
an argument, and the values in use are printed.

`make_pseudobulk()` makes the input from one split's raw counts, as
`Seurat::AggregateExpression(return.seurat = TRUE)` did in the pipeline: each
cluster's cells are pooled into a few samples, their counts summed, and each
sample log-normalized, `log1p(summed counts / sample total * 1e4)`. The cells
file holds log-normalized values only, so take the counts from the Seurat
object. By default (`method = "random"`) the cells are dealt into random parts
without replacement, which assumes each cluster is homogeneous and gives the
model several samples per cluster when there are few clusters. With
`method = "metadata"` each cluster is split by a cell metadata column
(`group`, e.g. `orig.ident`): each value is one sample and its replicate (to
pool libraries into two replicates as the pipeline did, map them to two
groups first).

```r
obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])     # v5 layers split by sample
counts <- SeuratObject::LayerData(obj, layer = "counts")   # one split's cells
pb <- make_pseudobulk(counts, obj$scmarco_cluster, n_samples = 2)
fit <- fit_mixture(pb$lognorm, pb$cluster, pb$replicate)
```

It adapts Fred P. Davis's mixture model and code
([opticlobe](https://github.com/fredpdavis/opticlobe), Davis et al. 2020) through
the mixture-model pipeline
([mixture_modeling_smk](https://github.com/chenyenchung/mixture_modeling_smk)),
and tables from that pipeline are read the same way.

### 2a. One table per split (preferred)

The mixture model runs separately for each split, and its output for one split
is used as is: a `data.frame` (`.rds`) with

- genes as row names;
- one probability column per **cluster ID** of that split (e.g. `151`,
  `ANm_12`): a number in [0, 1], or NA;
- the QC columns `diff.elpd`, `diff.elpd.se`, `elpd.z`, `elpd.q`,
  `ashman.D`, `ashman.D.lower`, `spread.share`, `bimodal.on_level`,
  `bimodal.off_level`, `bimodal.sd`, `bimodal.pi`, `unimodal.level`,
  `unimodal.sd`, `exprType` and `not.binarized`. Only `exprType` is required. A gene is
  bimodal **in this split** when its `exprType` is `"bimodal"`.

The model's raw output, a character matrix (every column is text because of
`exprType`), is accepted too: its numbers are parsed, and any other text in a
probability column is an error. In a `data.frame`, probability columns must be
numeric.

Every column that is not a QC column is a cluster. Since the split comes from
which table is read, the column names need no separator, and a cluster ID may
contain any character the build step allows.

Pass the tables in either of two ways:

```r
# one file per split
add_probability("data/<id>.clusters.h5ad",
                c(Male = "male_prob.rds", Female = "female_prob.rds"))

# one .rds holding a named list: list(Male = male_df, Female = female_df)
add_probability("data/<id>.clusters.h5ad", "prob_by_split.rds")

# tables already in memory
add_probability("data/<id>.clusters.h5ad", list(Male = male_df, Female = female_df))

# tables from outside call_mixture(), e.g. the mixture_modeling_smk pipeline
add_probability("data/<id>.clusters.h5ad",
                c(Male = "male_smk.rds", Female = "female_smk.rds"), external = TRUE)
```

Tables record the binarizer that made them (see *Versions* below). A table
without one, i.e. not made by `call_mixture()`, is refused unless
`external = TRUE` is given; every split is then recorded as `"external"`.

The split names must match the build step. Per-split bimodality is stored in
`varm/bimodal_by_split`. A gene counts as bimodal in `var/is_bimodal` when it
is bimodal in at least one split.

When every table has `bimodal.on_level`, `bimodal.off_level`, `bimodal.sd` and
`bimodal.pi`, they are stored per gene × split in `varm/mixture_on_level`,
`varm/mixture_off_level`, `varm/mixture_sd` and `varm/mixture_pi`, as given:
`call_mixture()` writes them unrounded, so the curves are exact; the older
pipeline rounded them to 2 decimals. The app's
Mixture model tab draws them: the distribution of cluster averages with the two
fitted Gaussians, and the P(On) curve with the clusters of interest marked. If
any table lacks them, a warning is given and the tab is hidden.

The `prob` argument also tells which splits were **modeled together**, recorded
in `uns/scmarco/mixture_models` (one model name per split, in `bimodal_splits`
order). Each named item (`Male = "male.rds"`) and each element of a named list held
in an unnamed `.rds` is its own model; a wide table (2b) is one model across all its splits. Items can be mixed,
e.g. `prob = c("mel_yak_wide.rds", San = "san.rds")` (Mel and Yak modeled together, San
alone). For a joint model the Mixture model tab pools the clusters of all its
splits, counting each split × cluster pair once. Datasets
binarized before this was added get the tab once `add_probability()` is
run again.

Tables from `call_mixture()` (binarizer 1.1) also carry, as attributes, the
flags that make a call ambiguous, each cluster's log2 fold over ambient RNA
and the parameters in use; `add_probability()` stores them:

| Where | What |
|---|---|
| `layers/call_flags` | integer per split × cluster × gene, the sum of 1 (chains disagree), 2 (below the detection floor: never ON) and 4 (above ambient RNA: never OFF); 0 where a table has no value |
| `layers/fold_ambient` | log2((UMIs + 1) / (expected ambient UMIs + 1)), for ranking pairs; written when at least one split was called with an ambient summary, NaN elsewhere |
| `varm/call_elpd_z`, `call_elpd_q`, `call_ashman_D`, `call_ashman_D_lower`, `call_spread_share` | genes × splits (`bimodal_splits` order): elpd z, its Benjamini-Hochberg q, Ashman's D (posterior median and 10th percentile), share of clusters whose chains disagree |
| `uns/scmarco/binarization/calls` | the parameters of `call_mixture()` (`summary`, `bh_q`, `min_D`, `spread_cutoff`, `reject_spread_frac`, `on_floor`, `ambient_log2fc`, `ambient_alpha`, `ambient` as 0/1, `rho_source`), one value per split in `bimodal_splits` order; contamination per library as `rho_split`, `rho_library`, `rho` |

They are stored only when every split's table has them (tables from outside
scMarco, or subset tables, have none), and replaced or removed with the
probabilities when `add_probability()` runs again. The validator checks their
shapes and values, and warns when a dataset of the current binarizer lacks
`layers/call_flags`.

### 2b. Legacy wide RDS (accepted)

This is the format `add_probability_to_h5.R` read, and the output of one model
fit across several splits (e.g. four species modeled together):

- a `data.frame` with genes as row names;
- one probability column per `"<split>@<cluster>"`, all splits in one table;
- the QC columns `diff.elpd`, `diff.elpd.se`, `bimodal.on_level`,
  `bimodal.off_level`, `bimodal.sd`, `bimodal.pi`, `unimodal.level`,
  `unimodal.sd` and `exprType`. Genes whose `exprType` is not `"bimodal"` are
  non-bimodal in every split.

### Rules, for both formats

| Rule | On violation |
|---|---|
| Every split given exists in `clusters.h5ad` | error |
| Every split in `clusters.h5ad` is given | error naming the missing split(s) |
| Every cluster column of a split has cells in that split in `clusters.h5ad` | error, listing the first few. A hyphen in a column name is matched to an underscore when that is the only difference (the legacy pipeline wrote hyphens). |
| No repeated cluster column within a split | error |
| Every gene exists in `clusters.h5ad` | error |
| `exprType` column present | error |
| `prob` in [0, 1] or missing | error |
| Genes in the dataset but missing from a split's table | stored as NA and non-bimodal in that split, with a warning giving the count |
| Clusters with cells in a split but missing from its table | stored as NA (never ON), with a warning |
| Genes that are non-bimodal in a split | their `prob` in that split is stored as NA, so they are never ON there |

---

## 3. Output files and states

The output is described fully in the README ("Data format"). In short:

| File | Contents |
|---|---|
| `<id>.clusters.h5ad` | obs = split × cluster pairs that exist (`split`, `cluster` ID, optional `cluster_name`, `n_cells`); layers `pct_pos`, `mean_lognorm` and, once binarized, `prob`; `var/is_bimodal`, `varm/bimodal_by_split`, `varm/mixture_*` and, from `call_mixture()` tables, `layers/call_flags`, `layers/fold_ambient` and `varm/call_*` once binarized; `uns/scmarco` metadata |
| `<id>.cells.h5ad` | X = log-normalized cells × genes (CSC); obs = `barcode`, `split`, `cluster`; `obsm/X_*` reductions |

| State | Has `layers/prob` | `uns/scmarco/binarization` | Validator | App |
|---|---|---|---|---|
| *built* | no | absent | `validate_h5ad(stage = "built")` passes | **refuses to start** |
| *binarized* | yes | `method`, `source` (the `prob` files, `<table>` for one passed in memory), `date`, `version` (binarizer; required) | `validate_h5ad()` passes | loads |

### Versions

Every object passed from fitting to the app records two versions:

| Version | Where | Meaning |
|---|---|---|
| scMarco schema (`2.0`) | `uns/scmarco/format_version` of both h5ad files; `schema_version` of a `fit_mixture()` object | the layout of the object. v1 was the legacy SQLite database (convert it with `sqlite_to_h5ad()`); v2 is still being consolidated, so fields it requires are added without a bump until it is released |
| binarizer (`1.1`) | `binarizer_version` of a `fit_mixture()` object; attribute `scmarco_binarizer` of a `call_mixture()` table; `uns/scmarco/binarization/version` | the binarizer (`fit_mixture()` + `call_mixture()`) that made the probabilities; bumped when a change alters them. `"external"`: tables from outside scMarco's binarizer (mixture_modeling_smk pipeline, hand-made, migrated legacy databases) |

| Step | Checks |
|---|---|
| `merge_mixture_fits()`, `call_mixture()` | refuse fits of another schema or binarizer version, or saved before fits recorded them (refit); merged chunks must share `n_iter`, `n_folds` and `seed` |
| `add_probability()`, `seurat_to_h5ad(prob = )` | record the tables' binarizer, not the package's; tables recording none are refused unless `external = TRUE`, which records every split as `"external"` (with a warning if some tables did record one); otherwise every split's table must come from the same binarizer; a newer one is refused, an older one warns |
| writer | refuses a `prob` layer without a binarization record |
| `validate_h5ad()` | binarized files need `binarization` with `method`, `source`, `date` and a valid `version`; newer binarizer is an error, older a warning; the cells file's `format_version` must match |
| app (`ds_check_config()`, `ds_load()`) | refuses another `format_version` (either file), a missing, invalid or newer binarizer version; warns on an older one; shows the binarizer in the sidebar |

A `call_mixture()` table keeps its attributes (binarizer, flags, folds,
parameters) through `saveRDS()` and through selecting genes (rows), but not
through selecting columns; a table without them counts as `"external"`. Files
binarized before the version was required can be given the record (as
`"external"`) with `deploy/misc/stamp_binarizer_version.R`.

Migrated legacy datasets are binarized from the start. Their `uns/scmarco/source`
records where they came from, e.g. `sqlite:db.sqlite:ozel_2021_new`.
