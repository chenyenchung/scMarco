# Generates inst/extdata/example.clusters.h5ad, the clusters file of a real
# dataset trimmed to 300 genes: the "example" table of the legacy scMarco
# SQLite database min_db.sqlite (fly optic lobe development, six stages
# from P15 to adult). The database is kept in the source repository under
# deploy/example/, not in the package.
#
#   Rscript inst/scripts/make_example_data.R   # from the package root
#
# It has no cells file and no bimodality information (every gene is used),
# and its clusters are bare IDs without names.

library(scMarco)
sqlite_to_h5ad(
    db = file.path("deploy", "example", "min_db.sqlite"),
    table = "example",
    stages = c("P15", "P30", "P40", "P50", "P70", "Adult"),
    split_label = "stage", ordered = TRUE, n_genes = 300,
    out = file.path("inst", "extdata", "example.clusters.h5ad")
)
