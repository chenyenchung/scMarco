# Load the scMarco dataset tools (seurat_to_h5ad(), add_probability(),
# set_cluster_names(), sqlite_to_h5ad(), find_pairs(), validate_h5ad(), ...)
# into the global environment:
#   source("script/load.R", chdir = TRUE)
for (f in c("../src/combo_core.R", "../src/backend_h5ad.R", "lib_write_h5ad.R",
            "lib_import.R", "lib_validate.R", "lib_build.R", "lib_migrate.R",
            "lib_pairs.R")) {
  source(f)
}
rm(f)
