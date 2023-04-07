# This script preprocesses the mixture modeling data to make it easy
# for Shiny app to access
library(readxl)

# Define input data
file_path <- "data/Mixture_modelling_all_stages.xlsx"
annot_path <- "data/20210504_annotation_clusters.xlsx"

# Get all tabs from the Excel sheet
sheet_list <- excel_sheets(file_path)

# Read tables and annotate with stage
expmtx <- lapply(
  sheet_list,
  function(sheet_name) {
    mat <- read_xlsx(
      file_path,
      sheet = sheet_name
    )
    colnames(mat)[[1]] <- "gene"

    # Strip the meta columns (last 11)
    mat <- mat[ , c(1:(ncol(mat) - 11))]

    # Annotate stage
    mat$stage <- strsplit(sheet_name, "_")[[1]][1]
    return(mat)
  }
)

# Format the data to long-form
explong <- lapply(
  expmtx,
  function(x) {
    longmat <- melt(
      x,
      id.vars = c("gene", "stage"),
      variable.name = "cluster",
      value.name = "prob"
    )
    return(longmat)
  }
)

# Rbind long-form tables
explong <- do.call("rbind", explong)

# Fix Excel name-to-date error
explong$gene[explong$gene == "43709"] <- "Sep1"
explong$gene[explong$gene == "43710"] <- "Sep2"
explong$gene[explong$gene == "43712"] <- "Sep4"
explong$gene[explong$gene == "43713"] <- "Sep5"
explong$gene[explong$gene == "43710"] <- "Dec1"

# Make sure cluster is cluster
explong$cluster <- as.character(explong$cluster)

# Order the stages
explong$stage <- factor(
  explong$stage,
  levels = c("P15", "P30", "P40", "P50", "P70", "Adult")
)

# Annotate the clusters
explong <- AnnotateCluster(explong, annot_path)

if (!dir.exists("int")) {dir.create("int")}
saveRDS(explong, "int/mixture_modeling.rds")
