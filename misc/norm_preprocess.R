# This script loads and preprocesses log-normalized expression
library(readxl)
library(reshape2)

file_path <- "data/log_normalized_average_expression_all_stages.xlsx"
annot_path <- "data/20210504_annotation_clusters.xlsx"

source("src/utils.R")

# Get all tabs
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
      value.name = "lognorm"
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

# Make sure cluster is character
explong$cluster <- as.character(explong$cluster)

# Order the stages
explong$stage <- factor(
  explong$stage,
  levels = c("P15", "P30", "P40", "P50", "P70", "Adult")
)

# # Scale the expression
# ## Split expression per cluster
# explist <- split(explong, explong$cluster)
#
# # Remove redundant variable to conserve RAM
# rm(expmtx)
# rm(explong)
#
# # Convert to wide form to scale and convert back
# # This step is pretty slow.
# explong <- lapply(
#   names(explist), function(matname) {
#     # Print current cluster name to learn how far the process has gone.
#     print(matname)
#     wmat <- dcast(
#       data = explist[[matname]], stage ~ gene, value.var = "lognorm"
#     )
#
#     # Fill NAs with 0
#     wmat[is.na(wmat)] <- 0
#
#     # Scale expression per gene
#     # The first column is stage
#     wmat[ , -1] <- scale(wmat[ , -1])
#
#     ## If all expression are 0, then scaling will return NaN! ##
#     lmat <- melt(
#       wmat, id.vars = "stage",
#       value.name = "scaled_exp", variable.name = "gene"
#     )
#
#     # Add cluster information back
#     lmat$cluster <- matname
#
#     # Convert gene names back to chr
#     lmat$gene <- as.character(lmat$gene)
#
#     return(lmat)
#   }
# )
#
# explong <- do.call("rbind", explong)

# Annotate the clusters
explong <- AnnotateCluster(explong, annot_path)

if (!dir.exists("int")) {dir.create("int")}
saveRDS(explong, "int/lognorm.rds")
