library(readxl)
library(tidyr)

tempf <- tempfile()

# Download mixture modeling from Ozel el al. (2021)
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE142787&format=file&file=GSE142787%5FMixture%5Fmodeling%2Exlsx",
  destfile = tempf
)

# Each stage is saved as a sheet, let's load all stages as a list.
mm_matrices <- lapply(
  excel_sheets(tempf), function(sheet_name) {
    sheet <- read_xlsx(
      tempf, sheet = sheet_name
    )
    colnames(sheet)[1] <- "gene"
    return(sheet)
  }
)

# Name the elements with their corresponding stage
names(mm_matrices) <- vapply(
  strsplit(excel_sheets(tempf), "_"), function(x) {
    x[1]
  },
  FUN.VALUE = character(1)
)

# To make the example slim, we are going to keep only 10 clusters for
# demonstration purpose.

to_keep <- c("1", "8", "9", "12", "14", "15", "19", "27", "31")

mm_matrices <- lapply(mm_matrices, function(sheet) sheet[ , c("gene", to_keep)])

# We then process log-normalized values similarly

temp_avg <- tempfile()

download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE142787&format=file&file=GSE142787%5FLog%5Fnormalized%5Faverage%5Fexpression%2Exlsx",
  destfile = temp_avg
)

# Each stage is saved as a sheet, let's load all stages as a list.
avg_matrices <- lapply(
  excel_sheets(temp_avg), function(sheet_name) {
    sheet <- read_xlsx(
      temp_avg, sheet = sheet_name
    )
    colnames(sheet)[1] <- "gene"
    return(sheet)
  }
)

# Name the elements with their corresponding stage
names(avg_matrices) <- vapply(
  strsplit(excel_sheets(temp_avg), "_"), function(x) {
    x[1]
  },
  FUN.VALUE = character(1)
)

avg_matrices <- lapply(avg_matrices, function(sheet) sheet[ , c("gene", to_keep)])

# Save for demonstration purpose
mm_path <- "example/data/mixture_model/"
dir.create(mm_path, recursive = TRUE)

for (sheet_name in names(mm_matrices)) {
  write.csv(
    mm_matrices[[sheet_name]],
    paste0(mm_path, sheet_name, ".csv"),
    row.names = FALSE
  )
}

avg_path <- "example/data/log_norm/"
dir.create(avg_path, recursive = TRUE)

for (sheet_name in names(avg_matrices)) {
  write.csv(
    mm_matrices[[sheet_name]],
    paste0(avg_path, sheet_name, ".csv"),
    row.names = FALSE
  )
}
