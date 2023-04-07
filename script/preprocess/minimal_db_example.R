library(tidyr)
library(DBI)
library(RSQLite)

db <- dbConnect(SQLite(), "../example/min_db.sqlite")

# Average expression
lognorm <- list.files("../example/data/log_norm/")

lognorm_t <- lapply(
  lognorm, function(i) {
    # Get stage
    stage <- sub("\\.csv$", "", i)

    # Read table
    x <- read.csv(paste0("../example/data/log_norm/",i), check.names = FALSE)

    # Store gene symbol
    x$gene <- row.names(x)

    # Pivot
    x <- pivot_longer(
      x, cols = -gene, names_to = "cluster", values_to = "value"
    )

    # Annotate
    x$stage <- stage
    x$type <- "lognorm"
    return(x)
  }
)
lognorm_t <- do.call(rbind, lognorm_t)
dbWriteTable(
  conn = db,
  name = "example",
  value = lognorm_t
)

# Mixture modeling probability
mm <- list.files("../example/data/mixture_model/")

mm_t <- lapply(
  mm, function(i) {
    # Get stage
    stage <- sub("\\.csv$", "", i)

    # Read table
    x <- read.csv(paste0("../example/data/mixture_model/",i), check.names = FALSE)

    # Store gene symbol
    x$gene <- row.names(x)

    # Pivot
    x <- pivot_longer(
      x, cols = -gene, names_to = "cluster", values_to = "value"
    )

    # Annotate
    x$stage <- stage
    x$type <- "prob"
    return(x)
  }
)

mm_t <- do.call(rbind, mm_t)

dbAppendTable(
  conn = db,
  name = "example",
  value = mm_t
)

dbDisconnect(db)
