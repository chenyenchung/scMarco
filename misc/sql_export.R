library(DBI)
library(RSQLite)

prob_tbl <- read.csv("int/probs.csv")
log_norm_tbl <- read.csv("int/lognorm.csv")
colnames(prob_tbl)[4] <- "value"
colnames(log_norm_tbl)[4] <- "value"
prob_tbl$type <- "prob"
log_norm_tbl$type <- "lognorm"

prob_db <- dbConnect(RSQLite::SQLite(), "int/db.sqlite")

dbWriteTable(prob_db, "ors", prob_tbl)
dbAppendTable(prob_db, "ors", log_norm_tbl)

dbSendStatement(prob_db, "CREATE UNIQUE INDEX ors_idx ON ors(gene, stage, cluster, type)")

prob_tbl <- read.csv("../ozel_et_al_2021/int/mixture_modeling.csv")
log_norm_tbl <- read.csv("../ozel_et_al_2021/int/lognorm.csv")

colnames(prob_tbl)[4] <- "value"
colnames(log_norm_tbl)[4] <- "value"
prob_tbl$type <- "prob"
log_norm_tbl$type <- "lognorm"

dbWriteTable(prob_db, "ozel_2021", prob_tbl)
dbAppendTable(prob_db, "ozel_2021", log_norm_tbl)


dbSendStatement(prob_db, "CREATE UNIQUE INDEX ozel_idx ON ozel_2021(gene, stage, cluster, type)")

dbDisconnect(prob_db)
