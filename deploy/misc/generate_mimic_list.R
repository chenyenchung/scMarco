library(readxl)

# Read MIMIC list from
mimic <- read_xlsx("data/MI-list-2016-06-30.xlsx")
code_int <- mimic[grepl("coding intron", mimic$Position), ]

genes <- sapply(
  strsplit(code_int$`Gene(s) Affected`, " "),
  function(x) {
    out <- unlist(x)
    out <- setdiff(out, c("[+]", "[-]"))
    return(out)
  }
)

genes <- unlist(genes)

writeLines(
  genes,
  "int/mimic_temp.txt"
)

