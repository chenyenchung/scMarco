# Cluster metadata tables for set_cluster_metadata(): the cell class and
# developmental origin of optic-lobe clusters, as levels Cell class > Origin >
# Temporal origin > Spatial origin > Notch status. Origin is "OPC" when all
# three are known, "Non-OPC" otherwise (including clusters missing from the
# source table).
#
# Cell class is Glia, Others (non-neuronal, non-glial clusters, and the P15
# developmental groups NB, GMC*, Im*, LPC, Apop) or Neuron (everything else,
# including unannotated clusters):
# - PNAS table: Glia when the comment mentions glia, Others when it says
#   non-neuronal or contaminating.
# - ozel_2021_simon_2026: from the cluster names (set from MOESM10): glia
#   names end in G (EnshG1, PerinGa, Fen_Sub_G) or are astrocyte (Astro*),
#   pseudocartridge (PseudoCG*) or LQ glia.
#
# - Reclustered datasets (ozel_2021, Kurmangaliyev_2020): the PNAS cluster table
#   (pnas.2614314123.sd01.xlsx, sheet 1), looked up by the PNAS cluster number
#   used as cluster ID (see pnas_cluster_ids.R). A merged ID such as 1+2 takes
#   the values its parts share.
# - ozel_2021_simon_2026: Simon et al. (2026) Supplementary Table
#   41593_2026_2256_MOESM8_ESM.xlsx (sheet "Developmental origin"), one 0/1
#   column per spatial domain (vDpp ... dDpp), per temporal window (Hth ...
#   D_BarHI) and N_ON, looked up by the Ozel et al. (2021) cluster number.
#
# Spatial domains are collapsed to Vsx, Optix and Dpp (dorsal and ventral
# alike) and joined with "/" in that order, e.g. Vsx/Optix, as in the PNAS table.
#
# Usage (needs readxl):
#   source("script/load.R", chdir = TRUE)
#   source("script/misc/origin_metadata.R")
#   meta <- pnas_origin_metadata("pnas.2614314123.sd01.xlsx", "data/x.clusters.h5ad")
#   meta <- moesm8_origin_metadata("41593_2026_2256_MOESM8_ESM.xlsx", "data/x.clusters.h5ad")
#   set_cluster_metadata("data/x.clusters.h5ad", meta)

ORIGIN_LEVELS <- c("Temporal origin", "Spatial origin", "Notch status")

# P15 developmental groups: neuroblasts, GMCs, immature neurons, lamina
# precursors, apoptotic cells
OTHER_GROUPS <- "^(NB|GMC[0-9]+|Im[0-9]+|LPC|Apop)$"
GLIA_NAMES <- "(G[0-9]*[ab]?$|^Astro|^PseudoCG|^LQ[0-9]*$)"

origin_text <- function(x) {
  #' Trimmed text with "", "NA" and NA as NA
  x <- trimws(as.character(x))
  x[!is.na(x) & x %in% c("", "NA")] <- NA
  x
}

origin_table <- function(clusters, by_id, source, out = NULL, cell_class = NULL) {
  #' data.frame(cluster_id, Cell class, Origin, <ORIGIN_LEVELS>) for every
  #' cluster of `clusters` from `by_id` (one row per source cluster number, row
  #' names = numbers). Merged IDs ("1+2") keep the values their parts share.
  #' Cell class comes from by_id$class, else from `cell_class` (id -> class);
  #' P15 groups are Others, and anything left is Neuron.
  ids <- levels(h5ad_open(clusters)$obs$cluster)
  levels <- c("class", ORIGIN_LEVELS)
  if (!"class" %in% names(by_id)) by_id$class <- NA_character_
  rows <- lapply(ids, function(id) {
    parts <- strsplit(id, "+", fixed = TRUE)[[1]]
    vals <- by_id[match(parts, rownames(by_id)), levels, drop = FALSE]
    vapply(levels, function(level) {
      v <- unique(vals[[level]])
      if (length(v) == 1) v else NA_character_
    }, "")
  })
  meta <- data.frame(cluster_id = ids, do.call(rbind, rows), check.names = FALSE)
  class <- meta$class
  if (!is.null(cell_class)) class[is.na(class)] <- cell_class[ids][is.na(class)]
  class[grepl(OTHER_GROUPS, ids)] <- "Others"
  class[is.na(class)] <- "Neuron"
  known <- rowSums(is.na(meta[ORIGIN_LEVELS])) == 0
  meta <- data.frame(cluster_id = ids, `Cell class` = class,
                     Origin = ifelse(known, "OPC", "Non-OPC"),
                     meta[ORIGIN_LEVELS], check.names = FALSE)
  if (!is.null(out)) {
    write.table(meta, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  }
  missing <- setdiff(unlist(strsplit(ids, "+", fixed = TRUE)), rownames(by_id))
  message(sprintf("%d clusters (%s), %d OPC; %d cluster number(s) not in %s (%s)",
                  length(ids),
                  paste(names(table(class)), table(class), sep = " ", collapse = ", "),
                  sum(known), length(missing), basename(source),
                  paste(head(missing, 10), collapse = ", ")))
  meta
}

pnas_origin_metadata <- function(xlsx, clusters, out = NULL) {
  #' Origin metadata of a reclustered dataset from the PNAS cluster table
  tbl <- as.data.frame(readxl::read_excel(xlsx, sheet = 1, col_types = "text",
                                         .name_repair = "minimal"))
  col <- function(pattern) {
    i <- grep(pattern, names(tbl))
    if (length(i) != 1) stop(xlsx, ": expected one column matching ", pattern)
    origin_text(tbl[[i]])
  }
  id <- origin_text(tbl[[1]])
  if (anyDuplicated(id)) stop("Cluster number(s) listed twice: ", paste(unique(id[duplicated(id)]), collapse = ", "))
  by_id <- data.frame(col("^Temporal Origin"), col("^Spatial Origin"), col("^Notch Status"),
                      row.names = id)
  names(by_id) <- ORIGIN_LEVELS
  # Spelled as in MOESM8 (D_BarHI), so both datasets share the value
  by_id[["Temporal origin"]][by_id[["Temporal origin"]] %in% "D/BarH"] <- "D/BarHI"
  note <- paste(col("^Annotation"), col("^Comments"))
  by_id$class <- ifelse(grepl("glia", note, ignore.case = TRUE), "Glia",
                        ifelse(grepl("non-neuronal|contaminating", note, ignore.case = TRUE),
                               "Others", NA))
  origin_table(clusters, by_id, xlsx, out)
}

moesm8_origin_metadata <- function(xlsx, clusters, out = NULL) {
  #' Origin metadata of ozel_2021_simon_2026 from Simon et al. (2026) MOESM8
  tbl <- as.data.frame(readxl::read_excel(xlsx, sheet = "Developmental origin",
                                         col_types = "text"))
  domains <- list(Vsx = c("vVsx", "dVsx"), Optix = c("vOptix", "dOptix"),
                  Dpp = c("vDpp", "dDpp"))
  windows <- c("Hth", "Hth_Opa", "Opa_Erm", "Erm_Ey", "Ey_Hbn", "Hbn_Opa_Slp",
               "Slp_D", "D_BarHI")
  need <- c("Cluster", unlist(domains), windows, "N_ON")
  if (!all(need %in% names(tbl))) {
    stop(xlsx, " lacks column(s) ", paste(setdiff(need, names(tbl)), collapse = ", "))
  }
  flag <- function(cols) {
    #' 0/1 matrix of `cols`; text such as "NA" becomes NA
    m <- sapply(cols, function(c) suppressWarnings(as.integer(origin_text(tbl[[c]]))))
    matrix(m, nrow(tbl), dimnames = list(NULL, cols))
  }
  pick <- function(m, labels) {
    #' "/"-joined labels of the columns set to 1; NA when none is or any is NA
    apply(m, 1, function(r) {
      if (anyNA(r) || !any(r == 1)) NA_character_ else paste(labels[r == 1], collapse = "/")
    })
  }
  sp <- sapply(domains, function(cols) {
    m <- flag(cols)
    ifelse(rowSums(is.na(m)) == ncol(m), NA, rowSums(m, na.rm = TRUE) > 0) * 1L
  })
  sp <- matrix(sp, nrow(tbl), dimnames = list(NULL, names(domains)))
  notch <- flag("N_ON")[, 1]

  id <- origin_text(tbl$Cluster)
  if (anyDuplicated(id)) stop("Cluster number(s) listed twice: ", paste(unique(id[duplicated(id)]), collapse = ", "))
  by_id <- data.frame(pick(flag(windows), gsub("_", "/", windows, fixed = TRUE)),
                      pick(sp, names(domains)),
                      ifelse(is.na(notch), NA, ifelse(notch == 1, "On", "Off")),
                      row.names = id)
  names(by_id) <- ORIGIN_LEVELS
  ch <- h5ad_open(clusters)
  name <- read_cluster_names(ch$obs, levels(ch$obs$cluster), clusters)
  cell_class <- ifelse(grepl(GLIA_NAMES, name), "Glia",
                       ifelse(grepl(OTHER_GROUPS, name), "Others", NA))
  origin_table(clusters, by_id, xlsx, out, setNames(cell_class, names(name)))
}
