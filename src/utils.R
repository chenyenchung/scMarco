ConvertToWideFormat <- function(
    mat_of_interest, formula = "gene ~ stage", value.var = "prob", stages
) {
  #' Convert a long-format expression matrix
  #' (with cluster, gene, stage, expression probability)
  #' To a list of wide format matrises (gene * stage)

  # Replace undetected genes with 0
  possible_cols <- c("gene", "stage", "cluster")
  # Get the cols used to dcast
  metaname <- setdiff(possible_cols, strsplit(formula, " ~ ")[[1]])

  # Split by cluster into a list
  mat_of_interest <- split(mat_of_interest, mat_of_interest[[metaname]])

  # Convert to wide format
  widemat <- lapply(
    names(mat_of_interest),
    function(x) {
      expwid <- reshape2::dcast(
        data = mat_of_interest[[as.character(x)]],
        formula = as.formula(formula),
        value.var = value.var,
        fun.aggregate = mean
      )

      # Check if there are missing stages
      missing <- setdiff(
        stages,
        colnames(expwid)
      )

      for (missing_col in missing) {
        expwid[[missing_col]] <- NA
      }


      expwid[is.na(expwid)] <- 0
      expwid[[metaname]] <- x
      return(expwid)
    }
  )
  names(widemat) <- names(mat_of_interest)
  return(widemat)
}

CheckFilterInputDb <- function(db, cut_off) {
  #' This function examines the anticipated input variable type
  #' for GetAllExpressedGenes() and GetAllExpressedClusters()

  # Input type-checking
  if (class(db) != "SQLiteConnection") {
    stop("The input data (db) is expected to be a SQLiteConnection.")
  }
  if (!is.numeric(cut_off) | cut_off > 1 | cut_off < 0) {
    stop("cut_off must be a number between 0 ~ 1.")
  }
}

CheckFilterInput <- function(exp_mat, cut_off) {
  #' This function examines the anticipated input variable type
  #' for GetAllExpressedGenes() and GetAllExpressedClusters()

  # Input type-checking
  if (!is.data.frame(exp_mat)) {
    stop("The input data (exp_mat) is expected to be a data.frame.")
  }
  if (!is.numeric(cut_off) | cut_off > 1 | cut_off < 0) {
    stop("cut_off must be a number between 0 ~ 1.")
  }
}

GetAllExpressedClusters <- function(exp_mat, cut_off) {
  #' This function gets all other clusters that at least express a gene
  #' once
  #' exp_mat: A data.frame (long form expression matrix)
  #' cut_off: A numeric value between 0 and 1

  # Input type-checking
  CheckFilterInput(exp_mat, cut_off)

  # Removed un-expressed genes
  exp_mat_exp <- exp_mat[exp_mat$prob >= cut_off, ]

  # Get the clusters each gene is expressed
  exp_cluster <- tapply(
    exp_mat_exp$cluster,
    exp_mat_exp$gene,
    unique
  )
  return(exp_cluster)
}

GetAllExpressedGenesDb <- function(
    db, db_tbl, cut_off, sql_where, stage = NULL, count = 1, stages
) {
  #' This function takes a long form expression matrix and returns a list
  #' of expressed genes per cluster based on the cut off provided.
  #' exp_mat: A data.frame (long form expression matrix)
  #' cut_off: A numeric value between 0 and 1
  #' stage: A Boolean vector of 6 items signifying selection of P15, 30, 40, 50,
  #' 70, and Adult
  #' count: The minimal number of stages for a gene to be expressed to be
  #' returned (only effective when stage = NULL)

  # Input type-checking
  CheckFilterInputDb(db, cut_off)

  # Check if WHERE statement exists already
  sql_where <- sql_where
  if (sql_where == "") {
    sql_where <- c()
  }

  # Removed un-expressed genes
  sql_where <- c(sql_where, paste("prob >=", cut_off))
  sql_where <- wherecat(sql_where)

  # exp_mat_exp <- exp_mat[exp_mat$prob >= cut_off, ]
  exp_mat_exp <- dbGetQuery(
    db,
    paste(
      "SELECT * FROM", db_tbl, sql_where)
  )

  # Get all genes expressed at certain stages
  exp_gene <- tapply(
    exp_mat_exp$gene,
    exp_mat_exp$cluster,
    unique
  )

  # If filter for specific stages
  if (any(stage)) {
    stage_choose <- stages[stage]
    exp_mat_exp <- subset(exp_mat_exp, stage %in% stage_choose)
    stage_num <- sum(stage)
  }

  # Otherwise filter for stage number
  stage_num <- count

  # Get all genes that passed the filter criteria for expression
  filtered_exp_gene <- tapply(
    exp_mat_exp$gene,
    exp_mat_exp$cluster,
    function(x) {
      gene_occurrence <- table(x)
      gene_pass <- names(gene_occurrence)[gene_occurrence >= stage_num]
      return(gene_pass)
    }
  )


  return(list(all = exp_gene, filtered = filtered_exp_gene))
}

FindDistinctiveGenes <- function(
    x,
    idents,
    max_presence = 1,
    group = NA,
    against_all = FALSE
) {
  #' Return genes that are is expressed in fewer than [max_presence] clusters (
  #' each stage is counted as 1 presence).
  #' in each group of selected clusters
  #' x: A gene list returned by GetAllExpressedGenes() in which each item
  #' represents a cluster
  #' idents: A character vector of cluster labels to analyze
  #' max_presence: A numeric value defining the maximum number of presence to
  #' be considered distinctive
  #' against_all: A logical value deciding if a gene must be NOT present in
  #' any other cluster to be considered as distinct. (TRUE: All clusters are
  #' considered; FALSE: Only consider clusters of interest)

  # Check input type
  if (!is.list(x)) {
    stop("x must be a list generated by GetAllExpressedGenes().")
  }
  if (!all(idents %in% names(x$filtered))) {
    label_not_found <- setdiff(idents, names(x$filtered))
    stop(
      paste0(
        "The following clusters are not found: ", label_not_found, "."
      )
    )
  }

  if (against_all) {
    # Get the clusters of interest
    exp_interest <- x$filtered[idents]

    # Get all genes expressed in other clusters at any stage
    other_clust <- x$all[setdiff(names(x$filtered), idents)]

    # Find genes that are expressed in more than [max_presence] cluster
    dup_list <- append(exp_interest, other_clust)
    dup_list <- unlist(dup_list)
    dup_freq <- table(dup_list)
    drop_list <- names(dup_freq)[dup_freq > max_presence]

    # Remove genes with multiple appearance
    distinct_list <- lapply(
      names(exp_interest),
      function(cluster) {

        # Find distinct genes by removing genes with multiple presence
        distinct_genes <- setdiff(exp_interest[[cluster]], drop_list)

        # If theres no distinct gene, go to the next cluster
        if (length(distinct_genes) == 0) {
          return(NULL)
        }

        # Return a table of distinct genes
        distinct_df <- data.frame(
          cluster = cluster,
          distinct_gene = distinct_genes,
          group = group
        )
        return(distinct_df)
      }
    )
  } else {
    distinct_list <- lapply(
      idents, function(ident) {
        # Get genes that passed filtering in this cluster
        this_cluster <- x$filtered[ident]

        # Get all genes that are expressed in other clusters
        expressed_in_others <- unlist(
          x$all[setdiff(idents, ident)]
        )

        # Filter by number of appearances
        if (max_presence > 1) {
          dup_freq <- table(expressed_in_others)
          drop_list <- names(dup_freq)[dup_freq > max_presence]
        } else {
          drop_list <- unique(expressed_in_others)
        }

        # Return genes that are "unique"
        # Find distinct genes by removing genes with multiple presence
        distinct_genes <- setdiff(x$filtered[[ident]], drop_list)

        # If theres no distinct gene, go to the next cluster
        if (length(distinct_genes) == 0) {
          return(NULL)
        }

        # Return a table of distinct genes
        distinct_df <- data.frame(
          cluster = ident,
          distinct_gene = distinct_genes,
          group = group
        )
        return(distinct_df)
      }
    )

  }

  # Return a message if there's no gene found
  if (all(sapply(distinct_list, is.null))) {
    return(
      data.frame(
        data.frame(
          cluster = "Not found",
          distinct_gene = "Not found",
          group = "Not found"
        )
      )
    )
  }

  # Return if there's no problem
  distinct_table <- do.call(rbind, distinct_list)
  return(distinct_table)
}

FindDistinctiveGenesPerGroup <- function(
    x, idents, max_presence = 1, group = NULL
) {
  # If there's no group
  if (is.null(group)) {
    return(FindDistinctiveGenes(x, idents, max_presence))
  }

  x_group <- split(x, group)

  per_group <- lapply(names(x_group), function(group_label) {
    cluster_of_interest <- x_group[[group_label]]
    return(
      FindDistinctiveGenes(x, idents, max_presence, group = group_label)
    )
  })

  per_group_table <- do.call(rbind, per_group)
}

vec2sqllist <- function(x) {
  # Escape single quotes by making them double (SQL)
  x <- gsub("'", "''", x)
  capture_output <- capture.output(cat(paste0("'", x, "'"), sep = ","))
  capture_output <- paste0("(", capture_output, ")")
  return(capture_output)
}

wherecat <- function(x) {
  out <- paste(x, collapse = " AND ")
  out <- paste("WHERE", out)
  return(out)
}
