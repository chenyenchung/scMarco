#!/usr/bin/env Rscript
source("./src/utils.R")

OnPlot <- function(db, db_tbl, sql_where, marker_tbl, idents, cut_off, stages) {
  #' This function takes a list of expression matrices
  #' and plots on off plot based on the cut.off values provided

  # Check if WHERE statement exists already
  sql_where <- sql_where
  if (sql_where == "") {
    sql_where <- c()
  }

  sql_where <- c(
    sql_where,
    paste(
      "gene IN", vec2sqllist(marker_tbl$distinct_gene)
    ),
    paste(
      "cluster IN", vec2sqllist(idents)
    )
  )
  to_plot <- dbGetQuery(
    db, paste(
      "SELECT * FROM", db_tbl, wherecat(sql_where)
    )
  )

  to_plot <- reshape2::melt(
    ConvertToWideFormat(to_plot, formula = "gene ~ stage", stages = stages),
    id.vars = c("cluster", "gene"),
    value.name = "prob",
    variable.name = "stage"
  )

  to_plot$exp <- ifelse(
    to_plot$prob >= cut_off,
    "Expressed",
    "Not Expressed"
  )


  box_size <- 18 - 0.75 * length(unique(to_plot$cluster))

  p <- ggplot(to_plot, aes(x = stage, y = gene, fill = exp)) +
    geom_point(pch = 22, size = box_size) +
    scale_fill_manual(
      values = c("Expressed" = "black", "Not Expressed" = "white"),
      na.value = "white"
    ) +
    scale_x_discrete(limits = stages) +
    labs(fill = "") +
    facet_wrap(~cluster) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.key = element_blank(),
      legend.text = element_text(size = 16)
    )
  return(p)
}

SplitPlot <- function(
    db, db_tbl, sql_where, genes, cut_off, count = FALSE, plot_all = TRUE, stages
) {
  #' This function takes a list of expression matrices
  #' and plots on off plot based on the cut.off values provided
  #' x: A data.frame that is a long form exrpession matrix
  #' genes: A character vector containing the genes of interest
  #' cut_off: A numeric value between 0 and 1 that defines the threshold to
  #' be seen as expressed
  #' count: A logical value determining whether the coexpressed plot should
  #' be shown binarized or the total value.
  #' plot_all: A logical value determining whether heatmap for coexpression
  #' and individual values should be shown together.

  # Floating point number inaccuracy tolerance
  tol <- 1e-6

  # Check if WHERE statement exists already
  sql_where <- sql_where
  if (sql_where == "") {
    sql_where <- c()
  }

  sql_where <- c(
    sql_where,
    paste(
      "gene IN", vec2sqllist(genes)
    )
  )
  to_plot <- dbGetQuery(
    db, paste(
      "SELECT * FROM", db_tbl, wherecat(sql_where)
    )
  )

  to_plot$exp <- to_plot$prob >= cut_off

  # Find all clusters that at least express the given genes once
  all_clust <- unique(unlist(GetAllExpressedClusters(to_plot, cut_off)))

  # Exit early if none of the selected genes are expressed
  if (length(all_clust) == 0) {
    p <- ggplot(
      data.frame(message = "Sorry.\nNone of your selected genes is expressed")
    ) +
      geom_text(aes(label = message), x = 0.5, y = 0.5, size = 10) +
      theme_void()
    return(p)
  }


  to_plot <- subset(
    to_plot,
    cluster %in% all_clust
  )

  # Split expression matrix by gene
  to_plot <- ConvertToWideFormat(
    to_plot, formula = "cluster ~ stage",
    value.var = "exp", stages = stages
  )

  # Get the cluster names to make sure every gene matrix
  # has the same order of clusters
  row_order <- all_clust

  # Iteratively count the number of times a gene is positively detected
  coexp <- Reduce(
    `+`,
    lapply(
      names(to_plot),
      function(gene_label) {
        # Get the matrix per label
        gene_mat <- to_plot[[gene_label]]

        # Deal with the case when a cluster is not present in all
        missing_clust <- setdiff(row_order, gene_mat$cluster)

        # By adding 0 for the missing clusters
        stage_pad_0 <- as.list(rep(0, length(stages)))
        names(stage_pad_0) <- stages

        if (length(missing_clust) != 0) {
          gene_mat <- rbind(
            gene_mat,
            data.frame(
              cluster = missing_clust,
              stage_pad_0,
              gene = gene_label,
              check.names = FALSE
            )
          )
        }
        # Set cluster names as row indices
        row.names(gene_mat) <- gene_mat$cluster

        # Binarize gene expression probability
        prob_cols <- stages
        return((1 - gene_mat[row_order, prob_cols]) < tol)
      }
    )
  )

  # Convert coexpression count to a data.frame
  coexp <- as.data.frame(coexp)

  # Order the clusters by number of co-expressed clusters
  ## Get number of co-expression per cluster
  exp_row <- rowSums(
    coexp[ , stages] == length(genes)
  )
  exp_row <- sort(exp_row, decreasing = FALSE)


  # Annotate the cluster labels
  coexp$cluster <- row_order
  coexp$gene <- "Coexpressed"

  to_plot[["Coexpression"]] <- coexp

  # Convert back to long form
  to_plot <- do.call(rbind, to_plot)
  to_plot <- reshape2::melt(
    to_plot,
    id.vars = c("cluster", "gene"),
    value.name = "exp",
    variable.name = "stage"
  )

  # Order stages
  to_plot$stage <- factor(to_plot$stage,
                          levels = stages)

  to_plot$cluster <- factor(
    to_plot$cluster,
    levels = names(exp_row)
  )

  # Determining output layout
  if (plot_all) {
    # Binarizing
    exp_val <- to_plot$exp
    exp_val[to_plot$gene != "Coexpressed"] <- ifelse(
      (1 - to_plot[to_plot$gene != "Coexpressed", "exp"]) < tol,
      "(Co)expressed",
      "Not (co)expressed"
    )
    exp_val[to_plot$gene == "Coexpressed"] <- ifelse(
      (length(genes) - to_plot[to_plot$gene == "Coexpressed", "exp"]) < tol ,
      "(Co)expressed",
      "Not (co)expressed"
    )
    to_plot$exp <- exp_val

    # Set panel order
    to_plot$gene <- factor(to_plot$gene)
    to_plot$gene <- relevel(to_plot$gene, "Coexpressed")

    box_size <- max(18 - 0.75 * length(unique(to_plot$cluster)), 8)
    p <- ggplot(
      to_plot,
      aes(x = stage, y = cluster, label = exp)
    ) +
      # geom_label(size = box_size) +
      geom_point(aes(fill = exp), pch = 22, size = box_size) +
      scale_fill_manual(
        values = c("(Co)expressed" = "black", "Not (co)expressed" = "white"),
        na.value = "white"
      ) +
      labs(fill = "") +
      facet_wrap(~gene) +
      theme(
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_rect(color = "black", fill = "white"),
        strip.text = element_text(size = 20, face = "bold"),
        legend.key = element_blank(),
        legend.text = element_text(size = 16)
      )
    return(p)
  }

  # If not plotting all
  co_plot <- subset(to_plot, gene == "Coexpressed")
  co_plot$count <- as.integer(co_plot$exp)
  co_plot$exp <- ifelse(
    (length(genes) - co_plot$exp) < tol ,
    "Co-expressed",
    "Not co-expressed"
  )


  box_size <- max(18 - 0.75 * length(unique(to_plot$cluster)), 8)

  p <- ggplot(
    co_plot,
    aes(x = stage, y = cluster, label = count)
  ) +
    scale_fill_manual(
      values = c("Co-expressed" = "black", "Not co-expressed" = "white"),
      na.value = "white"
    ) +
    labs(fill = "") +
    facet_wrap(~gene) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 20, face = "bold"),
      legend.key = element_blank(),
      legend.text = element_text(size = 16)
    )

  # Determine whether coexpression is presented as count or color
  if (count) {
    p <- p + geom_label(size = box_size)
  } else {
    p <- p + geom_point(aes(fill = exp), pch = 22, size = box_size)
  }
  return(p)
}

LinePlot <- function(
    db, db_tbl, sql_where, features, highlight_idents,
    idents = NULL,
    stages,
    highlight_col = "#57068c",
    lowlight_dim = 0.2,
    ggobj = FALSE
) {
  #' This functions plot log-normalized expression overtime
  #' x: A path to a SQLite database containing log-normalized average expression
  #' db_tbl: Name of table containing log-normalized average expression
  #' genes: A character vector containing the genes to plot
  #' highlight_idents: A character vector containing the clusters to highlight
  #' idents: A character vector containing the clusters to plot. All clusters
  #' will be plotted if NULL
  #' num_trend: A numeric value defining how many different gene expression
  #' trend to be used to split the line plot.
  #' highlight_col: A character value (hex color) defining the color of the
  #' lines corresponding to the clusters of interest
  #' lowlight_dim: A numeric value defining the alpha of all other clusters

  # Check if WHERE statement exists already
  sql_where <- sql_where
  if (sql_where == "") {
    sql_where <- c()
  }

  sql_where <- c(
    sql_where,
    paste(
      "gene IN", vec2sqllist(features)
    )
  )

  # Only keep the gene of interest
  gene_to_plot <- dbGetQuery(
    db, paste(
      "SELECT * FROM", db_tbl, wherecat(sql_where)
    )
  )

  if (!is.null(idents)) {
    gene_to_plot <- subset(gene_to_plot, cluster %in% idents)
  }

  # Set color panel for the clusters to highlight
  highlight_palette <- ifelse(
    unique(gene_to_plot$cluster) %in% highlight_idents, highlight_col, "black"
  )
  names(highlight_palette) <- unique(gene_to_plot$cluster)
  highlight_alpha <- ifelse(
    unique(gene_to_plot$cluster) %in% highlight_idents, 1, lowlight_dim
  )
  names(highlight_alpha) <- unique(gene_to_plot$cluster)

  if (ggobj) {
    gene_to_plot$hlcluster <- "Others"
    gene_to_plot$hlcluster[gene_to_plot$cluster %in% highlight_idents] <-
      gene_to_plot$cluster[gene_to_plot$cluster %in% highlight_idents]
    gene_to_plot$hlcluster <- factor(
      gene_to_plot$hlcluster,
      levels = c(highlight_idents, "Others")
    )
    p <- ggplot(
      gene_to_plot,
      aes(x = stage, y = lognorm, group = cluster)
    ) +
      geom_line(aes(color = hlcluster, alpha = cluster)) +
      scale_color_manual(
        values = c(
          viridisLite::viridis(n = length(highlight_idents)), "black"
        )
      ) +
      scale_alpha_manual(values = highlight_alpha) +
      scale_x_discrete(limits = stages,
                       labels = stages) +
      guides(alpha = "none") +
      labs(
        x = "",
        y = "Log-normalized expression",
        color = ""
      ) +
      theme_cowplot() +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "white")
      )

    return(p)
  }

  # Plot it
  p <- ggplot(
    gene_to_plot,
    aes(x = stage, y = lognorm, group = cluster)
  ) +
    geom_line(aes(color = cluster, alpha = cluster)) +
    scale_color_manual(values = highlight_palette) +
    scale_alpha_manual(values = highlight_alpha) +
    scale_x_discrete(limits = stages,
                     labels = stages) +
    guides(color = "none", alpha = "none") +
    labs(
      x = "",
      y = "Log-normalized expression"
    ) +
    theme_cowplot() +
    theme(
      legend.position = 'none',
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(color = "white")
    )


  pinteract <- ggplotly(p, tooltip = "alpha")

  return(pinteract)
}
