# Plot builders shared by the modules. All take a dataset from ds_load() and,
# where relevant, a binarized view from ds_binarize(); none touch Shiny.

library(ggplot2)

STATE_FILL <- c(
  "Expressed" = "#1B1B1B",
  "Not expressed" = "#FFFFFF",
  "Non-bimodal" = "#D3D3D3",
  "Absent" = "#FFFFFF"
)

theme_scmarco_tiles <- function() {
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 20, face = "bold"),
    legend.key = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom"
  )
}

relabel_clusters <- function(ds, f) {
  #' A factor of cluster IDs with its levels shown as display labels
  factor(f, levels = levels(f), labels = cluster_label(ds, levels(f)))
}

message_plot <- function(text) {
  ggplot(data.frame(message = text)) +
    geom_text(aes(label = message), x = 0.5, y = 0.5, size = 8,
              color = "#555555") +
    theme_void()
}

tile_layers <- function(df, box_size) {
  #' df needs columns x, y, state; absent cells are drawn as a light cross
  present <- df[df$state != "Absent", , drop = FALSE]
  absent <- df[df$state == "Absent", , drop = FALSE]
  list(
    geom_point(data = present, aes(x = x, y = y, fill = state),
               shape = 22, size = box_size),
    if (nrow(absent) > 0) {
      geom_point(data = absent, aes(x = x, y = y, shape = "Not in this split"),
                 size = box_size / 2.5, color = "#9A9A9A")
    },
    if (nrow(absent) > 0) {
      scale_shape_manual(values = c("Not in this split" = 4), name = NULL)
    },
    scale_fill_manual(values = STATE_FILL, drop = TRUE, name = NULL)
  )
}

pct_label_layer <- function(df, box_size) {
  #' df needs x, y, pct (0-1); labels are drawn only where given
  df <- df[!is.na(df$pct), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  df$label <- paste0(round(df$pct * 100), "%")
  geom_text(data = df, aes(x = x, y = y, label = label), color = "white",
            size = max(3, box_size / 4), fontface = "bold", inherit.aes = FALSE)
}

# ---------------------------------------------------------------------------
# On/Off plot (Find distinct genes within)
# ---------------------------------------------------------------------------

onoff_plot <- function(ds, b, genes, idents, show_pct = FALSE) {
  genes <- unique(genes)
  if (length(genes) == 0) {
    return(message_plot("No distinctive gene found."))
  }
  df <- state_table(b, genes, clusters = idents)
  df$x <- factor(df$split, levels = ds$splits)
  df$y <- factor(df$gene, levels = rev(genes))
  df$panel <- factor(cluster_label(ds, df$cluster), levels = cluster_label(ds, idents))
  box_size <- if (show_pct) 20 else 12

  p <- ggplot(df, aes(x = x, y = y)) +
    tile_layers(df, box_size) +
    facet_wrap(~panel) +
    theme_scmarco_tiles()

  if (show_pct && ds$has_pct) {
    pct <- ds_cluster_layer_long(ds, "pct_pos", genes, clusters = idents)
    lab <- merge(df[df$state == "Expressed", c("gene", "cluster", "split", "x", "y", "panel")],
                 pct, by = c("gene", "cluster", "split"))
    lab$pct <- lab$value
    p <- p + pct_label_layer(lab, box_size)
  }
  p
}

onoff_size <- function(n_genes, n_clusters, show_pct = FALSE) {
  rows <- max(1, ceiling(n_clusters / 3))
  per_gene <- if (show_pct) 100 else 80
  c(
    width = min(n_clusters, 3) * 400 + 200,
    height = max(per_gene * n_genes * rows + 150, 500)
  )
}

# ---------------------------------------------------------------------------
# Co-expression plot
# ---------------------------------------------------------------------------

coexp_data <- function(ds, b, genes) {
  #' Long table of per-gene states plus a "Co-expressed" panel for every
  #' cluster in which at least one of the genes is ON in some split.
  genes <- intersect(unique(genes), b$genes)
  if (length(genes) == 0) return(NULL)
  anyon <- any_on(b, genes = genes)
  clusters <- b$clusters[colSums(anyon) > 0]
  if (length(clusters) == 0) return(NULL)

  st <- state_table(b, genes, clusters = clusters)
  on <- b$on[genes, clusters, , drop = FALSE]
  co <- apply(on, c(2, 3), all)
  co_df <- expand.grid(
    cluster = clusters, split = b$splits, stringsAsFactors = FALSE
  )
  co_df$gene <- "Co-expressed"
  co_df$state <- ifelse(
    !b$present[cbind(match(co_df$cluster, b$clusters),
                     match(co_df$split, b$splits))],
    "Absent",
    ifelse(co[cbind(match(co_df$cluster, clusters),
                    match(co_df$split, b$splits))],
           "Expressed", "Not expressed")
  )
  df <- rbind(st[, c("gene", "cluster", "split", "state")], co_df)

  # Order clusters: most splits co-expressed at the top
  n_co <- tapply(co_df$state == "Expressed", co_df$cluster, sum)
  cluster_order <- names(sort(n_co[clusters], decreasing = FALSE))
  df$cluster <- factor(df$cluster, levels = cluster_order)
  df$gene <- factor(df$gene, levels = c("Co-expressed", genes))
  df
}

coexp_cell_fraction <- function(ds, genes) {
  #' True fraction of cells co-expressing ALL genes (UMI > 0), per
  #' cluster x split, from the cell-level file. NULL if unavailable.
  if (is.null(ds$cells_h)) return(NULL)
  obs <- ds_cells_obs(ds)
  all_pos <- rep(TRUE, nrow(obs))
  for (g in genes) {
    all_pos <- all_pos & ds_cell_expression(ds, g) > 0
  }
  frac <- aggregate(
    all_pos, by = list(cluster = as.character(obs$cluster),
                       split = as.character(obs$split)),
    FUN = mean
  )
  names(frac)[3] <- "pct"
  frac$gene <- "Co-expressed"
  frac
}

coexp_plot <- function(ds, b, genes, show_pct = FALSE) {
  df <- coexp_data(ds, b, genes)
  if (is.null(df)) {
    return(message_plot("None of the selected genes is expressed."))
  }
  df$x <- factor(df$split, levels = ds$splits)
  df$y <- relabel_clusters(ds, df$cluster)
  n_clust <- nlevels(df$cluster)
  box_size <- if (show_pct) {
    max(26 - 0.75 * n_clust, 12)
  } else {
    max(18 - 0.75 * n_clust, 8)
  }

  p <- ggplot(df, aes(x = x, y = y)) +
    tile_layers(df, box_size) +
    facet_wrap(~gene) +
    theme_scmarco_tiles()

  if (show_pct) {
    labels <- list()
    if (ds$has_pct) {
      pct <- ds_cluster_layer_long(ds, "pct_pos", genes,
                                   clusters = levels(df$cluster))
      names(pct)[names(pct) == "value"] <- "pct"
      labels[[1]] <- pct
    }
    frac <- coexp_cell_fraction(ds, intersect(genes, b$genes))
    if (!is.null(frac)) labels[[2]] <- frac[, c("gene", "cluster", "split", "pct")]
    if (length(labels) > 0) {
      lab <- do.call(rbind, labels)
      on_df <- df[df$state == "Expressed", c("gene", "cluster", "split", "x", "y")]
      on_df$gene <- as.character(on_df$gene)
      on_df$cluster <- as.character(on_df$cluster)
      lab <- merge(on_df, lab, by = c("gene", "cluster", "split"))
      lab$gene <- factor(lab$gene, levels = levels(df$gene))
      lab$y <- relabel_clusters(ds, factor(lab$cluster, levels = levels(df$cluster)))
      p <- p + pct_label_layer(lab, box_size)
    }
  }
  p
}

coexp_size <- function(p, n_genes, show_pct = FALSE) {
  n_clust <- if (is.null(p$data$cluster)) 1 else nlevels(p$data$cluster)
  panel_rows <- ceiling((n_genes + 1) / 3)
  per_clust <- if (show_pct) 70 else 50
  c(
    width = min(n_genes + 1, 3) * 500 + 200,
    height = max(per_clust * n_clust, 300) * panel_rows + 150
  )
}

# ---------------------------------------------------------------------------
# Expression trend (ordered splits only)
# ---------------------------------------------------------------------------

trend_plot <- function(ds, gene, highlight = character(0), lowlight_alpha = 0.2,
                       interactive = TRUE) {
  df <- ds_cluster_layer_long(ds, "mean_lognorm", gene)
  df <- df[!is.na(df$value), , drop = FALSE]
  if (nrow(df) == 0) return(message_plot("Gene not found."))
  df$split <- factor(df$split, levels = ds$splits)
  df$highlight <- df$cluster %in% highlight
  df$cluster_label <- cluster_label(ds, df$cluster)
  highlight <- cluster_label(ds, highlight)
  df <- df[order(df$highlight), ]

  if (interactive) {
    p <- ggplot(df, aes(x = split, y = value, group = cluster, text = cluster_label)) +
      geom_line(aes(color = highlight, alpha = highlight)) +
      scale_color_manual(values = c(`FALSE` = "#333333", `TRUE` = "#C0392B"),
                         guide = "none") +
      scale_alpha_manual(values = c(`FALSE` = lowlight_alpha, `TRUE` = 1),
                         guide = "none") +
      labs(x = NULL, y = "Log-normalized expression", title = gene) +
      cowplot::theme_cowplot()
    return(plotly::ggplotly(p, tooltip = c("text", "y")))
  }

  df$label <- ifelse(df$highlight, df$cluster_label, "Others")
  df$label <- factor(df$label, levels = c(highlight, "Others"))
  colors <- c(viridisLite::viridis(length(highlight)), "#333333")
  names(colors) <- c(highlight, "Others")
  ggplot(df, aes(x = split, y = value, group = cluster)) +
    geom_line(aes(color = label, alpha = highlight)) +
    scale_color_manual(values = colors, name = NULL) +
    scale_alpha_manual(values = c(`FALSE` = lowlight_alpha, `TRUE` = 1),
                       guide = "none") +
    labs(x = NULL, y = "Log-normalized expression", title = gene) +
    cowplot::theme_cowplot()
}

# ---------------------------------------------------------------------------
# UMAP / reductions
# ---------------------------------------------------------------------------

reduction_label <- function(key) {
  name <- sub("^X_", "", key)
  ifelse(grepl("^(umap|tsne)_full$", name), "All cells", name)
}

split_reductions <- function(keys, splits) {
  #' Reductions computed separately per split, stored as X_<split> (umap) or
  #' X_<split>_<reduction> (seurat_to_h5ad.R with one object per split).
  #' Returns list(<reduction> = c(<split> = key, ...)) with splits in order.
  groups <- list()
  for (s in splits[order(-nchar(splits))]) {
    stem <- paste0("X_", s)
    hit <- keys == stem | startsWith(keys, paste0(stem, "_"))
    for (key in keys[hit]) {
      red <- if (key == stem) "umap" else substring(key, nchar(stem) + 2)
      groups[[red]][[s]] <- key
    }
    # The longest split name wins, e.g. X_A_B_tsne is split "A_B", not "A"
    keys <- keys[!hit]
  }
  lapply(groups, function(g) unlist(g[intersect(splits, names(g))]))
}

reduction_choices <- function(ds) {
  #' Selector entries: shared reductions by obsm key, per-split groups as
  #' "split:<reduction>". Per-split keys are folded into their group.
  groups <- split_reductions(ds$reductions, ds$splits)
  shared <- setdiff(ds$reductions, unlist(groups))
  names(shared) <- reduction_label(shared)
  if (length(groups) == 0) return(shared)
  per_split <- paste0("split:", names(groups))
  names(per_split) <- paste0(
    "Each ", ds$split_label, " separately",
    if (length(groups) > 1) paste0(" (", names(groups), ")") else ""
  )
  c(shared, per_split)
}

resolve_reduction <- function(ds, choice) {
  #' A selector value -> the obsm key, or c(<split> = key, ...) for a
  #' per-split group
  if (!startsWith(choice, "split:")) return(choice)
  split_reductions(ds$reductions, ds$splits)[[sub("^split:", "", choice)]]
}

reduction_axis <- function(key) {
  #' Axis title stem: reductions named *tsne* (seurat_to_h5ad.R --reduction tsne) are tSNE
  if (grepl("tsne", key, ignore.case = TRUE)) "tSNE" else "UMAP"
}

square_extents <- function(df, by_split = TRUE) {
  #' Two corner points that stretch the shorter axis to match the longer, per
  #' split or (by_split = FALSE) over all cells for facets sharing their scales
  corners <- function(d) {
    half <- max(diff(range(d$dim1)), diff(range(d$dim2))) / 2
    mid1 <- mean(range(d$dim1))
    mid2 <- mean(range(d$dim2))
    data.frame(dim1 = mid1 + c(-half, half), dim2 = mid2 + c(-half, half))
  }
  if (!by_split) return(corners(df))
  do.call(rbind, lapply(split(df, df$split, drop = TRUE), function(d) {
    cbind(corners(d), split = d$split[1])
  }))
}

umap_plot <- function(ds, key, color_by = c("gene", "cluster", "split"),
                      gene = NULL, splits = ds$splits, facet = FALSE,
                      show_legend = FALSE, ncol = NULL) {
  #' key: one obsm key shared by all splits, or c(<split> = key, ...) when each
  #' split has its own reduction; then every split gets its own facet
  #' ncol: facets per row (NULL lets facet_wrap decide)
  color_by <- match.arg(color_by)
  obs <- ds_cells_obs(ds)
  per_split <- !is.null(names(key))
  if (per_split) {
    key <- key[intersect(names(key), splits)]
    splits <- names(key)
    facet <- TRUE
    coords <- matrix(NA_real_, nrow(obs), 2)
    for (s in splits) {
      in_split <- as.character(obs$split) == s
      coords[in_split, ] <- ds_reduction(ds, key[[s]])[in_split, 1:2]
    }
  } else {
    coords <- ds_reduction(ds, key)
  }
  df <- data.frame(
    dim1 = coords[, 1], dim2 = coords[, 2],
    cluster = cluster_label(ds, obs$cluster),
    split = factor(as.character(obs$split), levels = ds$splits)
  )
  if (color_by == "gene") {
    df$expression <- ds_cell_expression(ds, gene)
  }
  keep <- !is.na(df$dim1) & df$split %in% splits
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) {
    return(message_plot("No cells in this reduction for the selected splits."))
  }
  split_label <- tools::toTitleCase(ds$split_label)
  title <- if (per_split) {
    paste(reduction_axis(key[[1]]), "of each", ds$split_label)
  } else {
    reduction_label(key)
  }
  # Embedding coordinates carry no meaning, so drop tick labels, ticks and grid
  base <- theme_minimal() +
    theme(text = element_text(size = 18),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          aspect.ratio = 1)

  p <- if (color_by == "gene") {
    df <- df[order(df$expression), ]
    ggplot(df, aes(dim1, dim2, color = expression)) +
      geom_point(size = 0.8) +
      scale_color_viridis_c(name = "Log-normalized\nexpression") +
      labs(subtitle = gene)
  } else if (color_by == "split") {
    ggplot(df, aes(dim1, dim2, color = split)) +
      geom_point(size = 0.8, alpha = 0.7) +
      scale_color_brewer(palette = "Set1", name = split_label)
  } else {
    q <- ggplot(df, aes(dim1, dim2, color = cluster)) +
      geom_point(size = 0.8, alpha = 0.7)
    if (show_legend) {
      q + scale_color_discrete(name = "Cluster")
    } else {
      group_cols <- if (facet) c("split", "cluster") else "cluster"
      centers <- aggregate(
        df[, c("dim1", "dim2")], by = df[, group_cols, drop = FALSE],
        FUN = function(v) mean(v, trim = 0.1)
      )
      q + ggrepel::geom_text_repel(
        data = centers, aes(dim1, dim2, label = cluster), color = "black",
        size = 5, max.overlaps = Inf, inherit.aes = FALSE
      ) + scale_color_discrete(guide = "none")
    }
  }
  # Separate reductions share no coordinates, so each facet gets its own axes
  if (facet) {
    p <- p + facet_wrap(~split, ncol = ncol, scales = if (per_split) "free" else "fixed")
  }
  # Square panels (legend aside) with 1 unit on dim1 == 1 unit on dim2: pad each
  # panel's data range to a square; base sets the panel aspect ratio to 1.
  # coord_fixed() alone would size panels by the data range instead
  p <- p + geom_blank(data = square_extents(df, by_split = per_split),
                      aes(dim1, dim2), inherit.aes = FALSE)
  axis <- reduction_axis(key[[1]])
  p + labs(title = title, x = paste(axis, 1), y = paste(axis, 2)) + base
}

umap_size <- function(n_splits, facet, ncol = 2) {
  if (!facet) return(c(width = 1000, height = 800))
  ncol <- min(n_splits, ncol)
  c(width = ncol * 700 + 100, height = ceiling(n_splits / ncol) * 600 + 100)
}

# ---------------------------------------------------------------------------
# Mixture model: cluster averages and P(On)
# ---------------------------------------------------------------------------

MIXTURE_COLORS <- c("OFF" = "#2a78d6", "ON" = "#eb6834")

# The pipeline rounds parameters to 2 decimals, so a stored 0 means < 0.005
MIXTURE_FLOOR <- 0.005

mixture_p_on <- function(x, fit) {
  #' P(On | x) for two Gaussians with a shared SD: a logistic curve in x
  sd <- max(fit$sd, MIXTURE_FLOOR)
  pi_on <- min(max(fit$pi, MIXTURE_FLOOR), 1 - MIXTURE_FLOOR)
  slope <- (fit$on_level - fit$off_level) / sd^2
  plogis(qlogis(pi_on) + slope * (x - (fit$on_level + fit$off_level) / 2))
}

mixture_plot <- function(ds, gene, model, highlight = character(0), cut_off = 0.5) {
  #' model: a mixture model of the dataset (see mixture_model_choices()); when
  #' it was fit across several splits, their split x cluster averages are pooled
  fit <- ds_mixture_fit(ds, gene, model)
  if (is.null(fit)) return(message_plot("This dataset has no mixture-model fit."))
  splits <- mixture_model_splits(ds, model)
  where <- paste(splits, collapse = " + ")
  joint <- length(splits) > 1
  if (anyNA(unlist(fit[c("on_level", "off_level", "sd", "pi")]))) {
    return(message_plot(paste0(gene, " was not modeled in ", where, ".")))
  }
  avg <- ds_cluster_average(ds, gene, splits)
  if (is.null(avg)) return(message_plot("No cluster averages for this dataset."))
  scale_note <- attr(avg, "scale")
  avg <- avg[!is.na(avg$value), , drop = FALSE]
  sd <- max(fit$sd, MIXTURE_FLOOR)
  floored <- fit$sd < MIXTURE_FLOOR || fit$pi < MIXTURE_FLOOR ||
    fit$pi > 1 - MIXTURE_FLOOR

  centers <- c(fit$off_level, fit$on_level)
  lims <- range(c(avg$value, centers - 3 * sd, centers + 3 * sd, 0))
  grid <- seq(lims[1], lims[2], length.out = 500)
  binwidth <- diff(lims) / 40

  # Fitted components, scaled to the histogram's counts
  n <- nrow(avg)
  comp <- rbind(
    data.frame(x = grid, component = "OFF",
               y = n * binwidth * (1 - fit$pi) * dnorm(grid, fit$off_level, sd)),
    data.frame(x = grid, component = "ON",
               y = n * binwidth * fit$pi * dnorm(grid, fit$on_level, sd))
  )
  comp$component <- factor(comp$component, levels = names(MIXTURE_COLORS))

  # One arrow per split x cluster pair of a cluster of interest
  pts <- avg[avg$cluster %in% highlight, , drop = FALSE]
  pts <- pts[order(match(pts$cluster, highlight), match(pts$split, splits)), , drop = FALSE]
  names(pts)[names(pts) == "value"] <- "x"
  pts$y <- mixture_p_on(pts$x, fit)
  pts$label <- sprintf("%s%s\nP(On) = %.2f", cluster_label(ds, pts$cluster),
                       if (joint) paste0(", ", pts$split) else "", pts$y)
  absent <- setdiff(as.character(highlight), avg$cluster)
  avg <- data.frame(x = avg$value)

  # Shared x range; zooming (not limits) keeps histogram bins at the edges
  x_scale <- coord_cartesian(xlim = lims, expand = TRUE)
  base <- cowplot::theme_cowplot()

  top <- ggplot(avg, aes(x = x)) +
    geom_histogram(binwidth = binwidth, boundary = 0, fill = "#B5B5B5",
                   color = "white", linewidth = 0.3) +
    geom_line(data = comp, aes(x = x, y = y, color = component), linewidth = 1) +
    scale_color_manual(values = MIXTURE_COLORS, name = "Fitted component") +
    x_scale +
    labs(x = NULL, y = "Number of clusters",
         title = paste0(gene, " in ", where),
         subtitle = paste0(
           if (joint) sprintf("Modeled together; each %s's clusters counted separately\n",
                              ds$split_label),
           sprintf(
             "%s; OFF mean %.2f, ON mean %.2f, shared SD %.2f, ON weight %.2f",
             if (fit$bimodal) "Bimodal" else "Not bimodal (never ON here)",
             fit$off_level, fit$on_level, fit$sd, fit$pi
           )
         )) +
    base +
    theme(legend.position = "top", legend.justification = "left")
  if (nrow(pts) > 0) {
    top <- top + geom_vline(xintercept = pts$x, linetype = "dashed",
                            color = "#1B1B1B", linewidth = 0.4)
  }

  curve <- data.frame(x = grid, y = mixture_p_on(grid, fit))
  bottom <- ggplot(curve, aes(x = x, y = y)) +
    geom_hline(yintercept = cut_off, linetype = "dotted", color = "#7A7A7A") +
    geom_line(linewidth = 1, color = "#1B1B1B") +
    geom_rug(data = avg, aes(x = x), inherit.aes = FALSE,
             sides = "b", alpha = 0.4, color = "#555555") +
    coord_cartesian(xlim = lims, ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(x = "Log-normalized expression (cluster average)", y = "P(On)") +
    base
  if (nrow(pts) > 0) {
    bottom <- bottom +
      geom_point(data = pts, aes(x = x, y = y), size = 3.5, color = "#1B1B1B",
                 fill = "white", shape = 21, stroke = 1.5) +
      ggrepel::geom_label_repel(
        data = pts, aes(x = x, y = y, label = label),
        nudge_y = ifelse(pts$y > 0.5, -0.3, 0.3), direction = "both",
        min.segment.length = 0, box.padding = 0.6, size = 4.5,
        arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
        max.overlaps = Inf, seed = 1
      )
  }

  notes <- c(
    if (identical(scale_note, "log of mean")) {
      "Cluster average: log1p of the mean expression in linear space, the scale the model was fit on."
    } else {
      "Cluster average: mean of log-normalized values (no cell-level file); the model was fit on log of the mean, so clusters may sit off the curve."
    },
    sprintf("Dotted line: the probability threshold (%.2f).", cut_off),
    if (floored) sprintf("Parameters are stored to 2 decimals; values below %.3f were drawn as %.3f.",
                         MIXTURE_FLOOR, MIXTURE_FLOOR),
    if (length(absent) > 0) {
      lbl <- cluster_label(ds, absent)
      paste0("Not in ", where, ": ",
             paste(ifelse(is.na(lbl), absent, lbl), collapse = ", "), ".")
    }
  )
  bottom <- bottom + labs(caption = paste(notes, collapse = "\n")) +
    theme(plot.caption = element_text(hjust = 0, size = 11, color = "#555555"))

  cowplot::plot_grid(top, bottom, ncol = 1, align = "v", axis = "lr",
                     rel_heights = c(1, 1.15))
}
