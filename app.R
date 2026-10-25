# scMarco: select marker combinations from binarized single-cell expression.
#
# One codebase for every lab's instance: datasets, gene lists, title, theme
# and README all come from the config file. Set SCMARCO_CONFIG to use a
# config other than ./config.yaml.

library(shiny)
library(bslib)
library(yaml)
library(ggplot2)
library(markdown)
library(shinybusy)
options(bitmapType = "cairo")
options(ragg.max_dim = Inf)

config_path <- Sys.getenv("SCMARCO_CONFIG", "config.yaml")
config <- read_yaml(config_path)

source("src/combo_core.R")
source("src/backend_h5ad.R")
source("src/plot_core.R")
source("src/ui_helpers.R")
source("src/module_ggdownload.R")
source("src/module_find_exp.R")
source("src/module_find_positive.R")
source("src/module_find_distinct.R")
source("src/module_coexp.R")
source("src/module_lineplot.R")
source("src/module_umap_viz.R")
source("src/module_mixture.R")

if (length(config$datasets) == 0) {
  stop("No datasets configured in ", config_path)
}
dataset_labels <- vapply(config$datasets, function(x) x$label, character(1))
if (anyDuplicated(dataset_labels)) {
  stop("Dataset labels must be unique in ", config_path)
}
default_dataset <- names(config$datasets)[
  vapply(config$datasets, function(x) isTRUE(x$default), logical(1))
]
default_dataset <- if (length(default_dataset) > 0) default_dataset[1] else
  names(config$datasets)[1]

gene_list_labels <- vapply(
  config$gene_list, function(x) x$label, character(1)
)

# Datasets load on first use and are then shared across sessions (ds_load
# memoizes per process). Check them all now so a bad config still fails at
# startup, and warm the default so the landing page is ready.
ds_check_config(config$datasets)
message("Loading default dataset ", default_dataset, "...")
invisible(ds_load(default_dataset, config$datasets[[default_dataset]]))

# readme: one Markdown file or a list of them (e.g. lab intro + shared usage)
readme_paths <- if (is.null(config$readme)) "readme/usage.md" else unlist(config$readme)

ui <- fluidPage(
  theme = bs_theme(bootswatch = if (is.null(config$theme)) "minty" else config$theme),
  titlePanel(config$title),
  add_busy_bar(color = "#57068c", height = "12px"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "dataset",
        label = "Which dataset to use?",
        choices = setNames(names(dataset_labels), dataset_labels),
        selected = default_dataset
      ),
      selectInput(
        "gene_group_to_use",
        "Which genes to consider?",
        c("All" = "All", setNames(names(gene_list_labels), gene_list_labels)),
        selected = "All"
      ),
      uiOutput("meta_filter"),
      selectizeInput(
        'clusters_of_interest',
        'Cluster(s) of interest',
        choices = NULL,
        multiple = TRUE
      ),
      sliderInput(
        "margin", "Confidence margin (%)",
        min = 1, max = 50, step = 1, value = 10
      ),
      helpText(
        "ON if P(On) \u2265 100% \u2212 margin, OFF if P(On) \u2264 margin;",
        "in between is ambiguous (gray) and never counts as ON or OFF."
      ),
      uiOutput("splits_selector"),
      tabsetPanel(
        id = "controlTabs",
        type = "tabs",
        tabPanel("From cluster", find_exp_UI("fexp")),
        tabPanel("From gene", find_pos_UI("fpos")),
        tabPanel("Find distinct genes within", find_distinct_UI("fd_genes")),
        tabPanel("Co-expression plot", coexp_UI("coexp")),
        tabPanel("Plot expression trend", lineplot_UI("lplot")),
        tabPanel("UMAP", umap_viz_UI("umap")),
        tabPanel("Mixture model", mixture_UI("mixture"))
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "outputTabs",
        type = "tabs",
        tabPanel(title = "Read Me", lapply(readme_paths, includeMarkdown)),
        tabPanel(title = "Gene Table", find_exp_output("fexp")),
        tabPanel(title = "Cluster Table", find_pos_output("fpos")),
        tabPanel(title = "On/Off", find_distinct_output("fd_genes")),
        tabPanel(title = "Co-Expression", coexp_output("coexp")),
        tabPanel(title = "Expression Trend", lineplot_output("lplot")),
        tabPanel(title = "UMAP", umap_viz_output("umap")),
        tabPanel(title = "Mixture Model", mixture_output("mixture"))
      )
    )
  )
)

server <- function(input, output, session) {
  ds <- reactive({
    shiny::req(input$dataset)
    ds_load(input$dataset, config$datasets[[input$dataset]])
  })

  gene_filter <- reactive({
    if (input$gene_group_to_use == "All") return(NULL)
    readLines(config$gene_list[[input$gene_group_to_use]]$path)
  })

  # Genes used for combination finding (bimodal genes when known)
  genes <- reactive({
    g <- ds()$analysis_genes
    if (!is.null(gene_filter())) g <- intersect(g, gene_filter())
    g
  })

  # Genes for per-cell and trend views (every gene in the dataset)
  all_genes <- reactive({
    g <- ds()$all_genes
    if (!is.null(gene_filter())) g <- intersect(g, gene_filter())
    g
  })

  # Dropdown choices for views that also show unimodal genes, labeled as such
  all_gene_choices <- reactive(gene_choices(ds(), all_genes()))

  # Debounced: each distinct margin is a full binarize, so don't run one per
  # intermediate slider position
  margin <- debounce(reactive(input$margin / 100), 500)

  # Cluster metadata: the user picks which levels (config subset_by) to subset
  # by, and each picked level shows its own filter. Levels cascade coarse to
  # fine: each narrows the choices of the levels below it. Searches only see
  # clusters that match. An unpicked level is "all".
  output$meta_filter <- renderUI({
    meta <- ds()$cluster_meta
    if (is.null(meta)) return(NULL)
    tagList(
      selectizeInput("meta_levels", "Subset clusters by", choices = names(meta),
                     multiple = TRUE, options = list(placeholder = "none (all clusters)")),
      lapply(seq_along(meta), function(i) {
        conditionalPanel(
          sprintf("(input.meta_levels || []).indexOf(%s) > -1",
                  encodeString(names(meta)[i], quote = '"')),
          selectizeInput(
            paste0("meta_", i), paste0("Search within: ", names(meta)[i]),
            choices = sort(unique(meta[[i]][!is.na(meta[[i]])])), multiple = TRUE,
            options = list(placeholder = "all")
          )
        )
      }),
      helpText(textOutput("meta_scope", inline = TRUE))
    )
  })

  meta_filters <- reactive({
    meta <- ds()$cluster_meta
    if (is.null(meta)) return(list())
    active <- input$meta_levels
    # Values left over from another dataset are ignored, not matched
    setNames(lapply(seq_along(meta), function(i) {
      if (!names(meta)[i] %in% active) return(character(0))
      as.character(intersect(input[[paste0("meta_", i)]], meta[[i]]))
    }), names(meta))
  })

  # A level that is no longer picked goes back to "all"
  observeEvent(input$meta_levels, ignoreNULL = FALSE, {
    meta <- ds()$cluster_meta
    for (i in seq_along(meta)) {
      key <- paste0("meta_", i)
      if (!names(meta)[i] %in% input$meta_levels && length(input[[key]]) > 0) {
        updateSelectizeInput(session, key, selected = character(0))
      }
    }
  })

  universe <- reactive(meta_universe(ds(), meta_filters()))

  # Last choices sent per level, so a level is only updated when they change
  # (updating closes an open dropdown)
  meta_sent <- new.env(parent = emptyenv())
  observeEvent(ds(), rm(list = ls(meta_sent), envir = meta_sent))
  observeEvent(meta_filters(), {
    filters <- meta_filters()
    for (i in seq_along(filters)[-1]) {
      level <- names(filters)[i]
      choices <- meta_level_choices(ds(), filters, level)
      kept <- as.character(intersect(filters[[i]], choices))
      key <- paste0("meta_", i)
      if (!identical(meta_sent[[key]], choices) || !identical(kept, filters[[i]])) {
        updateSelectizeInput(session, key, choices = choices, selected = kept)
        meta_sent[[key]] <- choices
      }
      filters[i] <- list(kept)
    }
  })

  output$meta_scope <- renderText({
    sprintf("%d of %d clusters searched.", length(universe()), length(ds()$clusters))
  })

  b <- reactive(restrict_clusters(ds_binarize(ds(), margin()), universe()))

  # The cluster picker lists the clusters in scope; a new dataset starts empty
  shown_ds <- NULL
  observeEvent(universe(), {
    keep <- if (identical(shown_ds, ds()$id)) {
      intersect(isolate(input$clusters_of_interest), universe())
    } else {
      character(0)
    }
    shown_ds <<- ds()$id
    updateSelectizeInput(session, "clusters_of_interest",
                         choices = cluster_choices(ds(), universe()), selected = keep)
  })

  observeEvent(ds(), {
    first_mark(character(0))

    has_trend <- ds()$has_mean && ds()$ordered_splits
    has_umap <- length(ds()$reductions) > 0
    has_mixture <- ds()$has_mixture
    for (tab in list(
      list("controlTabs", "Plot expression trend", has_trend),
      list("outputTabs", "Expression Trend", has_trend),
      list("controlTabs", "UMAP", has_umap),
      list("outputTabs", "UMAP", has_umap),
      list("controlTabs", "Mixture model", has_mixture),
      list("outputTabs", "Mixture Model", has_mixture)
    )) {
      if (tab[[3]]) showTab(tab[[1]], tab[[2]]) else hideTab(tab[[1]], tab[[2]])
    }
  })

  output$splits_selector <- renderUI({
    splits <- ds()$splits
    label <- ds()$split_label
    tagList(
      checkboxGroupInput(
        "required_splits",
        paste0("Must be ON in these ", label, "s"),
        choices = splits,
        inline = TRUE
      ),
      sliderInput(
        "min_n",
        paste0("...and ON in at least this many ", label, "s"),
        min = 1, max = length(splits), step = 1, value = 1
      )
    )
  })

  required <- reactive({
    r <- intersect(input$required_splits, ds()$splits)
    if (is.null(r)) character(0) else r
  })

  min_n <- reactive({
    n <- if (is.null(input$min_n)) 1 else input$min_n
    # Requiring k splits implies at least k splits
    min(max(n, length(required()), 1), length(ds()$splits))
  })

  first_mark <- reactiveVal(character(0))

  ctx <- list(
    ds = ds,
    b = b,
    genes = genes,
    all_genes = all_genes,
    all_gene_choices = all_gene_choices,
    idents = reactive(input$clusters_of_interest),
    required = required,
    min_n = min_n,
    margin = margin,
    split_label = reactive(ds()$split_label),
    scope = reactive(describe_scope(ds(), meta_filters())),
    first_mark = first_mark,
    parent_session = session
  )

  find_exp_server("fexp", ctx)
  find_pos_server("fpos", ctx)
  coexp_session <- coexp_server("coexp", ctx)
  find_distinct_server("fd_genes", ctx, coexp_session)
  lineplot_server("lplot", ctx)
  umap_viz_server("umap", ctx)
  mixture_server("mixture", ctx)
}

shinyApp(ui = ui, server = server)
