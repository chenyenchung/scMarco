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

# Load every dataset once per process (shared across sessions)
datasets <- lapply(names(config$datasets), function(id) {
  message("Loading dataset ", id, "...")
  ds_load(id, config$datasets[[id]])
})
names(datasets) <- names(config$datasets)

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
      selectizeInput(
        'clusters_of_interest',
        'Cluster(s) of interest',
        choices = NULL,
        multiple = TRUE
      ),
      sliderInput(
        "pos_cut", "Probability threshold (ON if > threshold)",
        min = 0, max = 1, step = 0.05, value = 0.5
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
        tabPanel("UMAP", umap_viz_UI("umap"))
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
        tabPanel(title = "UMAP", umap_viz_output("umap"))
      )
    )
  )
)

server <- function(input, output, session) {
  ds <- reactive({
    shiny::req(input$dataset)
    datasets[[input$dataset]]
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

  b <- reactive(ds_binarize(ds(), input$pos_cut))

  observeEvent(ds(), {
    updateSelectizeInput(session, "clusters_of_interest",
                         choices = cluster_choices(ds()), selected = character(0))
    first_mark(character(0))

    has_trend <- ds()$has_mean && ds()$ordered_splits
    has_umap <- length(ds()$reductions) > 0
    for (tab in list(
      list("controlTabs", "Plot expression trend", has_trend),
      list("outputTabs", "Expression Trend", has_trend),
      list("controlTabs", "UMAP", has_umap),
      list("outputTabs", "UMAP", has_umap)
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
    idents = reactive(input$clusters_of_interest),
    required = required,
    min_n = min_n,
    cut = reactive(input$pos_cut),
    split_label = reactive(ds()$split_label),
    first_mark = first_mark,
    parent_session = session
  )

  find_exp_server("fexp", ctx)
  find_pos_server("fpos", ctx)
  coexp_session <- coexp_server("coexp", ctx)
  find_distinct_server("fd_genes", ctx, coexp_session)
  lineplot_server("lplot", ctx)
  umap_viz_server("umap", ctx)
}

shinyApp(ui = ui, server = server)
