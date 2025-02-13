# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(yaml)
library(ggplot2)
library(reshape2)
library(plotly)
library(markdown)
library(cowplot)
library(DBI)
library(RSQLite)
library(shinybusy)
library(bslib)
options(ragg.max_dim = Inf)
options(bitmapType='cairo')

config <- read_yaml("config.yaml")
if (is.null(config$cache_path)) {
  config$cache_path <- "int/cache.rds"
}

cache_last_mod <- ifelse(
  file.exists(config$cache_path),
  file.info(config$cache_path)$mtime,
  0
)

db <- dbConnect(RSQLite::SQLite(), config$database_path)

db_last_mod <- file.info(config$database_path)$mtime

if (!file.exists(config$cache_path) | cache_last_mod <= db_last_mod) {
  message("Generating new cache...")
  dbList <- names(config$databases)
  names(dbList) <- dbList

  stages <- lapply(
    dbList, function(tbl) {
      distinct_stages <- dbGetQuery(
        db,
        paste("SELECT DISTINCT stage FROM", tbl)
      )$stage
      return(distinct_stages)
    }
  )

  clusters <- lapply(
    dbList, function(tbl) {
      distinct_clusters <- dbGetQuery(
        db,
        paste("SELECT DISTINCT cluster FROM", tbl)
      )$cluster
      return(distinct_clusters)
    }
  )

  genes <- lapply(
    dbList, function(tbl) {
      distinct_genes <- dbGetQuery(
        db,
        paste("SELECT DISTINCT gene FROM", tbl)
      )$gene
      return(distinct_genes)
    }
  )

  saveRDS(
    list(stages = stages, clusters = clusters, genes = genes),
    config$cache_path
  )
}

cache <- readRDS(config$cache_path)

source("src/utils.R")
source("src/visualization.R")
source("src/module_find_distinct.R")
source("src/module_find_positive.R")
source("src/module_coexp.R")
source("src/module_lineplot.R")
source("src/module_ggdownload.R")
source("src/module_find_exp.R")

gene_list <- vapply(
  config$gene_list, function(x) {
    return(x[["label"]])
  },
  FUN.VALUE = character(1)
)

db_list <- vapply(
  config$databases, function(x) {
    return(x[["label"]])
  },
  FUN.VALUE = character(1)
)

ui <- fluidPage(
  theme = bs_theme(bootswatch = "journal"),
  # Application title
  titlePanel(config[["title"]]),
  add_busy_bar(color = "#57068c", height = "12px"),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      tags$div(style = "display: none",
        selectInput(
        inputId = "cache_last_g",
        label = "",
        choices = NULL,
        multiple = TRUE
      )),
      selectInput(
        inputId = "dataset",
        label = "Which dataset to use?",
        choices = unname(db_list),
        selected = unname(db_list)[1]
      ),
      selectInput(
        "gene_group_to_use",
        "Which genes to consider?",
        c("All", unname(gene_list)),
        selected = "All",
        multiple = FALSE
      ),
      selectizeInput(
        'clusters_of_interest',
        'Clusters of interest',
        choices = NULL,
        multiple = TRUE
      ),
      sliderInput(
        "pos_cut", "Probability threshold to be\nconsidered as expressed (0 - 1)",
        min = 0, max = 1, step = 0.05, value = 0.5
      ),
      tabsetPanel(
        id = "controlTabs",
        type = "tabs",
        tabPanel(
          "From cluster",
          uiOutput("fexp")
        ),
        tabPanel(
          "From gene",
          uiOutput("fpos")
        ),
        tabPanel(
          "Find distinct genes within",
          uiOutput("fd_genes")
        ),
        tabPanel(
          "Co-expression plot",
          coexp_UI("coexp")
        ),
        tabPanel(
          "Plot expression trend",
          lineplot_UI("lplot")
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "outputTabs",
        type = "tabs",
        tabPanel(title = "Read Me", includeMarkdown("README.md")),
        tabPanel(title = "Gene Table", find_exp_output("fexp")),
        tabPanel(title = "Cluster Table", find_exp_output("fpos")),
        tabPanel(title = "On/Off", find_distinct_output("fd_genes")),
        tabPanel(title = "Co-Expression", coexp_output("coexp")),
        tabPanel(title = "Expression Trend", lineplot_output("lplot"))
      )
    )
  )
)


server <- function(input, output, session) {

  # Filter genes of interest from config file if necessary
  sql_where <- eventReactive(
    input$gene_group_to_use, {
    if (input$gene_group_to_use == "All") {
      return("")
    } else {
      active_filter <- names(gene_list)[gene_list == input$gene_group_to_use]
      gene_of_interest <- readLines(config[["gene_list"]][[active_filter]][["path"]])
      query_sentence <- paste(
        "gene IN", vec2sqllist(gene_of_interest)
      )
      return(query_sentence)
    }
  }
  )

  # Retrieve cluster and stage information from the current dataset
  current_tbl <- eventReactive(
    input$dataset, {
      dataset_id <- names(db_list)[db_list == input$dataset]
      return(dataset_id)
    }
  )

  current_stages <- eventReactive(
    input$dataset, {
      return(config[['databases']][[current_tbl()]]$stages)
    }
  )

  update_genelist <- reactive({
    paste0(input$dataset, input$gene_group_to_use)
  })

  current_genes <- eventReactive(
    update_genelist(), {
      if (input$gene_group_to_use == "All") {
        return(cache[['genes']][[current_tbl()]])
      } else {
        active_filter <- names(gene_list)[gene_list == input$gene_group_to_use]
        genes <- intersect(
          # The list
          readLines(config[["gene_list"]][[active_filter]][["path"]]),
          # All genes
          cache[['genes']][[current_tbl()]]
        )
        return(genes)
      }
    }
  )

  current_clusters <- eventReactive(
    input$dataset, {return(cache$clusters[[current_tbl()]])}
  )

  prob_mat <- eventReactive(
    input$gene_group_to_use, {
      prob_mat <- dbGetQuery(
        db,
        paste(
          "SELECT * FROM",
          current_tbl(), wherecat(
            c(sql_where(), paste("value >", input$pos_cut), "type = 'prob'")
          )
        )
      )
      return(prob_mat)
    }
  )

  observe(
    {
      updateSelectizeInput(
        session, "clusters_of_interest",
        choices = current_clusters()
      )
      output$fd_genes <- renderUI(
        {find_distinct_UI("fd_genes", current_stages())}
      )
      output$fexp <- renderUI(
        {find_exp_UI("fexp", current_stages())}
      )

      output$fpos <- renderUI(
        {find_pos_UI("fpos", current_stages())}
      )
    }
  )


  ##### Find expressed genes for a cluster
  gene_exp <- find_exp_server(
    id = "fexp",
    parent_session = session,
    prob_mat = prob_mat,
    stages = current_stages,
    db = db,
    db_tbl = current_tbl,
    sql_where = sql_where,
    idents = reactive({input$clusters_of_interest}),
    cut_off = reactive({input$pos_cut})
  )

  ##### Find clusters that express a gene
  find_pos_server(
    id = "fpos",
    parent_session = session,
    stages = current_stages,
    db = db,
    db_tbl = current_tbl,
    sql_where = sql_where,
    idents = reactive({input$clusters_of_interest}),
    cut_off = reactive({input$pos_cut}),
    genes = current_genes
  )

  ##### Find distinct genes for given clusters
  find_distinct_server(
    id = "fd_genes",
    stages = current_stages,
    parent_session = session,
    db = db,
    db_tbl = current_tbl,
    sql_where = sql_where,
    prob_mat = prob_mat,
    idents = reactive({input$clusters_of_interest}),
    cut_off = reactive({input$pos_cut}),
    first_mark = reactive({input$cache_last_g}),
    coexp_session = coexp_session
  )

  ##### Plot coexpression plot for given genes
  coexp_session <- coexp_server(
    id = "coexp",
    parent_session = session,
    stages = current_stages,
    db = db,
    db_tbl = current_tbl,
    sql_where =  sql_where,
    cut_off = reactive(input$pos_cut),
    genes = current_genes
  )

  ##### Plot line plot for log-normalized expression for a selected marker
  lineplot_server(
    id =  "lplot",
    parent_session = session,
    stages = current_stages,
    db = db,
    db_tbl = current_tbl,
    sql_where,
    idents = reactive({input$clusters_of_interest}),
    genes = current_genes
  )
}

# Run the application
shinyApp(
  ui = ui,
  server = server,
  onStart = function() {
    # Close connection to SQLite db when the app is stopped
    onStop(
      function() {dbDisconnect(db)}
    )
  }
)
