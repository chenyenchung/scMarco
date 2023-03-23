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
options(ragg.max_dim = Inf)


config <- read_yaml("config.yaml")

current_db <- config$databases$ozel_2021
stages <- current_db$stages

db <- dbConnect(RSQLite::SQLite(), current_db$path)
mm_tbl <- "probs"
gex_tbl <- "lognorm"

source("src/utils.R")
source("src/visualization.R")
source("src/module_find_distinct.R")
source("src/module_coexp.R")
source("src/module_lineplot.R")
source("src/module_ggdownload.R")

gene_list <- vapply(
  config$gene_list, function(x) {
    return(x[["label"]])
  },
  FUN.VALUE = character(1)
)

ui <- fluidPage(
  # Application title
  titlePanel(config[["title"]]),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "gene_group_to_use",
        "Which genes to consider?",
        c("All", unname(gene_list)),
        selected = "All",
        multiple = FALSE
      ),
      selectInput(
        'clusters_of_interest',
        'Clusters of interest',
        dbGetQuery(db, paste("SELECT DISTINCT cluster FROM", mm_tbl)),
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
          "Find distinct genes for clusters of interest",
          find_distinct_UI("fd_genes", stages)
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
        tabPanel(title = "Read Me", includeMarkdown("doc/README.md")),
        tabPanel(title = "On/Off", find_distinct_output("fd_genes")),
        tabPanel(title = "Co-Expression", coexp_output("coexp")),
        tabPanel(title = "Expression Trend", lineplot_output("lplot"))
      )
    )
  )
)


server <- function(input, output, session) {

  # Filter genes of interest from config file if necessary
  sql_where <- reactive({
    req(input$gene_group_to_use)

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

  ##### Find distinct genes for given clusters
  find_distinct_server(
    "fd_genes", stages, session, db, mm_tbl,
    sql_where, reactive({input$clusters_of_interest}), reactive({input$pos_cut})
  )

  ##### Plot coexpression plot for given genes
  coexp_server(
    "coexp", session, stages, db, mm_tbl,
    sql_where, reactive(input$pos_cut)
  )

  ##### Plot line plot for log-normalized expression for a selected marker
  lineplot_server(
    "lplot", session, stages, db, gex_tbl,
    sql_where, reactive({input$clusters_of_interest})
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
