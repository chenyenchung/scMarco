# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(yaml)
library(data.table)
library(ggplot2)
library(reshape2)
library(plotly)
library(markdown)
library(cowplot)
options(ragg.max_dim = Inf)

# TODO: Config-dependent interface
# TODO: SQLite encapsulation of input

config <- read_yaml("config.yaml")
stages <- names(config$mixture_model_tbl)

norm <- fread("int/lognorm.csv", data.table = FALSE)
mm <- fread("int/probs.csv", data.table = FALSE)

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
  titlePanel("Marker Combo Selector"),

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
        unique(mm$cluster),
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
  mm_work <- reactive({
    req(input$gene_group_to_use)
    if (input$gene_group_to_use == "All") {
      return(mm)
    } else {
      active_filter <- names(gene_list)[gene_list == input$gene_group_to_use]
      gene_of_interest <- readLines(config[["gene_list"]][[active_filter]][["path"]])
      working <- subset(mm, gene %in% gene_of_interest)
      return(working)
    }
  }
  )

  ##### Find distinct genes for given clusters
  fd_genes <- find_distinct_server(
    "fd_genes", stages, session,
    mm_work, reactive({input$clusters_of_interest}), reactive({input$pos_cut})
  )

  ##### Plot line plot for log-normalized expression for a selected marker
  lineplot <- lineplot_server(
    "lplot", session, norm, stages,
    mm_work, reactive({input$clusters_of_interest})
  )

  ##### Plot coexpression plot for given genes
  splitplot <- coexp_server(
    "coexp", session, stages, mm,
    mm_work, reactive(input$pos_cut)
  )

}

# Run the application
shinyApp(ui = ui, server = server)
