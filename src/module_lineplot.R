lineplot_UI <- function(id) {
  tagList(
    selectizeInput(
      inputId = NS(id, 'gene_of_interest'),
      label = 'Which gene?',
      choices = NULL,
      multiple = FALSE
    ),
    sliderInput(
      NS(id, "line_alpha"),
      "How opaque other clusters should be?",
      min = 0.1, max = 1, step = 0.1, value = 0.2
    ),
    actionButton(
      NS(id, "do_lineplot"),
      "Plot log-normalized expression"
    ),
    uiOutput(NS(id, "cond_dl"))
  )
}

lineplot_output <- function(id) {
  tagList(
    plotly::plotlyOutput(NS(id, "lineplot"), height = "800px")
  )
}

lineplot_server <- function(id, ctx) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        updateSelectizeInput(
          session, "gene_of_interest", choices = ctx$all_genes(), server = TRUE
        )
      })

      request <- eventReactive(
        input$do_lineplot, {
          ok <- require_or_modal(
            nzchar(input$gene_of_interest),
            "Please select a gene", "Select a gene to plot its expression trend."
          )
          if (!ok) return(NULL)
          list(gene = input$gene_of_interest, highlight = ctx$idents(),
               alpha = input$line_alpha, ds = ctx$ds())
        }
      )

      observeEvent(
        input$do_lineplot, {
          updateTabsetPanel(
            ctx$parent_session,
            inputId = "outputTabs",
            selected = "Expression Trend"
          )
          output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_lplot"))})
        }
      )

      output$lineplot <- plotly::renderPlotly({
        req <- request()
        shiny::req(req)
        trend_plot(req$ds, req$gene, req$highlight, req$alpha, interactive = TRUE)
      })

      ggdownload_server(
        "dl_lplot",
        filename = reactive(c(request()$gene, cluster_label(request()$ds, request()$highlight))),
        payload = function() {
          req <- request()
          trend_plot(req$ds, req$gene, req$highlight, req$alpha, interactive = FALSE)
        },
        height = 1000,
        width = 1500
      )
    }
  )
}
