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
      "Plot log-normalized exp over stages"
    ),
    uiOutput(NS(id, "cond_dl"))
  )
}

lineplot_output <- function(id) {
  tagList(
    plotlyOutput(NS(id, "lineplot"), height = "800px")
  )
}

lineplot_server <- function(
    id, parent_session, stages, db, db_tbl,
    # Reactive params
    sql_where, idents, genes
) {
  moduleServer(
    id,
    function(input, output, session) {
      observe(
        {
          # Generate SQL statement based on global filter condition
          this_sql_where <- sql_where()
          if (this_sql_where != "") {
            # Grow WHERE clause if there is already elements
            this_sql_where <- c(this_sql_where, "type = 'lognorm'")
          } else {
            # Start from scratch if WHERE clause was empty
            this_sql_where <- "type = 'lognorm'"
          }

          input_q <- paste(
            "SELECT DISTINCT gene FROM", db_tbl(), wherecat(this_sql_where)
          )

          # Update the gene selection panel for line plot
          updateSelectizeInput(
            session,
            inputId = "gene_of_interest",
            choices = genes(),
            server = TRUE
          )
        }
      )

      lineplot <- eventReactive(
        input$do_lineplot,
        {
          LinePlot(
            db,
            db_tbl(),
            sql_where = sql_where(),
            features = input$gene_of_interest,
            highlight_idents = idents(),
            stages = stages(),
            lowlight_dim = input$line_alpha,
            highlight_col = "red"
          )
        }
      )

      # # Switch to the on/off plot panel when asked to plot
      observeEvent(
        input$do_lineplot,
        {
          updateTabsetPanel(
            parent_session,
            inputId = "outputTabs",
            selected = "Expression Trend"
          )
        }
      )

      output$lineplot <- renderPlotly(
        lineplot()
      )

      # Show download button after plot is generated
      observeEvent(input$do_lineplot, {
        output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_lplot"))})
      })

      ggdownload_server(
        "dl_lplot",
        filename = idents,
        payload = function() {
            LinePlot(
              db,
              db_tbl(),
              sql_where = sql_where(),
              features = input$gene_of_interest,
              highlight_idents = idents(),
              stages = stages(),
              lowlight_dim = input$line_alpha,
              ggobj = TRUE
            )
          },
        height = function() {return(1000)},
        width = 1500
      )
    }
  )
}
