coexp_UI <- function(id) {
  tagList(
    selectizeInput(
      inputId = NS(id, 'split_genes'),
      label = 'Which genes?',
      choices = NULL,
      multiple = TRUE
    ),
    actionButton(
      NS(id, "do_splitplot"),
      "Plot co-expression"
    ),
    uiOutput(NS(id, "cond_dl"))
  )
}

coexp_output <- function(id) {
  tagList(
    plotOutput(NS(id, "splitplot"), height = "auto")
  )
}

coexp_server <- function(
    id, parent_session, stages, db, db_tbl,
    # Reactive params
    sql_where, cut_off, genes
) {
  moduleServer(
    id,
    function(input, output, session) {
      observe(
        {
          # Generate SQL statement based on global filter condition
          sql_where <- sql_where()
          input_q <- ifelse(
            sql_where == "",
            paste("SELECT DISTINCT gene FROM", db_tbl()),
            paste("SELECT DISTINCT gene FROM", db_tbl(), wherecat(sql_where))
          )

          updateSelectizeInput(
            session,
            'split_genes',
            choices = genes(),
            server = TRUE
          )
        }
      )

      splitplot <- eventReactive(
        input$do_splitplot,
        {
          updateSelectizeInput(
            session,
            'split_genes',
            choices = genes(),
            server = TRUE
          )
          SplitPlot(
            db = db,
            db_tbl = db_tbl(),
            sql_where = sql_where(),
            genes = input$split_genes,
            cut_off = cut_off(),
            count = FALSE,
            plot_all = TRUE,
            stages = stages()
          )
        }
      )

      # Switch to the split-expression plot panel when asked to plot
      observeEvent(
        input$do_splitplot,
        {
          updateTabsetPanel(
            parent_session,
            inputId = "outputTabs",
            selected = "Co-Expression"
          )
        }
      )

      output$splitplot <- renderPlot(
        splitplot(),
        height = function() {
          num_clust <- length(unique(splitplot()$data$cluster))
          use_height <- max(
            50 * num_clust,
            300
          )
          return(use_height)
        }
      )

      # Show download button after plot is generated
      observeEvent(input$do_splitplot, {
        output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_coexp"))})
      })

      ggdownload_server(
        "dl_coexp",
        filename = reactive(input$split_genes),
        payload = splitplot,
        height = function() {
          num_clust <- length(unique(splitplot()$data$cluster))
          use_height <- 2 *
            (as.integer((length(input$split_genes) + 1) / 3) + 1) *
            max(
            50 * num_clust,
            300
          )
          return(use_height)
        },
        width = 2500
      )
      return(session)
    }
  )
}
