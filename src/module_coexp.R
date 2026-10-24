coexp_UI <- function(id) {
  tagList(
    selectizeInput(
      inputId = NS(id, 'split_genes'),
      label = 'Which genes?',
      choices = NULL,
      multiple = TRUE
    ),
    uiOutput(NS(id, "pct_toggle")),
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

coexp_server <- function(id, ctx) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        updateSelectizeInput(
          session, 'split_genes', choices = ctx$all_gene_choices(), server = TRUE
        )
      })

      output$pct_toggle <- renderUI({
        ds <- ctx$ds()
        if (!ds$has_pct && is.null(ds$cells_h)) return(NULL)
        label <- if (is.null(ds$cells_h)) {
          "Show % cells having >0 UMIs (per gene)"
        } else {
          "Show % cells having >0 UMIs (per gene, and co-expressing cells)"
        }
        checkboxInput(NS(id, "show_pct"), label, value = FALSE)
      })

      request <- eventReactive(
        input$do_splitplot, {
          ok <- require_or_modal(
            length(input$split_genes) > 0,
            "Please select genes",
            "Select one or more genes to plot their co-expression."
          )
          if (!ok) return(NULL)
          list(genes = input$split_genes, show_pct = isTRUE(input$show_pct),
               b = ctx$b(), ds = ctx$ds())
        }
      )

      splitplot <- reactive({
        req <- request()
        if (is.null(req)) return(message_plot("Select genes and click 'Plot co-expression'."))
        with_spinner("Plotting co-expression...", coexp_plot(
          req$ds, req$b, req$genes, show_pct = req$show_pct
        ))
      })

      plot_dims <- reactive({
        req <- request()
        if (is.null(req)) return(c(width = 800, height = 300))
        coexp_size(splitplot(), length(req$genes), req$show_pct)
      })

      observeEvent(
        input$do_splitplot, {
          updateTabsetPanel(
            ctx$parent_session,
            inputId = "outputTabs",
            selected = "Co-Expression"
          )
          output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_coexp"))})
        }
      )

      output$splitplot <- renderPlot(
        splitplot(),
        height = function() plot_dims()[["height"]]
      )

      ggdownload_server(
        "dl_coexp",
        filename = reactive(input$split_genes),
        payload = splitplot,
        height = function() plot_dims()[["height"]] * 2,
        width = function() plot_dims()[["width"]] * 2
      )
      return(session)
    }
  )
}
