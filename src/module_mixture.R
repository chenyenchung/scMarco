mixture_UI <- function(id) {
  tagList(
    selectizeInput(
      inputId = NS(id, "gene_of_interest"),
      label = "Which gene?",
      choices = NULL,
      multiple = FALSE
    ),
    uiOutput(NS(id, "model_ui")),
    actionButton(NS(id, "do_mixture"), "Plot mixture model"),
    uiOutput(NS(id, "cond_dl"))
  )
}

mixture_output <- function(id) {
  tagList(
    plotOutput(NS(id, "mixture_plot"), height = "850px")
  )
}

mixture_server <- function(id, ctx) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        updateSelectizeInput(
          session, "gene_of_interest", choices = ctx$all_genes(), server = TRUE
        )
      })

      models <- reactive(mixture_model_choices(ctx$ds()))

      # One choice per mixture model: a split, or splits modeled together
      output$model_ui <- renderUI({
        choices <- models()
        if (length(choices) < 2) return(NULL)
        selected <- isolate(input$model)
        if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
        joint <- any(grepl(" + ", names(choices), fixed = TRUE))
        label <- if (joint) "Which model?" else paste0("Which ", ctx$split_label(), "?")
        selectInput(NS(id, "model"), label, choices = choices, selected = selected)
      })

      request <- eventReactive(
        input$do_mixture, {
          ok <- require_or_modal(
            nzchar(input$gene_of_interest),
            "Please select a gene", "Select a gene to plot its mixture model."
          )
          if (!ok) return(NULL)
          choices <- models()
          model <- if (is.null(input$model) || !input$model %in% choices) {
            choices[[1]]
          } else {
            input$model
          }
          list(ds = ctx$ds(), gene = input$gene_of_interest, model = model,
               highlight = ctx$idents(), cut = ctx$cut())
        }
      )

      observeEvent(
        input$do_mixture, {
          updateTabsetPanel(ctx$parent_session, inputId = "outputTabs",
                            selected = "Mixture Model")
          output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_mixture"))})
        }
      )

      mixture_plot_obj <- reactive({
        r <- request()
        shiny::req(r)
        with_spinner("Plotting...", mixture_plot(
          r$ds, r$gene, r$model, highlight = r$highlight, cut_off = r$cut
        ))
      })

      output$mixture_plot <- renderPlot(mixture_plot_obj())

      ggdownload_server(
        "dl_mixture",
        filename = reactive({
          r <- request()
          c("mixture", r$gene, mixture_model_splits(r$ds, r$model),
            cluster_label(r$ds, r$highlight))
        }),
        payload = mixture_plot_obj,
        height = 1400,
        width = 1500
      )
    }
  )
}
