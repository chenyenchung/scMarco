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
        plotOutput(NS(id, "mixture_plot"), height = "auto")
    )
}

mixture_server <- function(id, ctx) {
    moduleServer(
        id,
        function(input, output, session) {
            observe({
                updateSelectizeInput(
                    session, "gene_of_interest", choices = ctx$all_gene_choices(), server = TRUE
                )
            })

            models <- reactive(mixture_model_choices(ctx$ds()))

            # One checkbox per mixture model (a split, or splits modeled
            # together); the checked ones are stacked top to bottom
            output$model_ui <- renderUI({
                choices <- models()
                if (length(choices) < 2) return(NULL)
                selected <- intersect(isolate(input$model), choices)
                if (length(selected) == 0) selected <- choices
                joint <- any(grepl(" + ", names(choices), fixed = TRUE))
                label <- if (joint) "Models to stack:" else
                    paste0(tools::toTitleCase(ctx$split_plural()), " to stack:")
                checkboxGroupInput(NS(id, "model"), label, choices = choices, selected = selected)
            })

            request <- eventReactive(
                input$do_mixture, {
                    ok <- require_or_modal(
                        nzchar(input$gene_of_interest),
                        "Please select a gene", "Select a gene to plot its mixture model."
                    )
                    if (!ok) return(NULL)
                    choices <- models()
                    # In dataset order; a single model has no checkboxes
                    model <- if (length(choices) < 2) choices else choices[choices %in% input$model]
                    ok <- require_or_modal(
                        length(model) > 0, paste("No", ctx$split_label(), "selected"),
                        paste0("Please check at least one ", ctx$split_label(), " to plot.")
                    )
                    if (!ok) return(NULL)
                    list(ds = ctx$ds(), gene = input$gene_of_interest, model = model,
                        highlight = ctx$idents(), band = ctx$band())
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
                    r$ds, r$gene, r$model, highlight = r$highlight, band = r$band
                ))
            })

            n_models <- reactive({
                r <- request()
                if (is.null(r)) 1 else length(r$model)
            })

            output$mixture_plot <- renderPlot(
                mixture_plot_obj(),
                height = function() 850 * n_models()
            )

            ggdownload_server(
                "dl_mixture",
                filename = reactive({
                    r <- request()
                    splits <- unlist(lapply(r$model, mixture_model_splits, ds = r$ds))
                    c("mixture", r$gene, splits, cluster_label(r$ds, r$highlight))
                }),
                payload = mixture_plot_obj,
                height = function() 1400 * n_models(),
                width = 1500
            )
        }
    )
}
