lineplot_UI <- function(id) {
    tagList(
        selectizeInput(
            inputId = NS(id, 'gene_of_interest'),
            label = 'Which gene?',
            choices = NULL,
            multiple = FALSE
        ),
        # The opacity of other clusters only applies to the line plot
        conditionalPanel(
            "output.trend_mode != 'box'", ns = NS(id),
            sliderInput(
                NS(id, "line_alpha"),
                "How opaque other clusters should be?",
                min = 0.1, max = 1, step = 0.1, value = 0.2
            )
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
        uiOutput(NS(id, "trend_ui"))
    )
}

lineplot_server <- function(id, ctx) {
    moduleServer(
        id,
        function(input, output, session) {
            output$trend_mode <- renderText({
                mode <- ctx$ds()$trend
                if (is.null(mode)) "" else mode
            })
            outputOptions(output, "trend_mode", suspendWhenHidden = FALSE)

            observe({
                updateSelectizeInput(
                    session, "gene_of_interest", choices = ctx$all_gene_choices(), server = TRUE
                )
            })

            request <- eventReactive(
                input$do_lineplot, {
                    ok <- require_or_modal(
                        nzchar(input$gene_of_interest),
                        "Please select a gene", "Select a gene to plot its expression."
                    )
                    if (!ok) return(NULL)
                    ds <- ctx$ds()
                    list(gene = input$gene_of_interest, highlight = ctx$idents(),
                        alpha = input$line_alpha, ds = ds, mode = ds$trend)
                }
            )

            # One panel per cluster of interest in box mode
            n_box_panels <- function(req) max(length(req$highlight), 1)

            draw <- function(req, interactive) {
                if (req$mode == "box") {
                    box_plot(req$ds, req$gene, req$highlight, interactive = interactive)
                } else {
                    trend_plot(req$ds, req$gene, req$highlight, req$alpha,
                        interactive = interactive, box = req$mode == "both")
                }
            }

            observeEvent(
                input$do_lineplot, {
                    updateTabsetPanel(
                        ctx$parent_session,
                        inputId = "outputTabs",
                        selected = "Expression Levels"
                    )
                    output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_lplot"))})
                }
            )

            output$trend_ui <- renderUI({
                req <- request()
                shiny::req(req)
                height <- if (req$mode == "box") ceiling(n_box_panels(req) / 3) * 400 + 100 else 800
                plotly::plotlyOutput(NS(id, "lineplot"), height = paste0(height, "px"))
            })

            output$lineplot <- plotly::renderPlotly({
                req <- request()
                shiny::req(req)
                draw(req, interactive = TRUE)
            })

            download_size <- reactive({
                req <- request()
                if (req$mode == "box") box_size(n_box_panels(req)) else c(width = 1500, height = 1000)
            })

            ggdownload_server(
                "dl_lplot",
                filename = reactive(c(request()$gene, cluster_label(request()$ds, request()$highlight))),
                payload = function() draw(request(), interactive = FALSE),
                height = function() download_size()[["height"]],
                width = function() download_size()[["width"]]
            )
        }
    )
}
