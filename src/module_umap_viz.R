umap_viz_UI <- function(id) {
  tagList(
    selectInput(NS(id, "umap_selection"), "Select UMAP:", choices = NULL),
    uiOutput(NS(id, "viz_type_ui")),
    conditionalPanel(
      condition = "input.viz_type == 'gene'",
      ns = NS(id),
      selectizeInput(NS(id, "gene_selection"), "Select gene:", choices = NULL)
    ),
    conditionalPanel(
      condition = "input.viz_type == 'cluster'",
      ns = NS(id),
      checkboxInput(NS(id, "show_legend"), "Show legend", value = FALSE)
    ),
    uiOutput(NS(id, "facet_ui")),
    helpText("Cells are limited to the checked splits above; with none checked, all are shown."),
    actionButton(NS(id, "plot_umap"), "Plot"),
    uiOutput(NS(id, "cond_dl"))
  )
}

umap_viz_output <- function(id) {
  tagList(
    plotOutput(NS(id, "umap_plot"), height = "auto")
  )
}

umap_viz_server <- function(id, ctx) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        ds <- ctx$ds()
        keys <- ds$reductions
        if (length(keys) == 0) {
          updateSelectInput(session, "umap_selection",
                            choices = list("No UMAPs available" = ""))
          return()
        }
        names(keys) <- reduction_label(keys)
        default <- if ("X_umap_full" %in% keys) "X_umap_full" else keys[[1]]
        updateSelectInput(session, "umap_selection", choices = keys,
                          selected = default)
      })

      label <- reactive(tools::toTitleCase(ctx$split_label()))

      output$viz_type_ui <- renderUI({
        radioButtons(
          NS(id, "viz_type"), "Color by:",
          choices = setNames(c("gene", "cluster", "split"),
                             c("Gene", "Cluster", label())),
          selected = isolate(if (is.null(input$viz_type)) "gene" else input$viz_type)
        )
      })

      output$facet_ui <- renderUI({
        checkboxInput(NS(id, "split_by_group"), paste("Split by", ctx$split_label()),
                      value = isolate(isTRUE(input$split_by_group)))
      })

      observe({
        updateSelectizeInput(session, "gene_selection",
                             choices = ctx$all_genes(), server = TRUE)
      })

      request <- eventReactive(
        input$plot_umap, {
          ok <- require_or_modal(
            nzchar(input$umap_selection),
            "No UMAP selected", "Please select a UMAP from the dropdown."
          )
          if (!ok) return(NULL)
          ok <- require_or_modal(
            input$viz_type != "gene" || nzchar(input$gene_selection),
            "No gene selected", "Please select a gene from the dropdown."
          )
          if (!ok) return(NULL)
          ds <- ctx$ds()
          splits <- ctx$required()
          if (length(splits) == 0) splits <- ds$splits
          list(ds = ds, key = input$umap_selection, color_by = input$viz_type,
               gene = input$gene_selection, splits = splits,
               facet = isTRUE(input$split_by_group),
               show_legend = isTRUE(input$show_legend))
        }
      )

      umap_plot_obj <- reactive({
        r <- request()
        shiny::req(r)
        with_spinner("Plotting...", umap_plot(
          r$ds, r$key, color_by = r$color_by, gene = r$gene, splits = r$splits,
          facet = r$facet, show_legend = r$show_legend
        ))
      })

      dims <- reactive({
        r <- request()
        if (is.null(r)) return(c(width = 1000, height = 600))
        umap_size(length(r$splits), r$facet)
      })

      observeEvent(
        input$plot_umap, {
          updateTabsetPanel(ctx$parent_session, inputId = "outputTabs",
                            selected = "UMAP")
          output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_umap"))})
        }
      )

      output$umap_plot <- renderPlot(
        umap_plot_obj(),
        height = function() dims()[["height"]]
      )

      ggdownload_server(
        "dl_umap",
        filename = reactive({
          r <- request()
          c("umap", sub("^X_", "", r$key), if (r$color_by == "gene") r$gene else r$color_by)
        }),
        payload = umap_plot_obj,
        height = function() dims()[["height"]] * 1.5,
        width = function() dims()[["width"]] * 1.5
      )
    }
  )
}
