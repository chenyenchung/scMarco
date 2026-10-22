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
        keys <- reduction_choices(ds)
        if (length(keys) == 0) {
          updateSelectInput(session, "umap_selection",
                            choices = list("No UMAPs available" = ""))
          return()
        }
        default <- if ("X_umap_full" %in% keys) "X_umap_full" else keys[[1]]
        updateSelectInput(session, "umap_selection", choices = keys,
                          selected = default)
      })

      label <- reactive(tools::toTitleCase(ctx$split_label()))

      # c(<split> = key, ...) when each split has its own reduction, else NULL
      split_keys <- reactive({
        choice <- input$umap_selection
        if (is.null(choice) || !startsWith(choice, "split:")) return(NULL)
        resolve_reduction(ctx$ds(), choice)
      })

      output$viz_type_ui <- renderUI({
        # Coloring by split is redundant when every split is its own panel
        choices <- if (is.null(split_keys())) c("gene", "cluster", "split") else c("gene", "cluster")
        selected <- isolate(input$viz_type)
        if (is.null(selected) || !selected %in% choices) selected <- "gene"
        radioButtons(
          NS(id, "viz_type"), "Color by:",
          choices = setNames(choices, c("Gene", "Cluster", label())[seq_along(choices)]),
          selected = selected
        )
      })

      ncol_slider <- function(n) {
        if (n < 2) return(NULL)
        value <- isolate(input$facet_ncol)
        if (is.null(value)) value <- 2
        sliderInput(NS(id, "facet_ncol"), "How many panels per row?",
                    min = 1, max = n, step = 1, value = min(value, n))
      }

      output$facet_ui <- renderUI({
        keys <- split_keys()
        if (is.null(keys)) {
          return(tagList(
            checkboxInput(NS(id, "split_by_group"), paste("Split by", ctx$split_label()),
                          value = isolate(isTRUE(input$split_by_group))),
            conditionalPanel(
              condition = "input.split_by_group",
              ns = NS(id),
              ncol_slider(length(ctx$ds()$splits))
            ),
            helpText("Cells are limited to the checked splits above; with none checked, all are shown.")
          ))
        }
        # No shared reduction: each chosen split is drawn as its own panel
        tagList(
          selectizeInput(
            NS(id, "shown_splits"), paste0(label(), "s to show:"),
            choices = names(keys), selected = names(keys), multiple = TRUE
          ),
          ncol_slider(length(keys))
        )
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
          keys <- split_keys()
          if (is.null(keys)) {
            key <- input$umap_selection
            name <- sub("^X_", "", key)
            splits <- ctx$required()
            if (length(splits) == 0) splits <- ds$splits
            facet <- isTRUE(input$split_by_group)
          } else {
            splits <- intersect(names(keys), input$shown_splits)
            ok <- require_or_modal(
              length(splits) > 0, paste("No", ctx$split_label(), "selected"),
              paste0("Please choose at least one ", ctx$split_label(), " to show.")
            )
            if (!ok) return(NULL)
            key <- keys[splits]
            name <- paste0(sub("^split:", "", input$umap_selection), "_by_", ds$split_label)
            facet <- TRUE
          }
          ncol <- input$facet_ncol
          if (is.null(ncol)) ncol <- 2
          list(ds = ds, key = key, name = name, color_by = input$viz_type,
               gene = input$gene_selection, splits = splits, facet = facet,
               show_legend = isTRUE(input$show_legend),
               ncol = min(ncol, length(splits)))
        }
      )

      umap_plot_obj <- reactive({
        r <- request()
        shiny::req(r)
        with_spinner("Plotting...", umap_plot(
          r$ds, r$key, color_by = r$color_by, gene = r$gene, splits = r$splits,
          facet = r$facet, show_legend = r$show_legend, ncol = r$ncol
        ))
      })

      dims <- reactive({
        r <- request()
        if (is.null(r)) return(c(width = 1000, height = 600))
        umap_size(length(r$splits), r$facet, r$ncol)
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
          c("umap", r$name, if (r$color_by == "gene") r$gene else r$color_by)
        }),
        payload = umap_plot_obj,
        height = function() dims()[["height"]] * 1.5,
        width = function() dims()[["width"]] * 1.5
      )
    }
  )
}
