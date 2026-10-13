find_distinct_UI <- function(id) {
  tagList(
    numericInput(
      NS(id, "k"),
      "A gene must be ON in exactly this many of the selected clusters",
      value = 1, min = 1, step = 1
    ),
    numericInput(
      NS(id, "n_show"),
      "Genes to plot per cluster (fewest outside off-targets first)",
      value = 10, min = 1, max = 50, step = 1
    ),
    uiOutput(NS(id, "pct_toggle")),
    actionButton(NS(id, "do_onplot"), "Find genes"),
    uiOutput(NS(id, "cond_dl"))
  )
}

find_distinct_output <- function(id) {
  tagList(
    actionButton(NS(id, "select"), "Select Genes"),
    plotOutput(NS(id, "onplot"), height = "auto"),
    tableOutput(NS(id, "hits"))
  )
}

find_distinct_server <- function(id, ctx, coexp_session) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pct_toggle <- renderUI({
        if (!ctx$ds()$has_pct) return(NULL)
        checkboxInput(NS(id, "show_pct"), "Show % cells having >0 UMIs", value = FALSE)
      })
      show_pct <- reactive(isTRUE(input$show_pct) && ctx$ds()$has_pct)

      found <- eventReactive(
        input$do_onplot, {
          idents <- ctx$idents()
          ok <- require_or_modal(
            length(idents) >= 2,
            "Please select multiple clusters",
            "To find distinguishing genes, select at least 2 clusters in 'Cluster(s) of interest'."
          )
          if (!ok) return(NULL)
          k <- as.integer(input$k)
          ok <- require_or_modal(
            !is.na(k) && k >= 1 && k < length(idents),
            "Invalid number of clusters",
            sprintf("The number must be between 1 and %d.", length(idents) - 1)
          )
          if (!ok) return(NULL)

          b <- ctx$b()
          tbl <- with_spinner("Finding distinct genes...", distinct_within(
            b, idents, ctx$required(), ctx$min_n(), k = k, genes = ctx$genes()
          ))
          if (is.null(tbl)) return(NULL)
          ok <- require_or_modal(
            nrow(tbl) > 0,
            "No distinctive genes found",
            HTML(paste0(
              "<p>No gene distinguishes the selected clusters. ",
              describe_splits(ctx$required(), ctx$min_n(), ctx$split_label()),
              "</p><p>Try another threshold, fewer required ",
              ctx$split_label(), "s, different clusters, or a broader gene list.</p>"
            ))
          )
          if (!ok) return(NULL)
          n_show <- max(1, as.integer(input$n_show), na.rm = TRUE)
          list(tbl = tbl, shown = top_per_cluster(tbl, n_show),
               idents = idents, b = b, ds = ctx$ds(), first = ctx$first_mark())
        }
      )

      onplot <- reactive({
        res <- found()
        if (is.null(res)) {
          return(message_plot("Select clusters and click 'Find genes'."))
        }
        onoff_plot(res$ds, res$b, unique(res$shown$gene), res$idents,
                   show_pct = show_pct())
      })

      plot_dims <- reactive({
        res <- found()
        if (is.null(res)) return(c(width = 800, height = 400))
        onoff_size(length(unique(res$shown$gene)), length(res$idents), show_pct())
      })

      output$onplot <- renderPlot(
        onplot(),
        height = function() plot_dims()[["height"]]
      )

      output$hits <- renderTable({
        # Every hit, not only the plotted ones. With a first marker, also check
        # the pair across the whole dataset: distinct_within() only guarantees
        # uniqueness among the selected clusters
        res <- found()
        if (is.null(res)) return(NULL)
        tbl <- data.frame(
          Cluster = cluster_label(res$ds, res$tbl$cluster),
          Gene = res$tbl$gene,
          `Clusters outside the selection also ON` = res$tbl$offtargets,
          check.names = FALSE
        )
        if (length(res$first) > 0) {
          pair_col <- sprintf("Other clusters with %s & gene both ON",
                              paste(res$first, collapse = " & "))
          tbl[[pair_col]] <- vapply(seq_len(nrow(res$tbl)), function(i) {
            combo <- unique(c(res$first, res$tbl$gene[i]))
            others <- combo_offtarget_clusters(res$b, combo, exclude = res$tbl$cluster[i])
            if (length(others) == 0) "none" else
              paste(cluster_label(res$ds, others), collapse = ", ")
          }, character(1))
        }
        tbl
      }, striped = TRUE)

      observeEvent(
        input$do_onplot, {
          updateTabsetPanel(
            ctx$parent_session,
            inputId = "outputTabs",
            selected = "On/Off"
          )
          output$cond_dl <- renderUI({ggdownload_UI(NS(id, "dl_fdist"))})
        }
      )

      ggdownload_server(
        "dl_fdist",
        filename = reactive(cluster_label(ctx$ds(), ctx$idents())),
        payload = onplot,
        height = function() plot_dims()[["height"]] * 2,
        width = function() plot_dims()[["width"]] * 2
      )

      observeEvent(
        input$select, {
          res <- found()
          if (is.null(res)) return()
          hits <- unique(res$shown$gene)
          if (length(res$first) == 0) {
            showModal(modalDialog(
              title = "Select genes to examine:",
              checkboxGroupInput(NS(id, "to_plot"), "", hits),
              actionButton(NS(id, "send"), "Send")
            ))
          } else {
            first <- paste(res$first, collapse = " & ")
            showModal(modalDialog(
              title = "Select a combination to examine:",
              radioButtons(NS(id, "to_plot"), "", paste(first, hits, sep = " & ")),
              actionButton(NS(id, "send"), "Send")
            ))
          }
        }
      )

      observeEvent(
        input$send, {
          removeModal(session)
          genes <- unique(unlist(strsplit(input$to_plot, " & ", fixed = TRUE)))
          updateSelectizeInput(
            coexp_session,
            inputId = "split_genes",
            choices = ctx$genes(),
            selected = genes,
            server = TRUE
          )
          updateTabsetPanel(
            ctx$parent_session,
            inputId = "controlTabs",
            selected = "Co-expression plot"
          )
        }
      )
    }
  )
}
