find_pos_UI <- function(id) {
  tagList(
    selectizeInput(
      NS(id, "gene"),
      label = 'Which gene(s)?',
      choices = NULL,
      multiple = TRUE
    ),
    radioButtons(
      NS(id, "mode"),
      "With several genes, a cluster must express",
      choices = c("all of them" = "all", "any of them" = "any"),
      selected = "all",
      inline = TRUE
    ),
    actionButton(NS(id, "do_fpos"), "Show clusters")
  )
}

find_pos_output <- function(id) {
  tagList(
    uiOutput(NS(id, "selectTbl"))
  )
}

find_pos_server <- function(id, ctx) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        updateSelectizeInput(
          session, 'gene', choices = ctx$genes(), server = TRUE
        )
      })

      found <- reactiveVal(NULL)

      observeEvent(
        input$do_fpos, {
          ok <- require_or_modal(
            length(input$gene) > 0,
            "Please select at least one gene",
            "To find clusters expressing genes, select one or more genes from the dropdown."
          )
          if (!ok) return()

          clusters <- clusters_expressing(
            ctx$b(), input$gene, ctx$required(), ctx$min_n(), mode = input$mode
          )
          ok <- require_or_modal(
            length(clusters) > 0,
            "No clusters found",
            paste0(
              "No cluster expresses ", paste(input$gene, collapse = ", "),
              " with probability > ", ctx$cut(), ". ",
              describe_splits(ctx$required(), ctx$min_n(), ctx$split_label())
            )
          )
          if (!ok) return()
          found(list(genes = input$gene, mode = input$mode, clusters = clusters,
                     ds = ctx$ds()))

          updateTabsetPanel(
            ctx$parent_session,
            inputId = "outputTabs",
            selected = "Cluster Table"
          )
        }
      )

      output$selectTbl <- renderUI({
        res <- found()
        if (is.null(res)) return(NULL)
        tagList(
          tags$p(style = "margin: 10px",
                 describe_splits(ctx$required(), ctx$min_n(), ctx$split_label())),
          tags$h4(
            style = "margin: 10px",
            checkboxGroupInput(
              NS(id, "select_clust"), "Clusters",
              cluster_choices(res$ds, res$clusters), selected = res$clusters
            )
          ),
          actionButton(NS(id, "pass"), "Select clusters")
        )
      })

      observeEvent(
        input$pass, {
          res <- found()
          updateSelectizeInput(
            ctx$parent_session,
            inputId = "clusters_of_interest",
            selected = input$select_clust
          )
          updateTabsetPanel(
            ctx$parent_session,
            inputId = "controlTabs",
            selected = "Find distinct genes within"
          )
          # Only an "all genes" combination is a marker to extend; "any" is
          # an exploratory union
          ctx$first_mark(if (res$mode == "all") res$genes else character(0))
        }
      )
    }
  )
}
