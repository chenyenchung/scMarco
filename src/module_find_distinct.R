find_distinct_UI <- function(id, stages) {
  tagList(
    sliderInput(
      NS(id, "pos_num"), "How many stages must the gene to be detected in?",
      min = 1, max = length(stages), step = 1, value = 1
    ),
    lapply(stages, function(x) checkboxInput(NS(id, x), x, FALSE)),
    actionButton(NS(id, "do_onplot"), "Find genes"),
    uiOutput(NS(id, "cond_dl"))
  )
}

find_distinct_output <- function(id) {
  tagList(
    plotOutput(NS(id, "onplot"), height = "auto")
  )
}

find_distinct_server <- function(
    id, stages, parent_session,
    # Reactive params
    db, idents, cut_off
    ) {
  moduleServer(
    id,
    function(input, output, session) {
      exp_genes <- reactive(
        {
          # Keep only genes that are expressed at provided stages
          # (If only selecting one: Any gene expressed at any stage)
          selected_stages <- vapply(
            stages,
            function(x) {input[[x]]},
            FUN.VALUE = logical(1)
          )

          exp_genes <- GetAllExpressedGenes(
            db(),
            cut_off = cut_off(),
            stage = selected_stages,
            count = input$pos_num,
            stages = stages
          )

          return(exp_genes)
        }
      )

      # # Filter expression table by clusters selected
      filmarkers <- eventReactive(
        input$do_onplot,
        {
          return_table <- FindDistinctiveGenesPerGroup(
            exp_genes(),
            idents = idents(),
            max_presence = 1
          )
          return(return_table)
        }
      )

      # Filter expression table by clusters selected
      # Plot On/Off plot for the markers
      onplot <- eventReactive(
        input$do_onplot,
        {
          OnPlot(
            x = db(),
            marker_tbl = filmarkers(),
            idents = idents(),
            cut_off = cut_off(),
            stages = stages
          )
        }
      )

      output$onplot <- renderPlot(
        onplot(),
        height = function() {
          num_gene <- length(unique(filmarkers()$distinct_gene))
          num_clust <- length(unique(filmarkers()$cluster))
          use_height <- max(
            80 * num_gene * max(1, round(num_clust/3)),
            800
          )
          return(use_height)
        }
      )

      # Switch to the on/off plot panel when asked to plot
      observeEvent(
        input$do_onplot,
        {
          updateTabsetPanel(
            parent_session,
            inputId = "outputTabs",
            selected = "On/Off"
          )
        }
      )

      # Show download button after plot is generated
      observeEvent(input$do_onplot, {
        output$cond_dl <- renderUI({
          tagList(
            ggdownload_UI(NS(id, "dl_fdist")),
            actionButton(NS(id, "pass"), "Select genes")
          )


        })
      })


      ggdownload_server(
        "dl_fdist",
        filename = idents,
        payload = onplot,
        height = function() {
          num_gene <- length(unique(filmarkers()$distinct_gene))
          num_clust <- length(unique(filmarkers()$cluster))
          use_height <- 3 * max(
            100 * num_gene * max(1, round(num_clust / 3)),
            800
          )
          return(use_height)
        },
        width = min(length(idents()), 3) * 1000 + 500
      )
    }
  )
}
