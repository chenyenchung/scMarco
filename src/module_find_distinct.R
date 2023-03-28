find_distinct_UI <- function(id, stages) {
  tagList(
    sliderInput(
      NS(id, "pos_num"), "How many stages must the gene to be detected in?",
      min = 1, max = length(stages), step = 1, value = 1
    ),
    checkboxGroupInput(NS(id, "stages"), "", stages),
    actionButton(NS(id, "do_onplot"), "Find genes"),
    uiOutput(NS(id, "cond_dl"))
  )
}

find_distinct_output <- function(id) {
  tagList(
    actionButton(NS(id, "select"), "Select Genes"),
    plotOutput(NS(id, "onplot"), height = "auto")
  )
}

find_distinct_server <- function(
    id, stages, parent_session, db, db_tbl, prob_mat,
    # Reactive params
    sql_where, idents, cut_off, first_mark, coexp_session
    ) {
  moduleServer(
    id,
    function(input, output, session) {
      exp_genes <- reactive(
        {
          # Keep only genes that are expressed at provided stages
          # (If only selecting one: Any gene expressed at any stage)
          selected_stages <- stages() %in% input$stages

          exp_genes <- GetAllExpressedGenesDb(
            prob_mat = prob_mat(),
            stage = selected_stages,
            count = input$pos_num,
            stages = stages()
          )

          return(exp_genes)
        }
      )

      # # Filter expression table by clusters selected
      filmarkers <- eventReactive(
        input$do_onplot,
        {
          return_table <- FindDistinctiveGenes(
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
            db = db,
            db_tbl = db_tbl(),
            sql_where = sql_where(),
            marker_tbl = filmarkers(),
            idents = idents(),
            cut_off = cut_off(),
            stages = stages()
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

          # Show download button after plot is generated
          output$cond_dl <- renderUI({
            tagList(
              ggdownload_UI(NS(id, "dl_fdist"))
            )
          })
        }
      )


      ggdownload_server(
        "dl_fdist",
        filename = idents,
        payload = onplot,
        height = function() {
          num_gene <- length(unique(filmarkers()$distinct_gene))
          num_clust <- length(unique(filmarkers()$cluster))
          use_height <- 3 * max(
            120 * num_gene * max(1, round(num_clust / 3)),
            800
          )
          return(use_height)
        },
        width = min(length(idents()), 3) * 1000 + 500
      )

      observeEvent(
        input$select, {
          prev_mark <- first_mark()

          if (any(prev_mark == "")) {
            showModal(
              modalDialog(
                title = "Select gene to examine:",
                checkboxGroupInput(
                  NS(id, "to_plot"), "",
                  filmarkers()$distinct_gene
                ),
                actionButton(NS(id, "send"), "Send")
              )
            )
          } else {
            showModal(
              modalDialog(
                title = "Select gene pair to examine:",
                radioButtons(
                  NS(id, "to_plot"), "",
                  paste(
                    prev_mark,
                    filmarkers()$distinct_gene,
                    sep = " & "
                  )
                ),
                actionButton(NS(id, "send"), "Send")
              )
            )
          }

        }
      )

      observeEvent(
        input$send, {
          removeModal(session)
          updateSelectizeInput(
            coexp_session,
            inputId = "split_genes",
            selected = strsplit(input$to_plot, " & ")[[1]],
            choices = strsplit(input$to_plot, " & ")[[1]]
          )
          updateTabsetPanel(
            parent_session,
            inputId = "controlTabs",
            selected = "Co-expression plot"
          )

        }
      )

    }
  )
}
