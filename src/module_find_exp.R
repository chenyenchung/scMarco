find_exp_UI <- function(id, stages) {
  tagList(
    sliderInput(
      NS(id, "pos_num"), "How many stages must the gene to be detected in?",
      min = 1, max = length(stages), step = 1, value = 1
    ),
    sliderInput(
      NS(id, "output_num"),
      "How many genes to show?",
      min = 1, max = 50, step = 1, value = 10
    ),
    checkboxGroupInput(
      NS(id, "stages"),
      "",
      stages
    ),
    actionButton(NS(id, "do_fexp"), "Show genes"),
  )
}

find_exp_output <- function(id) {
  tagList(
    uiOutput(NS(id, "selectTbl"))
  )
}

find_exp_server <- function(
    id, stages, db, db_tbl, prob_mat,
    # Reactive params
    sql_where, idents, cut_off, parent_session
) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(
        input$do_fexp, {
          # Keep only genes that are expressed at provided stages
          # (If only selecting one: Any gene expressed at any stage)
          selected_stages <- stages() %in% input$stages

          exp_genes <- GetAllExpressedGenesForThis(
            prob_mat = prob_mat(),
            stage = selected_stages,
            count = input$pos_num,
            stages = stages(),
            ident = idents(),
            output_limit = input$output_num
          )

          updateTabsetPanel(
            parent_session,
            inputId = "outputTabs",
            selected = "Gene Table"
          )

          # Present genes that are expressed in the cluster of interest
          # as a radio button menu
          gene_to_select <- names(exp_genes)

          names(gene_to_select) <- paste0(
            gene_to_select,
            " (",
            exp_genes,
            ")"
          )

          output$selectTbl <- renderUI({
            tagList(
              tags$h4(
                style = "margin: 10px",
                radioButtons(
                  NS(id, "gc"), "Gene (# of off-targets)",
                  gene_to_select
                )
              ),
              actionButton(NS(id, "pass"), "Select gene")
            )
          })
        }
      )
      other_clust <- eventReactive(
        input$pass, {
          others <- dbGetQuery(
            db,
            paste(
              "SELECT DISTINCT cluster FROM",
              db_tbl(), "WHERE gene =",
              paste0("'", input$gc, "'"),
              "AND type = 'prob'",
              "AND value >", cut_off()
            )
          )
          return(others$cluster)
        }
      )
      observeEvent(
        input$pass, {
          showModal(
            modalDialog(
              title = NULL,
              checkboxGroupInput(
                NS(id, "clust"),
                "Which clusters to focus on?",
                other_clust(),
                selected = other_clust(),
                inline = TRUE
              ),
              actionButton(NS(id, "submit"), "Send"),
              actionButton(NS(id, "get_them_all"), "Select All"),
              actionButton(NS(id, "not_at_all"), "Clear All")
            )
          )
        }
      )
      observeEvent(
        input$submit, {
          removeModal(session)
          updateSelectizeInput(
            parent_session,
            inputId = "clusters_of_interest",
            selected = input$clust
          )
          updateTabsetPanel(
            parent_session,
            inputId = "controlTabs",
            selected = "Find distinct genes within"
          )
          updateSelectInput(
            parent_session,
            inputId = "cache_last_g",
            choices = input$gc,
            selected = input$gc
          )
        }
      )
      observeEvent(
        input$get_them_all, {
          updateCheckboxGroupInput(
            session,
            "clust",
            selected = other_clust()
          )
        }
      )
      observeEvent(
        input$not_at_all, {
          updateCheckboxGroupInput(
            session,
            "clust",
            selected = FALSE
          )
        }
      )
      return(reactive({input$gc}))
    }
  )
}
