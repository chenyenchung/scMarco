find_pos_UI <- function(id, stages) {
  tagList(
    selectizeInput(
      NS(id, "gene"),
      label = 'Which gene?',
      choices = NULL,
      multiple = TRUE
    ),
    sliderInput(
      NS(id, "pos_num"), "How many stages must the gene to be detected in?",
      min = 1, max = length(stages), step = 1, value = 1
    ),
    checkboxGroupInput(
      NS(id, "stages"),
      "",
      stages
    ),
    actionButton(NS(id, "do_fpos"), "Show clusters"),
  )
}

find_pos_output <- function(id) {
  tagList(
    uiOutput(NS(id, "selectTbl"))
  )
}

find_pos_server <- function(
    id, stages, db, db_tbl,
    # Reactive params
    sql_where, idents, cut_off, parent_session, genes
) {
  moduleServer(
    id,
    function(input, output, session) {

      observe(
        {
          # Generate SQL statement based on global filter condition
          sql_where <- sql_where()
          input_q <- ifelse(
            sql_where == "",
            paste("SELECT DISTINCT gene FROM", db_tbl()),
            paste("SELECT DISTINCT gene FROM", db_tbl(), wherecat(sql_where))
          )

          updateSelectizeInput(
            session,
            'gene',
            choices = genes(),
            selected = NULL,
            server = TRUE
          )
        }
      )
      observeEvent(
        input$pos_num, {
          # This is a hack. Without this, the above observe event is triggered
          # before the input object is generated, and thus the options won't
          # be updated at startup.
          updateSelectizeInput(
            session,
            'gene',
            choices = genes(),
            selected = NULL,
            server = TRUE
          )
        }
      )

      observeEvent(
        input$do_fpos, {
          current_stages <- stages()
          selected_stages <- current_stages %in% input$stages
          stage_num <- max(
            sum(selected_stages), input$pos_num - sum(selected_stages)
          )

          clusters <- dbGetQuery(
          db,
          paste(
            "SELECT DISTINCT cluster FROM",
            db_tbl(), "WHERE gene =",
            paste0("'", input$gene, "'"),
            "AND type = 'prob'",
            "AND value >", cut_off(),
            "GROUP BY cluster",
            "HAVING",
            "COUNT(DISTINCT CASE WHEN stage IN", vec2sqllist(current_stages[selected_stages]),
            "THEN stage END) =", sum(selected_stages), "AND",
            "COUNT(DISTINCT CASE WHEN stage IN", vec2sqllist(current_stages[!selected_stages]),
            "THEN stage END) >=", stage_num
          )
        )

        # Switch to cluster selection table
        updateTabsetPanel(
          parent_session,
          inputId = "outputTabs",
          selected = "Cluster Table"
        )

        # Show clusters for selection
        output$selectTbl <- renderUI({
          tagList(
            tags$h4(
              style = "margin: 10px",
              checkboxGroupInput(
                NS(id, "select_clust"), "Clusters",
                clusters$cluster
              )
            ),
            actionButton(NS(id, "pass"), "Select cluster")
          )
        })

        }
      )
      observeEvent(
        input$pass, {
          updateSelectizeInput(
            parent_session,
            inputId = "clusters_of_interest",
            selected = input$select_clust
          )
          updateTabsetPanel(
            parent_session,
            inputId = "controlTabs",
            selected = "Find distinct genes within"
          )
          updateSelectInput(
            parent_session,
            inputId = "cache_last_g",
            choices = input$gene,
            selected = input$gene
          )
        }
      )
    }
  )
}
