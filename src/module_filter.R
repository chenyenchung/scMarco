lineplot_UI <- function(id) {
  tagList(
    selectizeInput(
      inputId = NS(id, 'gene_of_interest'),
      label = 'Which gene?',
      choices = NULL,
      multiple = FALSE
    ),
    sliderInput(
      NS(id, "line_alpha"),
      "How opaque other clusters should be?",
      min = 0.1, max = 1, step = 0.1, value = 0.2
    ),
    actionButton(
      NS(id, "do_lineplot"),
      "Plot log-normalized exp over stages"
    ),
    selectInput(
      NS(id, 'filetype'),
      "File type for download: ",
      choices = c("pdf", "png"),
      selected = "png",
      multiple = FALSE
    ),
    downloadButton(
      NS(id, 'download'),
      'Download'
    )

  )
}

lineplot_server <- function(
    id, norm, stages,
    # Reactive params
    db, idents
) {
  moduleServer(
    id,
    function(input, output, session) {
      observe(
        {
          # Update the gene selection panel for line plot
          updateSelectizeInput(
            session,
            inputId = "gene_of_interest",
            choices = unique(db()$gene),
            server = TRUE
          )
        }
      )

      lineplot <- eventReactive(
        input$do_lineplot,
        {
          LinePlot(
            norm,
            gene_to_plot = input$gene_of_interest,
            highlight_idents = idents(),
            stages = stages,
            lowlight_dim = input$line_alpha,
            highlight_col = "red"
          )
        }
      )


      output$download <- downloadHandler(
        filename = function(gene = input$gene_of_interest, ext = input$filetype) {
          paste0(gene, ".", ext)
        },
        content = function(file){
          ggsave(
            file,
            plot = LinePlot(
              norm,
              gene_to_plot = input$gene_of_interest,
              highlight_idents = idents(),
              stages = stages,
              lowlight_dim = input$line_alpha,
              highlight_col = "red",
              ggobj = TRUE
            )
          )
        }
      )

      # return(
      #   list(
      #     plot = lineplot,
      #     button = reactive(input$do_lineplot)
      #   )
      # )
    }
  )
}
