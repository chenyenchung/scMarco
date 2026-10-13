ggdownload_UI <- function(id) {
  tagList(
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

ggdownload_server <- function(id, filename, payload, height, width) {
  #' height, width: pixels, as numbers or functions returning numbers
  moduleServer(
    id, function(input, output, session) {
      output$download <- downloadHandler(
        filename = function() {
          stem <- filename()
          if (length(stem) == 0) {
            stem <- "download"
          }
          stem <- paste(stem, collapse = "_")
          stem <- gsub("[[:punct:] ]+", "_", stem)
          paste(stem, input$filetype, sep = ".")
        },
        content = function(file) {
          ggsave(
            file,
            plot = payload(),
            height = value_of(height),
            width = value_of(width),
            limitsize = FALSE,
            units = "px",
            dpi = 150
          )
        }
      )
    }
  )
}
