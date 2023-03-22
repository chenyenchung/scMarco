ggdownload_UI <- function(id) {
  uiOutput(NS(id, "dl_interface"))
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
  moduleServer(
    id, function(input, output, session) {
      output$download <- downloadHandler(
        filename = function(gene = filename(), ext = input$filetype) {
          if (length(gene) == 0) {
            gene <- "download"
          }
          if (length(gene) > 1) {
            gene <- paste(gene, collapse = "_")
          }
          gene <- gsub("[[:punct:] ]+", "_", gene)
          return(paste(gene, ext, sep = "."))
        },
        content = function(file){
          ggsave(
            file,
            plot = payload(),
            height = height(),
            width = width,
            limitsize = FALSE,
            units = "px"
          )
        }
      )
    }
  )
}
