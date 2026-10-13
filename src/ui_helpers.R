# Small Shiny helpers shared by the modules

require_or_modal <- function(ok, title, message) {
  #' Show an explanatory modal and return FALSE when `ok` is not TRUE
  if (isTRUE(ok)) return(TRUE)
  showModal(modalDialog(
    title = title, message, easyClose = TRUE, footer = modalButton("OK")
  ))
  FALSE
}

with_spinner <- function(text, expr) {
  shinybusy::show_modal_spinner(spin = "fading-circle", text = text)
  on.exit(shinybusy::remove_modal_spinner())
  tryCatch(expr, error = function(e) {
    shinybusy::remove_modal_spinner()
    showModal(modalDialog(
      title = "Something went wrong", conditionMessage(e),
      easyClose = TRUE, footer = modalButton("OK")
    ))
    NULL
  })
}

describe_splits <- function(required, min_n, split_label) {
  req <- if (length(required) == 0) {
    ""
  } else {
    paste0(" ON in ", paste(required, collapse = " + "), " and")
  }
  sprintf("Criteria:%s ON in at least %d %s(s).", req, min_n, split_label)
}

value_of <- function(x) {
  if (is.function(x)) x() else x
}
