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

describe_splits <- function(required, min_n, split_label, scope = NULL) {
  req <- if (length(required) == 0) {
    ""
  } else {
    paste0(" ON in ", paste(required, collapse = " + "), " and")
  }
  paste(c(sprintf("Criteria:%s ON in at least %d %s(s).", req, min_n, split_label),
          scope), collapse = " ")
}

describe_scope <- function(ds, filters) {
  #' "Searching within class = neuron > subclass = Mi, Tm (42 of 210
  #' clusters)." or NULL when the search covers every cluster
  active <- meta_active(filters)
  if (length(active) == 0) return(NULL)
  what <- paste(names(active), "=", vapply(active, paste, "", collapse = ", "),
                collapse = " > ")
  sprintf("Searching within %s (%d of %d clusters).", what,
          length(meta_universe(ds, filters)), length(ds$clusters))
}

value_of <- function(x) {
  if (is.function(x)) x() else x
}
