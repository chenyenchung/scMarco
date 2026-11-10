# The scMarco Shiny app: one codebase for every lab's instance. Datasets, gene
# lists, title, theme and README all come from the config file.

#' Read and check an scMarco app configuration
#'
#' Reads a YAML config and resolves every relative path in it (datasets, gene
#' lists, README files) against `root`.
#'
#' @param config Path of a YAML config file, or a list with the same fields.
#'   See the package README for the fields.
#' @param root Directory that relative paths in the config are relative to.
#'   Defaults to the working directory, so an app directory holding
#'   `config.yaml` and `data/` works as is.
#'
#' @return The config as a list, with absolute paths.
#' @examples
#' cfg <- list(datasets = list(example = list(
#'     label = "Example",
#'     clusters = system.file("extdata", "example.clusters.h5ad",
#'         package = "scMarco")
#' )))
#' str(read_scmarco_config(cfg)$datasets$example)
#' @export
read_scmarco_config <- function(config = "config.yaml", root = ".") {
    config_path <- if (is.character(config)) config else "the config"
    if (is.character(config)) config <- yaml::read_yaml(config)
    if (length(config$datasets) == 0) {
        stop("No datasets configured in ", config_path, call. = FALSE)
    }
    resolve <- function(p) {
        if (is.null(p)) return(NULL)
        p <- as.character(p)
        ifelse(grepl("^(/|~|[A-Za-z]:)", p), path.expand(p), file.path(root, p))
    }
    for (id in names(config$datasets)) {
        for (key in c("clusters", "cells")) {
            config$datasets[[id]][[key]] <- resolve(config$datasets[[id]][[key]])
        }
    }
    for (id in names(config$gene_list)) {
        config$gene_list[[id]]$path <- resolve(config$gene_list[[id]]$path)
    }
    config$readme <- resolve(unlist(config$readme))
    labels <- vapply(config$datasets, function(x) {
        if (is.null(x$label)) NA_character_ else as.character(x$label)
    }, character(1))
    if (anyNA(labels) || anyDuplicated(labels)) {
        stop("Every dataset needs a unique label in ", config_path, call. = FALSE)
    }
    config
}

#' The scMarco Shiny app
#'
#' `scmarco_app()` builds the app object (for `shiny::runApp()`, an `app.R`
#' file or deployment); `run_scmarco()` builds and runs it.
#'
#' Every configured dataset is checked at startup (see [ds_check_config()]),
#' the default dataset is loaded, and the others load on first use. Tabs
#' appear only when a dataset has the data they need.
#'
#' @inheritParams read_scmarco_config
#' @param ... Passed to [shiny::runApp()].
#'
#' @return `scmarco_app()`: a `shiny.appobj`. `run_scmarco()` does not
#'   return until the app stops.
#' @examples
#' cfg <- list(title = "scMarco", datasets = list(example = list(
#'     label = "Example",
#'     clusters = system.file("extdata", "example.clusters.h5ad",
#'         package = "scMarco")
#' )))
#' app <- scmarco_app(cfg)
#' if (interactive()) shiny::runApp(app)
#' @export
scmarco_app <- function(config = Sys.getenv("SCMARCO_CONFIG", "config.yaml"),
                        root = ".") {
    options(bitmapType = "cairo")
    options(ragg.max_dim = Inf)
    config <- read_scmarco_config(config, root)

    dataset_labels <- vapply(config$datasets, function(x) x$label, character(1))
    default_dataset <- names(config$datasets)[
        vapply(config$datasets, function(x) isTRUE(x$default), logical(1))
    ]
    default_dataset <- if (length(default_dataset) > 0) default_dataset[1] else
        names(config$datasets)[1]

    gene_list_labels <- vapply(
        config$gene_list, function(x) x$label, character(1)
    )

    # Datasets load on first use and are then shared across sessions (ds_load
    # memoizes per process). Check them all now so a bad config still fails at
    # startup, and warm the default so the landing page is ready.
    ds_check_config(config$datasets)
    message("Loading default dataset ", default_dataset, "...")
    invisible(ds_load(default_dataset, config$datasets[[default_dataset]]))

    # readme: the lab's Markdown files, then the packaged usage guide
    readme_paths <- config$readme
    if (!isFALSE(config$usage_guide)) {
        readme_paths <- c(readme_paths, system.file("app", "usage.md", package = "scMarco"))
    }

    ui <- fluidPage(
        theme = bslib::bs_theme(bootswatch = if (is.null(config$theme)) "minty" else config$theme),
        titlePanel(if (is.null(config$title)) "scMarco" else config$title),
        shinybusy::add_busy_bar(color = "#57068c", height = "12px"),
        sidebarLayout(
            sidebarPanel(
                selectInput(
                    inputId = "dataset",
                    label = "Which dataset to use?",
                    choices = stats::setNames(names(dataset_labels), dataset_labels),
                    selected = default_dataset
                ),
                uiOutput("dataset_versions"),
                selectInput(
                    "gene_group_to_use",
                    "Which genes to consider?",
                    c("All" = "All", stats::setNames(as.character(names(gene_list_labels)), gene_list_labels)),
                    selected = "All"
                ),
                uiOutput("meta_filter"),
                selectizeInput(
                    'clusters_of_interest',
                    'Cluster(s) of interest',
                    choices = NULL,
                    multiple = TRUE
                ),
                sliderInput(
                    "margin", "Confidence margin (%)",
                    min = 1, max = 50, step = 1, value = 10
                ),
                helpText(
                    "ON if P(On) \u2265 100% \u2212 margin, OFF if P(On) \u2264 margin;",
                    "in between is ambiguous (gray) and never counts as ON or OFF."
                ),
                uiOutput("splits_selector"),
                tabsetPanel(
                    id = "controlTabs",
                    type = "tabs",
                    tabPanel("From cluster", find_exp_UI("fexp")),
                    tabPanel("From gene", find_pos_UI("fpos")),
                    tabPanel("Find distinct genes within", find_distinct_UI("fd_genes")),
                    tabPanel("Co-expression plot", coexp_UI("coexp")),
                    tabPanel("Plot expression levels", lineplot_UI("lplot")),
                    tabPanel("UMAP", umap_viz_UI("umap")),
                    tabPanel("Mixture model", mixture_UI("mixture"))
                )
            ),
            mainPanel(
                tabsetPanel(
                    id = "outputTabs",
                    type = "tabs",
                    tabPanel(title = "Read Me", lapply(readme_paths, includeMarkdown)),
                    tabPanel(title = "Gene Table", find_exp_output("fexp")),
                    tabPanel(title = "Cluster Table", find_pos_output("fpos")),
                    tabPanel(title = "On/Off", find_distinct_output("fd_genes")),
                    tabPanel(title = "Co-Expression", coexp_output("coexp")),
                    tabPanel(title = "Expression Levels", lineplot_output("lplot")),
                    tabPanel(title = "UMAP", umap_viz_output("umap")),
                    tabPanel(title = "Mixture Model", mixture_output("mixture"))
                )
            )
        )
    )

    server <- function(input, output, session) {
        scmarco_server(input, output, session, config)
    }

    shinyApp(ui = ui, server = server)
}

#' @rdname scmarco_app
#' @export
run_scmarco <- function(config = Sys.getenv("SCMARCO_CONFIG", "config.yaml"),
                        root = ".", ...) {
    runApp(scmarco_app(config, root), ...)
}

scmarco_server <- function(input, output, session, config) {
    ds <- reactive({
        shiny::req(input$dataset)
        ds_load(input$dataset, config$datasets[[input$dataset]])
    })

    gene_filter <- reactive({
        if (input$gene_group_to_use == "All") return(NULL)
        readLines(config$gene_list[[input$gene_group_to_use]]$path)
    })

    # Genes used for combination finding (bimodal genes when known)
    genes <- reactive({
        g <- ds()$analysis_genes
        if (!is.null(gene_filter())) g <- intersect(g, gene_filter())
        g
    })

    # Genes for per-cell and trend views (every gene in the dataset)
    all_genes <- reactive({
        g <- ds()$all_genes
        if (!is.null(gene_filter())) g <- intersect(g, gene_filter())
        g
    })

    # Dropdown choices for views that also show unimodal genes, labeled as such
    all_gene_choices <- reactive(gene_choices(ds(), all_genes()))

    # Debounced: each distinct margin is a full binarize, so don't run one per
    # intermediate slider position
    margin <- debounce(reactive(input$margin / 100), 500)

    # Cluster metadata: the user picks which levels (config subset_by) to subset
    # by, and each picked level shows its own filter. Levels cascade coarse to
    # fine: each narrows the choices of the levels below it. Searches only see
    # clusters that match. An unpicked level is "all".
    output$dataset_versions <- renderUI({
        b <- ds()$binarization
        v <- b$version
        helpText(sprintf("Binarizer %s%s; schema %s", v,
            if (identical(v, BINARIZER_EXTERNAL)) " (tables not from scMarco's binarizer)" else "",
            SCMARCO_FORMAT_VERSION))
    })

    output$meta_filter <- renderUI({
        meta <- ds()$cluster_meta
        if (is.null(meta)) return(NULL)
        tagList(
            selectizeInput("meta_levels", "Subset clusters by", choices = names(meta),
                multiple = TRUE, options = list(placeholder = "none (all clusters)")),
            lapply(seq_along(meta), function(i) {
                conditionalPanel(
                    sprintf("(input.meta_levels || []).indexOf(%s) > -1",
                        encodeString(names(meta)[i], quote = '"')),
                    selectizeInput(
                        paste0("meta_", i), paste0("Search within: ", names(meta)[i]),
                        choices = sort(unique(meta[[i]][!is.na(meta[[i]])])), multiple = TRUE,
                        options = list(placeholder = "all")
                    )
                )
            }),
            helpText(textOutput("meta_scope", inline = TRUE))
        )
    })

    meta_filters <- reactive({
        meta <- ds()$cluster_meta
        if (is.null(meta)) return(list())
        active <- input$meta_levels
        # Values left over from another dataset are ignored, not matched
        stats::setNames(lapply(seq_along(meta), function(i) {
            if (!names(meta)[i] %in% active) return(character(0))
            as.character(intersect(input[[paste0("meta_", i)]], meta[[i]]))
        }), names(meta))
    })

    # A level that is no longer picked goes back to "all"
    observeEvent(input$meta_levels, ignoreNULL = FALSE, {
        meta <- ds()$cluster_meta
        for (i in seq_along(meta)) {
            key <- paste0("meta_", i)
            if (!names(meta)[i] %in% input$meta_levels && length(input[[key]]) > 0) {
                updateSelectizeInput(session, key, selected = character(0))
            }
        }
    })

    universe <- reactive(meta_universe(ds(), meta_filters()))

    # Last choices sent per level, so a level is only updated when they change
    # (updating closes an open dropdown)
    meta_sent <- new.env(parent = emptyenv())
    observeEvent(ds(), rm(list = ls(meta_sent), envir = meta_sent))
    observeEvent(meta_filters(), {
        filters <- meta_filters()
        for (i in seq_along(filters)[-1]) {
            level <- names(filters)[i]
            choices <- meta_level_choices(ds(), filters, level)
            kept <- as.character(intersect(filters[[i]], choices))
            key <- paste0("meta_", i)
            if (!identical(meta_sent[[key]], choices) || !identical(kept, filters[[i]])) {
                updateSelectizeInput(session, key, choices = choices, selected = kept)
                meta_sent[[key]] <- choices
            }
            filters[i] <- list(kept)
        }
    })

    output$meta_scope <- renderText({
        sprintf("%d of %d clusters searched.", length(universe()), length(ds()$clusters))
    })

    b <- reactive(restrict_clusters(ds_binarize(ds(), margin()), universe()))

    # The cluster picker lists the clusters in scope; a new dataset starts empty
    shown_ds <- NULL
    observeEvent(universe(), {
        keep <- if (identical(shown_ds, ds()$id)) {
            intersect(isolate(input$clusters_of_interest), universe())
        } else {
            character(0)
        }
        shown_ds <<- ds()$id
        updateSelectizeInput(session, "clusters_of_interest",
            choices = cluster_choices(ds(), universe()), selected = keep)
    })

    observeEvent(ds(), {
        first_mark(character(0))

        has_trend <- !is.null(ds()$trend)
        has_umap <- length(ds()$reductions) > 0
        has_mixture <- ds()$has_mixture
        for (tab in list(
            list("controlTabs", "Plot expression levels", has_trend),
            list("outputTabs", "Expression Levels", has_trend),
            list("controlTabs", "UMAP", has_umap),
            list("outputTabs", "UMAP", has_umap),
            list("controlTabs", "Mixture model", has_mixture),
            list("outputTabs", "Mixture Model", has_mixture)
        )) {
            if (tab[[3]]) showTab(tab[[1]], tab[[2]]) else hideTab(tab[[1]], tab[[2]])
        }
    })

    output$splits_selector <- renderUI({
        splits <- ds()$splits
        plural <- ds()$split_plural
        tagList(
            checkboxGroupInput(
                "required_splits",
                paste0("Must be ON in these ", plural),
                choices = splits,
                inline = TRUE
            ),
            sliderInput(
                "min_n",
                paste0("...and ON in at least this many ", plural),
                min = 1, max = length(splits), step = 1, value = 1
            )
        )
    })

    required <- reactive(intersect(input$required_splits, ds()$splits))

    min_n <- reactive({
        n <- if (is.null(input$min_n)) 1 else input$min_n
        # Requiring k splits implies at least k splits
        min(max(n, length(required()), 1), length(ds()$splits))
    })

    # Show the count actually used when the required splits raise it
    observeEvent(required(), {
        if (!is.null(input$min_n) && input$min_n < length(required())) {
            updateSliderInput(session, "min_n", value = length(required()))
        }
    })

    first_mark <- reactiveVal(character(0))

    ctx <- list(
        ds = ds,
        b = b,
        genes = genes,
        all_genes = all_genes,
        all_gene_choices = all_gene_choices,
        idents = reactive(input$clusters_of_interest),
        required = required,
        min_n = min_n,
        margin = margin,
        split_label = reactive(ds()$split_label),
        split_plural = reactive(ds()$split_plural),
        scope = reactive(describe_scope(ds(), meta_filters())),
        first_mark = first_mark,
        parent_session = session
    )

    find_exp_server("fexp", ctx)
    find_pos_server("fpos", ctx)
    coexp_session <- coexp_server("coexp", ctx)
    find_distinct_server("fd_genes", ctx, coexp_session)
    lineplot_server("lplot", ctx)
    umap_viz_server("umap", ctx)
    mixture_server("mixture", ctx)
}
