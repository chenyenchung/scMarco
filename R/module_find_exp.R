find_exp_UI <- function(id) {
    tagList(
        sliderInput(
            NS(id, "output_num"),
            "How many genes to show?",
            min = 1, max = 50, step = 1, value = 10
        ),
        checkboxInput(
            NS(id, "show_pairs"),
            "Also rank gene pairs by shared off-targets",
            value = FALSE
        ),
        actionButton(NS(id, "do_fexp"), "Show genes")
    )
}

find_exp_output <- function(id) {
    tagList(
        uiOutput(NS(id, "selectTbl")),
        uiOutput(NS(id, "pairs"))
    )
}

heterogeneity_note <- function(b, ds, ident, genes = b$genes, n_show = 10) {
    # A warning that `ident` may be a mix of cell types: the genes whose calls
    # diverge between its subclusters (ambiguous there), or NULL for none
    het <- b$heterogeneous[intersect(genes, b$genes), ident, , drop = FALSE]
    hit <- rowSums(het, dims = 1) > 0
    if (!any(hit)) return(NULL)
    flagged <- names(hit)[hit]
    shown <- paste(head(flagged, n_show), collapse = ", ")
    if (length(flagged) > n_show) shown <- paste0(shown, ", ...")
    tags$p(
        style = "margin: 10px", class = "text-warning",
        tags$b(cluster_label(ds, ident)), " may be heterogeneous: ",
        length(flagged), if (length(flagged) == 1) " gene's calls diverge" else " genes' calls diverge",
        " between its subclusters (", shown, "). They count as ambiguous in this cluster."
    )
}

find_exp_server <- function(id, ctx) {
    # ctx: list of shared reactives built in app.R
    moduleServer(
        id,
        function(input, output, session) {
            result <- reactiveVal(NULL)

            observeEvent(
                input$do_fexp, {
                    idents <- ctx$idents()
                    ok <- require_or_modal(
                        length(idents) == 1,
                        "Please select exactly one cluster",
                        "To find genes expressed in a cluster, select exactly one cluster in 'Cluster(s) of interest'."
                    )
                    if (!ok) return()

                    res <- with_spinner("Finding expressed genes...", {
                        b <- ctx$b()
                        genes <- from_cluster(
                            b, idents, ctx$required(), ctx$min_n(),
                            limit = input$output_num, genes = ctx$genes()
                        )
                        pairs <- if (isTRUE(input$show_pairs)) {
                            pair_offtargets(b, idents, ctx$required(), ctx$min_n(),
                                limit = input$output_num, genes = ctx$genes())
                        }
                        list(ident = idents, genes = genes, pairs = pairs, ds = ctx$ds(),
                            het_note = heterogeneity_note(b, ctx$ds(), idents, ctx$genes()))
                    })
                    if (is.null(res)) return()
                    result(res)

                    updateTabsetPanel(
                        ctx$parent_session,
                        inputId = "outputTabs",
                        selected = "Gene Table"
                    )
                }
            )

            output$selectTbl <- renderUI({
                res <- result()
                if (is.null(res)) return(NULL)
                if (length(res$genes) == 0) {
                    return(tagList(res$het_note, tags$p(
                        style = "margin: 10px",
                        "No gene passes the criteria in ", tags$b(cluster_label(res$ds, res$ident)), ". ",
                        describe_splits(ctx$required(), ctx$min_n(), ctx$split_label(), ctx$split_plural(), ctx$scope())
                    )))
                }
                gene_to_select <- names(res$genes)
                names(gene_to_select) <- paste0(gene_to_select, " (", res$genes, ")")
                tagList(
                    res$het_note,
                    tags$p(style = "margin: 10px",
                        describe_splits(ctx$required(), ctx$min_n(), ctx$split_label(), ctx$split_plural(), ctx$scope())),
                    tags$h4(
                        style = "margin: 10px",
                        radioButtons(
                            NS(id, "gc"), "Gene (# of off-target clusters)",
                            gene_to_select
                        )
                    ),
                    actionButton(NS(id, "pass"), "Select gene")
                )
            })

            output$pairs <- renderUI({
                res <- result()
                if (is.null(res) || is.null(res$pairs)) return(NULL)
                if (nrow(res$pairs) == 0) {
                    return(tags$p(style = "margin: 10px", "Fewer than two candidate genes; no pairs to rank."))
                }
                tbl <- res$pairs
                names(tbl) <- c("Gene 1", "Gene 2", "Clusters where both are ON",
                    "Off-targets of each gene, summed")
                tagList(
                    tags$h4(style = "margin: 10px", "Best pairs for ", cluster_label(res$ds, res$ident)),
                    renderTable(tbl, striped = TRUE)
                )
            })

            other_clust <- eventReactive(
                input$pass, {
                    # Every cluster where the gene is ON in any split, like the
                    # off-target count shown next to it
                    b <- ctx$b()
                    anyon <- any_on(b, genes = input$gc)
                    b$clusters[anyon[1, ]]
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
                                cluster_choices(ctx$ds(), other_clust()),
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
                        ctx$parent_session,
                        inputId = "clusters_of_interest",
                        selected = input$clust
                    )
                    updateTabsetPanel(
                        ctx$parent_session,
                        inputId = "controlTabs",
                        selected = "Find distinct genes within"
                    )
                    ctx$first_mark(input$gc)
                }
            )

            observeEvent(
                input$get_them_all, {
                    updateCheckboxGroupInput(session, "clust", selected = other_clust())
                }
            )

            observeEvent(
                input$not_at_all, {
                    updateCheckboxGroupInput(session, "clust", selected = character(0))
                }
            )
        }
    )
}
