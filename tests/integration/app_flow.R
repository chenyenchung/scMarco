# Drives the whole app server through the main workflow for every dataset in
# a config. Needs the converted data, so it is not part of tests/testthat.
#
#   SCMARCO_CONFIG=configs/all.yaml Rscript tests/integration/app_flow.R

library(testthat)
suppressMessages(shiny::shinyOptions(cache = NULL))

app_env <- new.env()
app <- suppressMessages(sys.source("app.R", envir = app_env, toplevel.env = app_env))
server <- app_env$server
config <- app_env$config

render_ok <- function(output_fun) {
  # Evaluate a render function and fail on errors
  res <- tryCatch(output_fun, error = function(e) e)
  if (inherits(res, "error") && !inherits(res, "shiny.silent.error")) stop(res)
  res
}

download_ok <- function(session, id, type) {
  do.call(session$setInputs, setNames(list(type), paste0(id, "-filetype")))
  path <- session$getOutput(paste0(id, "-download"))
  expect_true(file.exists(path), info = paste(id, type))
  expect_gt(file.info(path)$size, 1000)
}

for (ds_id in intersect(names(config$datasets), strsplit(Sys.getenv("ONLY", paste(names(config$datasets), collapse=",")), ",")[[1]])) {
  test_that(paste("main workflow runs on", ds_id), {
    shiny::testServer(server, {
      session$setInputs(dataset = ds_id, gene_group_to_use = "All", margin = 10)
      ds_obj <- ds()
      session$setInputs(required_splits = character(0), min_n = 1)

      # From cluster: pick a cluster with a candidate marker that has off-targets
      b_obj <- b()
      pick <- NULL
      # Try named clusters first so the labels are exercised
      named <- !is.na(ds_obj$cluster_names)
      for (cl in c(ds_obj$clusters[named], ds_obj$clusters[!named])) {
        fc <- from_cluster(b_obj, cl, character(0), 1, limit = 10,
                           genes = ds_obj$analysis_genes)
        if (length(fc) > 0 && any(fc >= 1 & fc <= 6)) {
          pick <- list(cluster = cl, gene = names(fc)[fc >= 1 & fc <= 6][1])
          break
        }
      }
      expect_false(is.null(pick))
      session$setInputs(clusters_of_interest = pick$cluster)
      session$setInputs(`fexp-output_num` = 10, `fexp-show_pairs` = TRUE)
      session$setInputs(`fexp-do_fexp` = 1)
      expect_match(as.character(output$`fexp-selectTbl`$html), pick$gene, fixed = TRUE)
      render_ok(output$`fexp-pairs`)
      expect_match(as.character(output$`fexp-pairs`$html),
                   htmltools::htmlEscape(cluster_label(ds_obj, pick$cluster)), fixed = TRUE)

      # Select gene -> clusters sent to Find distinct, gene becomes first mark
      session$setInputs(`fexp-gc` = pick$gene, `fexp-pass` = 1)
      clust <- b_obj$clusters[any_on(b_obj, genes = pick$gene)[1, ]]
      session$setInputs(`fexp-clust` = clust, `fexp-submit` = 1)
      expect_equal(first_mark(), pick$gene)

      session$setInputs(clusters_of_interest = clust)
      session$setInputs(`fd_genes-k` = 1, `fd_genes-do_onplot` = 1)
      render_ok(output$`fd_genes-onplot`)
      hits_html <- output$`fd_genes-hits`
      dt_hits <- distinct_within(b_obj, clust, character(0), 1,
                                 genes = ds_obj$analysis_genes)
      for (cl in unique(dt_hits$cluster)) {
        expect_match(hits_html, htmltools::htmlEscape(cluster_label(ds_obj, cl)), fixed = TRUE)
      }
      if (nrow(dt_hits) > 0) {
        expect_match(hits_html, paste(pick$gene, "&amp; gene both ON"), fixed = TRUE)
      }
      download_ok(session, "fd_genes-dl_fdist", "png")

      # Co-expression of the first mark plus one distinct gene
      dt <- distinct_within(b_obj, clust, character(0), 1,
                            genes = ds_obj$analysis_genes)
      genes <- unique(c(pick$gene, head(dt$gene, 1)))
      session$setInputs(`coexp-split_genes` = genes, `coexp-show_pct` = TRUE,
                        `coexp-do_splitplot` = 1)
      render_ok(output$`coexp-splitplot`)
      download_ok(session, "coexp-dl_coexp", "pdf")

      # From gene with two genes in "all" mode
      session$setInputs(`fpos-gene` = genes, `fpos-mode` = "all",
                        `fpos-do_fpos` = 1)
      render_ok(output$`fpos-selectTbl`)

      # Split criteria: require the first split
      session$setInputs(required_splits = ds_obj$splits[1])
      expect_equal(required(), ds_obj$splits[1])
      session$setInputs(`fexp-do_fexp` = 2)
      render_ok(output$`fexp-selectTbl`)
      session$setInputs(required_splits = character(0))

      if (ds_obj$has_mean && ds_obj$ordered_splits) {
        session$setInputs(`lplot-gene_of_interest` = pick$gene,
                          `lplot-line_alpha` = 0.2, `lplot-do_lineplot` = 1)
        render_ok(output$`lplot-lineplot`)
        download_ok(session, "lplot-dl_lplot", "png")
      }

      if (length(ds_obj$reductions) > 0) {
        choices <- reduction_choices(ds_obj)
        n_plot <- 0
        # The first shared reduction and every per-split group
        for (choice in unique(c(choices[!startsWith(choices, "split:")][1],
                                choices[startsWith(choices, "split:")]))) {
          if (is.na(choice)) next
          per_split <- startsWith(choice, "split:")
          shown <- head(names(resolve_reduction(ds_obj, choice)), 2)
          for (viz in if (per_split) c("gene", "cluster") else c("gene", "cluster", "split")) {
            n_plot <- n_plot + 1
            session$setInputs(
              `umap-umap_selection` = choice, `umap-viz_type` = viz,
              `umap-gene_selection` = pick$gene, `umap-split_by_group` = TRUE,
              `umap-shown_splits` = shown, `umap-show_legend` = FALSE,
              `umap-facet_ncol` = 1, `umap-plot_umap` = n_plot
            )
            render_ok(output$`umap-umap_plot`)
            if (per_split && viz == "cluster") {
              p <- umap_plot(ds_obj, resolve_reduction(ds_obj, choice)[shown],
                             color_by = "cluster", splits = shown, ncol = 1)
              expect_setequal(as.character(unique(p$data$split)), shown)
              expect_equal(p$facet$params$ncol, 1)
            }
          }
        }
        download_ok(session, "umap-dl_umap", "png")
      }
    })
  })
}
