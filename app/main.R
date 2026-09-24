# Phage Synergy Explorer
#
# A dataset-agnostic front end to a Bayesian Bliss interaction model.
# Two ways in, one engine:
#   Quick pair - type the effects of two phages and of their combination.
#   Data table - upload / paste / load a table of any size and layout.

box::use(
  DT[DTOutput, datatable, formatRound, renderDT],
  bslib,
  dplyr[left_join],
  ggplot2[ggsave],
  plotly[plotlyOutput, renderPlotly],
  readr[write_csv],
  shiny,
  utils[head],
)
box::use(
  app/logic/bliss[bliss_table, fit_bayes, jags_ready, summarise_treatments],
  app/logic/data_input,
  app/logic/plots[as_interactive, matrix_plot, pair_plot, posterior_density, posterior_forest],
)

number_row <- function(...) shiny$div(class = "d-flex gap-2", ...)

plot_help <- paste("Hover for the numbers, drag to zoom, double-click to reset.",
                   "The expand icon at the bottom right of the card fills the screen.")

sidebar_ui <- function(ns) {
  bslib$sidebar(
    width = 380, class = "bg-body-tertiary",
    bslib$accordion(
      id = ns("steps"), multiple = TRUE, open = c("input", "rule"),

      bslib$accordion_panel(
        "1 - Input", value = "input",
        shiny$radioButtons(ns("mode"), NULL, inline = TRUE,
                           choices = c("Quick pair" = "pair", "Data table" = "table")),

        shiny$conditionalPanel(
          "input.mode == 'pair'", ns = ns,
          shiny$textInput(ns("pair_strain"), "Bacterial strain", "Strain 1"),
          number_row(
            shiny$textInput(ns("pair_name_a"), "Phage A", "Phage A"),
            shiny$textInput(ns("pair_name_b"), "Phage B", "Phage B")
          ),
          shiny$textInput(ns("pair_a"), "Effect of phage A alone", "0.08"),
          shiny$textInput(ns("pair_b"), "Effect of phage B alone", "0.99"),
          shiny$textInput(ns("pair_ab"), "Effect of the combination", "0.39"),
          shiny$div(class = "hint",
                    "One value, or several replicates separated by commas or spaces",
                    "(same count in each field).")
        ),

        shiny$conditionalPanel(
          "input.mode == 'table'", ns = ns,
          shiny$radioButtons(ns("source"), "Source", inline = TRUE,
                             choices = c("Upload" = "upload", "Paste" = "paste",
                                         "Example" = "example")),
          shiny$conditionalPanel(
            "input.source == 'upload'", ns = ns,
            shiny$fileInput(ns("file"), "Table file", width = "100%",
                            accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls", ".rds",
                                       ".RData", ".rda")),
            shiny$uiOutput(ns("item_ui"))
          ),
          shiny$conditionalPanel(
            "input.source == 'paste'", ns = ns,
            shiny$textAreaInput(
              ns("pasted"), "Paste a table (with a header row)", rows = 6, width = "100%",
              placeholder = "Bacteria\tReplica\tP1\tP2\tP1+P2\nStrain1\t1\t0.08\t0.99\t0.39"
            )
          ),
          shiny$conditionalPanel(
            "input.source == 'example'", ns = ns,
            shiny$div(class = "hint pt-2",
                      "Simulated screen (not real data): 8 strains x 3 replicates, three phages,",
                      "their pairs and the triple. Strains 01-03 are built synergistic,",
                      "04-06 additive, 07-08 antagonistic."),
            # Served straight from app/static, no download handler needed.
            shiny$tags$a(href = "static/example_mock.csv", download = "phage_synergy_example.csv",
                         "Download it as a template (CSV)")
          )
        )
      ),

      bslib$accordion_panel(
        "2 - Columns", value = "columns",
        shiny$conditionalPanel(
          "input.mode == 'pair'", ns = ns,
          shiny$div(class = "hint", "Not needed for a typed-in pair.")
        ),
        shiny$conditionalPanel(
          "input.mode == 'table'", ns = ns,
          shiny$radioButtons(ns("layout"), "Table layout", inline = TRUE,
                             choices = c("Wide (one column per treatment)" = "wide",
                                         "Long (one row per measurement)" = "long")),
          shiny$selectInput(ns("strain_col"), "Strain column", choices = NULL, width = "100%"),
          shiny$selectInput(ns("rep_col"), "Replicate column", choices = NULL, width = "100%"),
          shiny$conditionalPanel(
            "input.layout == 'long'", ns = ns,
            shiny$selectInput(ns("treat_col"), "Treatment column", choices = NULL,
                              width = "100%"),
            shiny$selectInput(ns("value_col"), "Value column", choices = NULL, width = "100%")
          ),
          shiny$textInput(ns("sep"), "Combination separator", "", placeholder = "auto"),
          shiny$div(class = "hint", "Blank = detect from the column names."),
          shiny$selectizeInput(ns("ignore_cols"), "Columns to ignore", choices = NULL,
                               multiple = TRUE, width = "100%")
        )
      ),

      bslib$accordion_panel(
        "3 - Decision rule", value = "rule",
        shiny$radioButtons(ns("scale"), "Effect scale", inline = TRUE,
                           choices = c("Auto" = "auto", "0-1" = "unit", "0-100 %" = "percent")),
        shiny$checkboxInput(ns("auto_rope"), "Set the threshold from replicate noise", TRUE),
        shiny$numericInput(ns("rope"), "Practical threshold (ROPE)", 0.05,
                           min = 0, max = 1, step = 0.01),
        shiny$sliderInput(ns("certainty"), "Posterior certainty required for a call",
                          0.5, 0.99, 0.95, 0.01)
      ),

      bslib$accordion_panel(
        "4 - Bayesian model", value = "model",
        if (!jags_ready) shiny$div(class = "text-danger small",
                                   "JAGS not found: point estimates only."),
        shiny$selectizeInput(ns("fit_combos"), "Combinations to fit", choices = NULL,
                             multiple = TRUE, width = "100%"),
        number_row(
          shiny$numericInput(ns("chains"), "Chains", 3, 1, 6, 1),
          shiny$numericInput(ns("thin"), "Thin", 1, 1, 20, 1)
        ),
        number_row(
          shiny$numericInput(ns("burnin"), "Burn-in", 2000, 100, 50000, 500),
          shiny$numericInput(ns("iter"), "Iterations", 5000, 500, 100000, 500)
        ),
        shiny$actionButton(ns("run_bayes"), "Fit posterior", class = "btn-primary w-100",
                           disabled = if (jags_ready) NULL else NA),
        shiny$div(class = "hint", shiny$textOutput(ns("fit_hint"), inline = TRUE))
      ),

      bslib$accordion_panel(
        "5 - Export", value = "export",
        shiny$downloadButton(ns("dl_results"), "Results (CSV)",
                             class = "btn-outline-secondary w-100 mb-2"),
        shiny$downloadButton(ns("dl_plot"), "Current figure (PNG)",
                             class = "btn-outline-secondary w-100")
      )
    )
  )
}

#' @export
ui <- function(id) {
  ns <- shiny$NS(id)
  bslib$page_sidebar(
    title = shiny$tags$span(
      shiny$tags$img(src = "static/SynerPhage.png", height = "32px", class = "me-2"),
      shiny$tags$strong("Phage Synergy Explorer"),
      shiny$tags$span(class = "text-muted ms-2 fs-6", "Bayesian Bliss interaction analysis")
    ),
    theme = bslib$bs_theme(
      version = 5, preset = "shiny",
      primary = "#12626B", success = "#1F7A5A", danger = "#B5443A",
      "border-radius" = "0.6rem"
    ),
    fillable = FALSE,
    sidebar = sidebar_ui(ns),

    bslib$layout_columns(
      fill = FALSE,
      bslib$value_box("Strains", shiny$textOutput(ns("vb_strains")), theme = "primary"),
      bslib$value_box("Combinations", shiny$textOutput(ns("vb_combos")), theme = "secondary"),
      bslib$value_box("Threshold (ROPE)", shiny$textOutput(ns("vb_rope")), theme = "secondary"),
      bslib$value_box("Synergy calls", shiny$textOutput(ns("vb_calls")), theme = "success")
    ),

    bslib$navset_card_tab(
      id = ns("tabs"), full_screen = TRUE,

      bslib$nav_panel(
        "Interaction map",
        bslib$layout_columns(
          fill = FALSE, col_widths = c(4, 3, 5),
          shiny$selectInput(ns("map_fill"), "Colour by", choices = NULL, width = "100%"),
          shiny$selectInput(ns("map_sort"), "Strain order", width = "100%",
                            choices = c("Alphabetical" = "name", "By interaction" = "effect")),
          shiny$div(class = "hint pt-4", plot_help)
        ),
        shiny$uiOutput(ns("map_ui"))
      ),

      bslib$nav_panel("Results", DTOutput(ns("results_table"))),

      bslib$nav_panel(
        "Posterior",
        bslib$layout_columns(
          fill = FALSE, col_widths = c(5, 2, 2, 3),
          shiny$selectizeInput(ns("post_combos"), "Combinations shown", choices = NULL,
                               multiple = TRUE, width = "100%"),
          shiny$selectInput(ns("post_sort"), "Strain order", width = "100%",
                            choices = c("Alphabetical" = "name", "By effect" = "effect")),
          shiny$numericInput(ns("post_ncol"), "Panels per row", 3, 1, 6, 1),
          shiny$div(class = "hint pt-4", plot_help)
        ),
        shiny$uiOutput(ns("posterior_ui")),
        DTOutput(ns("convergence_table"))
      ),

      bslib$nav_panel("Data", shiny$htmlOutput(ns("detection")), DTOutput(ns("data_table"))),
      bslib$nav_panel("Method", bslib$card_body(shiny$uiOutput(ns("method"))))
    )
  )
}

#' @export
server <- function(id) {
  shiny$moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Turn an engine error into a message in the output rather than a stack trace.
    checked <- function(expr) {
      tryCatch(expr, error = function(e) shiny$validate(shiny$need(FALSE, conditionMessage(e))))
    }

    upload_ext <- shiny$reactive({
      shiny$req(input$file)
      data_input$file_type(input$file$name)
    })

    output$item_ui <- shiny$renderUI({
      shiny$req(input$file)
      items <- checked(data_input$source_items(input$file$datapath, upload_ext()))
      if (!length(items)) return(NULL)
      shiny$selectInput(ns("item"), "Sheet / object", choices = items, width = "100%")
    })

    table_df <- shiny$reactive({
      shiny$req(input$source)
      switch(input$source,
        upload = {
          shiny$req(input$file)
          checked(data_input$read_source(input$file$datapath, upload_ext(), input$item))
        },
        paste = {
          shiny$validate(shiny$need(nzchar(trimws(input$pasted %||% "")),
                                    "Paste a table with a header row."))
          checked(data_input$read_pasted(input$pasted))
        },
        example = checked(data_input$read_source(data_input$example_file, "csv"))
      )
    })

    # Column pickers follow the loaded table.
    shiny$observeEvent(table_df(), {
      cols <- names(table_df())
      strain <- data_input$guess_col(cols, "^(bacteri|strain|isolate|species|host|sample)")
      rep <- data_input$guess_col(cols, "^(replica|replicate|rep$|repeat|run)")
      treat <- data_input$guess_col(cols, "^(treatment|phage|condition|combination|association)")
      value <- data_input$guess_col(cols, "^(value|effect|auc|score|response|inhibition)")

      shiny$updateSelectInput(session, "strain_col", choices = c("(none)" = "", cols),
                              selected = strain)
      shiny$updateSelectInput(session, "rep_col", choices = c("(none)" = "", cols),
                              selected = rep)
      shiny$updateSelectInput(session, "treat_col", choices = cols, selected = treat)
      shiny$updateSelectInput(session, "value_col", choices = cols, selected = value)
      shiny$updateSelectizeInput(
        session, "ignore_cols", choices = cols,
        selected = cols[grepl(data_input$control_pattern, cols, ignore.case = TRUE)]
      )
      shiny$updateTextInput(session, "sep", value = data_input$detect_separator(
        setdiff(cols, c(strain, rep, treat, value))
      ))
    })

    # In a long table the combination names sit in the treatment column rather
    # than in the header, so the separator can only be read from those values.
    shiny$observeEvent(list(input$layout, input$treat_col), {
      shiny$req(identical(input$layout, "long"), input$treat_col %in% names(table_df()))
      shiny$updateTextInput(session, "sep", value = data_input$detect_separator(
        unique(trimws(as.character(table_df()[[input$treat_col]])))
      ))
    })

    prep <- shiny$reactive({
      if (identical(input$mode, "pair")) {
        df <- checked(data_input$pair_dataset(input$pair_strain, input$pair_name_a,
                                              input$pair_name_b, input$pair_a, input$pair_b,
                                              input$pair_ab))
        checked(data_input$prepare_dataset(df, "bacteria", "replica", "+", scale = input$scale))
      } else {
        df <- table_df()
        if (identical(input$layout, "long")) {
          df <- checked(data_input$widen_long(df, input$strain_col, input$rep_col,
                                              input$treat_col, input$value_col))
          checked(data_input$prepare_dataset(df, input$strain_col, input$rep_col, input$sep,
                                             scale = input$scale))
        } else {
          checked(data_input$prepare_dataset(df, input$strain_col, input$rep_col, input$sep,
                                             input$ignore_cols, input$scale))
        }
      }
    })

    bayes <- shiny$reactiveVal(NULL)

    treatment_summary <- shiny$reactive(summarise_treatments(prep()))

    point <- shiny$reactive({
      checked(bliss_table(prep(),
                          rope = if (isTRUE(input$auto_rope)) NULL else input$rope,
                          summary = treatment_summary()))
    })

    # Show the noise-derived threshold in the sidebar. No feedback loop: the
    # reactive above only reads input$rope when the automatic rule is off.
    shiny$observeEvent(point(), {
      if (isTRUE(input$auto_rope)) {
        shiny$updateNumericInput(session, "rope", value = round(point()$rope[1], 4))
      }
    })

    shiny$observeEvent(prep(), {
      combos <- names(prep()$combos)
      shiny$updateSelectizeInput(session, "fit_combos", choices = combos, selected = combos)
      bayes(NULL)
    })

    output$fit_hint <- shiny$renderText({
      n <- length(input$fit_combos)
      if (!n) return("Select at least one combination.")
      sprintf("%d model%s to fit, %d chains of %d iterations each.",
              n, if (n > 1) "s" else "", input$chains, input$iter)
    })

    shiny$observeEvent(input$run_bayes, {
      p <- prep()
      combos <- input$fit_combos
      res <- tryCatch(
        shiny$withProgress(message = "Fitting the Bliss model", value = 0, {
          fit_bayes(
            p, combos, rope = point()$rope[1],
            chains = input$chains, burnin = input$burnin, iter = input$iter,
            thin = input$thin, certainty = input$certainty,
            on_step = function(i, cmb) {
              shiny$setProgress(value = (i - 1) / length(combos),
                                detail = sprintf("%d/%d - %s", i, length(combos), cmb))
            }
          )
        }),
        error = function(e) {
          shiny$showNotification(conditionMessage(e), type = "error", duration = 12)
          NULL
        }
      )
      if (!is.null(res)) {
        bayes(res)
        shiny$showNotification("Posterior ready.", type = "message")
      }
    })

    results <- shiny$reactive({
      out <- point()
      b <- bayes()
      if (!is.null(b)) {
        out <- left_join(out, b$posterior, by = c("bacteria", "phage_combination"))
      }
      out
    })

    # --- headline numbers ---

    output$vb_strains <- shiny$renderText(length(prep()$strains))
    output$vb_combos <- shiny$renderText(length(prep()$combos))
    output$vb_rope <- shiny$renderText(sprintf("%.1f pts", 100 * point()$rope[1]))
    output$vb_calls <- shiny$renderText({
      res <- results()
      call <- if ("bayes_category" %in% names(res)) res$bayes_category else res$point_category
      sprintf("%d / %d", sum(call == "Synergy", na.rm = TRUE), nrow(res))
    })

    # --- figures ---

    # Colouring the map by the posterior only makes sense once it has been fitted.
    map_fill_choices <- shiny$reactive({
      ch <- c("Observed - expected" = "delta_observed")
      if ("synergy_delta" %in% names(results())) {
        ch <- c("Posterior synergy effect" = "synergy_delta", ch, "P(synergy)" = "p_synergy")
      }
      ch
    })

    shiny$observeEvent(map_fill_choices(), {
      ch <- map_fill_choices()
      keep <- shiny$isolate(input$map_fill)
      shiny$updateSelectInput(session, "map_fill", choices = ch,
                              selected = if (isTRUE(keep %in% ch)) keep else ch[[1]])
    })

    panels_per_row <- shiny$reactive({
      n <- input$post_ncol
      if (length(n) != 1L || is.na(n) || n < 1L) 3L else as.integer(n)
    })

    current_plot <- shiny$reactive({
      res <- results()
      if (nrow(res) == 1L) return(pair_plot(treatment_summary(), res[1, ]))
      ch <- map_fill_choices()
      fill <- if (isTRUE(input$map_fill %in% ch)) input$map_fill else ch[[1]]
      matrix_plot(res, fill, names(ch)[match(fill, ch)], "Bliss interaction map",
                  input$map_sort %||% "name")
    })

    # A 25-strain matrix squeezed into a 600 px box is unreadable. Give every
    # strain its own band of pixels and let the card scroll instead.
    output$map_ui <- shiny$renderUI({
      h <- max(520L, 26L * length(prep()$strains) + 230L)
      plotlyOutput(ns("main_plot"), height = sprintf("%dpx", h))
    })

    output$main_plot <- renderPlotly(as_interactive(current_plot()))

    shiny$observeEvent(bayes(), {
      shiny$req(bayes())
      ch <- unique(bayes()$posterior$phage_combination)
      shiny$updateSelectizeInput(session, "post_combos", choices = ch, selected = ch)
    })

    posterior_shown <- shiny$reactive({
      b <- bayes()
      shiny$req(b)
      if (!length(input$post_combos)) {
        b$posterior
      } else {
        b$posterior[b$posterior$phage_combination %in% input$post_combos, , drop = FALSE]
      }
    })

    output$posterior_ui <- shiny$renderUI({
      shiny$validate(shiny$need(!is.null(bayes()),
                                "No posterior yet: press 'Fit posterior' in the sidebar."))
      post <- posterior_shown()
      shiny$validate(shiny$need(nrow(post) > 0L, "Select at least one combination."))
      if (nrow(post) == 1L) return(plotlyOutput(ns("posterior_plot"), height = "460px"))
      rows <- ceiling(length(unique(post$phage_combination)) / panels_per_row())
      h <- rows * (18L * length(unique(post$bacteria)) + 130L) + 90L
      plotlyOutput(ns("posterior_plot"), height = sprintf("%dpx", min(h, 8000L)))
    })

    output$posterior_plot <- renderPlotly({
      b <- bayes()
      post <- posterior_shown()
      if (nrow(post) == 1L) {
        as_interactive(posterior_density(b$draws, post$phage_combination[[1]], 1L,
                                         post$bacteria[[1]], b$rope),
                       tooltip = c("x", "y"))
      } else {
        as_interactive(posterior_forest(post, b$rope, panels_per_row(),
                                        input$post_sort %||% "name"))
      }
    })

    output$convergence_table <- renderDT({
      b <- bayes()
      shiny$req(b)
      datatable(b$convergence, rownames = FALSE,
                caption = sprintf(paste("Convergence (target: R-hat < 1.1, ESS > 400).",
                                        "Prior odds of relevant synergy: %.4f"), b$prior_odds),
                options = list(dom = "tp", pageLength = 5))
    })

    # --- tables ---

    output$results_table <- renderDT({
      res <- results()
      keep <- intersect(c("bacteria", "phage_combination", "components", "n_phages",
                          "observed", "observed_sd", "bliss_expected_point", "delta_observed",
                          "point_category", "bliss_expected", "combo_effect", "synergy_delta",
                          "synergy_delta_lower", "synergy_delta_upper", "p_synergy",
                          "p_additivity", "p_antagonism", "bayes_factor_synergy",
                          "bayes_category"), names(res))
      res <- res[keep]
      num <- names(res)[vapply(res, is.numeric, logical(1))]
      datatable(res, rownames = FALSE, filter = "top",
                options = list(pageLength = 15, scrollX = TRUE)) |>
        formatRound(num, 3)
    })

    output$data_table <- renderDT({
      datatable(prep()$data, rownames = FALSE,
                options = list(pageLength = 10, scrollX = TRUE)) |>
        formatRound(setdiff(names(prep()$data), c("bacteria", "replica")), 3)
    })

    output$detection <- shiny$renderUI({
      p <- prep()
      cl <- p$clamped
      n_clamped <- cl$n_clamped_low + cl$n_clamped_high
      shiny$tags$div(
        class = "small text-muted mb-3",
        shiny$tags$p(
          shiny$tags$strong("Detected: "),
          sprintf(paste("separator '%s' | %d single phages (%s) | %d combinations |",
                        "%d strains | %d replicates"),
                  p$sep, length(p$singles), paste(head(p$singles, 8), collapse = ", "),
                  length(p$combos), length(p$strains), p$n_rep)
        ),
        if (p$rescaled) shiny$tags$p("Values were divided by 100 to reach the 0-1 scale."),
        if (n_clamped > 0) {
          shiny$tags$p(sprintf(
            "%d of %d values were clamped into [0.001, 0.999] for the Beta likelihood.",
            n_clamped, cl$n_values
          ))
        }
      )
    })

    output$method <- shiny$renderUI({
      shiny$markdown(sprintf("
### What the app computes

**Bliss independence** is the null model: two phages that do not interact leave
behind the product of what each leaves behind, so the expected combination
effect is `1 - prod(1 - p_k)`. For two phages this is `p1 + p2 - p1*p2`. The
generalised form is used, so a triple cocktail is handled the same way.

**Interaction** is the gap between the observed combination effect and that
expectation. A gap is only called when it exceeds the *region of practical
equivalence* (ROPE, currently %.3f), by default the median within-strain
replicate SD floored at 0.05 - a difference smaller than the assay's own noise
is not an interaction.

**The Bayesian fit** (`bliss_zhao_n`, JAGS) estimates latent single-phage
effects per strain on the logit scale, derives the Bliss expectation from them,
and adds a strain-specific deviation `delta` to the combination. Replicates
enter through Beta likelihoods with separate precisions for single and
combination measurements. The reported quantity is
`synergy_effect = p12 - pred`, with `P(synergy_effect > ROPE)` driving the call
at the selected certainty (%.2f). Anything less certain stays *Inconclusive*
rather than being forced into a category.

**Bayes factor** compares the posterior odds of a relevant synergy against the
prior odds obtained by simulating the same priors forward, so the evidence is
reported separately from the prior.

### What the app expects

Any table whose header names the treatments. A column is read as a combination
when its name splits into components that are themselves columns:
`P1+P2`, `p1_p2`, `phageA:phageB` and `P1+P2+P3` all work. Controls such as
`Nophage` or `all_phages` are ignored by default. Effects must be fractions of
growth suppressed; values above 1 are read as percentages.
", point()$rope[1], input$certainty))
    })

    # --- exports ---

    output$dl_results <- shiny$downloadHandler(
      filename = function() sprintf("phage_synergy_results_%s.csv", Sys.Date()),
      content = function(file) write_csv(results(), file)
    )

    # Export whichever figure is on screen, at a height that matches how many
    # strains and panels it actually has.
    on_posterior <- shiny$reactive(identical(input$tabs, "Posterior") && !is.null(bayes()))

    output$dl_plot <- shiny$downloadHandler(
      filename = function() {
        sprintf("phage_synergy_%s_%s.png", if (on_posterior()) "posterior" else "map",
                Sys.Date())
      },
      content = function(file) {
        n <- length(prep()$strains)
        if (on_posterior()) {
          post <- posterior_shown()
          p <- posterior_forest(post, bayes()$rope, panels_per_row(),
                                input$post_sort %||% "name")
          rows <- ceiling(length(unique(post$phage_combination)) / panels_per_row())
          h <- min(40, rows * (0.14 * n + 1.4))
        } else {
          p <- current_plot()
          h <- max(6, min(30, 0.3 * n + 3))
        }
        # The tooltip aesthetic is unknown to ggplot itself; see as_interactive().
        suppressWarnings(ggsave(file, p, width = 12, height = h, dpi = 200,
                                limitsize = FALSE, bg = "white"))
      }
    )
  })
}
