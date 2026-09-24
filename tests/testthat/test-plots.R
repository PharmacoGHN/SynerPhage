box::use(
  plotly[plotly_build],
  testthat[expect_s3_class, expect_setequal, skip_if_not, test_that],
)
box::use(
  app/logic/bliss,
  app/logic/data_input,
  app/logic/plots,
)

prep <- data_input$prepare_dataset(data_input$read_source(data_input$example_file, "csv"),
                                   "Bacteria", "Replica", ignore = "Nophage")
res <- bliss$bliss_table(prep)
prep_pair <- data_input$prepare_dataset(
  data_input$pair_dataset("Strain 1", "P1", "P2", "0.2 0.25", "0.4, 0.45", "0.9 0.88"),
  "bacteria", "replica", "+"
)

builds <- function(p, ...) expect_s3_class(plotly_build(plots$as_interactive(p, ...)), "plotly")

test_that("figures survive the conversion to interactive plots", {
  builds(plots$matrix_plot(res, "delta_observed", "Observed - expected", "map", "effect"))
  builds(plots$pair_plot(bliss$summarise_treatments(prep_pair), bliss$bliss_table(prep_pair)))
})

test_that("strain order by effect is a permutation of the strains", {
  expect_setequal(plots$strain_levels(res, "delta_observed", "effect"), prep$strains)
})

test_that("the forest plot converts to plotly", {
  skip_if_not(bliss$jags_ready, "JAGS not installed")
  fit <- bliss$fit_bayes(prep_pair, names(prep_pair$combos), rope = 0.05,
                         chains = 2, burnin = 500, iter = 1000)
  builds(plots$posterior_forest(fit$posterior, 0.05, 1))
})
