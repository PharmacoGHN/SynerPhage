box::use(
  testthat[expect_equal, expect_false, expect_identical, expect_true, test_that],
  tidyr[pivot_longer],
)
box::use(
  app/logic/bliss[bliss_table],
  app/logic/data_input,
)

sample_df <- data_input$read_source(data_input$example_file, "csv")
effects <- setdiff(names(sample_df), c("Bacteria", "Replica"))
wide <- bliss_table(data_input$prepare_dataset(sample_df, "Bacteria", "Replica",
                                               ignore = "Nophage"))

test_that("treatments are discovered from the header", {
  found <- data_input$split_treatments(effects, data_input$detect_separator(effects))
  expect_identical(data_input$detect_separator(effects), "+")
  expect_equal(length(setdiff(found$singles, "Nophage")), 3)
  expect_equal(length(found$combos), 4)
  expect_false("Nophage" %in% names(found$combos))
  expect_identical(data_input$detect_separator(c("pA", "pB", "pC", "pA_pB", "pB_pC")), "_")
})

test_that("long input reproduces the wide analysis", {
  long_df <- as.data.frame(pivot_longer(sample_df[names(sample_df) != "Nophage"],
                                        -c("Bacteria", "Replica"),
                                        names_to = "treatment", values_to = "effect"),
                           check.names = FALSE)
  res_long <- bliss_table(data_input$prepare_dataset(
    data_input$widen_long(long_df, "Bacteria", "Replica", "treatment", "effect"),
    "Bacteria", "Replica"
  ))
  expect_true(max(abs(res_long$delta_observed - wide$delta_observed)) < 1e-12)
})

test_that("0-100 input is rescaled to 0-1", {
  pct_df <- sample_df
  pct_df[effects] <- pct_df[effects] * 100
  res_pct <- bliss_table(data_input$prepare_dataset(pct_df, "Bacteria", "Replica",
                                                    ignore = "Nophage"))
  expect_true(max(abs(res_pct$delta_observed - wide$delta_observed)) < 1e-9)
})

test_that("a typed pair builds one combination and keeps its replicates", {
  prep_pair <- data_input$prepare_dataset(
    data_input$pair_dataset("Strain 1", "P1", "P2", "0.2 0.25", "0.4, 0.45", "0.9 0.88"),
    "bacteria", "replica", "+"
  )
  expect_identical(names(prep_pair$combos), "P1 + P2")
  expect_equal(prep_pair$n_rep, 2)
})
