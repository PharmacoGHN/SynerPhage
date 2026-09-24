box::use(
  testthat[expect_equal, expect_gt, expect_lt, expect_true, skip_if_not, test_that],
)
box::use(
  app/logic/bliss,
  app/logic/data_input,
)

prep <- data_input$prepare_dataset(data_input$read_source(data_input$example_file, "csv"),
                                   "Bacteria", "Replica", ignore = "Nophage")
res <- bliss$bliss_table(prep)
summ <- bliss$summarise_treatments(prep)

prep_pair <- data_input$prepare_dataset(
  data_input$pair_dataset("Strain 1", "P1", "P2", "0.2 0.25", "0.4, 0.45", "0.9 0.88"),
  "bacteria", "replica", "+"
)
res_pair <- bliss$bliss_table(prep_pair)

test_that("generalised Bliss matches the closed forms", {
  p1 <- c(0.2, 0.5)
  p2 <- c(0.4, 0.9)
  p3 <- c(0.1, 0.3)
  expect_lt(max(abs(bliss$bliss_expectation(list(p1, p2)) - (p1 + p2 - p1 * p2))), 1e-12)
  expect_lt(max(abs(bliss$bliss_expectation(list(p1, p2, p3)) -
                      (1 - (1 - p1) * (1 - p2) * (1 - p3)))), 1e-12)
})

test_that("the wide table is analysed end to end", {
  expect_equal(length(prep$combos), 4)
  expect_equal(nrow(res), length(prep$strains) * 4)

  row <- res[res$phage_combination == "PhiA+PhiB", ][1, ]
  mean_of <- function(treatment) {
    summ$mean_effect[summ$bacteria == row$bacteria & summ$treatment == treatment]
  }
  pa <- mean_of("PhiA")
  pb <- mean_of("PhiB")
  expect_lt(abs(row$bliss_expected_point - (pa + pb - pa * pb)), 1e-12)
  expect_lt(abs(row$delta_observed - (row$observed - row$bliss_expected_point)), 1e-12)
  expect_true(row$rope >= 0.05)
})

test_that("simulated interactions are recovered", {
  pair_calls <- res$point_category[res$phage_combination == "PhiA+PhiB"]
  expect_equal(pair_calls[1], "Synergy")
  expect_equal(pair_calls[8], "Antagonism")
  expect_equal(res_pair$point_category, "Synergy")
})

test_that("JAGS matrix rows stay with their strain", {
  mat <- matrix(bliss$rectangular(prep)[["PhiA+PhiB"]], nrow = length(prep$strains),
                byrow = TRUE)
  combo_means <- summ[summ$treatment == "PhiA+PhiB", ]
  observed_by_strain <- combo_means$mean_effect[match(prep$strains, combo_means$bacteria)]
  expect_lt(max(abs(rowMeans(mat) - observed_by_strain)), 1e-12)
})

test_that("the posterior agrees with the observed synergy", {
  skip_if_not(bliss$jags_ready, "JAGS not installed")
  fit <- bliss$fit_bayes(prep_pair, names(prep_pair$combos), rope = 0.05,
                         chains = 2, burnin = 500, iter = 1000)
  expect_lt(abs(fit$posterior$synergy_delta - res_pair$delta_observed), 0.25)
  expect_gt(fit$posterior$p_synergy, 0.9)
  expect_gt(fit$posterior$bayes_factor_synergy, 1)
})
