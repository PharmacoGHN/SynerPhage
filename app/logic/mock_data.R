# Simulated screen used as the app's built-in example. NOT experimental data.
# app/static/example_mock.csv is make_mock_data() written with write.csv(row.names = FALSE);
# tests/testthat/test-mock_data.R checks the two still agree.
#
# 8 strains x 3 replicates, phages PhiA / PhiB / PhiC, their three pairs and
# the triple, plus a no-phage control. Each strain carries a known deviation
# from Bliss on the logit scale, so the app should find synergy for the first
# strains, additivity in the middle and antagonism for the last ones.

box::use(
  stats[plogis, qlogis, rbeta, runif, setNames],
)

#' @export
make_mock_data <- function(seed = 42) {
  set.seed(seed)
  strains <- sprintf("Strain_%02d", 1:8)
  delta <- setNames(c(2.0, 1.5, 1.0, 0, 0, 0, -1.0, -1.5), strains)
  phi <- 60 # Beta precision: sets the replicate noise

  draw <- function(mu, n = 3) rbeta(n, mu * phi, (1 - mu) * phi)
  bliss <- function(p) 1 - prod(1 - p)

  rows <- lapply(strains, function(s) {
    p <- c(PhiA = runif(1, 0.1, 0.5), PhiB = runif(1, 0.2, 0.6), PhiC = runif(1, 0.05, 0.4))
    combo <- function(k) draw(plogis(qlogis(bliss(p[k])) + delta[[s]]))
    data.frame(
      Bacteria = s, Replica = 1:3,
      PhiA = draw(p[["PhiA"]]), PhiB = draw(p[["PhiB"]]), PhiC = draw(p[["PhiC"]]),
      "PhiA+PhiB" = combo(c("PhiA", "PhiB")),
      "PhiA+PhiC" = combo(c("PhiA", "PhiC")),
      "PhiB+PhiC" = combo(c("PhiB", "PhiC")),
      "PhiA+PhiB+PhiC" = combo(c("PhiA", "PhiB", "PhiC")),
      Nophage = draw(0.02),
      check.names = FALSE
    )
  })

  out <- do.call(rbind, rows)
  effects <- setdiff(names(out), c("Bacteria", "Replica"))
  out[effects] <- lapply(out[effects], round, 4)
  out
}
