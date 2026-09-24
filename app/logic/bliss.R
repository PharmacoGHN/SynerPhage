# Deterministic Bliss arithmetic and the Bayesian (JAGS) fit. Both take the
# list returned by data_input$prepare_dataset().

box::use(
  dplyr[.data, filter, group_by, left_join, mutate, rename, select, summarise],
  stats[median, sd, update],
  tidyr[pivot_longer, pivot_wider],
)
box::use(
  app/logic/bliss_model[bliss_zhao_n],
  app/logic/postprocessing[
    bayes_factor, build_posterior_long, categorise_bayes, categorise_point, convergence_table,
    prior_odds_synergy,
  ],
)

# JAGS is optional: without it the app still runs, it just cannot fit posteriors.
#' @export
jags_ready <- isTRUE(suppressWarnings(requireNamespace("rjags", quietly = TRUE)))

# --- deterministic Bliss ---------------------------------------------------

#' @export
summarise_treatments <- function(prep) {
  prep$data |>
    pivot_longer(-c("bacteria", "replica"), names_to = "treatment", values_to = "effect") |>
    group_by(.data$bacteria, .data$treatment) |>
    summarise(mean_effect = mean(.data$effect, na.rm = TRUE),
              sd_effect = sd(.data$effect, na.rm = TRUE),
              n_replicates = sum(!is.na(.data$effect)), .groups = "drop")
}

# Bliss independence for any number of components: the combination is expected
# to leave behind the product of what each phage alone leaves behind.
# For two phages, 1 - (1 - p1)(1 - p2) is the familiar p1 + p2 - p1 * p2.
#' @export
bliss_expectation <- function(p_list) 1 - Reduce(`*`, lapply(p_list, function(p) 1 - p))

#' @export
bliss_table <- function(prep, rope = NULL, summary = summarise_treatments(prep)) {
  means <- summary |>
    select("bacteria", "treatment", "mean_effect") |>
    pivot_wider(names_from = "treatment", values_from = "mean_effect")

  expected <- do.call(rbind, lapply(names(prep$combos), function(cmb) {
    comps <- prep$combos[[cmb]]
    data.frame(bacteria = means$bacteria,
               phage_combination = cmb,
               n_phages = length(comps),
               components = paste(comps, collapse = " + "),
               bliss_expected_point = bliss_expectation(means[comps]),
               stringsAsFactors = FALSE)
  }))

  observed <- summary |>
    filter(.data$treatment %in% names(prep$combos)) |>
    rename(phage_combination = "treatment", observed = "mean_effect",
           observed_sd = "sd_effect")

  out <- expected |>
    left_join(observed, by = c("bacteria", "phage_combination")) |>
    mutate(delta_observed = .data$observed - .data$bliss_expected_point)

  # Region of practical equivalence: the assay's own replicate noise, floored
  # at 5 points. A difference smaller than that noise is not called.
  if (is.null(rope) || !is.finite(rope)) {
    rope <- max(0.05, median(out$observed_sd, na.rm = TRUE))
  }
  if (!is.finite(rope)) rope <- 0.05

  out$rope <- rope
  out$point_category <- categorise_point(out$delta_observed, rope)
  out[order(out$bacteria, out$phage_combination), ]
}

# --- Bayesian fit ----------------------------------------------------------

# Rectangular (strain x replicate) layout, bacteria-major so that
# matrix(byrow = TRUE) puts one strain per row. Strains with fewer replicates
# are padded with NA, which JAGS treats as missing and imputes from the
# likelihood. Kept separate from the fit because this reshape is exactly where
# strain labels can drift away from their values (tests/testthat/test-bliss.R
# guards it).
#' @export
rectangular <- function(prep) {
  grid <- expand.grid(replica = seq_len(prep$n_rep), bacteria = prep$strains,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  left_join(grid[c("bacteria", "replica")], prep$data, by = c("bacteria", "replica"))
}

# Same three starting points and seeds as the published analysis, so a fit
# launched from the app reproduces it.
jags_inits <- function(seed_offset, chains) {
  base <- list(list(phi_single = 20, phi_combo = 20, sigma_delta = 0.15),
               list(phi_single = 10, phi_combo = 10, sigma_delta = 0.30),
               list(phi_single = 30, phi_combo = 30, sigma_delta = 0.05))
  lapply(seq_len(chains), function(k) {
    c(base[[(k - 1L) %% 3L + 1L]],
      list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 1000L * k + seed_offset))
  })
}

#' @export
fit_bayes <- function(prep, combos, rope, chains = 3, burnin = 2000, iter = 5000,
                      thin = 1, certainty = 0.95, on_step = NULL) {
  if (!jags_ready) stop("JAGS is not available: install JAGS and the rjags package.")
  if (!length(combos)) stop("Select at least one combination to fit.")
  # rjags needs the JAGS system library, so it is imported only once we know it
  # loads; a top-level import would stop the whole app on a server without JAGS.
  box::use(
    rjags[coda.samples, jags.model],
  )

  strains <- prep$strains
  n_bact <- length(strains)
  n_rep <- prep$n_rep

  full <- rectangular(prep)
  as_mat <- function(col) matrix(full[[col]], nrow = n_bact, byrow = TRUE)

  draws <- vector("list", length(combos))
  names(draws) <- combos

  for (i in seq_along(combos)) {
    cmb <- combos[[i]]
    comps <- prep$combos[[cmb]]
    if (is.null(comps)) stop("Unknown combination: ", cmb)
    if (!is.null(on_step)) on_step(i, cmb)

    p_single <- array(NA_real_, dim = c(n_bact, n_rep, length(comps)))
    for (k in seq_along(comps)) p_single[, , k] <- as_mat(comps[[k]])

    model <- jags.model(
      file = textConnection(bliss_zhao_n),
      data = list(n_bact = n_bact, n_rep = n_rep, n_p = length(comps),
                  P = p_single, Y_obs = as_mat(cmb)),
      inits = jags_inits(i, chains), n.chains = chains, quiet = TRUE
    )
    update(model, n.iter = burnin, progress.bar = "none")
    draws[[i]] <- coda.samples(
      model,
      variable.names = c("pred", "p12_mean", "synergy_effect",
                         "phi_single", "phi_combo", "sigma_delta"),
      n.iter = iter, thin = thin, progress.bar = "none"
    )
  }

  post <- build_posterior_long(draws, names(draws), strains, rope)
  prior_odds <- prior_odds_synergy(rope)
  post$bayes_factor_synergy <- bayes_factor(post$p_synergy, prior_odds)
  post$bayes_category <- categorise_bayes(post$p_synergy, post$p_additivity,
                                          post$p_antagonism, certainty)

  list(posterior = post, draws = draws, prior_odds = prior_odds,
       convergence = convergence_table(draws, names(draws)),
       rope = rope, certainty = certainty)
}
