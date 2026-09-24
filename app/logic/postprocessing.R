# Post-processing helpers for the Bayesian Bliss model.
#
# Everything here is keyed by NAME rather than by position: writing model
# output into result tables positionally is how strain and combination labels
# silently drift away from their values.

box::use(
  coda[effectiveSize, gelman.diag],
  stats[plogis, qlogis, quantile, rnorm, runif],
)

# --- data hygiene -----------------------------------------------------------

#' @export
clamp_probability <- function(x, lower = 0.001, upper = 0.999) {
  pmin(pmax(x, lower), upper)
}

# Report how many values a clamp would move, so the boundary correction is
# declared instead of happening silently.
#' @export
clamp_report <- function(x, lower = 0.001, upper = 0.999) {
  x <- as.matrix(x)
  data.frame(
    n_values = length(x),
    n_clamped_low = sum(x < lower, na.rm = TRUE),
    n_clamped_high = sum(x > upper, na.rm = TRUE)
  )
}

# --- MCMC extraction --------------------------------------------------------

# Pull a monitored vector node and return its columns in numeric index order.
# grep() alone returns lexicographic order, which puts [10] before [2].
#' @export
mcmc_node <- function(mcmc_matrix, node) {
  nm <- colnames(mcmc_matrix)
  hit <- which(startsWith(nm, paste0(node, "[")))
  # A vector node of length 1 (a single bacterium) loses its index in the coda
  # column names, so "pred[1]" arrives as plain "pred".
  if (length(hit) == 0 && sum(nm == node) == 1) {
    return(mcmc_matrix[, nm == node, drop = FALSE])
  }
  if (length(hit) == 0) stop("node not monitored: ", node)
  idx <- suppressWarnings(as.integer(sub("^[^\\[]+\\[([0-9]+)\\]$", "\\1", nm[hit])))
  if (anyNA(idx)) stop("node is not a plain vector: ", node)
  mcmc_matrix[, hit[order(idx)], drop = FALSE]
}

posterior_summary <- function(mcmc_matrix, node, strains) {
  m <- mcmc_node(mcmc_matrix, node)
  stopifnot(ncol(m) == length(strains))
  q <- apply(m, 2, quantile, probs = c(0.025, 0.5, 0.975))
  data.frame(
    bacteria = strains,
    median = q[2, ], lower = q[1, ], upper = q[3, ], mean = colMeans(m),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# --- tidy posterior table ---------------------------------------------------

# One tidy row per (strain, combination). Because the combination name travels
# with the model object, a reordering of either input cannot scramble labels.
#' @export
build_posterior_long <- function(output_rjags, combination_names, strains, threshold) {
  stopifnot(length(output_rjags) == length(combination_names))
  out <- lapply(seq_along(output_rjags), function(k) {
    m <- as.matrix(output_rjags[[k]])
    pred <- posterior_summary(m, "pred", strains)
    p12 <- posterior_summary(m, "p12_mean", strains)
    dlt <- posterior_summary(m, "synergy_effect", strains)
    ds <- mcmc_node(m, "synergy_effect")
    data.frame(
      bacteria = strains,
      phage_combination = combination_names[k],
      bliss_expected = pred$median,
      bliss_expected_lower = pred$lower,
      bliss_expected_upper = pred$upper,
      combo_effect = p12$median,
      combo_effect_lower = p12$lower,
      combo_effect_upper = p12$upper,
      synergy_delta = dlt$median,
      synergy_delta_mean = dlt$mean,
      synergy_delta_lower = dlt$lower,
      synergy_delta_upper = dlt$upper,
      p_synergy = colMeans(ds > threshold),
      p_additivity = colMeans(abs(ds) <= threshold),
      p_antagonism = colMeans(ds < -threshold),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# --- categorisation ---------------------------------------------------------

# Bayesian call: a category is declared only when the posterior is at least
# `certainty` sure of it, otherwise the strain-pair stays Inconclusive.
#' @export
categorise_bayes <- function(p_syn, p_add, p_ant, certainty = 0.95) {
  ifelse(is.na(p_syn) | is.na(p_add) | is.na(p_ant), NA_character_,
    ifelse(p_syn >= certainty, "Synergy",
      ifelse(p_ant >= certainty, "Antagonism",
        ifelse(p_add >= certainty, "Additivity", "Inconclusive")
      )
    )
  )
}

# Point-estimate call for the deterministic Bliss rule. NA stays NA rather than
# falling through to "Additivity".
#' @export
categorise_point <- function(score, threshold) {
  ifelse(is.na(score), NA_character_,
    ifelse(score > threshold, "Synergy",
      ifelse(score < -threshold, "Antagonism", "Additivity")
    )
  )
}

# --- Bayes factor -----------------------------------------------------------

# Prior odds that a strain shows a relevant synergy, obtained by simulating the
# model priors forward. Posterior odds / prior odds gives the Bayes factor;
# posterior odds alone would confound the evidence with the prior, which here
# is far from 50:50.
#' @export
prior_odds_synergy <- function(threshold, n_sim = 2e6, seed = 20251010) {
  set.seed(seed)
  p1 <- plogis(rnorm(n_sim, 0, 1 / sqrt(0.16)))
  p2 <- plogis(rnorm(n_sim, 0, 1 / sqrt(0.16)))
  pred <- p1 + p2 - p1 * p2
  delta <- rnorm(n_sim, 0, runif(n_sim, 0, 2))
  p12 <- plogis(qlogis(pred) + delta)
  p <- mean((p12 - pred) > threshold)
  p / (1 - p)
}

#' @export
bayes_factor <- function(p_posterior, prior_odds, eps = 1e-6) {
  post_odds <- pmin(pmax(p_posterior, eps), 1 - eps)
  post_odds <- post_odds / (1 - post_odds)
  post_odds / prior_odds
}

# --- convergence ------------------------------------------------------------

#' @export
convergence_table <- function(output_rjags, combination_names,
                              params = c("phi_single", "phi_combo", "sigma_delta")) {
  do.call(rbind, lapply(seq_along(output_rjags), function(k) {
    s <- output_rjags[[k]]
    gd <- try(gelman.diag(s, multivariate = FALSE, autoburnin = FALSE), silent = TRUE)
    psrf <- if (inherits(gd, "try-error")) NA_real_ else max(gd$psrf[, 1], na.rm = TRUE)
    ess <- effectiveSize(s)
    keep <- startsWith(names(ess), "synergy_effect") | names(ess) %in% params
    data.frame(
      combination = combination_names[k],
      max_rhat = round(psrf, 4),
      min_ess = round(min(ess[keep], na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
}
