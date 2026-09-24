# Replicate-level Bayesian Bliss model (JAGS), for a combination of n_p phages.
#
# For each strain the model estimates latent single-phage effects, derives the
# Bliss expectation pred = 1 - prod(1 - p_k) from them (p1 + p2 - p1*p2 for a
# pair), and adds a strain-specific deviation delta on the logit scale to give
# the combination effect. synergy_effect = p12_mean - pred.
#
# Data: P[n_bact, n_rep, n_p], Y_obs[n_bact, n_rep].

#' @export
bliss_zhao_n <- "
  model{
    for (b in 1:n_bact){
      for (k in 1:n_p){
        logit_p_mean[b, k] ~ dnorm(0, 0.16)
        p_mean[b, k] <- ilogit(logit_p_mean[b, k])

        for (r in 1:n_rep){
          P[b, r, k] ~ dbeta(p_mean[b, k] * phi_single, (1 - p_mean[b, k]) * phi_single)
        }
      }

      surv[b, 1] <- 1 - p_mean[b, 1]
      for (k in 2:n_p){
        surv[b, k] <- surv[b, k - 1] * (1 - p_mean[b, k])
      }
      pred[b] <- 1 - surv[b, n_p]

      delta[b] ~ dnorm(0, tau_delta)
      logit_p12_mean[b] <- logit(pred[b]) + delta[b]
      p12_mean[b] <- ilogit(logit_p12_mean[b])
      synergy_effect[b] <- p12_mean[b] - pred[b]

      for (r in 1:n_rep){
        Y_obs[b, r] ~ dbeta(p12_mean[b] * phi_combo, (1 - p12_mean[b]) * phi_combo)
      }
    }

    sigma_delta ~ dunif(0, 2)
    tau_delta <- pow(sigma_delta, -2)

    phi_single ~ dgamma(2, 0.1)
    phi_combo ~ dgamma(2, 0.1)
  }
"
