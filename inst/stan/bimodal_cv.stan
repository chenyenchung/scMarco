// bimodal.stan fit on the training samples (_t), with the log likelihood
// of the held-out samples (_h) for cross-validation.
// Adapted from level1_ordered_CV.stan by Fred P. Davis
// (https://github.com/fredpdavis/opticlobe; Davis et al. 2020, eLife
// 9:e50901), as used in mm_marker / mixture_modeling_smk.
//
// Binarizer 1.1 replaces the original prior on the means, normal(7, 5),
// which was set for bulk log1p(TPM) on 0-14, with a per-gene prior
// (fit_mixture()).
//
// Copyright (c) 2018-2020 Fred P. Davis; adaptation (c) 2023-2026
// Yen-Chung Chen. Released under the MIT license with Fred P. Davis's
// permission; see inst/COPYRIGHTS.

data {
  int<lower=2> n_t;
  int<lower=2> n_h;
  vector[n_t] log_e_t;
  vector[n_h] log_e_h;
  real prior_on_mean;                 // binarizer 1.1: per-gene prior on the
  real<lower=0> prior_on_sd;          //   means, N(midpoint, range / 3) of the
  real prior_off_mean;                //   gene's values (training fold for CV)
  real<lower=0> prior_off_sd;
}

parameters {
  positive_ordered[2] mu;             // OFF and ON means
  real<lower=0.001, upper=3> sd1;     // shared SD
  real<lower=0, upper=1> pi_on;       // ON weight
}

model {
  mu[1] ~ normal(prior_off_mean, prior_off_sd);
  mu[2] ~ normal(prior_on_mean, prior_on_sd);
  for (s in 1:n_t) {
    target += log_mix(pi_on, normal_lpdf(log_e_t[s] | mu[2], sd1),
                      normal_lpdf(log_e_t[s] | mu[1], sd1));
  }
}

generated quantities {
  vector[n_t] log_lik_t;
  vector[n_h] log_lik_h;
  for (s in 1:n_t) {
    log_lik_t[s] = log_mix(pi_on, normal_lpdf(log_e_t[s] | mu[2], sd1),
                           normal_lpdf(log_e_t[s] | mu[1], sd1));
  }
  for (s in 1:n_h) {
    log_lik_h[s] = log_mix(pi_on, normal_lpdf(log_e_h[s] | mu[2], sd1),
                           normal_lpdf(log_e_h[s] | mu[1], sd1));
  }
}
