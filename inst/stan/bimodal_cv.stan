// bimodal.stan fit on the training samples (_t), with the log likelihood
// of the held-out samples (_h) for cross-validation.
// Adapted from level1_ordered_CV.stan by Fred P. Davis
// (https://github.com/fredpdavis/opticlobe; Davis et al. 2020, eLife
// 9:e50901), as used in mm_marker / mixture_modeling_smk.
//
// Copyright (c) 2018-2020 Fred P. Davis; adaptation (c) 2023-2026
// Yen-Chung Chen. Released under the MIT license with Fred P. Davis's
// permission; see inst/COPYRIGHTS.

data {
  int<lower=2> n_t;
  int<lower=2> n_h;
  vector[n_t] log_e_t;
  vector[n_h] log_e_h;
}

parameters {
  positive_ordered[2] mu;             // OFF and ON means
  real<lower=0.001, upper=3> sd1;     // shared SD
  real<lower=0, upper=1> pi_on;       // ON weight
}

model {
  mu ~ normal(7, 5);
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
