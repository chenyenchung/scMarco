// unimodal.stan fit on the training samples (_t), with the log likelihood
// of the held-out samples (_h) for cross-validation.
// Adapted from level1_unimodal_CV.stan by Fred P. Davis
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
  real<lower=0, upper=14> mu1;
  real<lower=0.001, upper=3> sd1;
}

model {
  for (s in 1:n_t) {
    target += normal_lpdf(log_e_t[s] | mu1, sd1);
  }
}

generated quantities {
  vector[n_t] log_lik_t;
  vector[n_h] log_lik_h;
  for (s in 1:n_t) {
    log_lik_t[s] = normal_lpdf(log_e_t[s] | mu1, sd1);
  }
  for (s in 1:n_h) {
    log_lik_h[s] = normal_lpdf(log_e_h[s] | mu1, sd1);
  }
}
