// One gene's log-normalized expression across pseudobulk samples as a
// single Gaussian.
// Adapted from level1_unimodal.stan by Fred P. Davis
// (https://github.com/fredpdavis/opticlobe; Davis et al. 2020, eLife
// 9:e50901), as used in mm_marker / mixture_modeling_smk.
//
// Copyright (c) 2018-2020 Fred P. Davis; adaptation (c) 2023-2026
// Yen-Chung Chen. Released under the MIT license with Fred P. Davis's
// permission; see inst/COPYRIGHTS.

data {
  int<lower=2> n_samples;
  vector[n_samples] log_e;
}

parameters {
  real<lower=0, upper=14> mu1;
  real<lower=0.001, upper=3> sd1;
}

model {
  for (s in 1:n_samples) {
    target += normal_lpdf(log_e[s] | mu1, sd1);
  }
}
