// One gene's log-normalized expression across pseudobulk samples as a
// mixture of an OFF and an ON Gaussian with a shared SD. The means are
// ordered so the components cannot swap. Generated quantities: the log
// probability that each sample, replicate and cluster is ON.
// Adapted from level1_ordered.stan by Fred P. Davis
// (https://github.com/fredpdavis/opticlobe; Davis et al. 2020, eLife
// 9:e50901), as used in mm_marker / mixture_modeling_smk.
//
// Copyright (c) 2018-2020 Fred P. Davis; adaptation (c) 2023-2026
// Yen-Chung Chen. Released under the MIT license with Fred P. Davis's
// permission; see inst/COPYRIGHTS.

data {
  int<lower=2> n_samples;
  int<lower=1> n_replicates;
  int<lower=1> n_clusters;
  vector[n_samples] log_e;            // expression of each sample
  array[n_samples] int replicate;     // replicate of each sample
  array[n_samples] int cluster;       // cluster of each sample
}

parameters {
  positive_ordered[2] mu;             // OFF and ON means
  real<lower=0.001, upper=3> sd1;     // shared SD
  real<lower=0, upper=1> pi_on;       // ON weight
}

transformed parameters {
  real mu_on = mu[2];
  real mu_off = mu[1];
}

model {
  mu ~ normal(7, 5);
  for (s in 1:n_samples) {
    target += log_mix(pi_on, normal_lpdf(log_e[s] | mu_on, sd1),
                      normal_lpdf(log_e[s] | mu_off, sd1));
  }
}

generated quantities {
  vector[n_samples] log_p_on_sample;
  vector[n_replicates] log_p_on_replicate;
  vector[n_clusters] log_p_on_cluster;
  // summed log likelihood of each replicate and cluster under each component
  vector[n_replicates] replicate_on = rep_vector(0, n_replicates);
  vector[n_replicates] replicate_off = rep_vector(0, n_replicates);
  vector[n_clusters] cluster_on = rep_vector(0, n_clusters);
  vector[n_clusters] cluster_off = rep_vector(0, n_clusters);

  for (s in 1:n_samples) {
    real lp_on = normal_lpdf(log_e[s] | mu_on, sd1);
    real lp_off = normal_lpdf(log_e[s] | mu_off, sd1);
    log_p_on_sample[s] = log(pi_on) + lp_on - log_mix(pi_on, lp_on, lp_off);
    replicate_on[replicate[s]] += lp_on;
    replicate_off[replicate[s]] += lp_off;
    cluster_on[cluster[s]] += lp_on;
    cluster_off[cluster[s]] += lp_off;
  }
  for (c in 1:n_clusters) {
    log_p_on_cluster[c] = log(pi_on) + cluster_on[c] -
                          log_mix(pi_on, cluster_on[c], cluster_off[c]);
  }
  for (d in 1:n_replicates) {
    log_p_on_replicate[d] = log(pi_on) + replicate_on[d] -
                            log_mix(pi_on, replicate_on[d], replicate_off[d]);
  }
}
