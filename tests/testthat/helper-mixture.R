# One gene's fit as fit_mixture() returns it, without running Stan. The
# draws repeat the given parameters in 2 chains x 5 draws, so Ashman's D is
# (on_level - off_level) / sd at every draw; p_on_mean and chain_spread are
# per cluster (by default the median P(On) and no spread)
fake_gene_fit <- function(log_p_on, diff, diff_se, on_level = 1.23456789,
                          off_level = 0.0123456, sd = 0.00312, pi = 0.987654,
                          p_on_mean = exp(log_p_on),
                          chain_spread = rep(0, length(log_p_on))) {
  draws <- array(rep(c(on_level, off_level, sd, pi), each = 10), dim = c(5, 2, 4))
  list(
    bimodal = list(log_p_on_cluster = log_p_on, on_level = on_level,
                   off_level = off_level, sd = sd, pi = pi, draws = draws,
                   p_on_mean = p_on_mean, chain_spread = chain_spread),
    unimodal = list(level = 0.456, sd = 0.789),
    elpd = list(diff = diff, diff_se = diff_se)
  )
}

# Binarizer 1.0 record of the same fit: no draws or summaries over them
fake_gene_fit_v10 <- function(...) {
  g <- fake_gene_fit(...)
  g$bimodal[c("draws", "p_on_mean", "chain_spread")] <- NULL
  g
}
