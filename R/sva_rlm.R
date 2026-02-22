# stim_level includes levels for which RLM is computed and averaged
sva_rlm <- function(y, pars_vector, stim_level = c(1, 2)) {
  g_noise <- stats::dnorm(y, mean = pars_vector[2], sd = pars_vector[10])
  out <- matrix(nrow = 2, ncol = length(y))
  for (i in stim_level) {
    out[i, ] <- stats::dnorm(
      y,
      mean = pars_vector[6 + (i - 1) * 2],
      sd = pars_vector[16 + (i - 1) * 3]
    ) / g_noise
  }
  colMeans(out, na.rm = TRUE)
}

sva_safe_log2 <- function(x, use_log_scale) {
  if (!isTRUE(use_log_scale)) {
    return(x)
  }
  log2(x)
}
