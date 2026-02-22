# stim_level includes levels for which RLNA is computed and averaged
sva_rlna <- function(y, pars_vector, stim_level = c(1, 2)) {
  g_noise <- stats::dnorm(y, mean = pars_vector[2], sd = pars_vector[10])
  out <- matrix(nrow = 2, ncol = length(y))
  for (i in stim_level) {
    out[i, ] <- g_noise / stats::dnorm(
      y,
      mean = pars_vector[6 + (i - 1) * 2],
      sd = pars_vector[16 + (i - 1) * 3]
    )
  }
  colMeans(out, na.rm = TRUE)
}
