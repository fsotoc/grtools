# stim_level includes levels for which individual RLNA is computed and averaged
sva_rlna_ind <- function(y, pars_vector, sub, model_type, stim_level = c(1, 2)) {
  if (model_type == 1) {
    indpar_n <- 6
  } else if (model_type == 2) {
    indpar_n <- 8
  } else {
    stop("Unknown model type")
  }

  kap <- pars_vector[20 + (sub - 1) * indpar_n + 1]
  lam1 <- pars_vector[20 + (sub - 1) * indpar_n + 2]
  lam2 <- 1 - lam1

  g_noise <- stats::dnorm(y, mean = pars_vector[2], sd = pars_vector[10] / (kap * lam2))
  out <- matrix(nrow = 2, ncol = length(y))
  for (i in stim_level) {
    out[i, ] <- g_noise / stats::dnorm(
      y,
      mean = pars_vector[6 + (i - 1) * 2],
      sd = pars_vector[16 + (i - 1) * 3] / (kap * lam2)
    )
  }
  colMeans(out, na.rm = TRUE)
}
