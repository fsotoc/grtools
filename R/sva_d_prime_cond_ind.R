sva_d_prime_cond_ind <- function(y, pars_vector, sub, model_type) {
  if (model_type == 1) {
    indpar_n <- 6
  } else if (model_type == 2) {
    indpar_n <- 8
  } else {
    stop("Unknown model type")
  }

  kap <- pars_vector[20 + (sub - 1) * indpar_n + 1]
  lam1 <- pars_vector[20 + (sub - 1) * indpar_n + 2]

  mu_s1 <- pars_vector[5] + (sqrt(pars_vector[15]) / sqrt(pars_vector[16]) * pars_vector[17]) * (y - pars_vector[6])
  var_s1 <- (1 - pars_vector[17]^2) * pars_vector[15] / (kap * lam1)
  mu_s2 <- pars_vector[7] + (sqrt(pars_vector[18]) / sqrt(pars_vector[19]) * pars_vector[20]) * (y - pars_vector[8])
  var_s2 <- (1 - pars_vector[20]^2) * pars_vector[18] / (kap * lam1)
  (mu_s2 - mu_s1) / sqrt((var_s1 + var_s2) / 2)
}
