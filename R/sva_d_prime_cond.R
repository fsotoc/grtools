sva_d_prime_cond <- function(y, pars_vector) {
  mu_s1 <- pars_vector[5] + (sqrt(pars_vector[15]) / sqrt(pars_vector[16]) * pars_vector[17]) * (y - pars_vector[6])
  var_s1 <- (1 - pars_vector[17]^2) * pars_vector[15]
  mu_s2 <- pars_vector[7] + (sqrt(pars_vector[18]) / sqrt(pars_vector[19]) * pars_vector[20]) * (y - pars_vector[8])
  var_s2 <- (1 - pars_vector[20]^2) * pars_vector[18]
  (mu_s2 - mu_s1) / sqrt((var_s1 + var_s2) / 2)
}
