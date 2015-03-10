#' @export
print.sumstats_micro <- function(x){
  cat("Test of sampling independence:\n\n")
  x$sampling_independence$Expected_p <- round(as.numeric(x$sampling_independence$Expected_p),2)
  x$sampling_independence$Observed_p <- round(as.numeric(x$sampling_independence$Observed_p),2)
  x$sampling_independence$z <- round(as.numeric(x$sampling_independence$z),2)
  x$sampling_independence$p_value <- round(as.numeric(x$sampling_independence$p_value),5)
  names(x$sampling_independence) <- c("Stimulus", " Response", " Expected Prob", " Observed Prob",
                                      " Z stat", " P-value", " Pass?")
  print(x$sampling_independence, row.names=F)
  
  cat("\n\nTest of equal conditional d's:\n\n")
  x$equal_conditional_d_prime$dprime_hit <- round(as.numeric(x$equal_conditional_d_prime$dprime_hit),2)
  x$equal_conditional_d_prime$dprime_miss <- round(as.numeric(x$equal_conditional_d_prime$dprime_miss),2)
  x$equal_conditional_d_prime$z <- round(as.numeric(x$equal_conditional_d_prime$z),2)
  x$equal_conditional_d_prime$p_value <- round(as.numeric(x$equal_conditional_d_prime$p_value),5)
  names(x$equal_conditional_d_prime) <- c("Test", " d' Hits", " d' Miss", " Z stat", " P-value", " Pass?")
  print(x$equal_conditional_d_prime, row.names=F)
  
  cat("\n\nTest of equal conditional c:\n\n")
  x$equal_conditional_c$c_hit <- round(as.numeric(x$equal_conditional_c$c_hit),2)
  x$equal_conditional_c$c_miss <- round(as.numeric(x$equal_conditional_c$c_miss),2)
  x$equal_conditional_c$z <- round(as.numeric(x$equal_conditional_c$z),2)
  x$equal_conditional_c$p_value <- round(as.numeric(x$equal_conditional_c$p_value),5)
  names(x$equal_conditional_c) <- c("Test", " c Hits", " c Miss", " Z stat", " P-value", " Pass?")
  print(x$equal_conditional_c, row.names=F)
}