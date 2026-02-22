#' @export
print.sumstats_garner <- function(x, ...){
  
  cat("Garner interference test:\n\n")
  x$garner_interference$Baseline <- round(as.numeric(x$garner_interference$Baseline), 2)
  x$garner_interference$Filtering <- round(as.numeric(x$garner_interference$Filtering), 2)
  x$garner_interference$Difference <- round(as.numeric(x$garner_interference$Difference), 3)
  x$garner_interference$Stat <- round(as.numeric(x$garner_interference$Stat), 2)
  x$garner_interference$P_value <- round(as.numeric(x$garner_interference$P_value), 5)
  names(x$garner_interference) <- c("Test", " Baseline", " Filtering", " Difference",
                                    " Statistic", " P_value", " GI effect?")
  print(x$garner_interference, row.names=F)
  cat("_______________________________________________________________________\n")
  cat("  Note: Uses Wilcoxon test to compare RTs and a z-test to compare \n")
  cat("  proportions in the baseline and filtering conditions\n\n")
  
  cat("\nMarginal response invariance test:\n\n")
  x$marginal_response_invariance$z <- round(as.numeric(x$marginal_response_invariance$z), 2)
  x$marginal_response_invariance$p_value <- round(as.numeric(x$marginal_response_invariance$p_value), 5)
  names(x$marginal_response_invariance) <- c("Test", "  Z stat", "  P-value", "  Pass?")
  print(x$marginal_response_invariance, row.names=F)
  
  cat("\n\nMarginal response time invariance test:\n\n")
  x$marginal_rt_invariance$KS_stat <- round(as.numeric(x$marginal_rt_invariance$KS_stat), 4)
  x$marginal_rt_invariance$P_value <- round(as.numeric(x$marginal_rt_invariance$P_value), 5)
  names(x$marginal_rt_invariance) <- c("Distribution", "Test", " KS stat", " P-value", " Pass?")
  cat("          RT distributions from correct trials\n")
  print(x$marginal_rt_invariance[x$marginal_rt_invariance$Distribution=="Correct trials",2:5], row.names=F)
  if (  any(x$marginal_rt_invariance$Distribution=="Incorrect trials") ) {
    cat("\n          RT distributions from incorrect trials\n")
    print(x$marginal_rt_invariance[x$marginal_rt_invariance$Distribution=="Incorrect trials",2:5], row.names=F)
  }
  cat("_______________________________________________________________________\n")
  cat("  Note: Uses Kolmogorov-Smirnov test to compare RT distributions\n")
  cat("  across the two levels of the irrelevant dimension\n")
}
