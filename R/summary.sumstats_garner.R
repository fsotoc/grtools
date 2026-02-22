#' @export
summary.sumstats_garner <- function(object, ...){
  ssg_list <- object
  
  # summary
  sep <- c("NO", "NO")
  sep[ssg_list$garner_interference$GI_effect=="NO"] <- "yes?"
  sep <- c(sep, "NO", "NO")
  if (all(ssg_list$marginal_response_invariance$Pass=="YES")){
    sep[3] <- "yes?"
  }
  if (all(ssg_list$marginal_RT_invariance$Pass=="YES")){
    sep[4] <- "yes?"
  }
  
  ssg_list_summary <- data.frame(Test = c("Garner Interference - Accuracy",
                                         "Garner Interference - RT",
                                         "Marginal Response Invariance",
                                         "Marginal RT Invariance"),
                                Pass = sep)
  
  # print summary table
  names(ssg_list_summary) <- c("Separability Test", "  Pass?")
  print(ssg_list_summary, row.names=F)
  
}
