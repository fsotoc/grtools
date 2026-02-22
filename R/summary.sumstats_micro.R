#' @export
summary.sumstats_micro <- function(object, ...){
  micro <- object
  
  # summary for macroanalyses
  # code for conclusions
  conc <- c("yes", "NO", "?", "no?", "yes?", "NO", "YES")
  
  # check result for PS and DS (Table 1 in Kadlec, 1995)
  check_table = matrix(c(1, 1, 1, 1, 1,
                         0, 1, 1, 1, 1,
                         1, 1, 0, 5, 2,
                         1, 0, 1, 5, 2,
                         1, 0, 0, 5, 2,
                         0, 1, 0, 4, 4,
                         0, 0, 1, 4, 4,
                         0, 0, 0, 4, 4),
                       nrow=8, ncol=5, byrow=T)
  
  # A1B1 and A2B1 (A|B1)
  microvec <- c(prod(micro$sampling_independence$Pass[1:8]=="YES"),
                (micro$equal_conditional_d_prime$Pass[1]=="YES")*1,
                (micro$equal_conditional_c$Pass[1]=="YES")*1)
  microvec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(microvec)),]
  
  # put in table
  micro_summary <- data.frame(Stimulus_Pair="A|B1",
                              Sampling_Independence=conc[microvec[1]+6],
                              Conditional_d=conc[microvec[2]+6],
                              Conditional_c=conc[microvec[3]+6],
                              PI=conc[microvec[4]],
                              DS=conc[microvec[5]])
  
  # A1B2 and A2B2 (A|B2)
  microvec <- c(prod(micro$sampling_independence$Pass[9:16]=="YES"),
                (micro$equal_conditional_d_prime$Pass[2]=="YES")*1,
                (micro$equal_conditional_c$Pass[2]=="YES")*1)
  microvec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(microvec)),]
  
  # put in table
  micro_summary <- rbind(micro_summary, data.frame(Stimulus_Pair="A|B2",
                                                   Sampling_Independence=conc[microvec[1]+6],
                                                   Conditional_d=conc[microvec[2]+6],
                                                   Conditional_c=conc[microvec[3]+6],
                                                   PI=conc[microvec[4]],
                                                   DS=conc[microvec[5]]))
  
  # A1B1 and A1B2 (B|A1)
  microvec <- c(prod(micro$sampling_independence$Pass[c(1:4,9:12)]=="YES"),
                (micro$equal_conditional_d_prime$Pass[3]=="YES")*1,
                (micro$equal_conditional_c$Pass[3]=="YES")*1)
  microvec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(microvec)),]
  
  # put in table
  micro_summary <- rbind(micro_summary, data.frame(Stimulus_Pair="B|A1",
                                                   Sampling_Independence=conc[microvec[1]+6],
                                                   Conditional_d=conc[microvec[2]+6],
                                                   Conditional_c=conc[microvec[3]+6],
                                                   PI=conc[microvec[4]],
                                                   DS=conc[microvec[5]]))
  
  # A2B1 and A2B2 (B|A2)
  microvec <- c(prod(micro$sampling_independence$Pass[c(5:8,13:16)]=="YES"),
                (micro$equal_conditional_d_prime$Pass[4]=="YES")*1,
                (micro$equal_conditional_c$Pass[4]=="YES")*1)
  microvec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(microvec)),]
  
  # put in table
  micro_summary <- rbind(micro_summary, data.frame(Stimulus_Pair="B|A2",
                                                   Sampling_Independence=conc[microvec[1]+6],
                                                   Conditional_d=conc[microvec[2]+6],
                                                   Conditional_c=conc[microvec[3]+6],
                                                   PI=conc[microvec[4]],
                                                   DS=conc[microvec[5]]))
  
  # print results
  names(micro_summary) <- c("Stimulus pair", "  Sampling Independence", "  Equal Cond d'", "  Equal Cond c", "  PI", "  DS")
  print(micro_summary, row.names=F)
}
