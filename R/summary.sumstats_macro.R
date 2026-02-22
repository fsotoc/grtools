#' @export
summary.sumstats_macro <- function(object, ...){
  macro <- object
  # summary for macroanalyses
  # code for conclusions
  conc <- c("yes", "NO", "?", "no?", "yes?", "NO", "YES")
  
  # check result for PS and DS (Table 1 in Kadlec, 1995)
  check_table = matrix(c(1, 1, 1, 1, 1,
                         1, 1, 0, 1, 2,
                         1, 0, 1, 2, 1,
                         1, 0, 0, 2, 2,
                         0, 1, 1, 1, 4,
                         0, 1, 0, 1, 2,
                         0, 0, 1, 2, 3,
                         0 ,0 ,0 ,2 ,3),
                       nrow=8, ncol=5, byrow=T)
  
  # Dimension A
  macrovec <- c(prod(macro$marginal_response_invariance$Pass[1:2]=="YES"),
                (macro$equal_marginal_d_prime$Pass[1]=="YES")*1,
                (macro$equal_marginal_c$Pass[1]=="YES")*1)
  macrovec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(macrovec)),]
  
  # put in table
  macro_summary <- data.frame(Dimension="A", MRI=conc[macrovec[1]+6],
                              Marginal_d=conc[macrovec[2]+6],
                              Marginal_C=conc[macrovec[3]+6],
                              PS=conc[macrovec[4]],
                              DS=conc[macrovec[5]])
  
  # Dimension B
  macrovec <- c(prod(macro$marginal_response_invariance$Pass[3:4]=="YES"),
                (macro$equal_marginal_d_prime$Pass[2]=="YES")*1,
                (macro$equal_marginal_c$Pass[2]=="YES")*1)
  macrovec <- check_table[which(data.frame(t(check_table[,1:3])) %in% data.frame(macrovec)),]
  
  # put in table
  macro_summary <- rbind(macro_summary, data.frame(Dimension="B",MRI=conc[macrovec[1]+6],
                                                   Marginal_d=conc[macrovec[2]+6],
                                                   Marginal_C=conc[macrovec[3]+6],
                                                   PS=conc[macrovec[4]],
                                                   DS=conc[macrovec[5]]))
  
  # print results
  names(macro_summary) <- c("Dimension", "MRI", "Marginal d'", "Marginal c", "PS", "DS")
  print(macro_summary, row.names=F)
}
