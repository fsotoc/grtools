mrti_test <- function(dt, test_incorrect=F, test="ks"){
  # function to test marginal response time invariance, input is a dataframe 
  # with trials in rows and 4 columns: relevant_level, irrelevant_level, 
  # accuracy, rt by default, the function tests only rt distributions for 
  # correct trials, but if testincorrect=TRUE, then it also tests incorrect 
  # trials. The function uses the Kolmogorov-Smirnov test (adk test is not
  # supported anymore, so we are stuck with an insensitive test; argument "test"
  # is included to allow more powerful tests in the future)
  
  options(warn=-1)
  names(dt) <- c("relevant_level", "irrelevant_level", "accuracy", "rt")
  
  if (test=="ks") {
    
    # test rt distributions for correct trials
    # relevant level 1
    ksr <- ks.test(dt$rt[dt$relevant_level==1 & dt$irrelevant_level==1 & dt$accuracy==1],
                 dt$rt[dt$relevant_level==1 & dt$irrelevant_level==2 & dt$accuracy==1])
    
    if (ksr$p.value<.05){
      Pass="NO"
    } else {
      Pass="YES"
    }
    
    results <- data.frame(Distribution="Correct trials",
                        Test="Level 1 of Relevant Dimension",
                        KS_stat=ksr$statistic,
                        P_value=ksr$p.value,
                        Pass=Pass, stringsAsFactors=FALSE)
    
    # relevant level 2
    ksr <- ks.test(dt$rt[dt$relevant_level==2 & dt$irrelevant_level==1 & dt$accuracy==1],
                 dt$rt[dt$relevant_level==2 & dt$irrelevant_level==2 & dt$accuracy==1])
    
    if (ksr$p.value<.05){
      Pass="NO"
    } else {
      Pass="YES"
    }
    
    results<-rbind(results,c(Distribution="Correct trials",
                             Test="Level 2 of Relevant Dimension",
                             KS_stat=ksr$statistic,
                             P_value=ksr$p.value,
                             Pass=Pass))
    
    # if requested, test rt distributions for incorrect trials
    if (test_incorrect){
      # relevant level 1
      ksr <- ks.test(dt$rt[dt$relevant_level==1 & dt$irrelevant_level==1 & dt$accuracy==0],
                   dt$rt[dt$relevant_level==1 & dt$irrelevant_level==2 & dt$accuracy==0])
      
      if (ksr$p.value<.05){
        Pass="NO"
      } else {
        Pass="YES"
      }
      
      results<-rbind(results,c(Distribution="Incorrect trials",
                               Test="Level 1 of Relevant Dimension",
                               KS_stat=ksr$statistic,
                               P_value=ksr$p.value,
                               Pass=Pass))
      
      # relevant level 2
      ksr<-ks.test(dt$rt[dt$relevant_level==2 & dt$irrelevant_level==1 & dt$accuracy==0],
                   dt$rt[dt$relevant_level==2 & dt$irrelevant_level==2 & dt$accuracy==0])
      
      if (ksr$p.value<.05){
        Pass="NO"
      } else {
        Pass="YES"
      }
      
      results<-rbind(results,c(Distribution="Incorrect trials",
                               Test="Level 2 of Relevant Dimension",
                               KS_stat=ksr$statistic,
                               P_value=ksr$p.value,
                               Pass=Pass))
    }
    
    
    
    
  }

  
  return(results)
  
}