
# data should be a data frame with the following columns:
# col 1: block type, 1=baseline, 2=interference
# col 2: level of relevant dimension, 1 and 2
# col 3: level of irrelevant dimension, 1 and 2
# col 4: accuracy, 0=incorrect, 1=correct
# col 5: response time
#' @export
sumstats_garner <- function(data, test_incorrect=F){
  
  # garner interference test for accuracy
  p_base <- mean(data[data[,1]==1,4])
  n_base <- sum(data[,1]==1)
  p_filt <- mean(data[data[,1]==2,4])
  n_filt <- sum(data[,1]==2)
  ptest <- testp(p_base, p_filt, n_base, n_filt)
  if (ptest$p_value<.05){
    gi_effect <- "YES"
  } else {
    gi_effect <- "NO"
  }
  results <- list()
  results$garner_interference <- data.frame(Test="Proportion Correct", 
                                            Baseline=p_base,
                                            Filtering=p_filt,
                                            Difference=p_filt-p_base,
                                            Stat=ptest$z,
                                            P_value=ptest$p_value,
                                            GI_effect=gi_effect)
  
  # garner interference effect for correct RT
  m_base <- median(data[data[,1]==1 & data[,4]==1,5])
  m_filt <- median(data[data[,1]==2 & data[,4]==1,5])
  w_test <- wilcox.test(data[data[,1]==1 & data[,4]==1,5],
                        data[data[,1]==2 & data[,4]==1,5],
                        alternative="less")
  if (w_test$p.value<.05){
    gi_effect = "YES"
  } else {
    gi_effect = "NO"
  }
  results$garner_interference <- rbind(results$garner_interference,
                                       data.frame(Test="Median correct RT", 
                                            Baseline=m_base,
                                            Filtering=m_filt,
                                            Difference=m_filt-m_base,
                                            Stat=w_test$statistic,
                                            P_value=w_test$p.value,
                                            GI_effect=gi_effect)  )
  
  # we only use interference trials for invariance tests
  data <- data[data[1,]==2,2:5]
  
  # marginal response invariance
  # create an incomplete confusion matrix to use function mritest
  icm <- matrix(0, nrow=4, ncol=4)
  icm[1,1] <- sum(data[data[,1]==1 & data[,2]==1,4])    # stimulus is A1B1 and response is a1
  icm[1,2] <- sum(data[,1]==1 & data[,2]==1) - icm[1,1] # stimulus is A1B1 and response is a2
  icm[2,2] <- sum(data[data[,1]==2 & data[,2]==1,4])    # stimulus is A2B1 and response is a2
  icm[2,1] <- sum(data[,1]==2 & data[,2]==1) - icm[2,2] # stimulus is A2B1 and response is a1
  icm[3,3] <- sum(data[data[,1]==1 & data[,2]==2,4])    # stimulus is A1B2 and response is a1
  icm[3,4] <- sum(data[,1]==1 & data[,2]==2) - icm[3,3] # stimulus is A1B2 and response is a2
  icm[4,4] <- sum(data[data[,1]==2 & data[,2]==2,4])    # stimulus is A2B2 and response is a2
  icm[4,3] <- sum(data[,1]==2 & data[,2]==2) - icm[4,4] # stimulus is A2B2 and response is a1
  
  results$marginal_response_invariance <- mri_test(icm)[1:2,]
  results$marginal_response_invariance$Test <- c("Level 1 of relevant dimension",
                                                 "Level 2 of relevant dimension")
  
  # marginal RT invariance
  names(data) <- c("relevant_level", "irrelevant_level", "accuracy", "rt")
  results$marginal_rt_invariance <- mrti_test(data, test_incorrect=test_incorrect)
  
  # return object of class "sumstats_garner"
  class(results) <- "sumstats_garner"
  return(results)
}