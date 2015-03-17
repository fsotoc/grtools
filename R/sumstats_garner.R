#' Perform a summary statistics analysis of data from a Garner filtering task.
#' 
#' Performs an analysis of data from a 2x2 Garner filtering experiment, based on
#' summary statistics (see Ashby & Maddox, 1994).
#' 
#' @param trial_data Trial-by-trial data from a single participant in a Garner 
#'   filtering experiment. See "Details" for instructions on the correct format
#'   for this data frame.
#' @param test_incorrect If TRUE, the function runs tests on response times for
#'   both correct and incorrect trials. The default is FALSE.
#'   
#' @return An object of class "\code{sumstats_garner}"
#'   
#'   The function \code{summary} is used to obtain a summary of conclusions from
#'   the analysis about separability. Note that any reported violations of 
#'   separability can be due to either violations of perceptual or decisional 
#'   separability, which cannot be dissociated through a Garner filtering 
#'   experiment.
#'   
#' @details A 2x2 Garner filtering experiment involves stimuli that vary in two 
#'   dimensions, each with two levels, 1 and 2. The task of the participant is 
#'   to classify stimuli according to their level in one of these dimensions 
#'   (the relevant dimension), while ignoring variation in the other dimension 
#'   (the irrelevant dimension).
#'   
#'   There are two block types in the task. During baseline blocks, the 
#'   irrelevant dimension is fixed to a specific level. During interference 
#'   blocks, the irrelevant dimension is not fixed, but varies across trials. 
#'   Both accuracy and response times are gathered during the task.
#'   
#'   The data from a single participant in this task should be ordered in a data
#'   frame with rows representing individual trials and columns with the 
#'   following format:
#'   
#'   \itemize{ \item{Column 1: block type, with a value of 1 for baseline blocks
#'   and a value of 2 for interference blocks} \item{Column 2: Level of the 
#'   relevant dimension, with values 1 and 2} \item{Column 3: Level of the 
#'   irrelevant dimension, with values 1 and 2} \item{Column 4: Accuracy, with a
#'   value of 0 for incorrect trials and 1 for correct trials} \item{Column 5: 
#'   Response time} }
#'   
#'   To see an example data frame, type \code{data(garner_data)} in the R
#'   console. The data will be available as a \code{data.frame} named
#'   \code{garner_data}
#'   
#' @references Ashby, F. G., & Maddox, W. T. (1994). A response time theory of 
#'   separability and integrality in speeded classification. \emph{Journal of 
#'   Mathematical Psychology, 38}(4), 423-466.
#' 
#' @examples
#' # Load example data frame and see the first 10 rows
#' data(garner_data)
#' garner_data[1:10,]
#'   
#' # Run the analysis
#' garner_results <- sumstats_garner(garner_data)
#'   
#' # See a summary of results
#' summary(garner_results)
#'   
#' # Print to screen the details of each test
#' garner_results
#'   
#' @export
sumstats_garner <- function(trial_data, test_incorrect=F){
  
  # garner interference test for accuracy
  p_base <- mean(trial_data[trial_data[,1]==1,4])
  n_base <- sum(trial_data[,1]==1)
  p_filt <- mean(trial_data[trial_data[,1]==2,4])
  n_filt <- sum(trial_data[,1]==2)
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
  m_base <- median(trial_data[trial_data[,1]==1 & trial_data[,4]==1,5])
  m_filt <- median(trial_data[trial_data[,1]==2 & trial_data[,4]==1,5])
  w_test <- wilcox.test(trial_data[trial_data[,1]==1 & trial_data[,4]==1,5],
                        trial_data[trial_data[,1]==2 & trial_data[,4]==1,5],
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
  trial_data <- trial_data[trial_data[1,]==2,2:5]
  
  # marginal response invariance
  # create an incomplete confusion matrix to use function mritest
  icm <- matrix(0, nrow=4, ncol=4)
  icm[1,1] <- sum(trial_data[trial_data[,1]==1 & trial_data[,2]==1,4])    # stimulus is A1B1 and response is a1
  icm[1,2] <- sum(trial_data[,1]==1 & trial_data[,2]==1) - icm[1,1] # stimulus is A1B1 and response is a2
  icm[2,2] <- sum(trial_data[trial_data[,1]==2 & trial_data[,2]==1,4])    # stimulus is A2B1 and response is a2
  icm[2,1] <- sum(trial_data[,1]==2 & trial_data[,2]==1) - icm[2,2] # stimulus is A2B1 and response is a1
  icm[3,3] <- sum(trial_data[trial_data[,1]==1 & trial_data[,2]==2,4])    # stimulus is A1B2 and response is a1
  icm[3,4] <- sum(trial_data[,1]==1 & trial_data[,2]==2) - icm[3,3] # stimulus is A1B2 and response is a2
  icm[4,4] <- sum(trial_data[trial_data[,1]==2 & trial_data[,2]==2,4])    # stimulus is A2B2 and response is a2
  icm[4,3] <- sum(trial_data[,1]==2 & trial_data[,2]==2) - icm[4,4] # stimulus is A2B2 and response is a1
  
  results$marginal_response_invariance <- mri_test(icm)[1:2,]
  results$marginal_response_invariance$Test <- c("Level 1 of relevant dimension",
                                                 "Level 2 of relevant dimension")
  
  # marginal RT invariance
  names(trial_data) <- c("relevant_level", "irrelevant_level", "accuracy", "rt")
  results$marginal_rt_invariance <- mrti_test(trial_data, test_incorrect=test_incorrect)
  
  # return object of class "sumstats_garner"
  class(results) <- "sumstats_garner"
  return(results)
}