#NOTATION
#s_0_C: prior belief about successes on treatment C before Phase II.
#f_0_C: prior belief about failures on treatment C before Phase II.
#s_0_D: prior belief about successes on treatment D before Phase II.
#f_0_D: prior belief about failures on treatment D before Phase II.
#cost.per.patient.PhaseII: cost of recruiting a patient at Phase II.
#cost.per.patient.PhaseIII: cost of recruiting a patient at Phase III.
#Pop: Population benefiting from treatments.
#n: patients in Phase III.
#R: Revenue from the new treatment (per patient).
#alpha: significance level for Phase II.
#beta: type II error for Phase II
#k: minimal detectable difference
#n_1: Phase II sample size



library(hash)
library(Exact)
library(VGAM)
library(Matrix)
library(tidyverse)
library(tictoc)
##################################################################################################
#Creating every possible outcome for the chosen Phase II design
Outcomes_calculator <- function(n_1, alpha, beta, k, s_0_D, f_0_D, s_0_C, f_0_C){
  #Dividing the sample size with 1:1 allocation
  n1_divided = ceiling(n_1/2)
  
  #Creating every possible combination of successes and failures for treatment D
  s_1D = seq(0,n1_divided)
  f_1D = n1_divided - s_1D
  
  #Creating every possible combination of successes and failures for treatment C
  s_1C = rev(s_1D)
  f_1C = n1_divided - s_1C
  
  # Create the keys for both D and C
  key_D = paste(s_1D, ",", f_1D)
  key_C = paste(s_1C, ",", f_1C)
  
  PhaseII.succ_fail_sample.size.combinations = 
    expand.grid(key_D, key_C)%>% 
    select(Var1, Var2) %>%
    mutate(prop_successes_C = (as.numeric(Var2)-1 + s_0_C)/(n1_divided + s_0_C + f_0_C),
           prop_successes_D = prop_successes_C + k) %>%
    mutate(n = pmax(20, ceiling((qnorm(1-alpha) + qnorm(1-beta))^2 *
                                  (prop_successes_D  * (1 - prop_successes_D )
                                   + prop_successes_C * (1 - prop_successes_C)) / k^2)))
  return(PhaseII.succ_fail_sample.size.combinations)
}


#Calculating the term (theta_D - theta_C)
proportions_difference_patient<- function(s_2D, s_1_D, f_1_D, s_1_C, f_1_C, s_0_D, f_0_D, s_0_C, f_0_C, n){
  #Updating successes and failures for both treatments
  successes_D = (s_2D + s_1_D + s_0_D)
  successes_C = (s_2D + s_1_C + s_0_C)
  failures_D = (ceiling(n/2) - s_2D + f_0_D + f_1_D)
  failures_C = (ceiling(n/2) - s_2D + f_0_C + f_1_C)
  
  
  #Calculating proportion of success for both treatments
  proportions_D <- pmin(1,successes_D/(successes_D + failures_D))
  proportions_C <- pmin(1,successes_C/(successes_C + failures_C))
  
  #Outer product of subtraction of proportions for each treatment 
  diff_matrix <- outer(proportions_D, proportions_C, FUN = function(x, y) (x - y))
  
  return(diff_matrix)
}



Optimal_patient_final_version <- function(R,Pop,n_1,alpha,beta,k, cost.per.patient.PhaseII,cost.per.patient.PhaseIII, s_0_D, f_0_D,s_0_C, f_0_C, theta_D, theta_C){
  
  #PhaseII.succ_fail_sample.size.combinations has the following columns
  #(s_1D,f_1D), (s_1C,f_1C),  prop_successes_D,  prop_successes_C,  n
  PhaseII.succ_fail_sample.size.combinations = Outcomes_calculator(n_1, alpha, beta, k, s_0_D, f_0_D, s_0_C, f_0_C)
  
  #Here we populate Reward at t = 2 (only if we had go decision)
  dictionary_of_decision = hash()
  for (i in 1:nrow(PhaseII.succ_fail_sample.size.combinations)) {
    #Extract s_1D
    s_1D =  as.integer(unlist(strsplit(as.character(PhaseII.succ_fail_sample.size.combinations$Var1[i]), ","))[1])
    #Extract f_1D
    f_1D = as.integer(unlist(strsplit(as.character(PhaseII.succ_fail_sample.size.combinations$Var1[i]), ","))[2])
    
    #Extract s_1C
    s_1C =  as.integer(unlist(strsplit(as.character(PhaseII.succ_fail_sample.size.combinations$Var2[i]), ","))[1])
    #Extract f_1C
    f_1C = as.integer(unlist(strsplit(as.character(PhaseII.succ_fail_sample.size.combinations$Var2[i]), ","))[2])
    
    #Extract n
    n_onearm =  PhaseII.succ_fail_sample.size.combinations$n[i]
    
    n = 2*n_onearm
    #Calculate the reject H0 region using z unpooled test
    reject_region = exact.reject.region(ceiling(n/2),ceiling(n/2),
                                        "greater", alpha = alpha, method = "z-unpooled")
    
    #Find every possible combination of successes for both treatments for given Phase II result
    s_2D = 0:(ceiling(n/2))
    s_2C = 0:(ceiling(n/2))
    
    #Calculate probability of getting to s_2C and s_2D, to use in expected reward 
    prob_C_III = dbinom(s_2C, ceiling(n/2), theta_C)
    prob_D_III = dbinom(s_2D, ceiling(n/2), theta_D)
    
    #Putting everything together 
    reward = prob_D_III %*% t(prob_C_III) * (reject_region * Pop * (theta_D - theta_C))  
    #Summing everything gives the expected reward at t = 2
    #Subtracting the cost for recruiting n patients gives the expected reward at t = 1 for a go decision
    Cumulative_Reward = sum(reward) - ceiling(n/2)*(theta_D - theta_C)
    
    expected_prob_C_III = dbetabinom.ab(s_2C, ceiling(n/2), s_1C+s_0_C,f_1C+f_0_C)
    expected_prob_D_III = dbetabinom.ab(s_2D,ceiling(n/2), s_1D+s_0_D,f_1D+f_0_D)

    prop_mat_after = proportions_difference_patient(s_2D, s_1D, f_1D, s_1C, f_1C, s_0_D, f_0_D, s_0_C, f_0_C, n)
    
    Observed_reward = expected_prob_C_III %*% t(expected_prob_D_III) * (reject_region * Pop * prop_mat_after)  
    prop_mat_before = proportions_difference_patient(0, s_1D, f_1D, s_1C, f_1C, s_0_D, f_0_D, s_0_C, f_0_C, n)
    Observed_Cumulative_Reward = sum(Observed_reward) - ceiling(n/2)*prop_mat_before
    
    #If the reward is greater than 0, the optimal decision is a go decision, otherwise we have a no-go decision
    if (Observed_Cumulative_Reward >= 0) {
      denom = dim(reject_region)[1]*dim(reject_region)[2]
      numer = sum(reject_region)
      s1.f1.key = paste(s_1D,",",f_1D,",",s_1C,",",f_1C,",",n,",",numer/denom ,",",Cumulative_Reward)
      dictionary_of_decision[[s1.f1.key]] = 2
    } else {
      denom = dim(reject_region)[1]*dim(reject_region)[2]
      numer = sum(reject_region)
      s1.f1.key = paste(s_1D,",",f_1D,",",s_1C,",",f_1C,",",n,",",numer/denom ,",",0)
      dictionary_of_decision[[s1.f1.key]] = 1
    }
  } 
  #The output of dictionary_of_decision is a dictionary
  #where the key is (s_1D,f_1D,s_1C,f_1C, n, expected reward), 
  #for every possible state vector value(s_1D,f_1D,s_1C,f_1C)
  #and the value will be either 1 (no go) or 2 (go)
  
  #The goal here (t = 0) is to calculate the Phase II trial expected total reward
  
  Final.Expected.Reward = 0
  
  #There are used to find extra metrics
  prop_go = 0
  sample_sizes = 0
  #For every (s1,f1,reward)
  for (s1.f1.reward in ls(dictionary_of_decision)){
    keysplit = strsplit(s1.f1.reward, ",")
    
    #Extract s_1D
    s_1D = as.numeric(keysplit[[1]][1])
    #Extract f_1D
    f_1D = as.numeric(keysplit[[1]][2])
    
    #Extract s_1C
    s_1C = as.numeric(keysplit[[1]][3])
    #Extract f_1C
    f_1C = as.numeric(keysplit[[1]][4])
    
    #Extract n
    n = as.numeric(keysplit[[1]][5])
    
    #props = c(props, as.numeric(keysplit[[1]][6]))
    #Extract expected reward from t = 1
    reward = as.numeric(keysplit[[1]][7])
    
    #Find Phase II sample size
    T_a_D = s_1D + f_1D
    T_a_C = s_1C + f_1C
    n_a = T_a_D + T_a_C
    
    #Calculate probability of getting to s_1D and s_1C, to use in expected reward 
    prob_D_II = dbinom(s_1D, T_a_D, theta_D)
    prob_C_II = dbinom(s_1C, T_a_C, theta_C)
    
    
    decision = (dictionary_of_decision[[s1.f1.reward]] - 1)
    prop_go = prop_go + prob_D_II*prob_C_II*decision
    sample_sizes = sample_sizes +  prob_D_II*prob_C_II*n 
    
    #Summing every reward at t = 1 and multiplying by the corresponding probability to get expected reward
    Final.Expected.Reward = Final.Expected.Reward + prob_D_II*prob_C_II*reward 
  }
  #Subtracting the cost of recruiting n_a patients from the expected reward to get the final expected reward
  
  Final.Expected.Reward = Final.Expected.Reward - n_a/2*(theta_D - theta_C)
  
  prop_nogo = 1 - prop_go 
  return(list(dictionary_of_decision = dictionary_of_decision,
              prop_of_nogo = prop_nogo,
              prop_of_go = prop_go,
              average_size_III = sample_sizes,
              Final.Expected.Reward = Final.Expected.Reward))
}


