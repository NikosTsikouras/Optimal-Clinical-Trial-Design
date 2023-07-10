#NOTATION
# s_0_C: prior belief about successes on treatment C before Phase II.
# f_0_C: prior belief about failures on treatment C before Phase II.
# s_0_D: prior belief about successes on treatment D before Phase II.
# f_0_D: prior belief about failures on treatment D before Phase II.
# cost.per.patient.PhaseII: cost of recruiting a patient at Phase II.
# cost.per.patient.PhaseIII: cost of recruiting a patient at Phase III.
# Pop: Population benefiting from treatments.
# n: patients in Phase III.
# R: Revenue from the new treatment (per patient).
# alpha: significance level for Phase II.
# beta: type II error for Phase II
# k: minimal detectable difference
# n_1: Phase II sample size




library(hash)
library(Exact)
library(VGAM)
library(Matrix)
library(tidyverse)
library(tictoc)

##################################################################################################

#Calculating the term (theta_D - theta_C)
proportions_difference_patient<- function(s_2D, s_1_D, f_1_D, s_1_C, f_1_C, s_0_D, f_0_D, s_0_C, f_0_C, n){
  #Updating successes and failures for both treatments
  successes_D = (s_2D + s_1_D + s_0_D)
  successes_C = (s_2D + s_1_C + s_0_C)
  failures_D = (ceiling(n/2) - s_2D + f_0_D)
  failures_C = (ceiling(n/2) - s_2D + f_0_C)
  
  
  #Calculating proportion of success for both treatments
  proportions_D <- pmin(1,successes_D/(successes_D + failures_D))
  proportions_C <- pmin(1,successes_C/(successes_C + failures_C))
  
  #Outer product of subtraction of proportions for each treatment 
  diff_matrix <- outer(proportions_D, proportions_C, FUN = function(x, y) (x - y))
  
  return(diff_matrix)
}



#Estimating the probability of adoption
prop_adop <- function(s_2D, s_2C, s_0_D,f_0_D, s_1D, f_1D, s_0_C, f_0_C,s_1C,f_1C, n){
  #Updating the parameters for the two beta distributions
  alpha1 = s_2D  + s_0_D + s_1D
  beta1 = ceiling(n/2) - s_2D + f_0_D + f_1D
  
  alpha2 = s_2C+ s_0_C + s_1C
  beta2 = ceiling(n/2) - s_2C + f_0_C + f_1C
  
  #Creating the parameters for the normal approximation
  mux = apply(cbind(alpha1,beta1),1,FUN= function(i)i[1]/(i[1]+i[2]))
  sigmax = apply(cbind(alpha1,beta1),1,FUN = function(i)i[1]*i[2]/((i[1]+i[2])^2*(i[1]+i[2]+1)))
  
  #Adding 0.2 for the probability of adoption
  muy = apply(cbind(alpha2,beta2),1,FUN= function(i)i[1]/(i[1]+i[2])) + 0.2
  sigmay = apply(cbind(alpha2,beta2),1,FUN = function(i)i[1]*i[2]/((i[1]+i[2])^2*(i[1]+i[2]+1)))
  
  mat1 <- cbind(mux, sigmax)
  mat2 <- cbind(muy, sigmay)
  
  
  # Compute the numerator and denominator terms for each element
  numerator <- outer(mat1[, 1], mat2[, 1], FUN = "-")
  denominator <- outer(mat1[, 2], mat2[, 2], FUN = "+")^(1/2)
  
  # Compute the final matrix by applying pnorm to the element-wise division
  mat <- pnorm(numerator / denominator)
  return(mat)
}

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




Optimal_pharma_frequentist <- function(R,Pop,n_1,alpha,beta,k, cost.per.patient.PhaseII,cost.per.patient.PhaseIII, s_0_D, f_0_D,s_0_C, f_0_C, theta_C, theta_D){
  
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
    
    #Extract n for one arm
    n_onearm =  PhaseII.succ_fail_sample.size.combinations$n[i]
    
    #Double the value to get the total Phase III sample size
    n = 2*n_onearm
    
    #Calculate the reject H0 region using z unpooled test
    reject_region = exact.reject.region(ceiling(n/2),ceiling(n/2),
                                        "greater", alpha = alpha, method = "z-unpooled")
    
    #Find every possible combination of successes for both treatments for given Phase II result
    s_2D = 0:(ceiling(n/2))
    s_2C = 0:(ceiling(n/2))
    
    #Find Phase II sample size
    T_a_D = s_1D + f_1D
    T_a_C = s_1C + f_1C
    n_a = T_a_D + T_a_C
    
    #Calculate frequentist probability of getting to s_1D and s_1C, to use in expected frequentist reward 
    prob_D_II = dbinom(s_1D, T_a_D, theta_D)
    prob_C_II = dbinom(s_1C, T_a_C, theta_C)
    
    #Calculate frequentist probability of getting to s_2C and s_2D, to use in expected frequentist reward 
    prob_C_III = dbinom(s_2C, ceiling(n/2), theta_C)
    prob_D_III = dbinom(s_2D, ceiling(n/2), theta_D)
    
    #(theta_D - theta_C) term 
    #Notice how we are using every possible piece of information (Prior, Phase II and Phase III)
    #for this as we are in the end of Phase III
    prop_mat = proportions_difference_patient(s_2D, s_1D, f_1D, s_1C, f_1C, s_0_D, f_0_D, s_0_C, f_0_C, n)
    
    #probability of adoption term
    #Notice how we are using every possible piece of information (Prior, Phase II and Phase III)
    #for this as we are in the end of Phase III
    prop_adoption = prop_adop(s_2D, s_2C, s_0_D, f_0_D, s_1D, f_1D, s_0_C, f_0_C,s_1C,f_1C,n)
    
    #Putting everything together to get the reward at t = 2
    reward = prob_D_III %*% t(prob_C_III) * (reject_region  * R * Pop * prop_mat* prop_adoption) 
    
    #Summing everything gives the expected reward at t = 2
    #Subtracting the cost for recruiting n patients gives the expected reward at t = 1 for a go decision
    #Here we use the probability of getting to the specific state of t = 1 
    #(i.e. the corresponding s_1D, f_1D, s_1C,f_1C) given prior values
    #Andrew and Peter I am emphasizing this as I want your opinion here
    Expected_Frequentist_Cumulative_Reward = prob_D_II * prob_C_II* (sum(reward) - n * cost.per.patient.PhaseIII)
    
    #Below we are doing the same but for the Observed expected reward, which will be used for taking the decision
    #That is we use the Bayesian probabilities
    expected_prob_C_III = dbetabinom.ab(s_2C, ceiling(n/2), s_1C+s_0_C,f_1C+f_0_C)
    expected_prob_D_III = dbetabinom.ab(s_2D,ceiling(n/2), s_1D+s_0_D,f_1D+f_0_D)
    
    expected_prob_C_II = dbetabinom.ab(s_1C, s_1C + f_1C, s_0_C,f_0_C)
    expected_prob_D_II = dbetabinom.ab(s_1D, s_1D + f_1D, s_0_D,f_0_D)
    
    Observed_reward = expected_prob_D_III %*% t(expected_prob_C_III) * (reject_region * R * Pop * prop_mat* prop_adoption) 
    
    Expected_Observed_Cumulative_Reward = expected_prob_D_II %*% t(expected_prob_C_II)*(sum(Observed_reward) -  n * cost.per.patient.PhaseIII)
    
    #Putting everything together 
    #If the Expected Observed Reward is greater than 0, the optimal decision is a go decision
    #with reward being the reward based on the true probabilities
    #otherwise we have a no-go decision and the reward is 0
    if (Expected_Observed_Cumulative_Reward >= 0) {
      #denom = dim(reject_region)[1]*dim(reject_region)[2]
      #numer = sum(reject_region)
      s1.f1.key = paste(s_1D,",",f_1D,",",s_1C,",",f_1C,",",n,",",",",Expected_Observed_Cumulative_Reward,",",Expected_Frequentist_Cumulative_Reward)
      dictionary_of_decision[[s1.f1.key]] = 2
    } else {
      #denom = dim(reject_region)[1]*dim(reject_region)[2]
      #numer = sum(reject_region)
      s1.f1.key = paste(s_1D,",",f_1D,",",s_1C,",",f_1C,",",n,",",",",0,",",0)
      dictionary_of_decision[[s1.f1.key]] = 1
    }
  } 
  
  
  #The output of dictionary_of_decision is a dictionary
  #where the key is (s_1D,f_1D,s_1C,f_1C, n, Reward based on Bayesian probability, Reward based on true probability), 
  #and the value will be either 1 (no go) or 2 (go)
  
  #We go to t = 0  to calculate the Phase II trial expected total reward
  
  #We will calculate both True Reward and Bayesian Reward
  #I am unsure which one is useful. Probably the Observed one is more important as if its below 0 we are 
  #probably less inclined to go through with the trial
  #Andrew and Peter please see this
  Frequentist_Final.Expected.Reward = 0
  Observed_Final.Expected.Reward = 0
  
  #Here we also calculate the proportion of no go decisions (and implicitly the go decisions) and the sample size
  Observed_prop_go = 0
  Observed_sample_sizes = 0
  
  Frequentist_prop_go = 0
  Frequentist_sample_sizes = 0
  
  #For every (s1,f1,reward) for the chose Phase II trial
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
    
    #Extract Expected_Observed_Cumulative_Reward
    Expected_Observed_Cumulative_Reward = as.numeric(keysplit[[1]][7])

    #Extract Expected_Frequentist_Cumulative_Reward
    Expected_Frequentist_Cumulative_Reward = as.numeric(keysplit[[1]][8])
    
    #Find Phase II sample size
    T_a_D = s_1D + f_1D
    T_a_C = s_1C + f_1C
    n_a = T_a_D + T_a_C
    
    #Calculate probability (true and bayesian) of getting to s_1D and s_1C, to use in expected reward 
    prob_D_II = dbinom(s_1D, T_a_D, theta_D)
    prob_C_II = dbinom(s_1C, T_a_C, theta_C)
    
    expected_prob_C_II = dbetabinom.ab(s_1C, T_a_C, s_0_C,f_0_C)
    expected_prob_D_II = dbetabinom.ab(s_1D, T_a_D, s_0_D,f_0_D)
    
    #This is used for the proportion of no go decisions, we need to subtract 1 as I have used 2 for go and 1 for no go
    decision = (dictionary_of_decision[[s1.f1.reward]] - 1)

    
    #Below we use the probabilities to reach each state in order to get the expectation 
    Observed_prop_go = Observed_prop_go + expected_prob_C_II*expected_prob_D_II*decision
    Observed_sample_sizes = Observed_sample_sizes +  expected_prob_C_II*expected_prob_D_II*n
    
    #Here we only sum the rewards we have calculated before
    #as we have already weighted them with the corresponding probabilities
    #Andrew and Peter please see this
    Observed_Final.Expected.Reward= Observed_Final.Expected.Reward + Expected_Observed_Cumulative_Reward 
    
    
    
    #Same as before
    Frequentist_prop_go = Frequentist_prop_go + prob_D_II*prob_C_II*decision
    Frequentist_sample_sizes = Frequentist_sample_sizes +  prob_D_II*prob_C_II*n 
    Frequentist_Final.Expected.Reward = Frequentist_Final.Expected.Reward + Expected_Frequentist_Cumulative_Reward 
    
  }
  #Subtracting the cost of recruiting n_a (Phase II) patients from the expected reward to get the final expected reward
  Frequentist_Final.Expected.Reward = Frequentist_Final.Expected.Reward - n_a*(cost.per.patient.PhaseII)
  Observed_Final.Expected.Reward = Observed_Final.Expected.Reward - n_a*(cost.per.patient.PhaseII)

  #Here we get the proportion of no go
  Frequentist_prop_nogo = 1 - Frequentist_prop_go 
  Observed_prop_nogo = 1 - Observed_prop_go 
  
  return(list(dictionary_of_decision = dictionary_of_decision,
              Frequentist_prop_of_nogo = Frequentist_prop_nogo,
              Frequentist_prop_of_go = Frequentist_prop_go,
              Frequentist_average_size_III = Frequentist_sample_sizes,
              Frequentist_Final.Expected.Reward = Frequentist_Final.Expected.Reward,
              Observed_prop_of_nogo = Observed_prop_nogo,
              Observed_prop_of_go = Observed_prop_go,
              Observed_average_size_III = Observed_sample_sizes,
              Observed_Final.Expected.Reward = Observed_Final.Expected.Reward
              ))
}



tic()
Opt1 = Optimal_pharma_frequentist(R = 0.01,
                                  Pop = 1000000,
                                  n_1 = 10,
                                  alpha = 0.025,
                                  beta = 0.2,
                                  k = 0.2,
                                  cost.per.patient.PhaseII = 1,
                                  cost.per.patient.PhaseIII = 1,
                                  s_0_D = 1,
                                  f_0_D = 1,
                                  s_0_C = 4,
                                  f_0_C= 6,
                                  theta_D = 0.6,
                                  theta_C = 0.4)
Opt1
toc()


dbetabinom.ab(s_2C, ceiling(n/2), s_1C+s_0_C,f_1C+f_0_C)
Optimal_patient_final_version(R = 0.01,
                           Pop = 1000000,
                           n_1 = 50,
                           alpha = 0.025,
                           beta = 0.2,
                           k = 0.2,
                           cost.per.patient.PhaseII = 1,
                           cost.per.patient.PhaseIII = 1,
                           s_0_D = 1,
                           f_0_D = 1,
                           s_0_C = 1,
                           f_0_C= 4,
                           theta_D = 0.8,
                           theta_C = 0.1)

#################################################
#trexw
tic()
Frequentist_Revenues_D_0.6_C_0.4 = c()
Frequentist_prop_no_go_D_0.6_C_0.4 = c()
Frequentist_average_sample_D_0.6_C_0.4 = c()
s_0D = 1.2 
f_0D = 0.8
for (i in seq(1,20,1)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = 14,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 0,
                                   f_0_D = 0,
                                   s_0_C = 0,
                                   f_0_C= 0,
                                   theta_D = 0.6,
                                   theta_C = 0.4)
  Frequentist_Revenues_D_0.6_C_0.4 = c(Frequentist_Revenues_D_0.6_C_0.4, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_D_0.6_C_0.4 = c(Frequentist_prop_no_go_D_0.6_C_0.4, Opt$prop_of_nogo)
  Frequentist_prop_succ_D_0.6_C_0.4 = c(Frequentist_prop_succ_D_0.6_C_0.4, Opt$prop_succ)
  Frequentist_average_sample_D_0.6_C_0.4 = c(Frequentist_average_sample_D_0.6_C_0.4, Opt$average_size_III)
  s_0D = s_0D * 1.1
  f_0D = f_0D * 1.1
}
toc()

plot(c(seq(1,20,1)),Frequentist_Revenues_D_0.6_C_0.4, type = "l")
plot(c(seq(1,20,1)),Frequentist_prop_no_go_D_0.6_C_0.4, type = "l")

#################################################
tic()
Frequentist_Revenues_D_0.5_C_0.5 = c()
Frequentist_prop_no_go_D_0.5_C_0.5 = c()
Frequentist_prop_succ_D_0.5_C_0.5 = c()
Frequentist_average_sample_D_0.5_C_0.5 = c()
s_0D = 1
f_0D = 1
for (i in seq(1,20,1)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = 100,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = s_0D,
                                   f_0_D = f_0D,
                                   s_0_C = 25,
                                   f_0_C= 25,
                                   theta_D = 0.5,
                                   theta_C = 0.5)
  Frequentist_Revenues_D_0.5_C_0.5 = c(Frequentist_Revenues_D_0.5_C_0.5, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_D_0.5_C_0.5 = c(Frequentist_prop_no_go_D_0.5_C_0.5, Opt$prop_of_nogo)
  Frequentist_prop_succ_D_0.5_C_0.5 = c(Frequentist_prop_succ_D_0.5_C_0.5, Opt$prop_succ)
  Frequentist_average_sample_D_0.5_C_0.5 = c(Frequentist_average_sample_D_0.5_C_0.5, Opt$average_size_III)
  s_0D = s_0D * 1.1
  f_0D = f_0D * 1.1
}
toc()
plot(c(seq(1,20,1)),Frequentist_Revenues_D_0.5_C_0.5, type = "l")
plot(c(seq(1,20,1)),Frequentist_prop_no_go_D_0.5_C_0.5, type = "l")






##################### run this ##################### 


tic()
Frequentist_Revenues_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c()
Frequentist_prop_no_go_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c()
Frequentist_prop_succ_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c()
Frequentist_average_sample_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c()
for (i in seq(2,70,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 0.6,
                                   f_0_D = 1.4,
                                   s_0_C = 15,
                                   f_0_C= 35,
                                   theta_D = 0.3,
                                   theta_C = 0.3)
  Frequentist_Revenues_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c(Frequentist_Revenues_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c(Frequentist_prop_no_go_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c(Frequentist_prop_succ_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, Opt$prop_succ)
  Frequentist_average_sample_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3 = c(Frequentist_average_sample_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, Opt$average_size_III)
}
toc()
plot(seq(2,70,2),Frequentist_Revenues_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, type = "l")
plot(seq(2,70,2),Frequentist_prop_no_go_sD0.6_fD1.4_sC15_sD35_D_0.3_C_0.3, type = "l")


#############################################################





tic()
Frequentist_Revenues_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c()
Frequentist_prop_no_go_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c()
Frequentist_prop_succ_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c()
Frequentist_average_sample_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c()
for (i in seq(2,70,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 25,
                                   f_0_C= 25,
                                   theta_D = 0.5,
                                   theta_C = 0.5)
  Frequentist_Revenues_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c(Frequentist_Revenues_sD1_fD1_sC25_sD25_D_0.5_C_0.5, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c(Frequentist_prop_no_go_sD1_fD1_sC25_sD25_D_0.5_C_0.5, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c(Frequentist_prop_succ_sD1_fD1_sC25_sD25_D_0.5_C_0.5, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC25_sD25_D_0.5_C_0.5 = c(Frequentist_average_sample_sD1_fD1_sC25_sD25_D_0.5_C_0.5, Opt$average_size_III)
}
toc()

plot(seq(2,70,2),Frequentist_Revenues_sD1_fD1_sC25_sD25_D_0.5_C_0.5, type = "l")
plot(seq(2,70,2),Frequentist_prop_no_go_sD1_fD1_sC25_sD25_D_0.5_C_0.5, type = "l")
#############################################################


tic()
Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c()
Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c()
Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c()
Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c()
for (i in c(2,10,16,20,24,30,36,40,50,70,100)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 100000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1.5,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 3,
                                   f_0_C= 2,
                                   theta_D = 0.8,
                                   theta_C = 0.6)
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c(Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.2, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c(Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.2, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c(Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.2, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.2 = c(Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.2, Opt$average_size_III)
}
toc()

Opt
plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.2, type = "l")

plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.2, type = "l")
toc()
876-849


tic()
Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c()
Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c()
Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c()
Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c()
for (i in c(2,10,16,20,24,30,36,40,50,70,100)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1.5,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 3.9,
                                   f_0_C= 2.1,
                                   theta_D = 0.8,
                                   theta_C = 0.65)
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c(Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.21, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c(Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.21, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c(Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.21, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c(Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.21, Opt$average_size_III)
}
toc()

plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.21, type = "l")

plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.21, type = "l")

tic()
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c()
Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c()
Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c()
Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c()
for (i in c(2,10,16,20,24,30,36,40,50,70,100)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1.5,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 4.2,
                                   f_0_C= 1.8,
                                   theta_D = 0.8,
                                   theta_C = 0.7)
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c(Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.211, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c(Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.211, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.211 = c(Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.211, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.21 = c(Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.211, Opt$average_size_III)
}
toc()
plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.211, type = "l")

plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.211, type = "l")
tic()
Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c()
Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c()
Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c()
Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c()
for (i in c(2,10,16,20,24,30,36,40,50,70,100)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1.5,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 4.5,
                                   f_0_C= 1.5,
                                   theta_D = 0.8,
                                   theta_C = 0.75)
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c(Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.23, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c(Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.23, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c(Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.3_C_0.23, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.23 = c(Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.3_C_0.23, Opt$average_size_III)
}
toc()
plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.3_C_0.2, type = "l")

plot(c(2,10,16,20,24,30,36,40,50,70,100),Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.3_C_0.2, type = "l")

#############################################################
#etrksa

tic()
Frequentist_Revenues_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_prop_no_go_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_prop_succ_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_average_sample_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c()
for (i in seq(2,70,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1.1,
                                   f_0_D = 0.9,
                                   s_0_C = 3,
                                   f_0_C= 2,
                                   theta_D = 0.3,
                                   theta_C = 0.2)
  Frequentist_Revenues_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_Revenues_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_prop_no_go_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_prop_succ_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, Opt$prop_succ)
  Frequentist_average_sample_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_average_sample_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, Opt$average_size_III)
}
toc()


plot(seq(2,38,2),Frequentist_Revenues_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, type = "l")
plot(seq(2,38,2),Frequentist_average_sample_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, type = "l")

plot(seq(2,38,2),Frequentist_prop_no_go_sD0.8_fD1.2_sC10_sD40_D_0.4_C_0.2, type = "l")
toc()
876-849

27/10
#############################################################


tic()
Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c()
Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c()
for (i in seq(2,70,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 10,
                                   f_0_C= 40,
                                   theta_D = 0.4,
                                   theta_C = 0.2)
  Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.4_C_0.2, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.4_C_0.2, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_prop_succ_sD1_fD1_sC10_sD40_D_0.4_C_0.2, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.4_C_0.2 = c(Frequentist_average_sample_sD1_fD1_sC10_sD40_D_0.4_C_0.2, Opt$average_size_III)
}
toc()

plot(seq(2,70,2),Frequentist_Revenues_sD1_fD1_sC10_sD40_D_0.4_C_0.2, type = "l")

plot(seq(2,70,2),Frequentist_prop_no_go_sD1_fD1_sC10_sD40_D_0.4_C_0.2, type = "l")



#############################################################


tic()
Frequentist_Revenues_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_prop_no_go_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_prop_succ_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_average_sample_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c()
for (i in seq(2,100,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 20,
                                   f_0_C= 30,
                                   theta_D = 0.6,
                                   theta_C = 0.4)
  Frequentist_Revenues_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_Revenues_sD1_fD1_sC20_sD30_D_0.6_C_0.4, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_prop_no_go_sD1_fD1_sC20_sD30_D_0.6_C_0.4, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_prop_succ_sD1_fD1_sC20_sD30_D_0.6_C_0.4, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_average_sample_sD1_fD1_sC20_sD30_D_0.6_C_0.4, Opt$average_size_III)
}
toc()

plot(seq(2,100,2),Frequentist_Revenues_sD1_fD1_sC20_sD30_D_0.6_C_0.4, type = "l")
plot(seq(2,100,2),Frequentist_prop_no_go_sD1_fD1_sC20_sD30_D_0.6_C_0.4, type = "l")





#############################################################
Opt$Final.Expected.Reward

tic()
Frequentist_Revenues_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_prop_no_go_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_prop_succ_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c()
Frequentist_average_sample_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c()
for (i in seq(2,100,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1.2,
                                   f_0_D = 0.8,
                                   s_0_C = 20,
                                   f_0_C= 30,
                                   theta_D = 0.6,
                                   theta_C = 0.4)
  Frequentist_Revenues_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_Revenues_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_prop_no_go_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_prop_succ_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, Opt$prop_succ)
  Frequentist_average_sample_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4 = c(Frequentist_average_sample_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, Opt$average_size_III)
}
toc()

plot(seq(2,100,2),Frequentist_Revenues_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, type = "l")
plot(seq(2,100,2),Frequentist_prop_no_go_sD1.2_fD0.8_sC20_sD30_D_0.6_C_0.4, type = "l")








#############################################################


tic()
Frequentist_Revenues_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_prop_no_go_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_prop_succ_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_average_sample_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c()
for (i in seq(2,100,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1.6,
                                   f_0_D = 0.4,
                                   s_0_C = 15,
                                   f_0_C= 35,
                                   theta_D = 0.8,
                                   theta_C = 0.3)
  Frequentist_Revenues_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_Revenues_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_prop_no_go_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_prop_succ_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, Opt$prop_succ)
  Frequentist_average_sample_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_average_sample_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, Opt$average_size_III)
}
toc()

plot(seq(2,100,2),Frequentist_Revenues_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, type = "l")




#############################################################

tic()
Frequentist_Revenues_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c()
Frequentist_prop_no_go_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c()
Frequentist_prop_succ_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c()
Frequentist_average_sample_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c()
for (i in seq(2,100,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 30,
                                   f_0_C= 20,
                                   theta_D = 0.8,
                                   theta_C = 0.6)
  Frequentist_Revenues_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c(Frequentist_Revenues_sD1_fD1_sC30_sD20_D_0.8_C_0.6, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c(Frequentist_prop_no_go_sD1_fD1_sC30_sD20_D_0.8_C_0.6, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c(Frequentist_prop_succ_sD1_fD1_sC30_sD20_D_0.8_C_0.6, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC30_sD20_D_0.8_C_0.6 = c(Frequentist_average_sample_sD1_fD1_sC30_sD20_D_0.8_C_0.6, Opt$average_size_III)
}
toc()

plot(seq(2,100,2),Frequentist_Revenues_sD1.6_fD0.4_sC15_sD35_D_0.8_C_0.3, type = "l")






#############################################################
#trexw ayto

tic()
Frequentist_Revenues_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_prop_no_go_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_prop_succ_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c()
Frequentist_average_sample_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c()
for (i in seq(2,100,2)){
  print(i)
  Opt = Optimal_pharma_frequentist(R = 0.01,
                                   Pop = 1000000,
                                   n_1 = i,
                                   alpha = 0.025,
                                   beta = 0.2,
                                   k = 0.2,
                                   cost.per.patient.PhaseII = 1,
                                   cost.per.patient.PhaseIII = 1,
                                   s_0_D = 1,
                                   f_0_D = 1,
                                   s_0_C = 15,
                                   f_0_C= 35,
                                   theta_D = 0.8,
                                   theta_C = 0.3)
  Frequentist_Revenues_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_Revenues_sD1_fD1_sC15_sD35_D_0.8_C_0.3, Opt$Final.Expected.Reward)
  Frequentist_prop_no_go_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_prop_no_go_sD1_fD1_sC15_sD35_D_0.8_C_0.3, Opt$prop_of_nogo)
  Frequentist_prop_succ_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_prop_succ_sD1_fD1_sC15_sD35_D_0.8_C_0.3, Opt$prop_succ)
  Frequentist_average_sample_sD1_fD1_sC15_sD35_D_0.8_C_0.3 = c(Frequentist_average_sample_sD1_fD1_sC15_sD35_D_0.8_C_0.3, Opt$average_size_III)
}
toc()

plot(seq(2,100,2),Frequentist_prop_no_go_sD2.5_fD1_sC15_sD35_D_0.8_C_0.3, type = "l")


save.image(file='teliko_tomorrow_OPENTHIS.RData')
