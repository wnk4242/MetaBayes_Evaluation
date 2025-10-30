################################
##Fixed-Effect Meta-Anlytic BF##
################################
#library(BayesFactor)
#bf = meta.ttestBF(t = df$t, n1 = df$n1, n2 = df$n2, nullInterval = c(0, Inf), rscale = sqrt(2)/2); bf


######################################
##Evidence Updating BF (Complete BF)##
######################################
#####
# #Genearting some data sets
# set.seed(01)
# #sample size for original and replication studies
# n=10
# #Original study
# #Sample size for experimental and control groups is 10 (each group has ten)
# df0 <- data.frame(
#   group = rep(letters[1:2], each = 10),
#   result = c(runif(10, min=0, max=5), runif(10, min=0, max=5))
# )
# df0
# #Replication study 1
# df1 <- data.frame(
#   group = rep(letters[1:2], each = 10),
#   result = c(runif(10, min=0, max=5), runif(10, min=0, max=5))
# )
# df1
# #Replication study 2
# df2 <- data.frame(
#   group = rep(letters[1:2], each = 10),
#   result = c(runif(10, min=0, max=5), runif(10, min=0, max=5))
# )
# df2
# #Use data.table library to find mean and sd of each group
# library(data.table)
# setDT(df0) #convert data frame to data table 
# setDT(df1)
# setDT(df2)
# df0[,list(mean=mean(result)), by=group]
# df0[,list(sd=sd(result)), by=group]
# df1[,list(mean=mean(result)), by=group]
# df1[,list(sd=sd(result)), by=group]
# df2[,list(mean=mean(result)), by=group]
# df2[,list(sd=sd(result)), by=group]
# #If you want to use the above simulated data in the combineStudies() function, 
# #you need to convert their structure using the following code
# #n_x <- rep(10,3)
# #n_y <- rep(10,3)
# #mean_x <- cbind(xbar_orig,xbar_rep1,xbar_rep2);mean_x <- unlist(unname(mean_x))
# #mean_y <- cbind(ybar_orig,ybar_rep1,ybar_rep2);mean_y <- unlist(unname(mean_y))
# #var_x <- cbind(s2_origx,s2_rep1x,s2_rep2x);var_x <- unlist(unname(var_x))
# #var_y <- cbind(s2_origy,s2_rep2y,s2_rep2y);var_y <- unlist(unname(var_y))
#####
#Manually compute t_all and EU Bayes factor
#If there were just the original study and one replication study
#original study stats
# n_origx=n
# xbar_orig=df0[,list(mean=mean(result)), by=group][1,2]
# s2_origx=(df0[,list(sd=sd(result)), by=group][1,2])^2
# n_origy=n
# ybar_orig=df0[,list(mean=mean(result)), by=group][2,2]
# s2_origy=(df0[,list(sd=sd(result)), by=group][2,2])^2
# #rep 1 study stats
# n_repx=n
# xbar_rep=df1[,list(mean=mean(result)), by=group][1,2]
# s2_repx=(df1[,list(sd=sd(result)), by=group][1,2])^2
# n_repy=n
# ybar_rep=df1[,list(mean=mean(result)), by=group][2,2]
# s2_repy=(df1[,list(sd=sd(result)), by=group][2,2])^2
# #rep 2 study stats
# n_repx=n
# xbar_rep=df2[,list(mean=mean(result)), by=group][1,2]
# s2_repx=(df2[,list(sd=sd(result)), by=group][1,2])^2
# n_repy=n
# ybar_rep=df2[,list(mean=mean(result)), by=group][2,2]
# s2_repy=(df2[,list(sd=sd(result)), by=group][2,2])^2
# # n_origx=10
# # xbar_orig=3.18
# # s2_origx=1.4^2
# # n_origy=10
# # ybar_orig=2.44
# # s2_origy=1.5^2
# # #rep 1 study stats
# # n_repx=10
# # xbar_rep=2.31
# # s2_repx=1.63^2
# # n_repy=10
# # ybar_rep=7.61
# # s2_repy=1.24^2
# # #rep 2 study stats
# # n_repx=10
# # xbar_rep=2.88
# # s2_repx=1.21^2
# # n_repy=10
# # ybar_rep=7.19
# # s2_repy=1.3^2
# n_allx=n_origx+n_repx
# n_ally=n_origy+n_repy
# xbar_all=((n_origx*xbar_orig)+(n_repx*xbar_rep))/n_allx
# ybar_all=((n_origy*ybar_orig)+(n_repy*ybar_rep))/n_ally
# v_origx=n_origx-1
# v_origy=n_origy-1
# v_repx=n_repx-1
# v_repy=n_repy-1
# ssq_x=(v_origx*s2_origx)+(n_origx*(xbar_orig)^2)+(v_repx*s2_repx)+(n_repx*(xbar_rep)^2)-(n_allx*(xbar_all)^2)
# ssq_y=(v_origy*s2_origy)+(n_origy*(ybar_orig)^2)+(v_repy*s2_repy)+(n_repy*(ybar_rep)^2)-(n_ally*(ybar_all)^2)
# s2=(ssq_x+ssq_y)/(n_allx+n_ally-2)
# t_all=(xbar_all-ybar_all)/sqrt(s2*((1/n_allx)+(1/n_ally)));t_all
# #t_all=-3.0448
# #In JSAP, use t_all, n_allx, n_ally in Bayesian Independent Samples T-Test, select Group 1 != Group 2
# #yields a BF10 = -3.046

#####
# If there were multiple groups of replications
# n_origx=n
# xbar_orig=df0[,list(mean=mean(result)), by=group][1,2]
# s2_origx=(df0[,list(sd=sd(result)), by=group][1,2])^2
# n_origy=n
# ybar_orig=df0[,list(mean=mean(result)), by=group][2,2]
# s2_origy=(df0[,list(sd=sd(result)), by=group][2,2])^2
# #rep 1 study stats
# n_rep1x=n
# xbar_rep1=df1[,list(mean=mean(result)), by=group][1,2]
# s2_rep1x=(df1[,list(sd=sd(result)), by=group][1,2])^2
# n_rep1y=n
# ybar_rep1=df1[,list(mean=mean(result)), by=group][2,2]
# s2_rep1y=(df1[,list(sd=sd(result)), by=group][2,2])^2
# #rep 2 study stats
# n_rep2x=n
# xbar_rep2=df2[,list(mean=mean(result)), by=group][1,2]
# s2_rep2x=(df2[,list(sd=sd(result)), by=group][1,2])^2
# n_rep2y=n
# ybar_rep2=df2[,list(mean=mean(result)), by=group][2,2]
# s2_rep2y=(df2[,list(sd=sd(result)), by=group][2,2])^2
# n_allx=n_origx+n_rep1x+n_rep2x
# n_ally=n_origy+n_rep1y+n_rep2y
# xbar_all=((n_origx*xbar_orig)+(n_rep1x*xbar_rep1)+(n_rep2x*xbar_rep2))/n_allx
# ybar_all=((n_origy*ybar_orig)+(n_rep1y*ybar_rep1)+(n_rep2y*ybar_rep2))/n_ally
# v_origx=n_origx-1
# v_origy=n_origy-1
# v_rep1x=n_rep1x-1
# v_rep1y=n_rep1y-1
# v_rep2x=n_rep2x-1
# v_rep2y=n_rep2y-1
# ssq_x=(v_origx*s2_origx)+(n_origx*(xbar_orig)^2)+(v_rep1x*s2_rep1x)+(n_rep1x*(xbar_rep1)^2)+(v_rep2x*s2_rep2x)+(n_rep2x*(xbar_rep2)^2)-(n_allx*(xbar_all)^2)
# ssq_y=(v_origy*s2_origy)+(n_origy*(ybar_orig)^2)+(v_rep1y*s2_rep1y)+(n_rep1y*(ybar_rep1)^2)+(v_rep2y*s2_rep2y)+(n_rep2y*(ybar_rep2)^2)-(n_ally*(ybar_all)^2)
# s2=(ssq_x+ssq_y)/(n_allx+n_ally-2)
# #Get overall t value
# t_all=(xbar_all-ybar_all)/sqrt(s2*((1/n_allx)+(1/n_ally)));t_all
# #Get complete BF
# bf012 <- ttest.tstat(t = unlist(t_all), n1 = n_allx, n2 = n_ally, simple = TRUE);bf012
# 


# Use a custom function to compute the combined t value and EU Bayes factor
#######
# Use combineStudies() to calculate overall t value and EU Bayes factor
combStudies <- function(n_x, mean_x, var_x, n_y, mean_y, var_y) {
  #n_x is sample size for treatment group in a replication
  #n_y is sample size for control group in a replication
  #mean_x is mean of treatment group of a replication
  #mean_y is mean of control group of a replication
  #var_x is variance of treatment group in a replication
  #var_y is variance of control group in a replication
  # Initialize combined statistics
  n_all_x <- 0; mean_all_x <- 0; var_all_x <- 0; ss_x <- 0
  n_all_y <- 0; mean_all_y <- 0; var_all_y <- 0; ss_y <- 0
  # Calculate combined stats for the treatment groups (x)
  for (i in 1:length(n_x)) {
    n_all_x <- n_all_x + n_x[i]
    mean_all_x <- mean_all_x + n_x[i] * mean_x[i]
    ss_x <-  ss_x + (n_x[i]-1)*var_x[i] + n_x[i]*((mean_x[i])^2)
  }
  mean_all_x <- mean_all_x / n_all_x
  ss_x <- ss_x - n_all_x * (mean_all_x)^2
  # Calculate combined stats for the control groups (y)
  for (i in 1:length(n_y)) {
    n_all_y <- n_all_y + n_y[i]
    mean_all_y <- mean_all_y + n_y[i] * mean_y[i]
    ss_y <-  ss_y + (n_y[i]-1) * var_y[i] + n_y[i] * ((mean_y[i])^2)
  }
  mean_all_y <- mean_all_y / n_all_y
  ss_y <- ss_y - n_all_y * (mean_all_y)^2
  s2 <- (ss_x + ss_y) / (n_all_x + n_all_y-2)
  # Calculate the combined t value
  t_all <- (mean_all_x - mean_all_y) / sqrt(s2 * ((1/n_all_x) + (1/n_all_y)))
  return(c(t_all,n_all_x, n_all_y))
}
# combStats <- combStudies(n_x = df$n1, mean_x = df$m1, var_x = df$v1, n_y = df$n2, mean_y = df$m2, var_y = df$v2)
# EUBF <- ttest.tstat(t = combStats[1], n1 = combStats[2], n2 = combStats[3], simple = TRUE); EUBF


#######################################
##Bayes factor based on meta-analysis##
#######################################
# BFroMeta() works with BFbMACalc() in helper function.R
# Two-sided BF returning BF10 = H1/H0
# BFroMeta <- function(p_value, var_theta, theta_A = 0.1){
#   mn0 <-  pi*theta_A^2/(2*var_theta)
#   z_m <- qnorm(1-p_value/2) #if p-value is already devided by 2, then remove /2
#   BF <- sqrt(1+mn0)*exp(-(z_m)^2/(2+(2/mn0))) #Spiegelhalter 2004 p 131
#   #return(c(1/BF, mn0, z_m))  
#   return(c(1/BF))  
# }

# One-sided BF (H1: theta > 0) returning BF10 = H1/H0
BFroMeta <- function(p_value, var_theta, theta_A = 0.1){
  mn0   <- pi * theta_A^2 / (2 * var_theta)            # m/n0
  z_abs <- qnorm(1 - p_value/2)                         # from TWO-SIDED p (magnitude only)
  BF_two <- sqrt(1 + mn0) * exp( - (z_abs^2) * mn0 / (2 * (1 + mn0)) )
  adj      <- 2 * pnorm( z_abs * sqrt( mn0 / (1 + mn0) ) ) #one-sided adjustment
  BF_one <- BF_two / adj #one-sided BF
  return(c(1 / BF_one))
}


# df <- simMA(20,0.5,0.03,"none","low",fixed.n = 20)
# ma <- rma(yi=df$d, vi=df$v, data=df, method = "DL");ma
# BF <- BFroMeta(p_value=ma$pval, var_theta = (ma$se)^2, theta_A = 0.1);BF

# There are two ways to obtain z_m. First approach is to use the combined ES divided by its SE.
# The second approach is to use the p value of the combined ES and use qnorm(1-p_value/2). 
# The results are almost identical.

# Here's a small simulation to prove that two approaches obtain same results:
# Also, no matter you choose one-sided or two-sided tests for the studies to be included into a meta-analysis,
# the meta-analysis's combined effect size's p value will be the same. In other words, there is no such thing as
# "one-sided" or "two-sided" BFbMA, there is only BFbMA (unless you change its 'spike-smear' prior).
# Therefore, when you analyze one-sided or two-sided MABF, just use the same BFbMA compared with other MABFs.
# library(metafor)
# 
# set.seed(12323)  # reproducibility
# 
# # Step 1: Simulate raw data for 10 experiments
# n_studies <- 10
# n_treat <- 30
# n_control <- 30
# true_effect <- 0.5  # true difference in means
# 
# raw_data <- lapply(1:n_studies, function(i) {
#   control <- rnorm(n_control, mean = 0, sd = 1)
#   treatment <- rnorm(n_treat, mean = true_effect, sd = 1)
#   data.frame(study = i,
#              group = rep(c("control", "treatment"), each = n_control),
#              value = c(control, treatment))
# })
# 
# # Step 2: Functions for two-sided and one-sided t-tests
# t_test_two_sided <- function(x, y) {
#   t.test(x, y, alternative = "two.sided", var.equal = TRUE)
# }
# 
# t_test_one_sided <- function(x, y) {
#   t.test(x, y, alternative = "greater", var.equal = TRUE)
# }
# 
# # Step 3: Compute effect sizes and store results
# results <- data.frame(study = 1:n_studies,
#                       m1i = NA, sd1i = NA, n1i = NA,
#                       m2i = NA, sd2i = NA, n2i = NA,
#                       p_two_sided = NA, p_one_sided = NA)
# 
# for (i in 1:n_studies) {
#   dat <- raw_data[[i]]
#   x <- dat$value[dat$group == "treatment"]
#   y <- dat$value[dat$group == "control"]
#   
#   res_two <- t_test_two_sided(x, y)
#   res_one <- t_test_one_sided(x, y)
#   
#   results[i, ] <- c(i,
#                     mean(x), sd(x), length(x),
#                     mean(y), sd(y), length(y),
#                     res_two$p.value,
#                     res_one$p.value)
# }
# 
# # Step 4: Compute SMDs (Hedges' g) using metafor::escalc
# es <- escalc(measure = "SMD",
#              m1i = m1i, sd1i = sd1i, n1i = n1i,
#              m2i = m2i, sd2i = sd2i, n2i = n2i,
#              data = results)
# 
# # Step 5: Random-effects meta-analysis (results are same regardless of p-values)
# meta_two_sided <- rma(yi, vi, data = es, method = "REML")
# meta_one_sided <- rma(yi, vi, data = es, method = "REML")
# 
# # Print results
# print(meta_two_sided$pval)
# print(meta_one_sided$pval)
# 
# # use p value to obtain BF
# BFroMeta <- function(p_value, var_theta, theta_A = 0.1){
#   mn0 <-  pi*theta_A^2/(2*var_theta)
#   z_m <- qnorm(1-p_value/2) #if p-value is already devided by 2, then remove /2
#   BF <- sqrt(1+mn0)*exp(-(z_m)^2/(2+(2/mn0))) #Spiegelhalter 2004 p 131
#   return(c(BF))  
# }
# BFroMeta(meta_two_sided$pval, (meta_two_sided$se)^2, 1.44)
# 
# # use es and se to obtain BF
# BFroMeta <- function(es, se, theta_A = 0.1){
#   mn0 <-  pi*theta_A^2/(2*(se^2))
#   z_m <- es/se 
#   BF <- sqrt(1+mn0)*exp(-(z_m)^2/(2+(2/mn0))) #Spiegelhalter 2004 p 131
#   return(c(BF))  
# }
# 
# BFroMeta(meta_two_sided$b, meta_two_sided$se, 1.44)
############################################################
##Inclusion Bayes factor (Bayesian Averaged Meta-Analysis)##
############################################################
#to run the bf10_t function, you need to source("./baMA functions.R")
#Initiate a variable BFplus0 to contain Bayes factor for each primary study in the simulated meta-analysis
#You don't need the indivisual BFplus0 to calculate the inclusion BF for a meta-analysis
# BFplus0 <- numeric(nrow(d))
# #calculate Bayes factor for each primary study in the simulated meta-analysis
# for (i in seq_len(nrow(d))) {
#   BFplus0[i] <- bf10_t(t = d$t[i], ny = d$n1[i], nx = d$n2[i],
#                        independentSamples = TRUE, prior.location = 0.34999,
#                        prior.scale = 0.1021, prior.df = 3)$BFplus0
# }
# library(metaBMA)
# #Create a data set
# d <- simMA(20,0,0.03,"none","low",fixed.n = 60)
# #ma <- rma(yi=d$d, vi=d$v, data=d, method = "DL");ma
# # prior for true effect size
# priorEStesting <- prior(family = "custom",
#                         function(x) LaplacesDemon::dst(x, mu = 0.34999, sigma = 0.1021, nu = 3)/
#                         pst(0, mu = 0.34999, sigma = 0.1021, nu = 3, lower.tail = FALSE),
#                         label = "d", lower = 0, upper = Inf)
# # prior for heterogeneity
# priorTau <- prior(family = "custom",
#                   function(x) extraDistr::dinvgamma(x, alpha = 1, beta = .15),
#                   label = "tau", lower = 0, upper = Inf)
# # conduct Bayesian model averaged meta-analysis to calculate inclusion BF
# time_taken <- system.time({
#   baMA <- meta_bma(y = d$d, SE = d$se, d = priorEStesting,
#                         tau = priorTau, summarize = "integrate")#[[9]]$incl.BF
# })
# print(time_taken)
# baMA
