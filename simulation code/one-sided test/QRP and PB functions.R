#0409 Update: changed expDataB MaxN from 100 to 1000
#I added m1, v1, m2, v2 to the output
#Added a participant-level data generation process (include participant-level standard deviation and raw data of treatment and control groups)
#Added a simMA_pbar
##############################################################################################################################################
# library(truncdist)
# getN <- function(k=1, min.n = 5, max.n = 1905, shape=1.15326986, scale=0.04622745) {
#   library(invgamma)
#   ns <- round(rtrunc(n=k, spec="invgamma", a=min.n, b=max.n, shape=shape, scale=scale))
# }


#==============
#   Outlier   #
#==============
outlier=function(x,mean,sd){
  out=if(abs((x-mean)/sd)<2){0}else{1} #if is outlier, assign 1, if not, assign 0
}

#############
#============
# expDataB  #
#============
# Removed g3 and g4 because moderator is not used for my dissertation
# If the maxN is small like 30, the observed correlation may not remain to be cbdv
# g1 is responses of subjects in treatment group on two correlated DVs
# g1 is responses of subjects in control group on two correlated DVs
# expDataB <- function(delta, tau, cbdv, maxN = 100){
#   delta_i = delta + rnorm(1, 0, tau)
#   D = matrix(delta_i, 2, maxN) # 2 rows, maxN columns
#   g1 = MASS::mvrnorm(maxN, rep(0,2), matrix(c(1,cbdv,cbdv,1),2,2)) + delta_i 
#   g2 = MASS::mvrnorm(maxN, rep(0,2), matrix(c(1,cbdv,cbdv,1),2,2))
#   G = array(c(g1,g2,D),dim=c(maxN, 2, 3)) #in this array named G: ,,1 is g1, ,,2 is g2, ,,3 is D
#   return(G) 
# }

expDataB <- function(delta, tau, cbdv, maxN = 1000){
  sigma <- runif(n = 1, min = 0.5, max = 2.5)
#  sigma <- 1
  delta_i = delta + rnorm(1, 0, tau)
  D = matrix(delta_i, 2, maxN) # 2 rows, maxN columns
  g22 = rnorm(n = 1, mean = 0, sd = 1)
  g11 = g22 + delta_i*sigma
  g2 = MASS::mvrnorm(maxN, rep(g22,2), matrix(c(1,cbdv,cbdv,1),2,2))
  g1 = MASS::mvrnorm(maxN, rep(g11,2), matrix(c(1,cbdv,cbdv,1),2,2))
  G = array(c(g1,g2,D),dim=c(maxN, 2, 3)) #in this array named G: ,,1 is g1, ,,2 is g2, ,,3 is D
  return(G)
}


##############
#==========
# testIt  #
#==========
testIt=function(DV, out){
  # a set of conditionals that determine the data to be analyzed.
  # no exclusion of outliers
  if(out==0){  
    Y = DV[,1]
    X = DV[,2]  
  }
  # exclusion of outliers
  # outliers are labeled 1 in DV[,3], nonoutliers are labeled 0
  if(out==1){
    Y = subset(DV[,1], DV[,3] < 1)
    X = subset(DV[,2], DV[,3] < 1)   
  }
  test = t.test(Y~X,var.equal=T) #see R syntax: https://statistics.laerd.com/r-tutorials/independent-samples-t-test-using-r-excel-and-rstudio-3.php
  n1 = length(subset(Y,X==1))
  n2 = length(subset(Y,X==2))
  v1 = var(subset(Y,X==1))
  v2 = var(subset(Y,X==2))
  N  = n1+n2
  t  = as.numeric(test[1])             
  p  = test$p.value
  m1 = as.numeric(test$estimate[1])
  m2 = as.numeric(test$estimate[2])
  df = N - 2
  S = sqrt( ((n1 - 1)*v1 + (n2 - 1)*v2) / df )
  d  = (m1-m2)/S
  out <- c(d,p,t,N,n1,n2,m1,v1,m2,v2)
  return(out)
}
##############
#===============
#    analyB    # 
#===============
#added an else condition for when no sig. results are found
analyB <- function(g1, g2, D, multDV, out){
  #Create combo groups  
  G1=rbind(g1); G2=rbind(g2)
  #create X codes
  X1.1=replicate(length(G1[,1]),1); X1.2=replicate(length(G1[,1]),1) #maxN number of 1
  X2.1=replicate(length(G2[,1]),2); X2.2=replicate(length(G2[,1]),2) #maxN number of 2
  #Create outlier codes
  o1.1=mapply(outlier,G1[,1],mean(G1[,1]),sd(G1[,1]))
  o1.2=mapply(outlier,G1[,2],mean(G1[,2]),sd(G1[,2]))
  o2.1=mapply(outlier,G2[,1],mean(G2[,1]),sd(G2[,1]))
  o2.2=mapply(outlier,G2[,2],mean(G2[,2]),sd(G2[,2]))
  #combine codes with outcome values
  c1=cbind(G1[,1],X1.1,o1.1); c2=cbind(G1[,2],X1.2,o1.2)
  c3=cbind(G2[,1],X2.1,o2.1); c4=cbind(G2[,2],X2.2,o2.2)
  #make "datasets"
  DV1=rbind(c1,c3) #responses of subjects in treatment and control groups on the first DV
  DV2=rbind(c2,c4) #responses of subjects in treatment and control groups on the second DV
  #test in each possible way
  #  first DV
  t100=testIt(DV1,0) #conduct t test on DV1 including outliers
  t101=testIt(DV1,1) #conduct t test on DV1 excluding outliers
  #  second DV
  t200=testIt(DV2,0) #conduct t test on DV2 including outliers
  t201=testIt(DV2,1) #conduct t test on DV2 excluding outliers
  #pull the best result given options
  #no outlier exclusion is always first priority, for both DVs
  #even if you set multDV=1, meaning there are two DVs, if the first DV is sig., R won't select DV2 even if it is sig.
  #add an else condition for when no conditions are met
  # Initialize an empty list to store candidate values
  candidates = list(t100 = t100, t101 = t101, t200 = t200, t201 = t201)
  
  best = if(t100[2]<.05 & t100[1]>0){t100}                     #DV1 and no outlier removal
  else if(out==1 & t101[2]<.05 & t101[1]>0){t101}              #DV1 with outlier removal
  else if(multDV==1 & t200[2]<.05 & t200[1]>0){t200}           #DV2 and no outlier removal
  else if(multDV==1 & out==1 & t201[2]<.05 & t201[1]>0){t201}  #DV2 with outlier removal
  else {
    #message("None of the conditions were met!")
    t100
  }
  
  # Print all candidate values
  #print("All candidate values:")
  #print(candidates)
  
  #get additional info for the best results
  d = best[1]
  p = best[2]
  t = best[3]
  N = best[4]
  n1 = best[5]
  n2 = best[6]
  m1 = best[7]
  v1 = best[8]
  m2 = best[9]
  v2 = best[10]
  df = N - 2
  #d_v= (n1 + n2)/(n1 * n2) + (d^2 / (2 *df)) * (n1 + n2) / df
  d_v= (n1 + n2)/(n1 * n2) + d^2 / (2 *df)
  d_se = sqrt(d_v)
  pow = pwr::pwr.t2n.test(d, n1=n1, n2=n2)
  #pow = pwr.t2n.test(d, n1=n1, n2=n2)
  pwr = pow$power
  #return the best result
  out=c(d, p, t, N, d_v, d_se, pwr, n1, n2, D, m1, v1, m2, v2)
  #all.out <- list(out, t100, t101, t200, t201) #DV1 is combined treatment group
  return(out)
}

#==============
#   simData.noQRP   #
#==============

# generate the results from an unbiased experiment
# delta is the true effect
# tau indicated heterogeneity
# minN and meanN are fed to a negative binomial for 
# sample size

# results from an unbiased experiment 
simData.noQRP <- function(delta, tau, fixed.n=NULL){
  
  # get the per-group sample size - 
  # either a fixed n, or sampled from the gamma distribution 
  #  if (!is.null(fixed.n)) {
  n <- fixed.n
  #  } else {
  #    n <- getN(k=1)
  #  }
  
  #calculate the treatement effect as a function of the 
  #true effect, delta, and tau
  delta_i = delta + rnorm(1, 0, tau)
  
  #generate two independent vectors of raw data 
  #the mean is zero and error is randomly distributed
  #and equal between groups
  sigma <- runif(n = 1, min = 0.5, max = 2.5) #standard deviation sigma
#  sigma <- 1
  Ycc = rnorm(1, 0, 1)
  Yee = Ycc + delta_i*sigma
  Yc = rnorm(n = n, mean = Ycc, sd = sigma)
  Ye = rnorm(n = n, mean = Yee, sd = sigma)
  #get the summary stats
  m1 = mean(Ye)
  v1 = var(Ye)
  m2 = mean(Yc)
  v2 = var(Yc)
  n1 = n
  n2 = n
  df = 2*n-2
  
  #get the pooled variance (formula from Hedge's g as)
  S = sqrt( ((n1 - 1)*v1 + (n2 - 1)*v2) / df )
  
  #compare the two distributions
  test = t.test(Ye, Yc)
  
  #calculate d, the variance of d, the p-value, the t-stat, and n.
  d = (m1 - m2)/S
  
  d_v = (n1 + n2)/(n1 * n2) + (d^2 / (2*(n1+n2)))
  d_se = sqrt(d_v)
  
  p = test$p.value
  t = as.numeric(test$statistic)
  
  #get power
  pow = pwr.t2n.test(d, n1 = n1, n2 = n2)
  pwr = pow$power 
  
  #output 
  out = c(d, p, t, n1+n2, d_v, d_se, pwr,  n1, n2, delta_i, m1, v1, m2, v2)  
}

# result <- as.data.frame(t(result))
# 
# names(result) = c(  'd',       # effect size, d
#                     'p',       # p value for the two group comparison
#                     't',       # t value for the two group comparison
#                     'N',       # total N
#                     'v',       # variance for the effect size
#                     'se',      # standard error for the effect size
#                     'pow',     # power given the true effect for the two group comparison
#                     'n1',      # experimental group sample size
#                     'n2',      # control group sample size
#                     'delta_i', # the study-level true effect
#                     'm1',
#                     'v1',
#                     'm2',
#                     'v2')     
# 
# result

#==================
#     simData.QRP     #     
#==================
# Produces results, a, from a p-hacked experiment.
simData.QRP <- function(delta, tau, QRP.strategy, maxN = 3000, fixed.n=NULL){
  #correlation between multiple DVs is set to 0.20 as default
  cbdv = 0.2
  # if QRP strategy is NONE
  if (QRP.strategy=='none'){
    a = simData.noQRP(delta, tau, fixed.n=fixed.n)
  }
  
  #if QRP strategy is LIGHT
  else if (QRP.strategy=='lit'){
    #get data for a study using QRPs
    G <- expDataB(delta, tau, cbdv)
    
    #determine the starting per-group sample size
    # get the per-group sample size - 
    # either a fixed n, or sampled from the gamma distribution 
    #    if (!is.null(fixed.n)) {
    s <- fixed.n
    #    } else {
    #      s <- getN(k=1)
    #    }
    
    # Divide sample size by 2: the idea is that the main factor of interest defined the two group sizes. A moderator factor is then added to create a 2*2, but because the moderator is not the main focus, the empirical sample sizes should only be used for the two groups formed by the main factor--not the four groups formed by the 2*2 split.
    #s <- round(s/2) #the starting per group sample size, which will be added upon per peek, smaller than the maxN
    #run the first analysis with some QRPs applied
    #all the data has been generated by expDataB(), use [1:s,] to reveal a step size of existing data after each peek
    a = analyB(g1 = G[,,1][1:s,],  #treatment group, 1:the current sample size
               g2 = G[,,2][1:s,],  #control group, 1:the current sample size
               D = G[,,3][1,1],   # the study-lvl true effect
               multDV=0, out=1) # LIGHT 
    
  }
  
  #if QRP strategy is MODERATE
  else if (QRP.strategy=='mod'){
    #get data for a study using QRPs
    G <- expDataB(delta, tau, cbdv)
    
    #determine the starting per-group sample size
    # get the per-group sample size - 
    # either a fixed n, or sampled from the gamma distribution 
    #    if (!is.null(fixed.n)) {
    s <- fixed.n
    #    } else {
    #      s <- getN(k=1)
    #    }
    
    # Divide sample size by 2: the idea is that the main factor of interest defined the two group sizes. A moderator factor is then added to create a 2*2, but because the moderator is not the main focus, the empirical sample sizes should only be used for the two groups formed by the main factor--not the four groups formed by the 2*2 split.
    #s <- round(s/2) #the starting per group sample size, which will be added upon per peek, smaller than the maxN
    #run the first analysis with some QRPs applied
    #all the data has been generated by expDataB(), use [1:s,] to reveal a step size of existing data after each peek
    a = analyB(g1 = G[,,1][1:s,],  #treatment group, 1:the current sample size
               g2 = G[,,2][1:s,],  #control group, 1:the current sample size
               D = G[,,3][1,1],   # the study-lvl true effect
               multDV=1, out=1) # MODERATE 
    
  }
  #if QRP strategy is AGGRESSIVE
  # ("aggressive" now is called "strong" in the paper)
  else if (QRP.strategy=='agg'){
    G = expDataB(delta,tau,cbdv,maxN)
    #determine the starting per-group sample size
    # get the per-group sample size - 
    # either a fixed n, or sampled from the gamma distribution 
    #    if (!is.null(fixed.n)) {
    s <- fixed.n
    #    } else {
    #      s <- getN(k=1)
    #    }
    
    # Divide sample size by 2: the idea is that the main factor of interest defined the two group sizes. A moderator factor is then added to create a 2*2, but because the moderator is not the main focus, the empirical sample sizes should only be used for the two groups formed by the main factor--not the four groups formed by the 2*2 split.
    #s <- round(s/2)
    #run the first analysis with some QRPs applied
    a = analyB(g1 = G[,,1][1:s,], 
               g2 = G[,,2][1:s,],
               D = G[,,3][1,1],      
               multDV=1,out=1) # AGGRESIVE 
    #define optional stopping parameters for AGGRESSIVE strategy
    colLim = 3 
    add = s*0.1
    #add = 3
    #see if you can benefit from optional stopping
    for (i in 1:colLim){
      if(a[1] > 0 & a[2] < .05){break}
      s=s+add
      a = analyB(g1 = G[,,1][1:s,],
                 g2 = G[,,2][1:s,],
                 D = G[,,3][1,1],     
                 multDV=1,out=1) # AGGRESSIVE
    }
  }
  else{print('ERROR: define QRP strategy')}
  #return the result
  return(a)    
}    

#======================
#   publish_decision   #     
#======================
# First publication bias function: Fixed publication bias probability
# Publication bias function to decide if a study will be published
# If p <.05, then 90% this study will be published
# If p > .05, then there is varying degrees of pub bias dictating a study will be published depending on which option you choose
publish_decision <- function(p_value, d, prob) {
  # Ensure that prob is one of the specified values
  if (!(prob %in% c("low", "medium", "high"))) {
    stop("prob must be 'low', 'medium', or 'high'")
  }
  
  # Define probabilities based on 'prob' argument
  #0.025<=p<0.05 and correct direction
  prob_values <- list(
    low = 1,
#   low = 0.9,    # 90% probability to assign 1
    medium = 0.6, # 80% probability to assign 1
    high = 0.5    # 60% probability to assign 1
  )
  
  #p>=0.05 and correct direction
  prob_values2 <- list(
    low = 1,
#   low = 0.9,    # 90% probability to assign 1
    medium = 0.2, # 20% probability to assign 1
#    high = 0
    high = 0.1    # 0% probability to assign 1
  )
  
  #p<0.05 and incorrect direction or p>=0.05 and incorrect direction
  prob_values3 <- list(
    # low = 0,    # 90% probability to assign 1
    # medium = 0, # 10% probability to assign 1
    # high = 0    # 5% probability to assign 1
    low = 1,    # 90% probability to assign 1
    medium = 0.1, # 10% probability to assign 1
    high = 0.05    # 5% probability to assign 1
  )
  
# #   # Function to assign 'publish' based on p-value and probability
#   assign_publish <- function(p_value, d, prob_values, prob) {
#     #0.025<=p<0.05 and correct direction
#     if (p_value >= 0.025 & p_value < 0.05 & d > 0) {
#       prob_assign <- prob_values[[prob]]
#       return(rbinom(1, 1, prob_assign))
#       #p<0.025 and correct direction
#     }else if (p_value < 0.025 & d > 0) {
#       return(rbinom(1, 1, 1)) # 100% probability to assign 1
#       #p>=0.05 and correct direction
#     }else if (p_value >= 0.05 & d > 0) {
#       prob_assign <- prob_values2[[prob]]
#       return(rbinom(1, 1, prob_assign))
#       #p<0.05 and incorrect direction or p>=0.05 and incorrect direction
#     } else if (d < 0) {
#       prob_assign <- prob_values3[[prob]]
#       return(rbinom(1, 1, prob_assign))
#     }
#   }
#   # Assign and return 'publish' value
#   publish <- assign_publish(p_value, d, prob_values, prob)
#   return(list(p_value = p_value, publish = publish))
# }

## results <- replicate(1000, publish_decision(0.01,-0.1, "high"), simplify = FALSE)
## results_df <- do.call(rbind, lapply(results, as.data.frame))
## sum(results_df$publish)


# Second publication bias function: random publication bias probability 
# This is an adjusted inverse sigmoid function
inverse_sigmoid <- function(x, startX = 0.05, endX = 0.025, x0 = (endX + startX)/2, maxY = 0.9, k = -1000, L = maxY - prob_assign, prob_assign=prob_assign) {
  return((maxY - prob_assign) / (1 + exp(-k * (x - x0)))+prob_assign)
}
# Function to assign 'publish' based on p-value and probability
# 1=publish, 0=not publish
  assign_publish <- function(p_value, d, prob_values, prob) {
    #0.025<=p<0.05 and correct direction
    if (p_value >= 0.025 & p_value < 0.05 & d > 0) {
      prob_assign <- prob_values[[prob]]
      prob_assigned <- inverse_sigmoid(p_value, startX = 0.05, endX = 0.025, maxY = 0.9, k = -1000, L = maxY - prob_assign, prob_assign=prob_assign)
      return(rbinom(1, 1, prob_assigned))
      #p<0.025 and correct direction
    }else if (p_value < 0.025 & d > 0) {
      return(rbinom(1, 1, 1)) # 100% probability to assign 1
      #p>=0.05 and correct direction
    }else if (p_value >= 0.05 & d > 0) {
      prob_assign <- prob_values2[[prob]]
      return(rbinom(1, 1, prob_assign))
      #p<0.05 and incorrect direction or p>=0.05 and incorrect direction
    } else if (d < 0) {
      prob_assign <- prob_values3[[prob]]
      return(rbinom(1, 1, prob_assign))
    }
  }
  # Assign and return 'publish' value
  publish <- assign_publish(p_value, d, prob_values, prob)
  return(list(p_value = p_value, publish = publish))
}
## results <- replicate(1000, publish_decision(0.05,0.1, "high"), simplify = FALSE)
## results_df <- do.call(rbind, lapply(results, as.data.frame))
## sum(results_df$publish)


#==================
#     simMA         #     
#==================
#The current simMA don't inflict any pub bias to generated data 
#simMA does not contain a progress bar because run_func has a progress bar
# simMA <- function(k, delta, tau, qrpEnv= c("none", "medium", "high"),  censorFunc = c("low", "medium", "high"), verbose=FALSE, fixed.n=NULL) {
# 
#   # validate parameters
#   if (length(censorFunc) == 1) {
#     censorFunc <- match.arg(censorFunc, c("low", "medium", "high"))
#   }
#   qrpEnv <- match.arg(qrpEnv, c("none",  "medium", "high"))
#   # Define the QRP environments:
#   # get the proportions of studies produced under each QRP strategy
#   if (qrpEnv == 'none'){
#     noneP = 1; litP = 0; modP = 0; aggP = 0
#   } else if (qrpEnv == 'medium'){
#     noneP = 0.3; litP = 0.5; modP = 0.2; aggP = 0
#   } else if (qrpEnv == 'high'){
#     noneP = 0; litP = 0.3; modP = 0.5; aggP = 0.2
#   } else {
#     print('ERROR: qrpEnv must be none, low, medium, or high')
#   }
#   # if (qrpEnv == 'none'){
#   #   noneP = 1; litP = 0; modP = 0; aggP = 0
#   # } else if (qrpEnv == 'medium'){
#   #   noneP = 0.2; litP = 0.3; modP = 0.4; aggP = 0
#   # } else if (qrpEnv == 'high'){
#   #   noneP = 0; litP = 0.1; modP = 0.4; aggP = 0.5
#   # } else {
#   #   print('ERROR: qrpEnv must be none, low, medium, or high')
#   # }
#   datMA <- data.frame()
# 
#   # repeatedly add a new study from that environment until the desired number of k is achieved
#   repeat {
#     thisStudiesHackingStyle <- sample(x = c("none", "lit", "mod", "agg"), size=1, replace=TRUE, prob = c(noneP, litP, modP, aggP))
# 
#     if (thisStudiesHackingStyle == "none") {
#       res <- simData.noQRP(delta=delta, tau=tau, fixed.n=fixed.n)
#       res[15] = 0 #QRP style
#     } else if (thisStudiesHackingStyle == "lit") {
#       res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="lit", fixed.n=fixed.n)
#       res[15] = 1 #QRP style
#     } else if (thisStudiesHackingStyle == "mod") {
#       res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="mod", fixed.n=fixed.n)
#       res[15] = 2 #QRP style
#     } else if (thisStudiesHackingStyle == "agg") {
#       res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="agg", fixed.n=fixed.n)
#       res[15] = 3 #QRP style
#     }
# 
# 
#     # inflict publication bias via the censoring function
#     if (is.character(censorFunc) && censorFunc == "low") {
#       publish <- publish_decision(res[2], res[1], "low")
#     } else if (is.character(censorFunc) && censorFunc == "medium") {
#       publish <- publish_decision(res[2], res[1], "medium")
#     } else if (is.vector(censorFunc) && censorFunc =="high") {
#       publish <- publish_decision(res[2], res[1], "high")
#     } else {
#       stop("Wrong specification of censor function!")
#     }
# 
#     if (publish[2] == 1) {
#       datMA <- rbind(datMA, res)
#       if (verbose==TRUE) print(nrow(datMA))
#     }
# 
#     if (nrow(datMA) >= k) {break}
# 
#   } # of repeat
# 
#   #name columnes
#   colnames(datMA) = c( 'd',       # effect size, d
#                        'p',       # p value for the two group comparison
#                        't',       # t value for the two group comparison
#                        'N',       # total N
#                        'v',       # variance for the effect size
#                        'se',      # standard error for the effect size
#                        'pow',     # power given the true effect for the two group comparison
#                        'n1',      # experimental group sample size
#                        'n2',      # control group sample size
#                        'delta_i', # the study-level true effect
#                        'm1',
#                        'v1',
#                        'm2',
#                        'v2',
#                        'qrp')     # 0 = 'none', 1 = 'mod', 2 = 'agg'
# 
# 
#   # Add Hedge's correction factor
#   df = datMA$n1 + datMA$n2 - 2
#   J = 1- 3/(4*df - 1)
#   datMA$g = datMA$d*J
#   datMA$g_v = datMA$v*J^2
#   datMA$g_se = sqrt(datMA$g_v)
# 
#   return(datMA)
# }
#The following function is the same function as simMA except it has a progress bar
#old name: simMA_pbar
#I use this function to generate original studies
simORIGs <- function(num_deltai, delta, tau, qrpEnv= c("none", "medium", "high"),  censorFunc = c("low", "medium", "high"), verbose=FALSE, fixed.n=NULL) {  
  
  # Display a message with the current combination of parameters at the start
  message(sprintf("Running simORIGs with num_deltai = %s, delta = %s, tau = %s, qrpEnv = '%s', censorFunc = '%s', fixed.n = %s", 
                  num_deltai, delta, tau, qrpEnv, censorFunc, ifelse(is.null(fixed.n), "NULL", fixed.n)))
  
  
  # validate parameters
  if (length(censorFunc) == 1) {
    censorFunc <- match.arg(censorFunc, c("low", "medium", "high"))
  }
  qrpEnv <- match.arg(qrpEnv, c("none",  "medium", "high"))
  # Define the QRP environments:
  # get the proportions of studies produced under each QRP strategy
  if (qrpEnv == 'none'){
    noneP = 1; litP = 0; modP = 0; aggP = 0
  } else if (qrpEnv == 'medium'){
    noneP = 0.3; litP = 0.5; modP = 0.2; aggP = 0
  } else if (qrpEnv == 'high'){
    noneP = 0; litP = 0.3; modP = 0.5; aggP = 0.2
  } else {
    print('ERROR: qrpEnv must be none, low, medium, or high')
  }
  
  datMA <- data.frame()
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = num_deltai, style = 3)
  
  # repeatedly add a new study from that environment until the desired number of num_deltai is achieved
  repeat {		
    thisStudiesHackingStyle <- sample(x = c("none", "lit", "mod", "agg"), size=1, replace=TRUE, prob = c(noneP, litP, modP, aggP))
    
    if (thisStudiesHackingStyle == "none") {
      res <- simData.noQRP(delta=delta, tau=tau, fixed.n=fixed.n)
      res[15] = 0 #QRP style
    } else if (thisStudiesHackingStyle == "lit") {
      res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="lit", fixed.n=fixed.n)
      res[15] = 1 #QRP style
    } else if (thisStudiesHackingStyle == "mod") {
      res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="mod", fixed.n=fixed.n)
      res[15] = 2 #QRP style
    } else if (thisStudiesHackingStyle == "agg") {
      res <- simData.QRP(delta=delta, tau=tau, QRP.strategy="agg", fixed.n=fixed.n)
      res[15] = 3 #QRP style
    }
    
    
    # inflict publication bias via the censoring function
    if (is.character(censorFunc) && censorFunc == "low") {
      publish <- publish_decision(res[2], res[1], "low")
    } else if (is.character(censorFunc) && censorFunc == "medium") {
      publish <- publish_decision(res[2], res[1], "medium")
    } else if (is.vector(censorFunc) && censorFunc =="high") {
      publish <- publish_decision(res[2], res[1], "high")
    } else {
      stop("Wrong specification of censor function!")
    }
    
    if (publish[2] == 1) {
      datMA <- rbind(datMA, res)
      if (verbose==TRUE) print(nrow(datMA))
    }
    
    # Update progress bar
    setTxtProgressBar(pb, nrow(datMA))
    
    if (nrow(datMA) >= num_deltai) {break}
    
  } # of repeat
  
  # Close progress bar
  close(pb)
  
  #name columnes
  colnames(datMA) = c( 'd',       # effect size, d
                       'p',       # p value for the two group comparison
                       't',       # t value for the two group comparison
                       'N',       # total N
                       'v',       # variance for the effect size
                       'se',      # standard error for the effect size
                       'pow',     # power given the true effect for the two group comparison
                       'n1',      # experimental group sample size
                       'n2',      # control group sample size
                       'delta_i', # the study-level true effect
                       'm1',
                       'v1',
                       'm2',
                       'v2',
                       'qrp')     # 0 = 'none', 1 = 'mod', 2 = 'agg'
  
  
  # Add Hedge's correction factor
  df = datMA$n1 + datMA$n2 - 2
  J = 1- 3/(4*df - 1)
  datMA$g = datMA$d*J
  datMA$g_v = datMA$v*J^2
  datMA$g_se = sqrt(datMA$g_v)											 
  
  return(datMA)
}

