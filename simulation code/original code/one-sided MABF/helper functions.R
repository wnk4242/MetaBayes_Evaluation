# 08022025 Update: Limited prior in EUBFcalc to only positive (i.e., changed H1: mu != 0 to H1: mu > 0)
# 0705 Update: Added a new function combine_matrix for FEMA outcome combinations
# 0524 Update: Changed the True positive criteria in ratesCalc_TPFN: added ES must be bigger than 0
# 0523 Update: Added some new functions to create bar plot; added some new functions to calculate rates for original studies
# 0516 Update: Added a new function MAcalc_betap for conducting MA in Step 1.2
# 0515 Update: Updated ratesCalc2, now can calculate TSR, FFR, TFR, FSR when true effect is null or not
# 0511 Update: Corrected mistakes in rateCalc2
# 0502 Update: Added two naming functions for MABF outcomes
# 0421 Update: Corrected a serious bug in run_func()
# 0417 Update: run_func() is updated with foreach
# 0415 Update: progress bar in run_func() is removed
# 0408 Update: a new iBFcalc() function
# 0409 Update: added explaination to iBFcalc() function change
# Note I changed k ot num_deltai
# Function used to run a function for a pre-specified number of times
# This function is only used in beta testing
# run_and_record_df <- function(func, times, ...) {
#   # Initialize a list to store results
#   results_list <- vector("list", times)
#   # Create a progress bar
#   pb <- txtProgressBar(min = 0, max = times, style = 3)
#   # Loop through the specified number of times
#   for (i in 1:times) {
#     # Run the provided function and store the result
#     results_list[[i]] <- func(...)
#     # Update the progress bar
#     setTxtProgressBar(pb, i)
#   }
#   # Close the progress bar
#   close(pb)
#   # Return the list of results
#   return(results_list)
# }

# Function used to run a function for a pre-specified number of times
# run_func and run_combs are used to generate data
# For example, together they can generate meta-analyses
# run_func <- function(func, times,...) {
#   # Initialize a list to store results
#   results_list <- vector("list", times)
#   # Loop through the specified number of times
#   for (i in 1:times) {
#     # Run the provided function and store the result
#     results_list[[i]] <- func(...)
#   }
#   # Return the list of results
#   return(results_list)
# }
#Following run_func uses foreach instead of for
#wrap func(...) in a list() so that the results
#of each iteration are stored in their own individual lists
run_func <- function(func, times, ...) {
  # Use foreach with %dopar% for parallel execution
  results_list <- foreach(i = 1:times, .combine = 'c') %dopar% {
    list(func(...))
  }
  # Return the list of results
  return(results_list)
}
# Function used to receive different combinations of parameters and deliver them to run_func
# old name: run_combs
# This function will not be used in dissertation
# pass_params_oirg <- function(times, k, delta, tau, fixed.n, qrpEnv, censorFunc) {
#   # Create a string describing the current parameter combination
#   # current_combination <- paste("k:", k, "|", "delta:", delta, "|", "tau:", tau, "|", "fixed.n:", fixed.n, "|", "qrpEnv:", qrpEnv, "|", "censorFunc:", censorFunc, sep="")
#   current_combination <- paste("censorFunc:", censorFunc, "|", "qrpEnv:", qrpEnv, "|", "fixed.n:", fixed.n, "|", "tau:", tau, "|", "delta:", delta, "|", "k:", k, sep="")
#   # Run simMA for the current combination of parameters, passing the combination string
#   df_lists <- run_func(simMA, current_combination = current_combination, times = times, k = k, delta = delta, tau = tau, qrpEnv = qrpEnv, censorFunc = censorFunc, verbose = FALSE, fixed.n = fixed.n)
# }
# Function used to receive different combinations of parameters and deliver them to run_func
pass_params <- function(times, num_deltai, delta_i, fixed.n) {
  # Run simMA for the current combination of parameters, passing the combination string
  df_lists <- run_func(simREPs, times = times, num_deltai = num_deltai, delta_i = delta_i, verbose = FALSE, fixed.n = fixed.n)
}
# Function used to name the lists created by run_combs() based on parameter combinations
# old name: name_Lists
# deleted k in params_id
naming_orig <- function(times, num_deltai, delta, tau, fixed.n, qrpEnv, censorFunc) {
  params_id <- paste(delta, tau, fixed.n, qrpEnv, censorFunc, sep = "_")
}
# Function used to name replications
naming_rep <- function(times, delta_i, num_deltai, fixed.n) {
  params_id <- paste(num_deltai, fixed.n, round(delta_i, digits=6), sep = "_")
}
# Function used to name MABF outcomes
naming_MABFlists <- function(delta, orig.n, qrpEnv, censorFunc) {
  params_id <- paste(delta, orig.n, qrpEnv, censorFunc, sep = "_")
}
naming_MABFsubls <- function(rep.number, rep.n) {
  params_id <- paste(rep.number, rep.n, sep = "_")
}
# Function used to generate a single replication
genREP <- function(delta_i, fixed.n=NULL){
  n <- fixed.n
  sigma <- 1
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
  out = c(d, p, t, n1+n2, d_v, d_se, pwr, n1, n2, delta_i, m1, v1, m2, v2)  
}

# Function used to simulate replications
simREPs<- function(num_deltai, delta_i, verbose=FALSE, fixed.n=NULL) {  
  # Run the genREP function in parallel k times
  results <- replicate(num_deltai, genREP(delta_i=delta_i, fixed.n=fixed.n), simplify = FALSE)
  # Combine all results into a single data frame
  datMA <- do.call(rbind, results)
  datMA <- as.data.frame(datMA)
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
                       'v2')    
  # Add Hedge's correction factor
  df = datMA$n1 + datMA$n2 - 2
  J = 1- 3/(4*df - 1)
  datMA$g = datMA$d*J
  datMA$g_v = datMA$v*J^2
  datMA$g_se = sqrt(datMA$g_v)											 
  return(datMA)
}

# Function used to categorize BF
# Count how many values fall within the interval (inclusive)
fallbtw <- function(vector, lower_bound, upper_bound) {
  count_within_interval <- sum(vector >= lower_bound & vector < upper_bound)
  proportion_within_interval <- count_within_interval / length(vector)
  return(proportion_within_interval)
}
# # Count how many values equal to 1
# countones <- function(vector) {
#   count_ones <- sum(vector = 1)
#   # Calculate the proportion of values within the interval
#   proportion_ones <- count_ones / length(vector)
#   # Return the proportion
#   return(proportion_ones)
# }
# Turn BFs into an evidence strength table
EviStrTab <- function(EUBFlist){
  # Proportions of BF according to categories of evidence strength 
  BF10Anecdotal <- fallbtw(EUBFlist, 1, 3)
  BF10Moderate <- fallbtw(EUBFlist, 3, 10)
  BF10Strong <- fallbtw(EUBFlist, 10, 30)
  BF10VStrong <- fallbtw(EUBFlist, 30, 100)
  BF10Etreme <- fallbtw(EUBFlist, 100, Inf)
  BF01Anecdotal <- fallbtw(EUBFlist, 0.33, 1)
  BF01Moderate <- fallbtw(EUBFlist, 0.1, 0.33)
  BF01Strong <- fallbtw(EUBFlist, 0.03, 0.1)
  BF01VStrong <- fallbtw(EUBFlist, 0.01, 0.03)
  BF01Etreme <- fallbtw(EUBFlist, 0, 0.01)
  # Combine the results into a data frame
  evidence_strength_table <- data.frame(
    Category = c("Anecdotal", "Moderate", "Strong", "Very Strong", "Xtreme"),
    BF10 = c(BF10Anecdotal, BF10Moderate, BF10Strong, BF10VStrong, BF10Etreme),
    BF01 = c(BF01Anecdotal, BF01Moderate, BF01Strong, BF01VStrong, BF01Etreme)
  )
  return(evidence_strength_table)
}

# Function used to generate bar plots for each table in a list
# genBarplots0 is used for MABFs except for Inclusion Bayes factor of Bayesian Averaged MA
genBarplots0 <- function(results_list, color_BF10 = "navy", color_BF01 = "orange") {
  tables <- lapply(results_list, EviStrTab)
  for (key in names(tables)) {
    # Reshape the current data frame from wide to long format
    df_long <- pivot_longer(tables[[key]], cols = c(BF10, BF01), names_to = "Evidence", values_to = "Proportion")
    # Create bar plot for the current data frame
    plot <- ggplot(df_long, aes(fill = Evidence, y = Proportion, x = Category)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal() +
      labs(title = paste("Proportion of Evidence Strength (", key, ")", sep = ""),
           x = "Strength of Evidence",
           y = "Proportion") +
      scale_fill_manual(values = c("BF10" = color_BF10, "BF01" = color_BF01)) +
      theme(text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5))
    # Print the plot
    print(plot)
  }
}
#genBarplots is used for generate bar plots for original studies, which will not be used for dissertation
genBarplots <- function(results_list, params_df, color_BF10 = "navy", color_BF01 = "orange") {
  for (i in seq_along(results_list)) {
    tables <- EviStrTab(results_list[[i]])
    # Reshape the current data frame from wide to long format
    df_long <- pivot_longer(tables, cols = c(BF10, BF01), names_to = "Evidence", values_to = "Proportion")
    # Extract parameters for the current combination
    current_params <- params_df[i, ]
    # Create a string representing parameter legend
    params_legend <- paste("num_deltai = ", current_params$num_deltai, "\n", "delta = ", current_params$delta, "\n", "tau = ", current_params$tau, "\n", "fixed.n = ", current_params$fixed.n, "\n", "QRP level = ", current_params$qrpEnv, "\n", "PB level = ", current_params$censorFunc, sep="")
    # Add parameter combination as a legend
    plot_title <- paste("Proportion of Evidence Strength", sep = "\n")
    # Create bar plot for the current data frame
    plot <- ggplot(df_long, aes(fill = Evidence, y = Proportion, x = Category)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal() +
      labs(title = plot_title,
           x = "Strength of Evidence",
           y = "Proportion",
           fill = NULL) +
      scale_fill_manual(values = c("BF10" = color_BF10, "BF01" = color_BF01)) +
      theme(text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            legend.box.background = element_rect(color = "black", size = 1)) +
      guides(fill = guide_legend(title = params_legend, label.position = "right", title.position = "top", title.theme = element_text(size = 12), label.theme = element_text(angle = 0, vjust = 1)))
    # Print the plot
    print(plot)
  }
}
# genRepBarplots is used for replications, which will be used for dissertation
genRepBarplots <- function(results_list, params_df, color_BF10 = "navy", color_BF01 = "orange") {
  for (i in seq_along(results_list)) {
    tables <- EviStrTab(results_list[[i]])
    # Reshape the current data frame from wide to long format
    df_long <- pivot_longer(tables, cols = c(BF10, BF01), names_to = "Evidence", values_to = "Proportion")
    # Extract parameters for the current combination
    current_params <- params_df[i, ]
    # Create a string representing parameter legend
    params_legend <- paste("num_deltai = ", current_params$num_deltai, "\n", "delta = ", current_params$delta, "\n", "fixed.n = ", current_params$fixed.n, "\n", sep="")
    # Add parameter combination as a legend
    plot_title <- paste("Proportion of Evidence Strength", sep = "\n")
    # Create bar plot for the current data frame
    plot <- ggplot(df_long, aes(fill = Evidence, y = Proportion, x = Category)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal() +
      labs(title = plot_title,
           x = "Strength of Evidence",
           y = "Proportion",
           fill = NULL) +
      scale_fill_manual(values = c("BF10" = color_BF10, "BF01" = color_BF01)) +
      theme(text = element_text(size = 12),
            plot.title = element_text(hjust = 0.5),
            legend.box.background = element_rect(color = "black", size = 1)) +
      guides(fill = guide_legend(title = params_legend, label.position = "right", title.position = "top", title.theme = element_text(size = 12), label.theme = element_text(angle = 0, vjust = 1)))
    # Print the plot
    print(plot)
  }
}

# beta testing function used to test the inclusion BF based on Bayesian Averaged Meta-analysis
# This function is used in betTesting.R. It will not be used for final simulations.
# This simBAMA comes with a progress bar
simBAMA <- function(n_iterations, baMA_args, qrpEnv_levels, censorFunc_levels, simMA_args) {
  # baMA_args is used to load priors, simMA_args is used to specify arugments in simMA() other than QRP and publicatoin bias
  # Calculate total number of simulations for progress bar (combinations * iterations)
  total_sims <- length(qrpEnv_levels) * length(censorFunc_levels) * n_iterations
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = total_sims, style = 3)
  progress_counter <- 0 # To track progress across all loops
  # Initialize a list to store results, structured by qrpEnv and censorFunc levels
  results_list <- list()
  for(qrpEnv in qrpEnv_levels) {
    for(censorFunc in censorFunc_levels) {
      combination_key <- paste(qrpEnv, censorFunc, sep = "_")
      results_list[[combination_key]] <- vector("list", n_iterations)
      for(i in 1:n_iterations) {
        # Update simMA_args with current qrpEnv and censorFunc without altering the original list
        temp_simMA_args <- c(simMA_args, list(qrpEnv = qrpEnv, censorFunc = censorFunc))
        # Generate a new dataset with dynamic simMA arguments
        d <- do.call("simMA", temp_simMA_args)
        # Conduct Bayesian model averaged meta-analysis for each simulation
        baMA <- do.call("meta_bma", c(list(y = d$d, SE = d$se), baMA_args))[[9]]$incl.BF
        # Store each baMA result
        results_list[[combination_key]][[i]] <- baMA
        # Update progress counter and progress bar
        progress_counter <- progress_counter + 1
        setTxtProgressBar(pb, progress_counter)
      }
    }
  }
  # Close the progress bar after the loop is done
  close(pb)
  # Return the list of results
  return(results_list)
}


# Calculate Inclusion BF of Bayesian Averaged MA 
# This version conduct complete analysis which is very time-consuming
# iBFcalc <- function(df_list){
#    #p <- progressor(along = df_list)
#     future_lapply(df_list, function(x) {
#    #p()  # Signal progress
#     meta_bma(y = x$d, SE = x$se, d = priorEStesting, tau = priorTau, summarize = "integrate")[[9]]$incl.BF
#   })
# }
# This version only calculates inclusion Bayes Factor but use both fixed and random effects models
# iBFcalc <- function(df_list){
#   future_lapply(df_list, function(x) {
#     fixed=meta_fixed(y = x$d, SE = x$se, d=priorEStesting)
#     random=meta_random(y = x$d, SE = x$se, d=priorEStesting, tau=priorTau, summarize = 'integrate')
#     inclusion(list(fixed,random),include="H1")$incl.BF
#   })
# }
# This version only calculates inclusion Bayes Factor but only use fixed-effect model
# This decision is made based on two rationale: 
# 1.the studies to meta-analyze are replications with very small or null tau, so FE model can be considered appropriate
# 2.this paper:  A tutorial on Bayesian model-averaged meta-analysis in JASP
iBFcalc <- function(df_list){
  future_lapply(df_list, function(x) {
    fixed=meta_fixed(y = x$d, SE = x$se, d=priorEStesting)
    fixed$BF[2]
  })
}


# Calculate Fixed effect MA Bayes factor 
# Added: FEMABFs <- do.call(c,result)
# Changed some variable names
FEMABFcalc <- function(df_list){
  p <- progressor(along = df_list)
  FEMABF <- future_lapply(df_list, function(x) {
    p()  # Signal progress
    t <-  x$t
    n1 <-  x$n1
    n2 <-  x$n2
    # Use meta.ttestBF function
    # meta.ttestBF produces many extra outcomes
    all_outcomes <- meta.ttestBF(t = t, n1 = n1, n2 = n2, nullInterval = c(0, Inf), rscale = sqrt(2)/2)
    # Extract FEMA Bayes factor
    FEMABF_only <- exp(all_outcomes[1]@bayesFactor[1][1,])
  })
  FEMABF <- do.call(c, FEMABF)
}
# Calculate Evidence Updating Bayes factor 
EUBFcalc <- function(df_list){
  p <- progressor(along = df_list)
  future_lapply(df_list, function(x) {
    p()  # Signal progress
    n1 <-  x$n1
    n2 <-  x$n2
    m1 <-  x$m1
    m2 <-  x$m2
    v1 <-  x$v1
    v2 <-  x$v2
    # Use combStudies function
    combStats <- combStudies(n_x = n1, mean_x = m1, var_x = v1, n_y = n2, mean_y = m2, var_y = v2)
    # Use ttest.tstat function
    EUBF <- ttest.tstat(t = combStats[1], n1 = combStats[2], n2 = combStats[3], nullInterval = c(0, Inf), rscale = sqrt(2)/2, simple = TRUE)
    unname(EUBF)
  })
}

# Calculate meta-analysis statistics for BFbMAcalc()
MAcalc <- function(df_list){
  p <- progressor(along = df_list)
  future_lapply(df_list, function(x) {
    #p()
    d <-  x$d
    vi <-  x$v
    out <- rma(yi=d, vi=v, data=x, method = "DL")
    #  return(list(out$pval, out$se))
  })
}

# Only obtain observed effect size and its p
MAcalc_betap <- function(df_list){
  p <- progressor(along = df_list)
  future_lapply(df_list, function(x) {
    #p()
    d <-  x$d
    vi <-  x$v
    out <- rma(yi=d, vi=v, data=x, method = "DL")
    return(list(beta=out$beta,pval=out$pval))
    #return(out$pval)
  })
}


# Calculate Bayes factor based on meta-analysis
BFbMAcalc<- function(ma_list){
  p <- progressor(along = ma_list)
  future_lapply(ma_list, function(x) {
    p()
    p_value <-  x$pval
    var_theta <- (x$se)^2
    BFbMA <- BFroMeta(p_value=p_value, var_theta = var_theta, theta_A = 0.1)
    #  return(list(out$pval, out$se))
  })
}


# Function to combine vectors from meta-analytic BF lists (e.g., FEMABF_lists)
# The MABFs are saved as lists of subgroups (subgroup contain MABF obtained under a unique parameter combination). 
# I want to group every two sublists as a group. 
# So, sublist1 and sublist2 is now a group, sublist3 and sublist4 is a group, and so on. 
# In each group, I want to combine the rows of the first num_deltai (set to 1000) vectors in the first sublist with 
# the first num_deltai vectors in the second sublist, and then combine the second num_deltai vectors in subgroup1 with 
# the second num_deltai vectors in the subgroup2, and so on.
# Here, the chunk_size is equal to num_deltai
combine_vectors <- function(sublist1, sublist2, chunk_size) {
  combined_list = list()  # Initialize an empty list to store combined vectors
  # Loop through the vectors according to the specified chunk size
  for(i in seq(1, length(sublist1), by = chunk_size)) {
    # Use rbind to combine vectors i to i+chunk_size-1 from both sublists
    combined_vectors = do.call(rbind, c(sublist1[i:(i+chunk_size-1)], sublist2[i:(i+chunk_size-1)]))
    combined_list[[length(combined_list) + 1]] = combined_vectors
  }
  return(combined_list)
}

# Function to create bar plots for origin studies showing their classification rates
# You can select the parameter to be QRP level (qrpEnv) or original sample size (orig.n.level)
barplot_XXRorig_overall <- function(trueES, ratetype, param, dataset) {
  dataset <- dataset %>%
    #filter(orig.n.level=="large")%>%
    group_by(.data[[param]]) %>%
    summarise(Average_RATE = mean(.data[[ratetype]], na.rm = TRUE))
  # Create the plot
  plot <- ggplot(dataset, aes(x = .data[[param]], y = Average_RATE, fill = .data[[param]])) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percent(Average_RATE, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), color = "white", size = 5) +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title = paste0("Overall ", ratetype, " of the significance method when the true effect size is ", trueES),
         x = "QRP Environment",
         y = paste0("Overall ", ratetype, " (%)"),
         fill = "QRP Level") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(plot)
}

#Functions to create bar plots for classification accuracy rates by MABF method#
# barplot_XXRbyMABF_overall calculate all types of overall rates and creates bar plots by MABF methods
barplot_XXRbyMABF_overall <- function(trueES, ratetype, dataset) {
  # Group by method and calculate the average ADR
  dataset <- dataset %>%
    group_by(method) %>%
    summarise(Average_RATE = mean(.data[[ratetype]], na.rm = TRUE))
  # Create the plot
  plot <- ggplot(dataset, aes(x = method, y = Average_RATE, fill = method)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percent(Average_RATE, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), color = "white", size = 5) +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title = paste0("Overall ", ratetype, " by each MABF when the true effect size is ", trueES),
         x = "Method",
         y = paste0("Overall ", ratetype, " (%)"),
         fill = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(plot)
}

# barplot_XXRbyMABF_cond calculate all types of conditional rates, creates bar plots by MABF methods, conditioned on various factors
barplot_XXRbyMABF_cond <- function(trueES, ratetype, dataset, paramcomb = list()) {
  # Apply filters conditionally based on non-NULL parameters in the list
  if (!is.null(paramcomb$repsampsize)) {
    dataset <- dataset %>% filter(rep.n.level == paramcomb$repsampsize)
  }
  if (!is.null(paramcomb$repnumber)) {
    dataset <- dataset %>% filter(rep.number.level == paramcomb$repnumber)
  }
  if (!is.null(paramcomb$origsampsize)) {
    dataset <- dataset %>% filter(orig.n.level == paramcomb$origsampsize)
  }
  if (!is.null(paramcomb$qrplevel)) {
    dataset <- dataset %>% filter(qrpEnv == paramcomb$qrplevel)
  }
  # Group by method and calculate the average ADR
  dataset <- dataset %>%
    group_by(method) %>%
    summarise(Average_RATE = mean(.data[[ratetype]], na.rm = TRUE))
  # Create the plot
  plot <- ggplot(dataset, aes(x = method, y = Average_RATE, fill = method)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = percent(Average_RATE, accuracy = 0.1)), 
              position = position_stack(vjust = 0.5), color = "white", size = 5) +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title = paste0("Average ", ratetype, " by each MABF when the true effect size is ", trueES,
                        if (!is.null(paramcomb$repsampsize)) paste0("\nReplication sample size: ", paramcomb$repsampsize, ""),
                        if (!is.null(paramcomb$repnumber)) paste0("\nNumber of replications: ", paramcomb$repnumber, ""),
                        if (!is.null(paramcomb$origsampsize)) paste0("\nOriginal sample size: ", paramcomb$origsampsize, ""),
                        if (!is.null(paramcomb$qrplevel)) paste0("\nQRP and PB level: ", paramcomb$qrplevel, "")),
         x = "Method",
         y = paste0("Average ", ratetype, " (%)"),
         fill = "Method") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
}

# Function to generate and save bar plots for each parameter combination
savebarplot_XXRbyMABF_cond <- function(...) {
  paramcomb <- list(...)
  directory <- dir2save
  dataset <- data2plot
  ratetype <- ratetype2calc
  trueES <- TrueES
  plot <- barplot_XXRbyMABF_cond(trueES,ratetype,dataset,paramcomb)
  # Dynamically construct the filename based on non-NULL parameters
  filename_parts <- lapply(names(paramcomb), function(name) {
    if (!is.null(paramcomb[[name]])) {
      paste0(name, "_", paramcomb[[name]])
    }
  })
  filename <- paste0(directory,"/",ratetype,"/",ratetype,"_",paste(na.omit(unlist(filename_parts)), collapse = "_"), ".jpg")
  ggsave(filename = filename, plot = plot, width = 8, height = 6)
}

#Combine beta and p matricies in sublists for FEMA outcomes
combine_matrices <- function(matrices1, matrices2) {
  mapply(cbind, matrices1, matrices2, SIMPLIFY = FALSE)
}