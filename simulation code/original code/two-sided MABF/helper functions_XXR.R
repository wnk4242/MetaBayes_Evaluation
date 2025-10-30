# Following helper functions are dedicated to calculating classification rates with different criteria
# Yof2025
# 0327 Update: Added processMA_sublist_AiO, processMA_matrix_AiO, ratesCalcMA_ROC_AiO
# Yof2024
# 0727 Update: Added functions for original study ES
# 0720 Update: Updated row-wise ratesCalc function
# 0715 Update: Added plot_line()
# 0709 Update: Added ratesCalc functions for ROC plot for fixed-effects meta-analysis (FEMA)
# 0704 Update: Updated the ratesCalc function family according to updated confusion matrix (the original study direction part)
# 0615 Update: Modified all functions to allow for change of criteria: added BFcutoff, pcutoff, EScutoff arguments
# 0604 Update: Added function: ratesCalc_ROC
# 0526 Update: In ratesCalc_TPFN, added the direction of the observed ES as a new criteria for true positives  and false negatives
#              In ratesCalc2, added the direction of the observed ES as a new criteria for TS1, FF1, TF1, FS1 in 

#################################################################
#Function to calculate rates for creating ROC plot (matrix-wise)#
#################################################################
# ratesCalc_ROC function
# This function is used to calculate true/false positive/negative rates used fro ROC plot
# It will run though different threshold/cutoff for every condition in which MABF is calculated
ratesCalc_ROC <- function(bf_matrix, BFcutoff) {
  BFs_no_effect <- bf_matrix[1:num_deltai, ]
  BFs_true_effect <- bf_matrix[(num_deltai + 1):(num_deltai * 2), ]
  TN <- sum(BFs_no_effect < BFcutoff)
  FP <- sum(BFs_no_effect > BFcutoff)
  TP <- sum(BFs_true_effect > BFcutoff)
  FN <- sum(BFs_true_effect < BFcutoff)
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  TNR <- TN / (TN + FP)
  FNR <- FN / (FN + TP)
  data.frame(TPR = TPR, FPR = FPR, TNR = TNR, FNR = FNR, TP = TP, FN = FN, FP = FP, TN = TN)
}

# Function to process each matrix
# nested in process_sublist
process_matrix <- function(cutoff, bf_matrix, sublist_name, matrix_name) {
  rates <- ratesCalc_ROC(bf_matrix, cutoff)
  rates$Sublist <- sublist_name
  rates$Matrix <- matrix_name
  rates$threshold <- cutoff
  return(rates)
}

# Function to process each sublist in XXBF_lists_0.xnull_regrouped
# e.g. EUBF_lists_0.2null_regrouped
process_sublist <- function(sublist, sublist_name) {
  map_dfr(names(sublist), function(matrix_name) {
    message("  Processing matrix: ", matrix_name)
    bf_matrix <- sublist[[matrix_name]]
    map_dfr(BFcutoff, function(cutoff) {
      process_matrix(cutoff, bf_matrix, sublist_name, matrix_name)
    })
  })
}

# ratesCalc_ROC_AiO function, which calculate all kinds of rates except ADR
# AiO is All in One
# Can specify cutoffs
ratesCalc_ROC_AiO <- function(bf_matrix, BFcutoff, pcutoff, EScutoff) {
  # Split the matrix into no effect and true effect
  BFs_no_effect_nop <- bf_matrix[1:num_deltai, 4:c(rep_times+3)] #the BFs in no effect rows
  BFs_true_effect_nop <- bf_matrix[c(num_deltai+1):c(num_deltai*2), 4:c(rep_times+3)] #the BFs in true effect rows
  BFs_no_effect_p <- bf_matrix[1:num_deltai, 1:c(rep_times+3)] #the BFs in no effect rows
  BFs_true_effect_p <- bf_matrix[c(num_deltai+1):c(num_deltai*2), 1:c(rep_times+3)] #the BFs in true effect rows

  # Calculate TPR, TNR, FPR, FNR
  TN <- sum(BFs_no_effect_nop < BFcutoff)
  FP <- sum(BFs_no_effect_nop > BFcutoff)
  TP <- sum(BFs_true_effect_nop > BFcutoff)
  FN <- sum(BFs_true_effect_nop < BFcutoff)
  # TPR <- TP / (TP + FN)
  # FPR <- FP / (FP + TN)
  # TNR <- TN / (TN + FP)
  # FNR <- FN / (FN + TP)

  # Calculate counts for true success and false failure in true effect
  TS1 <- sum(BFs_true_effect_p[,2] < pcutoff & BFs_true_effect_p[,3] > EScutoff & BFs_true_effect_nop > BFcutoff) #Added ES > EScutoff
  FF1 <- sum(BFs_true_effect_p[,2] > pcutoff & BFs_true_effect_nop > BFcutoff) 

  # Calculate counts for true failure and false success in true effect
  TF1 <- sum(BFs_true_effect_p[,2] < pcutoff & BFs_true_effect_p[,3] > EScutoff & BFs_true_effect_nop < BFcutoff) #Added ES > EScutoff
  FS1 <- sum(BFs_true_effect_p[,2] > pcutoff & BFs_true_effect_nop < BFcutoff) 

  # Calculate counts for false success and true failure in no effect
  FS2 <- sum(BFs_no_effect_p[,2] < pcutoff & BFs_no_effect_p[,3] > EScutoff & BFs_no_effect_nop > BFcutoff) #Added ES > EScutoff
  TF2 <- sum(BFs_no_effect_p[,2] > pcutoff & BFs_no_effect_nop > BFcutoff)

  # Calculate counts for false failure and true success in no effect
  FF2 <- sum(BFs_no_effect_p[,2] < pcutoff & BFs_no_effect_p[,3] > EScutoff & BFs_no_effect_nop < BFcutoff) #Added ES > EScutoff
  TS2 <- sum(BFs_no_effect_p[,2] > pcutoff & BFs_no_effect_nop < BFcutoff)

  # Calculate successful/failed affirmation/correction rates
  # SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate
  # FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate
  # FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate
  # SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate
  # SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
  # FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
  # FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
  # SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
  # SAR = (TS1 + TS2)/(TS1 + TS2 + TF1 + TF2)
  # FAR = (TF1 + TF2)/(TS1 + TS2 + TF1 + TF2)
  # FCR = (FS1 + FS2)/(FS1 + FS2 + FF1 + FF2)
  # SCR = (FF1 + FF2)/(FS1 + FS2 + FF1 + FF2)

  # Calculate overall successes and failures
  Successes <- TS1 + TS2 + FS1 + FS2
  Failures <-  TF1 + TF2 + FF1 + FF2

  # Split successes and failures by true and false effect
  Successes1 <- TS1 + FS1
  Failures1 <-  TF1 + FF1
  Successes2 <- TS2 + FS2
  Failures2 <-  TF2 + FF2
  
  # Calcualte True/False Success/Failure rates
  # TSR = (TS1 + TS2)/Successes
  # FSR = (FS1 + FS2)/Successes
  # TFR = (TF1 + TF2)/Failures
  # FFR = (FF1 + FF2)/Failures
  # 
  # TSR1 = TS1/Successes1
  # FSR1 = FS1/Successes1
  # TFR1 = TF1/Failures1
  # FFR1 = FF1/Failures1
  # 
  # TSR2 = TS2/Successes2
  # FSR2 = FS2/Successes2
  # TFR2 = TF2/Failures2
  # FFR2 = FF2/Failures2
  
  # TFSF = true/false/success/failure
  TFSFrates_df <- data.frame(
    # TPR,FPR,TNR,FNR,
    # 
    # SAR, FAR, SCR, FCR,
    # SAR1, FAR1, SCR1, FCR1,
    # SAR2, FAR2, SCR2, FCR2,
    # 
    # TSR,FSR,TFR,FFR,
    # TSR1,FSR1,TFR1,FFR1,
    # TSR2,FSR2,TFR2,FFR2,
    
    TP,FN,FP,TN,
    Successes,Failures,Successes1,Failures1,Successes2,Failures2,
    TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2
  )
}
# Function to process each matrix
# nested in process_sublist
process_matrix_AiO <- function(BFcutoff, pcutoff, EScutoff, bf_matrix, sublist_name, matrix_name) {
  rates <- ratesCalc_ROC_AiO(bf_matrix, BFcutoff, pcutoff, EScutoff)
  rates$Sublist <- sublist_name
  rates$Matrix <- matrix_name
  rates$threshold_BF <- BFcutoff
  rates$threshold_p <- pcutoff
  rates$threshold_ES <- EScutoff
  return(rates)
}

# Function to process each sublist in XXBF_lists_0.xnull_regrouped
# e.g. EUBF_lists_0.2null_regrouped
process_sublist_AiO <- function(sublist, sublist_name) {
  map_dfr(names(sublist), function(matrix_name) {
    message("  Processing matrix: ", matrix_name)
    bf_matrix <- sublist[[matrix_name]]
    map_dfr(cutoff, function(cutoff) {
      process_matrix_AiO(cutoff$BFcutoff, cutoff$pcutoff, cutoff$EScutoff, bf_matrix, sublist_name, matrix_name)
    })
  })
}

# BFs_true_effect_p[,2] represents original p value; BFs_true_effect_p[,3] represents original observed ES; BFs_true_effect_nop represents BF
# Function to calculate TNR and FPR for the original study based on p value (contained in a single sublist)
ratesCalc_TNFP <- function(sublist) {
  # Ensure p-values are numeric
  p_values <- as.numeric(sublist[["p"]])
  es_values <- as.numeric(sublist[["d"]])
  # Determine TN and FP
  TN <- sum(p_values > 0.05, na.rm = TRUE)
  FP <- sum(p_values < 0.05 & es_values > 0, na.rm = TRUE) #Added ES > 0
  # Calculate TNR and FPR
  TNR <- TN / (FP + TN)
  FPR <- FP / (FP + TN)
  # Return as a dataframe with numeric values
  TNFP_orig <- data.frame(TNR = as.numeric(TNR), FPR = as.numeric(FPR), TN = as.numeric(TN), FP = as.numeric(FP))
  return(TNFP_orig)
}
# Function to calculate TPR and FNR for the original study based on p value (contained in a single sublist)
ratesCalc_TPFN <- function(sublist) {
  # Ensure p-values are numeric
  p_values <- as.numeric(sublist[["p"]])
  es_values <- as.numeric(sublist[["d"]])
  # Determine TP and FN
  TP <- sum(p_values < 0.05 & es_values > 0, na.rm = TRUE) #Added ES > 0
  FN <- sum(p_values > 0.05, na.rm = TRUE) 
  # Calculate TPR and FNR
  TPR <- TP / (FN + TP)
  FNR <- FN / (FN + TP)
  # Return as a dataframe with numeric values
  TPFN_orig <- data.frame(TPR = as.numeric(TPR), FNR = as.numeric(FNR), TP = as.numeric(TP), FN = as.numeric(FN))
  return(TPFN_orig)
}
# Function to calculate the average effect size for the original study based on null effect
esCalc_null <- function(sublist) {
  p_values <- as.numeric(sublist[["p"]])
  es_values <- as.numeric(sublist[["d"]])
  # Filter rows based on the criteria
  filtered_sublist <- sublist[!(p_values < 0 & es_values < 0), ]
  # Calculate the average and median of the retained 'd' values
  avg_d <- mean(filtered_sublist[["d"]])
  med_d <- median(filtered_sublist[["d"]])
  # Return as a dataframe
  result <- data.frame(avg_d_null = avg_d, med_d_null = med_d)
  return(result)
}
# Function to calculate the average effect size for the original study based on true effect = 0.2 and 0.5
esCalc_true <- function(sublist) {
  p_values <- as.numeric(sublist[["p"]])
  es_values <- as.numeric(sublist[["d"]])
  # Filter rows based on the criteria
  filtered_sublist <- sublist[!(p_values < 0 & es_values < 0), ]
  # Calculate the average and median of the retained 'd' values
  avg_d <- mean(filtered_sublist[["d"]])
  med_d <- median(filtered_sublist[["d"]])
  # Return as a dataframe
  result <- data.frame(avg_d_true = avg_d, med_d_true = med_d)
  return(result)
}
# Calculate ADR, TPR, FPR, TNR, FNR for a MABF method
# Function to calculate ADR, TPR, FPR, TNR, FNR for a MABF method
ratesCalc <- function(bf_matrix, BFcutoff = 3, num_deltai = 500) {
  #Split the matrix into no effect and true effect
  BFs_no_effect <- bf_matrix[1:num_deltai, ] #These BFs are supposed to fall within the no-effect range but some may not
  BFs_true_effect <- bf_matrix[c(num_deltai+1):c(num_deltai*2), ] #These BFs are supposed to fall within the true-effect range but some may not
  
  # Calculate counts for true negative and false positive in no effect
  # You can change the 0.33 and 3 to 0.1 and 10 for more strict threshold
  TN <- sum(BFs_no_effect < round(1/BFcutoff, 2))
  FP <- sum(BFs_no_effect > BFcutoff)
  
  # Calculate counts for true positive and false negative in true effect
  TP <- sum(BFs_true_effect > BFcutoff)
  FN <- sum(BFs_true_effect < round(1/BFcutoff, 2))
  
  # Calculate rates
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  TNR <- TN / (TN + FP)
  FNR <- FN / (FN + TP)
  
  # Calculate proportions of indeterminate values
  anecdote_no_effect <- sum(BFs_no_effect >= round(1/BFcutoff, 2) & BFs_no_effect <= BFcutoff)
  anecdote_true_effect <- sum(BFs_true_effect >= round(1/BFcutoff, 2) & BFs_true_effect <= BFcutoff)
  anecdote_no_effect_prop <- sum(BFs_no_effect >= round(1/BFcutoff, 2) & BFs_no_effect <= BFcutoff) / (nrow(BFs_no_effect)*ncol(BFs_no_effect))
  anecdote_true_effect_prop <- sum(BFs_true_effect >= round(1/BFcutoff, 2) & BFs_true_effect <= BFcutoff) / (nrow(BFs_no_effect)*ncol(BFs_no_effect)) #BFs_no_effect and BF_true_effect have the same dimensions
  
  TFPNrates_df <- data.frame(
    TPR = TPR,
    FPR = FPR,
    TNR = TNR,
    FNR = FNR,
    ADRNE = anecdote_no_effect_prop,
    ADRTE = anecdote_true_effect_prop,
    ADR = (anecdote_no_effect + anecdote_true_effect)/(nrow(BFs_no_effect)*ncol(BFs_no_effect) + nrow(BFs_no_effect)*ncol(BFs_no_effect)),
    TP, FN, FP, TN,
    AD = anecdote_no_effect + anecdote_true_effect,
    ADNE = anecdote_no_effect,
    ADTE = anecdote_true_effect
  )
}
# Function to calculate TSR,FSR,TFR,FFR for a MABF method
ratesCalc2 <- function(bf_matrix,  pcutoff = 0.05, EScutoff = 0, BFcutoff = 3, num_deltai = 500) {
  # Split the matrix into no effect and true effect
  BFs_no_effect_nop <- bf_matrix[1:num_deltai, 4:c(rep_times+3)] #the BFs in no effect rows
  BFs_true_effect_nop <- bf_matrix[c(num_deltai+1):c(num_deltai*2), 4:c(rep_times+3)] #the BFs in true effect rows
  BFs_no_effect_p <- bf_matrix[1:num_deltai, 1:c(rep_times+3)] #the BFs in no effect rows
  BFs_true_effect_p <- bf_matrix[c(num_deltai+1):c(num_deltai*2), 1:c(rep_times+3)] #the BFs in true effect rows
  
  # Calculate counts for true success and false failure in true effect
  TS1 <- sum(BFs_true_effect_p[,2] < pcutoff & BFs_true_effect_p[,3] > EScutoff & BFs_true_effect_nop > BFcutoff) #Added ES > 0
  FF1 <- sum(BFs_true_effect_p[,2] > pcutoff & BFs_true_effect_nop > BFcutoff) 
  
  # Calculate counts for true failure and false success in true effect
  TF1 <- sum(BFs_true_effect_p[,2] < pcutoff & BFs_true_effect_p[,3] > EScutoff & BFs_true_effect_nop < round(1/BFcutoff, 2)) #Added ES > 0
  FS1 <- sum(BFs_true_effect_p[,2] > pcutoff & BFs_true_effect_nop < round(1/BFcutoff, 2)) 
  
  # Calculate counts for false success and true failure in no effect
  FS2 <- sum(BFs_no_effect_p[,2] < pcutoff & BFs_no_effect_p[,3] > EScutoff & BFs_no_effect_nop > BFcutoff) #Added ES > 0
  TF2 <- sum(BFs_no_effect_p[,2] > pcutoff & BFs_no_effect_nop > BFcutoff)
  
  # Calculate counts for false failure and true success in no effect
  FF2 <- sum(BFs_no_effect_p[,2] < pcutoff & BFs_no_effect_p[,3] > EScutoff & BFs_no_effect_nop < round(1/BFcutoff, 2)) #Added ES > 0
  TS2 <- sum(BFs_no_effect_p[,2] > pcutoff & BFs_no_effect_nop < round(1/BFcutoff, 2))
  
  # Calculate successful/failed affirmation/correction rates
  SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate; When original p-value is correct, what's the probability that MABF is also correct
  FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate; When original p-value is correct, what's the probability that MABF is not correct
  FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate; When original p-value is not correct, what's the probability that MABF is not correct, either
  SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate; When original p-value is not correct, what's the probability that MABF is correct
  SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
  FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
  FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
  SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
  SAR = (TS1 + TS2)/(TS1 + TS2 + TF1 + TF2)
  FAR = (TF1 + TF2)/(TS1 + TS2 + TF1 + TF2)
  FCR = (FS1 + FS2)/(FS1 + FS2 + FF1 + FF2)
  SCR = (FF1 + FF2)/(FS1 + FS2 + FF1 + FF2)
  
  # Calculate overall successes and failures
  Successes <- TS1 + TS2 + FS1 + FS2
  Failures <-  TF1 + TF2 + FF1 + FF2
  
  # Split successes and failures by true and false effect
  Successes1 <- TS1 + FS1
  Failures1 <-  TF1 + FF1
  Successes2 <- TS2 + FS2
  Failures2 <-  TF2 + FF2
  
  TSR = (TS1 + TS2)/Successes
  FSR = (FS1 + FS2)/Successes
  TFR = (TF1 + TF2)/Failures
  FFR = (FF1 + FF2)/Failures
  
  TSR1 = TS1/Successes1 #True success rate; Of all the claimed successful replications (i.e., original study and replications are consistent), what proportion of them are true successes (i.e., both original and replications reflects the true underlying effect)
  FSR1 = FS1/Successes1 #False success rate; Of all the claimed successful replications (i.e., original study and replications are consistent),  what proportion of them are false successes (i.e., neither original and replications reflects the true underlying effect)
  TFR1 = TF1/Failures1 #True failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are true failures (i.e., only the original study reflects the true underlying effect)
  FFR1 = FF1/Failures1 #False failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are false failures (i.e., only the replications reflect the true underlying effect)
  
  TSR2 = TS2/Successes2
  FSR2 = FS2/Successes2
  TFR2 = TF2/Failures2
  FFR2 = FF2/Failures2
  
  # Calculate rates
  TFSFrates_df <- data.frame(
    SAR, FAR, SCR, FCR,
    SAR1, FAR1, SCR1, FCR1,
    SAR2, FAR2, SCR2, FCR2,
    
    TSR,FSR,TFR,FFR,
    TSR1,FSR1,TFR1,FFR1,
    TSR2,FSR2,TFR2,FFR2,
    
    Successes,Failures,
    Successes1,Failures1,
    Successes2,Failures2,
    
    TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2
  )
}

########################################################################
# Function to calculate all kinds of row-wise rates (ADR, TPFPR, TSFSR)#
########################################################################
# Calculate ADR, TPFPR, TSFSR based on 500 null effects and 500 true effects, and repeat this process 500 times for each scenario
# For true effect = 0.2
ratesCalc_rowise <- function(dataset, pcutoff = 0.05, EScutoff = 0, BFcutoff = 3, num_deltai = 500) {
  # Create empty lists to store results
  results_null <- list()
  results_true <- list()
  
  # Get the names of sublists and matrices
  sublists_names <- names(dataset)
  
  # Iterate through each sublist and matrix
  for (i in 1:length(dataset)) {
    sublist_name <- sublists_names[i]
    message("Processing sublist ", sublist_name, " (", i, " of ", length(dataset), ")")
    
    matrices_names <- names(dataset[[i]])
    
    for (j in 1:length(dataset[[i]])) {
      matrix_name <- matrices_names[j]
     # message("  Processing matrix ", matrix_name, " (", j, " of ", length(dataset[[i]]), ")")
      
      current_matrix <- dataset[[i]][[j]]
      
      # Iterate through each row
      for (k in 1:nrow(current_matrix)) {
        if (k %% 100 == 0) {
          #message("    Processing row ", k, " of ", nrow(current_matrix))
        }
        
        if (k <= num_deltai) {
          # Null effect part
          p_null <- current_matrix[k, 2]
          es_null <- current_matrix[k, 3]
          BFs_null <- current_matrix[k, 4:(3 + num_deltai)]
          
          TN <- sum(BFs_null < round(1/BFcutoff, 2)) 
          FP <- sum(BFs_null > BFcutoff) 
          FS2 <- sum(p_null < pcutoff & es_null > EScutoff & BFs_null > BFcutoff) #Added es_null
          TF2 <- sum(p_null > pcutoff & BFs_null > BFcutoff)
          FF2 <- sum(p_null < pcutoff & es_null > EScutoff & BFs_null < round(1/BFcutoff, 2)) #Added es_null
          TS2 <- sum(p_null > pcutoff & BFs_null < round(1/BFcutoff, 2))
          AD_null <- sum(BFs_null >= round(1/BFcutoff, 2) & BFs_null <= BFcutoff) 
          AD_null_1 <- sum(BFs_null >= 1 & BFs_null <= 1)
          AD_null_anecdotal <- sum(BFs_null >= 0.33 & BFs_null <= 3) - AD_null_1
          AD_null_moderate <- sum(BFs_null >= 0.1 & BFs_null <= 10) - AD_null_anecdotal - AD_null_1
          AD_null_strong <- sum(BFs_null >= 0.03 & BFs_null <= 30) - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_vstrong <- sum(BFs_null >= 0.01 & BFs_null <= 100) - AD_null_strong - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_estrong <- num_deltai - AD_null_vstrong - AD_null_strong - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_anecdotal_N <- sum(BFs_null >= 0.33 & BFs_null < 1) #null direction
          AD_null_anecdotal_A <- sum(BFs_null > 1 & BFs_null <= 3) #alternative direction
          AD_null_moderate_N <- sum(BFs_null >= 0.1 & BFs_null < 0.33)
          AD_null_moderate_A <- sum(BFs_null > 3 & BFs_null <= 10)
          AD_null_strong_N <- sum(BFs_null >= 0.03 & BFs_null < 0.1)
          AD_null_strong_A <- sum(BFs_null > 10 & BFs_null <= 30)
          AD_null_vstrong_N <- sum(BFs_null >= 0.01 & BFs_null < 0.03)
          AD_null_vstrong_A <- sum(BFs_null > 30 & BFs_null <= 100)
          AD_null_estrong_N <- sum(BFs_null < 0.01)
          AD_null_estrong_A <- sum(BFs_null > 100)
          
          Successes2 <- TS2 + FS2
          Failures2 <- TF2 + FF2
          
          TNR = TN/(TN + FP)
          FPR = FP/(FP + TN)
          
          TSR2 = TS2/Successes2
          FSR2 = FS2/Successes2
          TFR2 = TF2/Failures2
          FFR2 = FF2/Failures2
          
          SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
          FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
          FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
          SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
          
          ADR_null <- AD_null/num_deltai
          ADR_null_1 <- AD_null_1/num_deltai
          ADR_null_anecdotal <- AD_null_anecdotal/num_deltai
          ADR_null_moderate <- AD_null_moderate/num_deltai
          ADR_null_strong <- AD_null_strong/num_deltai
          ADR_null_vstrong <- AD_null_vstrong/num_deltai
          ADR_null_estrong <- AD_null_estrong/num_deltai
          
          ADOdds_null_anecdotal <- AD_null_anecdotal_N/AD_null_anecdotal_A
          ADOdds_null_moderate <- AD_null_moderate_N/AD_null_moderate_A
          ADOdds_null_strong <- AD_null_strong_N/AD_null_strong_A
          ADOdds_null_vstrong <- AD_null_vstrong_N/AD_null_vstrong_A
          ADOdds_null_estrong <- AD_null_estrong_N/AD_null_estrong_A
          
          results_null <- append(results_null, list(data.frame(sublist_name = sublist_name, matrix_name = matrix_name, 
                                                               ADR_null, ADR_null_1, ADR_null_anecdotal, ADR_null_moderate, 
                                                               ADR_null_strong, ADR_null_vstrong, ADR_null_estrong,
                                                               ADOdds_null_anecdotal, ADOdds_null_moderate, 
                                                               ADOdds_null_strong, ADOdds_null_vstrong, ADOdds_null_estrong,
                                                               TNR = TNR, FPR = FPR,
                                                               SAR2, FAR2, FCR2, SCR2,
                                                               TSR2 = TSR2, FSR2 = FSR2, TFR2 = TFR2, FFR2 = FFR2,
                                                               AD_null = AD_null, AD_null_1 = AD_null_1, AD_null_anecdotal = AD_null_anecdotal,
                                                               AD_null_moderate = AD_null_moderate, AD_null_strong = AD_null_strong,
                                                               AD_null_vstrong = AD_null_vstrong, AD_null_estrong = AD_null_estrong,
                                                               AD_null_anecdotal_N = AD_null_anecdotal_N, AD_null_anecdotal_A = AD_null_anecdotal_A,
                                                               AD_null_moderate_N = AD_null_moderate_N, AD_null_moderate_A = AD_null_moderate_A,
                                                               AD_null_strong_N = AD_null_strong_N, AD_null_strong_A = AD_null_strong_A,
                                                               AD_null_vstrong_N = AD_null_vstrong_N, AD_null_vstrong_A = AD_null_vstrong_A,
                                                               AD_null_estrong_N = AD_null_estrong_N, AD_null_estrong_A = AD_null_estrong_A,
                                                               TN = TN, FP = FP, 
                                                               TS2 = TS2, FF2 = FF2, TF2 = TF2, FS2 = FS2)))
        } else {
          # True effect part
          p_true <- current_matrix[k, 2]
          es_true <- current_matrix[k, 3]
          BFs_true <- current_matrix[k, 4:(3 + num_deltai)]
          
          TP <- sum(BFs_true > BFcutoff) 
          FN <- sum(BFs_true < round(1/BFcutoff, 2)) 
          TS1 <- sum(p_true < pcutoff & es_true > EScutoff & BFs_true > BFcutoff) #Added es_true
          FF1 <- sum(p_true > pcutoff & BFs_true > BFcutoff)
          TF1 <- sum(p_true < pcutoff & es_true > EScutoff & BFs_true < round(1/BFcutoff, 2)) #Added es_true
          FS1 <- sum(p_true > pcutoff & BFs_true < round(1/BFcutoff, 2))
          AD_true <- sum(BFs_true >= round(1/BFcutoff, 2) & BFs_true <= BFcutoff)
          AD_true_1 <- sum(BFs_true >= 1 & BFs_true <= 1)
          AD_true_anecdotal <- sum(BFs_true >= 0.33 & BFs_true <= 3) - AD_true_1
          AD_true_moderate <- sum(BFs_true >= 0.1 & BFs_true <= 10) - AD_true_anecdotal - AD_true_1
          AD_true_strong <- sum(BFs_true >= 0.03 & BFs_true <= 30) - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_vstrong <- sum(BFs_true >= 0.01 & BFs_true <= 100) - AD_true_strong - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_estrong <- num_deltai - AD_true_vstrong - AD_true_strong - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_anecdotal_N <- sum(BFs_true >= 0.33 & BFs_true < 1) #true direction
          AD_true_anecdotal_A <- sum(BFs_true > 1 & BFs_true <= 3) #alternative direction
          AD_true_moderate_N <- sum(BFs_true >= 0.1 & BFs_true < 0.33)
          AD_true_moderate_A <- sum(BFs_true > 3 & BFs_true <= 10)
          AD_true_strong_N <- sum(BFs_true >= 0.03 & BFs_true < 0.1)
          AD_true_strong_A <- sum(BFs_true > 10 & BFs_true <= 30)
          AD_true_vstrong_N <- sum(BFs_true >= 0.01 & BFs_true < 0.03)
          AD_true_vstrong_A <- sum(BFs_true > 30 & BFs_true <= 100)
          AD_true_estrong_N <- sum(BFs_true < 0.01)
          AD_true_estrong_A <- sum(BFs_true > 100)
          
          Successes1 <- TS1 + FS1
          Failures1 <- TF1 + FF1
          
          TPR = TP/(TP + FN)
          FNR = FN/(FN + TP)
          
          TSR1 = TS1/Successes1
          FSR1 = FS1/Successes1
          TFR1 = TF1/Failures1
          FFR1 = FF1/Failures1
          
          SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate
          FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate
          FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate
          SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate
          
          ADR_true <- AD_true/num_deltai
          ADR_true_1 <- AD_true_1/num_deltai
          ADR_true_anecdotal <- AD_true_anecdotal/num_deltai
          ADR_true_moderate <- AD_true_moderate/num_deltai
          ADR_true_strong <- AD_true_strong/num_deltai
          ADR_true_vstrong <- AD_true_vstrong/num_deltai
          ADR_true_estrong <- AD_true_estrong/num_deltai
          
          ADOdds_true_anecdotal <- AD_true_anecdotal_A/AD_true_anecdotal_N
          ADOdds_true_moderate <- AD_true_moderate_A/AD_true_moderate_N
          ADOdds_true_strong <- AD_true_strong_A/AD_true_strong_N
          ADOdds_true_vstrong <- AD_true_vstrong_A/AD_true_vstrong_N
          ADOdds_true_estrong <- AD_true_estrong_A/AD_true_estrong_N
          
          results_true <- append(results_true, list(data.frame(sublist_name = sublist_name, matrix_name = matrix_name,  
                                                               ADR_true, ADR_true_1, ADR_true_anecdotal, ADR_true_moderate, 
                                                               ADR_true_strong, ADR_true_vstrong, ADR_true_estrong, 
                                                               ADOdds_true_anecdotal, ADOdds_true_moderate, 
                                                               ADOdds_true_strong, ADOdds_true_vstrong, ADOdds_true_estrong,
                                                               TPR = TPR, FNR = FNR,
                                                               SAR1, FAR1, FCR1, SCR1,
                                                               TSR1 = TSR1, FSR1 = FSR1, TFR1 = TFR1, FFR1 = FFR1,
                                                               AD_true = AD_true, AD_true_1 = AD_true_1, AD_true_anecdotal = AD_true_anecdotal,
                                                               AD_true_moderate = AD_true_moderate, AD_true_strong = AD_true_strong,
                                                               AD_true_vstrong = AD_true_vstrong, AD_true_estrong = AD_true_estrong,
                                                               AD_true_anecdotal_N = AD_true_anecdotal_N, AD_true_anecdotal_A = AD_true_anecdotal_A,
                                                               AD_true_moderate_N = AD_true_moderate_N, AD_true_moderate_A = AD_true_moderate_A,
                                                               AD_true_strong_N = AD_true_strong_N, AD_true_strong_A = AD_true_strong_A,
                                                               AD_true_vstrong_N = AD_true_vstrong_N, AD_true_vstrong_A = AD_true_vstrong_A,
                                                               AD_true_estrong_N = AD_true_estrong_N, AD_true_estrong_A = AD_true_estrong_A,
                                                               TP = TP, FN = FN, 
                                                               TS1 = TS1, FF1 = FF1, TF1 = TF1, FS1 = FS1)))
        }
      }
    }
  }
  
  # Convert lists to data frames
  results_null_df <- do.call(rbind, results_null)
  results_true_df <- do.call(rbind, results_true)
  
  # Concatenate sublist_name and matrix_name into one column for results_null_df
  results_null_df <- results_null_df %>%
    unite("sublist_matrix", sublist_name, matrix_name, sep = "_", remove = FALSE)
  
  # Concatenate sublist_name and matrix_name into one column for results_true_df
  results_true_df <- results_true_df %>%
    unite("sublist_matrix", sublist_name, matrix_name, sep = "_", remove = FALSE)
  
  # Break column names into factors
  results_null_df <- results_null_df %>% 
    separate(col = "sublist_name", into = c("true.effect", "orig.n", "QRP.level", "PB.level"), sep = "_") %>% 
    separate(col = "matrix_name", into = c("rep.number", "rep.n"), sep = "_") 
  
  results_true_df <- results_true_df %>% 
    separate(col = "sublist_name", into = c("true.effect", "orig.n", "QRP.level", "PB.level"), sep = "_") %>% 
    separate(col = "matrix_name", into = c("rep.number", "rep.n"), sep = "_") 
  
  # Change values in true.effect column
  results_null_df <- results_null_df %>%
    mutate(true.effect = replace(true.effect, true.effect == "0,0.2", 0))
  
  results_true_df <- results_true_df %>%
    mutate(true.effect = replace(true.effect, true.effect == "0,0.2", 0.2))
  
  return(list(results_null_df = results_null_df, results_true_df = results_true_df))
}

# For true effect = 0.5
ratesCalc_rowise2 <- function(dataset, pcutoff = 0.05, EScutoff = 0, BFcutoff = 3, num_deltai = 500) {
  # Create empty lists to store results
  results_null <- list()
  results_true <- list()
  
  # Get the names of sublists and matrices
  sublists_names <- names(dataset)
  
  # Iterate through each sublist and matrix
  for (i in 1:length(dataset)) {
    sublist_name <- sublists_names[i]
    message("Processing sublist ", sublist_name, " (", i, " of ", length(dataset), ")")
    
    matrices_names <- names(dataset[[i]])
    
    for (j in 1:length(dataset[[i]])) {
      matrix_name <- matrices_names[j]
      # message("  Processing matrix ", matrix_name, " (", j, " of ", length(dataset[[i]]), ")")
      
      current_matrix <- dataset[[i]][[j]]
      
      # Iterate through each row
      for (k in 1:nrow(current_matrix)) {
        if (k %% 100 == 0) {
          #message("    Processing row ", k, " of ", nrow(current_matrix))
        }
        
        if (k <= num_deltai) {
          # Null effect part
          p_null <- current_matrix[k, 2]
          es_null <- current_matrix[k, 3]
          BFs_null <- current_matrix[k, 4:(3 + num_deltai)]
          
          TN <- sum(BFs_null < round(1/BFcutoff, 2)) 
          FP <- sum(BFs_null > BFcutoff) 
          FS2 <- sum(p_null < pcutoff & es_null > EScutoff & BFs_null > BFcutoff) #Added es_null
          TF2 <- sum(p_null > pcutoff & BFs_null > BFcutoff)
          FF2 <- sum(p_null < pcutoff & es_null > EScutoff & BFs_null < round(1/BFcutoff, 2)) #Added es_null
          TS2 <- sum(p_null > pcutoff & BFs_null < round(1/BFcutoff, 2))
          AD_null <- sum(BFs_null >= round(1/BFcutoff, 2) & BFs_null <= BFcutoff) 
          AD_null_1 <- sum(BFs_null >= 1 & BFs_null <= 1)
          AD_null_anecdotal <- sum(BFs_null >= 0.33 & BFs_null <= 3) - AD_null_1
          AD_null_moderate <- sum(BFs_null >= 0.1 & BFs_null <= 10) - AD_null_anecdotal - AD_null_1
          AD_null_strong <- sum(BFs_null >= 0.03 & BFs_null <= 30) - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_vstrong <- sum(BFs_null >= 0.01 & BFs_null <= 100) - AD_null_strong - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_estrong <- num_deltai - AD_null_vstrong - AD_null_strong - AD_null_moderate - AD_null_anecdotal - AD_null_1
          AD_null_anecdotal_N <- sum(BFs_null >= 0.33 & BFs_null < 1) #null direction
          AD_null_anecdotal_A <- sum(BFs_null > 1 & BFs_null <= 3) #alternative direction
          AD_null_moderate_N <- sum(BFs_null >= 0.1 & BFs_null < 0.33)
          AD_null_moderate_A <- sum(BFs_null > 3 & BFs_null <= 10)
          AD_null_strong_N <- sum(BFs_null >= 0.03 & BFs_null < 0.1)
          AD_null_strong_A <- sum(BFs_null > 10 & BFs_null <= 30)
          AD_null_vstrong_N <- sum(BFs_null >= 0.01 & BFs_null < 0.03)
          AD_null_vstrong_A <- sum(BFs_null > 30 & BFs_null <= 100)
          AD_null_estrong_N <- sum(BFs_null < 0.01)
          AD_null_estrong_A <- sum(BFs_null > 100)
          
          Successes2 <- TS2 + FS2
          Failures2 <- TF2 + FF2
          
          TNR = TN/(TN + FP)
          FPR = FP/(FP + TN)
          
          TSR2 = TS2/Successes2
          FSR2 = FS2/Successes2
          TFR2 = TF2/Failures2
          FFR2 = FF2/Failures2
          
          SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
          FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
          FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
          SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
          
          ADR_null <- AD_null/num_deltai
          ADR_null_1 <- AD_null_1/num_deltai
          ADR_null_anecdotal <- AD_null_anecdotal/num_deltai
          ADR_null_moderate <- AD_null_moderate/num_deltai
          ADR_null_strong <- AD_null_strong/num_deltai
          ADR_null_vstrong <- AD_null_vstrong/num_deltai
          ADR_null_estrong <- AD_null_estrong/num_deltai
          
          ADOdds_null_anecdotal <- AD_null_anecdotal_N/AD_null_anecdotal_A
          ADOdds_null_moderate <- AD_null_moderate_N/AD_null_moderate_A
          ADOdds_null_strong <- AD_null_strong_N/AD_null_strong_A
          ADOdds_null_vstrong <- AD_null_vstrong_N/AD_null_vstrong_A
          ADOdds_null_estrong <- AD_null_estrong_N/AD_null_estrong_A
          
          results_null <- append(results_null, list(data.frame(sublist_name = sublist_name, matrix_name = matrix_name, 
                                                               ADR_null, ADR_null_1, ADR_null_anecdotal, ADR_null_moderate, 
                                                               ADR_null_strong, ADR_null_vstrong, ADR_null_estrong,
                                                               ADOdds_null_anecdotal, ADOdds_null_moderate, 
                                                               ADOdds_null_strong, ADOdds_null_vstrong, ADOdds_null_estrong,
                                                               TNR = TNR, FPR = FPR,
                                                               SAR2, FAR2, FCR2, SCR2,
                                                               TSR2 = TSR2, FSR2 = FSR2, TFR2 = TFR2, FFR2 = FFR2,
                                                               AD_null = AD_null, AD_null_1 = AD_null_1, AD_null_anecdotal = AD_null_anecdotal,
                                                               AD_null_moderate = AD_null_moderate, AD_null_strong = AD_null_strong,
                                                               AD_null_vstrong = AD_null_vstrong, AD_null_estrong = AD_null_estrong,
                                                               AD_null_anecdotal_N = AD_null_anecdotal_N, AD_null_anecdotal_A = AD_null_anecdotal_A,
                                                               AD_null_moderate_N = AD_null_moderate_N, AD_null_moderate_A = AD_null_moderate_A,
                                                               AD_null_strong_N = AD_null_strong_N, AD_null_strong_A = AD_null_strong_A,
                                                               AD_null_vstrong_N = AD_null_vstrong_N, AD_null_vstrong_A = AD_null_vstrong_A,
                                                               AD_null_estrong_N = AD_null_estrong_N, AD_null_estrong_A = AD_null_estrong_A,
                                                               TN = TN, FP = FP, 
                                                               TS2 = TS2, FF2 = FF2, TF2 = TF2, FS2 = FS2)))
        } else {
          # True effect part
          p_true <- current_matrix[k, 2]
          es_true <- current_matrix[k, 3]
          BFs_true <- current_matrix[k, 4:(3 + num_deltai)]
          
          TP <- sum(BFs_true > BFcutoff) 
          FN <- sum(BFs_true < round(1/BFcutoff, 2)) 
          TS1 <- sum(p_true < pcutoff & es_true > EScutoff & BFs_true > BFcutoff) #Added es_true
          FF1 <- sum(p_true > pcutoff & BFs_true > BFcutoff)
          TF1 <- sum(p_true < pcutoff & es_true > EScutoff & BFs_true < round(1/BFcutoff, 2)) #Added es_true
          FS1 <- sum(p_true > pcutoff & BFs_true < round(1/BFcutoff, 2))
          AD_true <- sum(BFs_true >= round(1/BFcutoff, 2) & BFs_true <= BFcutoff)
          AD_true_1 <- sum(BFs_true >= 1 & BFs_true <= 1)
          AD_true_anecdotal <- sum(BFs_true >= 0.33 & BFs_true <= 3) - AD_true_1
          AD_true_moderate <- sum(BFs_true >= 0.1 & BFs_true <= 10) - AD_true_anecdotal - AD_true_1
          AD_true_strong <- sum(BFs_true >= 0.03 & BFs_true <= 30) - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_vstrong <- sum(BFs_true >= 0.01 & BFs_true <= 100) - AD_true_strong - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_estrong <- num_deltai - AD_true_vstrong - AD_true_strong - AD_true_moderate - AD_true_anecdotal - AD_true_1
          AD_true_anecdotal_N <- sum(BFs_true >= 0.33 & BFs_true < 1) #true direction
          AD_true_anecdotal_A <- sum(BFs_true > 1 & BFs_true <= 3) #alternative direction
          AD_true_moderate_N <- sum(BFs_true >= 0.1 & BFs_true < 0.33)
          AD_true_moderate_A <- sum(BFs_true > 3 & BFs_true <= 10)
          AD_true_strong_N <- sum(BFs_true >= 0.03 & BFs_true < 0.1)
          AD_true_strong_A <- sum(BFs_true > 10 & BFs_true <= 30)
          AD_true_vstrong_N <- sum(BFs_true >= 0.01 & BFs_true < 0.03)
          AD_true_vstrong_A <- sum(BFs_true > 30 & BFs_true <= 100)
          AD_true_estrong_N <- sum(BFs_true < 0.01)
          AD_true_estrong_A <- sum(BFs_true > 100)
          
          Successes1 <- TS1 + FS1
          Failures1 <- TF1 + FF1
          
          TPR = TP/(TP + FN)
          FNR = FN/(FN + TP)
          
          TSR1 = TS1/Successes1
          FSR1 = FS1/Successes1
          TFR1 = TF1/Failures1
          FFR1 = FF1/Failures1
          
          SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate
          FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate
          FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate
          SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate
          
          ADR_true <- AD_true/num_deltai
          ADR_true_1 <- AD_true_1/num_deltai
          ADR_true_anecdotal <- AD_true_anecdotal/num_deltai
          ADR_true_moderate <- AD_true_moderate/num_deltai
          ADR_true_strong <- AD_true_strong/num_deltai
          ADR_true_vstrong <- AD_true_vstrong/num_deltai
          ADR_true_estrong <- AD_true_estrong/num_deltai
          
          ADOdds_true_anecdotal <- AD_true_anecdotal_A/AD_true_anecdotal_N
          ADOdds_true_moderate <- AD_true_moderate_A/AD_true_moderate_N
          ADOdds_true_strong <- AD_true_strong_A/AD_true_strong_N
          ADOdds_true_vstrong <- AD_true_vstrong_A/AD_true_vstrong_N
          ADOdds_true_estrong <- AD_true_estrong_A/AD_true_estrong_N
          
          results_true <- append(results_true, list(data.frame(sublist_name = sublist_name, matrix_name = matrix_name,  
                                                               ADR_true, ADR_true_1, ADR_true_anecdotal, ADR_true_moderate, 
                                                               ADR_true_strong, ADR_true_vstrong, ADR_true_estrong, 
                                                               ADOdds_true_anecdotal, ADOdds_true_moderate, 
                                                               ADOdds_true_strong, ADOdds_true_vstrong, ADOdds_true_estrong,
                                                               TPR = TPR, FNR = FNR,
                                                               SAR1, FAR1, FCR1, SCR1,
                                                               TSR1 = TSR1, FSR1 = FSR1, TFR1 = TFR1, FFR1 = FFR1,
                                                               AD_true = AD_true, AD_true_1 = AD_true_1, AD_true_anecdotal = AD_true_anecdotal,
                                                               AD_true_moderate = AD_true_moderate, AD_true_strong = AD_true_strong,
                                                               AD_true_vstrong = AD_true_vstrong, AD_true_estrong = AD_true_estrong,
                                                               AD_true_anecdotal_N = AD_true_anecdotal_N, AD_true_anecdotal_A = AD_true_anecdotal_A,
                                                               AD_true_moderate_N = AD_true_moderate_N, AD_true_moderate_A = AD_true_moderate_A,
                                                               AD_true_strong_N = AD_true_strong_N, AD_true_strong_A = AD_true_strong_A,
                                                               AD_true_vstrong_N = AD_true_vstrong_N, AD_true_vstrong_A = AD_true_vstrong_A,
                                                               AD_true_estrong_N = AD_true_estrong_N, AD_true_estrong_A = AD_true_estrong_A,
                                                               TP = TP, FN = FN, 
                                                               TS1 = TS1, FF1 = FF1, TF1 = TF1, FS1 = FS1)))
        }
      }
    }
  }
  
  # Convert lists to data frames
  results_null_df <- do.call(rbind, results_null)
  results_true_df <- do.call(rbind, results_true)
  
  # Concatenate sublist_name and matrix_name into one column for results_null_df
  results_null_df <- results_null_df %>%
    unite("sublist_matrix", sublist_name, matrix_name, sep = "_", remove = FALSE)
  
  # Concatenate sublist_name and matrix_name into one column for results_true_df
  results_true_df <- results_true_df %>%
    unite("sublist_matrix", sublist_name, matrix_name, sep = "_", remove = FALSE)
  
  # Break column names into factors
  results_null_df <- results_null_df %>% 
    separate(col = "sublist_name", into = c("true.effect", "orig.n", "QRP.level", "PB.level"), sep = "_") %>% 
    separate(col = "matrix_name", into = c("rep.number", "rep.n"), sep = "_") 
  
  results_true_df <- results_true_df %>% 
    separate(col = "sublist_name", into = c("true.effect", "orig.n", "QRP.level", "PB.level"), sep = "_") %>% 
    separate(col = "matrix_name", into = c("rep.number", "rep.n"), sep = "_") 
  
  # Change values in true.effect column
  results_null_df <- results_null_df %>%
    mutate(true.effect = replace(true.effect, true.effect == "0,0.5", 0))
  
  results_true_df <- results_true_df %>%
    mutate(true.effect = replace(true.effect, true.effect == "0,0.5", 0.5))
  
  return(list(results_null_df = results_null_df, results_true_df = results_true_df))
}




#########################################################
#Function to calculate rates for fixed-effects MA method#
#########################################################

# Calculate TPR, FPR, TNR, FNR for Fixed-Effects MA
ratesCalcMA <- function(MA_matrix, pcutoff_r = 0.05, EScutoff_r = 0, num_deltai = 500) {
  #Split the matrix into no effect and true effect
  beta_no_effect <- MA_matrix[1:num_deltai, 1:num_deltai] #observed beta of null effect
  beta_true_effect <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 1:num_deltai] #observed beta of true effects
  p_no_effect <- MA_matrix[1:num_deltai, c(num_deltai+1):c(num_deltai*2)] #p value of null effect
  p_true_effect <- MA_matrix[c(num_deltai+1):c(num_deltai*2), c(num_deltai+1):c(num_deltai*2)] #p value of true effect
  
  # Calculate counts for true negative and false positive in no effect
  TN <- sum(p_no_effect > pcutoff_r)
  FP <- sum(beta_no_effect > EScutoff_r & p_no_effect < pcutoff_r)
  
  # Calculate counts for true positive and false negative in true effect
  TP <- sum(beta_true_effect > EScutoff_r & p_true_effect < pcutoff_r)
  FN <- sum(p_true_effect > pcutoff_r)
  
  # Calculate rates
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  TNR <- TN / (TN + FP)
  FNR <- FN / (FN + TP)
  
  TFPNrates_df <- data.frame(
    TPR = TPR,
    FPR = FPR,
    TNR = TNR,
    FNR = FNR,
    TP, FN, FP, TN
  )
}
# Function to calculate TSR,FSR,TFR,FFR for Fixed-Effects MA
ratesCalcMA2 <- function(MA_matrix,  pcutoff_o = 0.05, EScutoff_o = 0, pcutoff_r = 0.05, EScutoff_r = 0, num_deltai = 500) {
  # Split the matrix into no effect and true effect
  beta_no_effect_nop <- MA_matrix[1:num_deltai, 4:c(rep_times+3)] #observed beta of null effect, not original ES or p
  beta_true_effect_nop <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 4:c(rep_times+3)] #observed beta of true effects,not original ES or p
  p_no_effect_nop <- MA_matrix[1:num_deltai, c(rep_times+4):c(rep_times*2+3)] #p value of null effect, not original ES or p
  p_true_effect_nop <- MA_matrix[c(num_deltai+1):c(num_deltai*2), c(rep_times+4):c(rep_times*2+3)] #p value of true effect, not original ES or p
  
  BFs_no_effect_p <- MA_matrix[1:num_deltai, 1:c(rep_times*2+3)] # contains original p values and observed ES of null effect
  BFs_true_effect_p <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 1:c(rep_times*2+3)] # contains p values and observed ES of true effect
  
  # Calculate counts for true success and false failure in true effect
  TS1 <- sum(BFs_true_effect_p[,2] < pcutoff_o & BFs_true_effect_p[,3] > EScutoff_o & beta_true_effect_nop > EScutoff_r & p_true_effect_nop < pcutoff_r) 
  FF1 <- sum(BFs_true_effect_p[,2] > pcutoff_o & beta_true_effect_nop > EScutoff_r & p_true_effect_nop < pcutoff_r) 
  
  # Calculate counts for true failure and false success in true effect
  TF1 <- sum(BFs_true_effect_p[,2] < pcutoff_o & BFs_true_effect_p[,3] > EScutoff_o & p_true_effect_nop > pcutoff_r) 
  FS1 <- sum(BFs_true_effect_p[,2] > pcutoff_o & p_true_effect_nop > pcutoff_r) 
  
  # Calculate counts for false success and true failure in no effect
  FS2 <- sum(BFs_no_effect_p[,2] < pcutoff_o & BFs_no_effect_p[,3] > EScutoff_o & beta_no_effect_nop > EScutoff_r & p_no_effect_nop < pcutoff_r) 
  TF2 <- sum(BFs_no_effect_p[,2] > pcutoff_o & beta_no_effect_nop > EScutoff_r & p_no_effect_nop < pcutoff_r)
  
  # Calculate counts for false failure and true success in no effect
  FF2 <- sum(BFs_no_effect_p[,2] < pcutoff_o & BFs_no_effect_p[,3] > EScutoff_o & p_no_effect_nop > pcutoff_r) 
  TS2 <- sum(BFs_no_effect_p[,2] > pcutoff_o & p_no_effect_nop > pcutoff_r)
  
  # Calculate successful/failed affirmation/correction rates
  SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate; When original p-value is correct, what's the probability that MABF is also correct
  FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate;When original p-value is correct, what's the probability that MABF is not correct
  FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate; When original p-value is not correct, what's the probability that MABF is not correct, either
  SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate; When original p-value is not correct, what's the probability that MABF is correct
  SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
  FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
  FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
  SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
  SAR = (TS1 + TS2)/(TS1 + TS2 + TF1 + TF2)
  FAR = (TF1 + TF2)/(TS1 + TS2 + TF1 + TF2)
  FCR = (FS1 + FS2)/(FS1 + FS2 + FF1 + FF2)
  SCR = (FF1 + FF2)/(FS1 + FS2 + FF1 + FF2)
  
  # Calculate overall successes and failures
  Successes <- TS1 + TS2 + FS1 + FS2
  Failures <-  TF1 + TF2 + FF1 + FF2
  
  # Split successes and failures by true and false effect
  Successes1 <- TS1 + FS1
  Failures1 <-  TF1 + FF1
  Successes2 <- TS2 + FS2
  Failures2 <-  TF2 + FF2
  
  TSR = (TS1 + TS2)/Successes #True success rate; Of all the claimed successful replications (i.e., original study and replications are consistent), what proportion of them are true successes (i.e., both original and replications reflects the true underlying effect)
  FSR = (FS1 + FS2)/Successes #False success rate; Of all the claimed successful replications (i.e., original study and replications are consistent),  what proportion of them are false successes (i.e., neither original and replications reflects the true underlying effect)
  TFR = (TF1 + TF2)/Failures #True failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are true failures (i.e., only the original study reflects the true underlying effect)
  FFR = (FF1 + FF2)/Failures #False failure rate; Of all the claimed failed replications (i.e., original and replications are not consistent), what proportion of them are false failures (i.e., only the replications reflect the true underlying effect)
  
  TSR1 = TS1/Successes1
  FSR1 = FS1/Successes1
  TFR1 = TF1/Failures1
  FFR1 = FF1/Failures1
  
  TSR2 = TS2/Successes2
  FSR2 = FS2/Successes2
  TFR2 = TF2/Failures2
  FFR2 = FF2/Failures2
  
  # Calculate rates
  TFSFrates_df <- data.frame(
    SAR, FAR, SCR, FCR,
    SAR1, FAR1, SCR1, FCR1,
    SAR2, FAR2, SCR2, FCR2,
    
    TSR,FSR,TFR,FFR,
    TSR1,FSR1,TFR1,FFR1,
    TSR2,FSR2,TFR2,FFR2,
    
    Successes,Failures,
    Successes1,Failures1,
    Successes2,Failures2,
    
    TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2
  )
}

##########################################################################
#Function to calculate rates for fixed-effects MA method for creating ROC#
##########################################################################
ratesCalcMA_ROC <- function(MA_matrix, pcutoff_r, EScutoff_r = 0) {
  #Split the matrix into no effect and true effect
  beta_no_effect <- MA_matrix[1:num_deltai, 1:num_deltai] #observed beta of null effect
  beta_true_effect <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 1:num_deltai] #observed beta of true effects
  p_no_effect <- MA_matrix[1:num_deltai, c(num_deltai+1):c(num_deltai*2)] #p value of null effect
  p_true_effect <- MA_matrix[c(num_deltai+1):c(num_deltai*2), c(num_deltai+1):c(num_deltai*2)] #p value of true effect
  
  # Calculate counts for true negative and false positive in no effect
  TN <- sum(p_no_effect > pcutoff_r)
  FP <- sum(beta_no_effect > EScutoff_r & p_no_effect < pcutoff_r)
  
  # Calculate counts for true positive and false negative in true effect
  TP <- sum(beta_true_effect > EScutoff_r & p_true_effect < pcutoff_r)
  FN <- sum(p_true_effect > pcutoff_r)
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  TNR <- TN / (TN + FP)
  FNR <- FN / (FN + TP)
  data.frame(TPR = TPR, FPR = FPR, TNR = TNR, FNR = FNR, TP = TP, FN = FN, FP = FP, TN = TN)
}

# Function to process each matrix
# nested in process_sublist
processMA_matrix <- function(cutoff, MA_matrix, sublist_name, matrix_name) {
  rates <- ratesCalcMA_ROC(MA_matrix, cutoff, EScutoff_r =0)
  rates$Sublist <- sublist_name
  rates$Matrix <- matrix_name
  rates$threshold <- cutoff
  return(rates)
}

# Function to process each sublist in XXBF_lists_0.xnull_regrouped
# e.g. EUBF_lists_0.2null_regrouped
processMA_sublist <- function(sublist, sublist_name) {
  map_dfr(names(sublist), function(matrix_name) {
    message("  Processing matrix: ", matrix_name)
    MA_matrix <- sublist[[matrix_name]]
    map_dfr(pcutoff_r, function(cutoff) {
      processMA_matrix(cutoff, MA_matrix, sublist_name, matrix_name)
    })
  })
}

# ratesCalcMA_ROC_AiO function, which calculate all kinds of rates except ADR
# AiO is All in One
# Can specify cutoffs
# For fixed-effect MA
ratesCalcMA_ROC_AiO <- function(MA_matrix, pcutoff_r, EScutoff_r, pcutoff_o, EScutoff_o) {
  #Split the matrix into no effect and true effect
  beta_no_effect_nop <- MA_matrix[1:num_deltai, 4:c(rep_times+3)] #observed beta of null effect, not original ES or p
  beta_true_effect_nop <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 4:c(rep_times+3)] #observed beta of true effects,not original ES or p
  p_no_effect_nop <- MA_matrix[1:num_deltai, c(rep_times+4):c(rep_times*2+3)] #p value of null effect, not original ES or p
  p_true_effect_nop <- MA_matrix[c(num_deltai+1):c(num_deltai*2), c(rep_times+4):c(rep_times*2+3)] #p value of true effect, not original ES or p
  
  ps_no_effect_p <- MA_matrix[1:num_deltai, 1:c(rep_times*2+3)] # contains original p values and observed ES of null effect
  ps_true_effect_p <- MA_matrix[c(num_deltai+1):c(num_deltai*2), 1:c(rep_times*2+3)] # contains p values and observed ES of true effect
  
  # Calculate counts for true negative and false positive in no effect
  TN <- sum(p_no_effect_nop > pcutoff_r)
  FP <- sum(beta_no_effect_nop > EScutoff_r & p_no_effect_nop < pcutoff_r)
  
  # Calculate counts for true positive and false negative in true effect
  TP <- sum(beta_true_effect_nop > EScutoff_r & p_true_effect_nop < pcutoff_r)
  FN <- sum(p_true_effect_nop > pcutoff_r)
  
  # TPR <- TP / (TP + FN)
  # FPR <- FP / (FP + TN)
  # TNR <- TN / (TN + FP)
  # FNR <- FN / (FN + TP)
  
  # Calculate counts for true success and false failure in true effect
  # Calculate counts for true success and false failure in true effect
  TS1 <- sum(ps_true_effect_p[,2] < pcutoff_o & ps_true_effect_p[,3] > EScutoff_o & beta_true_effect_nop > EScutoff_r & p_true_effect_nop < pcutoff_r) 
  FF1 <- sum(ps_true_effect_p[,2] > pcutoff_o & beta_true_effect_nop > EScutoff_r & p_true_effect_nop < pcutoff_r) 
  
  # Calculate counts for true failure and false success in true effect
  TF1 <- sum(ps_true_effect_p[,2] < pcutoff_o & ps_true_effect_p[,3] > EScutoff_o & p_true_effect_nop > pcutoff_r) 
  FS1 <- sum(ps_true_effect_p[,2] > pcutoff_o & p_true_effect_nop > pcutoff_r) 
  
  # Calculate counts for false success and true failure in no effect
  FS2 <- sum(ps_no_effect_p[,2] < pcutoff_o & ps_no_effect_p[,3] > EScutoff_o & beta_no_effect_nop > EScutoff_r & p_no_effect_nop < pcutoff_r) 
  TF2 <- sum(ps_no_effect_p[,2] > pcutoff_o & beta_no_effect_nop > EScutoff_r & p_no_effect_nop < pcutoff_r)
  
  # Calculate counts for false failure and true success in no effect
  FF2 <- sum(ps_no_effect_p[,2] < pcutoff_o & ps_no_effect_p[,3] > EScutoff_o & p_no_effect_nop > pcutoff_r) 
  TS2 <- sum(ps_no_effect_p[,2] > pcutoff_o & p_no_effect_nop > pcutoff_r)
  
  # Calculate successful/failed affirmation/correction rates
  # SAR1 = TS1/(TS1 + TF1) #Successful Affirmation Rate
  # FAR1 = TF1/(TS1 + TF1) #Failed Affirmation Rate
  # FCR1 = FS1/(FS1 + FF1) #Failed Correction Rate
  # SCR1 = FF1/(FS1 + FF1) #Successful Correction Rate
  # SAR2 = TS2/(TS2 + TF2) #Successful Affirmation Rate
  # FAR2 = TF2/(TS2 + TF2) #Failed Affirmation Rate
  # FCR2 = FS2/(FS2 + FF2) #Failed Correction Rate
  # SCR2 = FF2/(FS2 + FF2) #Successful Correction Rate
  # SAR = (TS1 + TS2)/(TS1 + TS2 + TF1 + TF2)
  # FAR = (TF1 + TF2)/(TS1 + TS2 + TF1 + TF2)
  # FCR = (FS1 + FS2)/(FS1 + FS2 + FF1 + FF2)
  # SCR = (FF1 + FF2)/(FS1 + FS2 + FF1 + FF2)
  
  # Calculate overall successes and failures
  Successes <- TS1 + TS2 + FS1 + FS2
  Failures <-  TF1 + TF2 + FF1 + FF2
  
  # Split successes and failures by true and false effect
  Successes1 <- TS1 + FS1
  Failures1 <-  TF1 + FF1
  Successes2 <- TS2 + FS2
  Failures2 <-  TF2 + FF2
  
  # Calcualte True/False Success/Failure rates
  # TSR = (TS1 + TS2)/Successes
  # FSR = (FS1 + FS2)/Successes
  # TFR = (TF1 + TF2)/Failures
  # FFR = (FF1 + FF2)/Failures
  # 
  # TSR1 = TS1/Successes1
  # FSR1 = FS1/Successes1
  # TFR1 = TF1/Failures1
  # FFR1 = FF1/Failures1
  # 
  # TSR2 = TS2/Successes2
  # FSR2 = FS2/Successes2
  # TFR2 = TF2/Failures2
  # FFR2 = FF2/Failures2
  
  # TFSF = true/false/success/failure
  TFSFrates_df <- data.frame(
    # TPR,FPR,TNR,FNR,
    # 
    # SAR, FAR, SCR, FCR,
    # SAR1, FAR1, SCR1, FCR1,
    # SAR2, FAR2, SCR2, FCR2,
    # 
    # TSR,FSR,TFR,FFR,
    # TSR1,FSR1,TFR1,FFR1,
    # TSR2,FSR2,TFR2,FFR2,
    
    TP,FN,FP,TN,
    Successes,Failures,Successes1,Failures1,Successes2,Failures2,
    TS1, TS2, FF1, FF2, TF1, TF2, FS1, FS2
  )
}
# Function to process each matrix
# nested in process_sublist
processMA_matrix_AiO <- function(pcutoff_r, EScutoff_r, pcutoff_o, EScutoff_o, MA_matrix, sublist_name, matrix_name) {
  rates <- ratesCalcMA_ROC_AiO(MA_matrix, pcutoff_r, EScutoff_r, pcutoff_o, EScutoff_o)
  rates$threshold_p_r <- pcutoff_r
  rates$threshold_p <- pcutoff_o
  rates$threshold_ES <- EScutoff_o
  rates$Sublist <- sublist_name
  rates$Matrix <- matrix_name
  
  return(rates)
}

# Function to process each sublist in XXBF_lists_0.xnull_regrouped
# e.g. EUBF_lists_0.2null_regrouped
processMA_sublist_AiO <- function(sublist, sublist_name) {
  map_dfr(names(sublist), function(matrix_name) {
    message("  Processing matrix: ", matrix_name)
    MA_matrix <- sublist[[matrix_name]]
    map_dfr(cutoffs, function(cutoff) {
      processMA_matrix_AiO(cutoff$pcutoff_r, EScutoff_r=0, cutoff$pcutoff_o, cutoff$EScutoff_o, MA_matrix, sublist_name, matrix_name)
    })
  })
}




# Plot functions for rates
# Function to create box plot
extract_cutoff_value <- function(folder_path) {
  match <- str_match(folder_path, "pcutoff_o=([0-9.]+)")
  return(match[2])
}
plot_box <- function(data, trueES, ratetype, ratetype_char, use_third_var = FALSE, third_var = NULL, thirdvar_labeller = NULL, apply_filter = FALSE, filter_var = NULL, filter_level = NULL) {
  # Extract cutoff value from the base_folder_path
  cutoff_value <- extract_cutoff_value(base_folder_path)
  # Convert third_var and filter_var to symbols if they are provided
  if (!is.null(third_var)) {
    third_var_sym <- sym(third_var)
  }
  if (!is.null(filter_var)) {
    filter_var_sym <- sym(filter_var)
  }
  # Create a labeller function for the third variable if provided
  if (!is.null(thirdvar_labeller)) {
    thirdvar_labeller_func <- as_labeller(thirdvar_labeller, label_parsed)
  }
  # Dynamically construct the facet formula
  if (use_third_var && !is.null(third_var)) {
    facet_formula <- as.formula(paste(third_var, "~ rep.number.level"))
  } else {
    facet_formula <- as.formula(". ~ rep.number.level")
  }
  # Create the plot title
  title_text <- bquote("True effect = " * .(trueES) * ", " * italic(p[orig]) == .(cutoff_value))
  # Apply filter if required
  if (apply_filter && !is.null(filter_var) && !is.null(filter_level)) {
    data <- data %>%
      filter(!!filter_var_sym == filter_level)
  }
  # Create the plot
  plot <- data %>%
    {
      if (use_third_var && !is.null(third_var)) {
        group_by(., Method, rep.number.level, !!third_var_sym)
      } else {
        group_by(., Method, rep.number.level)
      }
    } %>%
    ggplot(aes(x = rep.n.level, y = !!sym(ratetype), color = Method)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3, fatten = 0.3) +  # Exclude outliers
    {
      if (use_third_var && !is.null(third_var)) {
        facet_grid(facet_formula, labeller = labeller(rep.number.level = rep.number_labeller, !!sym(third_var) := thirdvar_labeller_func))
      } else {
        facet_grid(facet_formula, labeller = labeller(rep.number.level = rep.number_labeller))
      }
    } +
    scale_y_continuous(labels = percent, breaks = seq(0, 1, by = 0.1)) +
    labs(x = 'Replication Sample Size', y = ratetype_char, title = title_text) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text = element_text(face = "bold")
    )
  return(plot)
}

# Using lines to connect the median of a method across different conditions
# Define the function to extract cutoff value from the folder path
plot_line <- function(data, trueES, ratetype, ratetype_char, use_third_var = FALSE, third_var = NULL, thirdvar_labeller = NULL, apply_filter = FALSE, filter_var = NULL, filter_level = NULL) {
  # Extract cutoff value from the base_folder_path
  cutoff_value <- extract_cutoff_value(base_folder_path)
  # Convert third_var and filter_var to symbols if they are provided
  if (!is.null(third_var)) {
    third_var_sym <- sym(third_var)
  }
  if (!is.null(filter_var)) {
    filter_var_sym <- sym(filter_var)
  }
  # Create a labeller function for the third variable if provided
  if (!is.null(thirdvar_labeller)) {
    thirdvar_labeller_func <- as_labeller(thirdvar_labeller, label_parsed)
  }
  # Dynamically construct the facet formula
  if (use_third_var && !is.null(third_var)) {
    facet_formula <- as.formula(paste(third_var, "~ rep.number.level"))
  } else {
    facet_formula <- as.formula(". ~ rep.number.level")
  }
  # Create the plot title
  title_text <- bquote("True effect = " * .(trueES) * ", " * italic(p[orig]) == .(cutoff_value))
  # Apply filter if required
  if (apply_filter && !is.null(filter_var) && !is.null(filter_level)) {
    data <- data %>%
      filter(!!filter_var_sym == filter_level)
  }
  # Calculate the median for each group
  data_summary <- data %>%
    {
      if (use_third_var && !is.null(third_var)) {
        group_by(., Method, rep.number.level, rep.n.level, !!third_var_sym)
      } else {
        group_by(., Method, rep.number.level, rep.n.level)
      }
    } %>%
    summarise(median_rate = median(!!sym(ratetype), na.rm = TRUE), .groups = 'drop')
  # Create the plot
  plot <- data_summary %>%
    ggplot(aes(x = rep.n.level, y = median_rate, color = Method, group = Method)) +
    geom_line() +
    geom_point() +
    {
      if (use_third_var && !is.null(third_var)) {
        facet_grid(facet_formula, labeller = labeller(rep.number.level = rep.number_labeller, !!sym(third_var) := thirdvar_labeller_func))
      } else {
        facet_grid(facet_formula, labeller = labeller(rep.number.level = rep.number_labeller))
      }
    } +
    scale_y_continuous(labels = percent, breaks = seq(0, 1, by = 0.1)) +
    labs(x = 'Replication Sample Size', y = ratetype_char, title = title_text) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(colour = "black", fill = "white"),
      strip.text = element_text(face = "bold")
    )
  return(plot)
}

#####
# Function to calculate precision and recall
# APR = Accuracy, Precision, Recall
ratesCalc_APR <- function(bf_matrix, BFcutoff) {
  BFs_no_effect <- bf_matrix[1:num_deltai, ]
  BFs_true_effect <- bf_matrix[(num_deltai + 1):(num_deltai * 2), ]
  TN <- sum(BFs_no_effect < BFcutoff)
  FP <- sum(BFs_no_effect > BFcutoff)
  TP <- sum(BFs_true_effect > BFcutoff)
  FN <- sum(BFs_true_effect < BFcutoff)
  accuracy = (TP + TN)/(TP + TN + FP + FN)
  precision = TP/(TP + FP)
  recall = TP/(TP + FN)
  data.frame(Accuracy = accuracy, Precision = precision, Recall = recall, TP = TP, FN = FN, FP = FP, TN = TN)
}

# Function to process each matrix for precision and recall calculation
# nested in process_sublist_pr
process_matrix_apr <- function(cutoff, bf_matrix, sublist_name, matrix_name) {
  rates <- ratesCalc_APR(bf_matrix, cutoff)
  rates$Sublist <- sublist_name
  rates$Matrix <- matrix_name
  rates$threshold <- cutoff
  return(rates)
}

# Function to process each sublist in XXBF_lists_0.xnull_regrouped
# e.g. EUBF_lists_0.2null_regrouped
process_sublist_apr <- function(sublist, sublist_name) {
  map_dfr(names(sublist), function(matrix_name) {
    message("  Processing matrix: ", matrix_name)
    bf_matrix <- sublist[[matrix_name]]
    map_dfr(BFcutoff, function(cutoff) {
      process_matrix_apr(cutoff, bf_matrix, sublist_name, matrix_name)
    })
  })
}
