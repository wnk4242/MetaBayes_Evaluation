# This script is used to verify results of count numbers to see if FS2P + FS2N + TF2 + FF2P + FF2N + TS2 + ADNE = 250000 (i.e., 500 original study x500 replications)

# We need the following datasets:
# rates_EUBF_0.2null_c1 in MABFanalyses/matrix-wise/rates4Plot/fixed original cutoff/pcutoff_o=0.05, EScutoff_o=0/c1 
x <- EUBF_lists_0.2null_regrouped_deltap[["0,0.2_200_high_high"]][["10_400"]] #in hprc/MABF4ROC/4TSFS
ADNE <- 5611
# Fasle Success 2
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] < .05 & x[i, 3] > 0 & x[i, 4:503] > 3, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))
FS2P <- total_cases


# False Success 2 (negative direction)
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] <.05 & x[i, 3] < 0 & x[i, 4:503] > 3, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    #print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))

FS2N <- total_cases

# True Failure 2
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] >.05 & x[i, 4:503] > 3, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    #print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))

TF2 <- total_cases


# False Failure 2
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] <.05 & x[i, 3]>0 & x[i, 4:503] <.33, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    #print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))

FF2P <- total_cases

# False Failure 2 (negative direction)
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] <.05 & x[i, 3] < 0 & x[i, 4:503] < .33, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    #print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))

FF2N <- total_cases

# True Success 2
# Assuming your matrix is named 'x'
total_cases <- 0

# Loop through each row to apply the criteria
for (i in 1:500) {
    count_in_row <- sum(x[i, 2] >.05 & x[i, 4:503] <.33, na.rm = TRUE)
    total_cases <- total_cases + count_in_row
    #print(paste("Row", i, "cases:", count_in_row))  # Print the count for each row
}

# Print the final total number of cases
print(paste("Final total cases:", total_cases))

TS2 <- total_cases

FS2P + FS2N + TF2 + FF2P + FF2N + TS2 + ADNE