# changed end from 108 to 243
# Initialize variables
start <- 1  # starting number
end <- 243  # ending limit because params_o has 243 conditions
rows_to_keep <- c()  # vector to store the rows_to_keep

while (start <= end) {
  # Add the current start and the next increment by 4
  for (i in 1:2) {
    if (start <= end) {
      rows_to_keep <- c(rows_to_keep, start)
      start <- start + 4
    }
  }
  # Add the switch value (next start point)
  if (start <= end) {
    rows_to_keep <- c(rows_to_keep, start)
    start <- start + 1  # Move to the next start point
  }
}
#matrix(rows_to_keep, ncol = 3, byrow = TRUE)
