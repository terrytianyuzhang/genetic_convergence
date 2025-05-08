# Simulation under null hypothesis
simulate_null_data <- function(n_samples=500, n_features=100) {
  mu_con <- rep(0, n_features)
  mu_tr1 <- c(rep(1, 5), rep(0, n_features - 5))
  mu_tr2 <- c(rep(0, 5), rep(0.5, 5), rep(0, n_features - 10))
  
  control <- mvrnorm(n_samples, mu_con, diag(n_features))
  treatment1 <- mvrnorm(n_samples, mu_tr1, diag(n_features))
  treatment2 <- mvrnorm(n_samples, mu_tr2, diag(n_features))
  
  rownames(control) <- paste("Sample", 1:n_samples, sep = "_")
  rownames(treatment1) <- paste("Sample", 1:n_samples, sep = "_")
  rownames(treatment2) <- paste("Sample", 1:n_samples, sep = "_")
  
  colnames(control) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment1) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment2) <- paste("Feature", 1:n_features, sep = "_")
  
  return(list(control = control, treatment1 = treatment1, treatment2 = treatment2))
}

# Set the number of samples
n_samples <- 50

# Run the simulation
data <- simulate_null_data(n_samples = n_samples)
control <- data$control
treatment1 <- data$treatment1
treatment2 <- data$treatment2


test_result <- convergence_testing(control, treatment1, treatment2)
test_result$test_statistic
test_result$p_value

# Dynamically create the split indices based on n_samples
split_index <- list(1:floor(n_samples/2), (floor(n_samples/2) + 1):n_samples)

# Assign the split indices to each group
control_split_index <- split_index
tr1_split_index <- split_index
tr2_split_index <- split_index

result <- process_fold(1, control, treatment1, treatment2, control_split_index, tr1_split_index, tr2_split_index, "dense_pca", "lasso", TRUE)
print(result)

control <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
estimate_leading_pc(control)

# Run simulation
n_sims <- 100
test_statistics <- c()

for (i in 1:n_sims) {
  if (i %% 10 == 0) {
    cat("Simulation", i, "\n")
  }
  data <- simulate_null_data()
  control <- data$control
  treatment1 <- data$treatment1
  treatment2 <- data$treatment2
  T <- calculate_test_statistic(data$control, data$treatment1, data$treatment2)
  
  test_statistics <- c(test_statistics, T)
}



# Plot results
par(mfrow = c(1, 2))

# Histogram
hist(test_statistics, breaks = 30, prob = TRUE, main = "Distribution of Test Statistics")
curve(dnorm(x, mean = 0, sd = 1), add = TRUE, col = "red", lwd = 2)

# Q-Q plot
qqnorm(test_statistics)
qqline(test_statistics, col = "red")