# Number of simulations
rm()
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
source(glue(work_directory, 'R/convergence.R'))

simulate_null_data <- function(n_samples=500, n_features=100) {
  mu_con <- rep(0, n_features)
  mu_tr1 <- c(rep(5,5), rep(0, n_features - 5))
  mu_tr2 <- c(rep(0.2, 5),  rep(0, n_features - 10), rep(10, 5))
  
  ratio <- 3
  control <- mvrnorm(n_samples*ratio, mu_con, diag(c(1:5, rep(1, n_features - 5))))
  treatment1 <- mvrnorm(n_samples, mu_tr1, diag(n_features))
  treatment2 <- mvrnorm(n_samples, mu_tr2, diag(n_features))
  
  rownames(control) <- paste("Sample", 1:(n_samples*ratio), sep = "_")
  rownames(treatment1) <- paste("Sample", 1:n_samples, sep = "_")
  rownames(treatment2) <- paste("Sample", 1:n_samples, sep = "_")
  
  colnames(control) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment1) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment2) <- paste("Feature", 1:n_features, sep = "_")
  
  return(list(control = control, treatment1 = treatment1, treatment2 = treatment2))
}

n_simulations <- 30

# Store test statistics
n_samples <- 50
n_features <- 200
test_statistics <- numeric(n_simulations)
beta_final_values <- beta_values <- vector("list", n_simulations)
n_folds <- 10

for (i in 1:n_simulations) {
  set.seed(i)
  # Simulate null data
  data <- simulate_null_data(n_samples = n_samples, n_features = n_features)
  control <- data$control
  treatment1 <- data$treatment1
  treatment2 <- data$treatment2
  
  # Run the test
  test_result <- convergence_testing(control, treatment1, treatment2, 
                                     lambda_type = 'lambda.min', n_folds = n_folds)
  
  # Store the test statistic
  test_statistics[i] <- test_result$test_statistic
  
  # Collect first_beta for each fold
  beta_values[[i]] <- lapply(1:n_folds, function(M) test_result$fold_data[[M]]$first_beta)
  # beta_final_values[[i]] <- lapply(1:5, function(M) test_result$fold_data[[M]]$final_beta)
  
}
# test_result$fold_data[[2]]$first_beta
# Plot histogram
ggplot(data.frame(test_statistics), aes(x = test_statistics)) +
  geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.6, color = "black") +
  labs(title = "Distribution of Test Statistic Under Null Hypothesis", x = "Test Statistic", y = "Frequency")

# Q-Q plot
qqnorm(test_statistics)
abline(c(0,1), col = "red")

shapiro_test <- shapiro.test(test_statistics)
print(shapiro_test)

# Plot first_beta values across folds and simulations
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert beta_values into a tidy data frame
beta_df <- do.call(rbind, lapply(seq_along(beta_values), function(i) {
  data.frame(
    simulation = i,
    fold = rep(1:n_folds, each = length(beta_values[[i]][[1]])),  # Adjust based on beta dimensions
    beta_index = rep(seq_along(beta_values[[i]][[1]]), times = n_folds),
    first_beta = unlist(beta_values[[i]])
  )
}))

# Create a separate data frame for test statistics (one per simulation)
test_stats_df <- data.frame(
  simulation = 1:n_simulations,
  test_statistic = test_statistics
)

# Merge beta_df with test statistics for annotation
beta_df <- beta_df %>%
  left_join(test_stats_df, by = "simulation")

# Plot first_beta with test statistics as annotations
ggplot(beta_df, aes(x = beta_index, y = first_beta, color = as.factor(fold))) +
  geom_line(alpha = 0.6) +
  facet_wrap(~simulation, scales = "free_y") +
  labs(title = "First Beta Values Across Simulations and Folds",
       x = "Beta Index",
       y = "First Beta",
       color = "Fold") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10, face = "bold")) +  # Make facet labels bold
  geom_text(data = test_stats_df, 
            aes(x = Inf, y = Inf, 
                label = paste("Test Stat:", round(test_statistic, 2))),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1, size = 4)

