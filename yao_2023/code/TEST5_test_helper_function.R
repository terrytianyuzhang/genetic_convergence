# Number of simulations
rm()
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
source(glue(work_directory, 'R/convergence.R'))

simulate_null_data <- function(n_samples=500, n_features=100) {
  mu_con <- rep(0, n_features)
  mu_tr1 <- c(rep(5, 5), rep(0, n_features - 5))
  mu_tr2 <- c(rep(0, 10), rep(0, n_features - 20), rep(5, 10))
  
  ratio <- 1
  control <- mvrnorm(n_samples * ratio, mu_con, diag(c(1:5, rep(1, n_features - 5))))
  treatment1 <- mvrnorm(n_samples, mu_tr1, diag(n_features))
  treatment2 <- mvrnorm(n_samples, mu_tr2, diag(n_features))
  
  rownames(control) <- paste("Sample", 1:(n_samples * ratio), sep = "_")
  rownames(treatment1) <- paste("Sample", 1:n_samples, sep = "_")
  rownames(treatment2) <- paste("Sample", 1:n_samples, sep = "_")
  
  colnames(control) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment1) <- paste("Feature", 1:n_features, sep = "_")
  colnames(treatment2) <- paste("Feature", 1:n_features, sep = "_")
  
  return(list(control = control, treatment1 = treatment1, treatment2 = treatment2))
}

# Define grouping structure for Group Lasso
generate_group_structure <- function(n_features, group_size) {
  n_groups <- ceiling(n_features / group_size)
  return(rep(1:n_groups, each = group_size, length.out = n_features))
}

n_simulations <- 1

# Store test statistics
n_samples <- 200
n_features <- 100
group_size <- 5  # Example: Group features into 20 groups

test_statistics <- numeric(n_simulations)
beta_final_values <- beta_values <- vector("list", n_simulations)
n_folds <- 5

for (i in 1:n_simulations) {

  # Simulate null data
  data <- simulate_null_data(n_samples = n_samples, n_features = n_features)
  control <- data$control
  treatment1 <- data$treatment1
  treatment2 <- data$treatment2
  
  # Generate grouping vector
  group_vector <- generate_group_structure(n_features, group_size)
  
  # Run the test using Group Lasso
  test_result <- convergence_testing(control, treatment1, treatment2, 
                                     classifier_method = "lasso",
                                     lambda_type = 'lambda.min', 
                                     n_folds = n_folds,
                                     group = group_vector)
  
  # Store the test statistic
  test_statistics[i] <- test_result$test_statistic
  
  # Collect first_beta for each fold
  beta_values[[i]] <- lapply(1:n_folds, function(M) test_result$fold_data[[M]]$first_beta)
  # beta_final_values[[i]] <- lapply(1:n_folds, function(M) test_result$fold_data[[M]]$final_beta)
}

active_features <- collect_active_features(test_result)
collect_active_features(test_result, group = group_vector)
