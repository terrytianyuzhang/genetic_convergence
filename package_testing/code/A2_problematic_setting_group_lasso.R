library(HMC)
library(MASS)
library(ggplot2)

set.seed(123)
library(HMC)
library(MASS)
library(ggplot2)

set.seed(123)

# Helper: expand module-level specs into a per-feature vector
.expand_module_vector <- function(x, n_modules, module_size, n_features) {
  if (length(x) == 1L) {
    rep(x, n_features)
  } else if (length(x) == n_modules) {
    rep(x, each = module_size)
  } else if (length(x) == n_features) {
    x
  } else {
    stop("Length of vector must be 1, n_modules, or n_features.")
  }
}

# Helper: build block-diagonal covariance from per-module SDs or a single SD
.build_block_cov <- function(sd_within, n_modules, module_size, n_features) {
  sds <- .expand_module_vector(sd_within, n_modules, module_size, n_features)
  diag(sds^2)
}

# Helper: feature names as feature_(moduleIndex)_(orderWithinModule)
.feature_names <- function(n_modules, module_size, prefix = "feature") {
  mods <- rep(seq_len(n_modules), each = module_size)
  ords <- rep(seq_len(module_size), times = n_modules)
  paste0(prefix, "_", mods, "_", ords)
}

# -------- Modular simulator --------
simulate_modular_data <- function(
    n_samples   = 500,
    n_modules   = 7,            # number of modules
    module_size = 30,           # dimensions per module
    ratio       = 3,            # control:treatment sample ratio
    # Per-module mean shifts (can be length 1, n_modules, or n_features)
    mu_control  = 0,
    mu_tr1      = 0,
    mu_tr2      = 0,
    # Per-module within-feature SDs (same length rules as above)
    sd_within   = 1,
    # Alternatively supply full covariance matrices (overrides sd_within if given)
    Sigma_control = NULL,
    Sigma_tr1     = NULL,
    Sigma_tr2     = NULL,
    feature_prefix = "feature"
) {
  n_features <- n_modules * module_size
  
  # Expand means
  mu_con <- .expand_module_vector(mu_control, n_modules, module_size, n_features)
  mu_t1  <- .expand_module_vector(mu_tr1,     n_modules, module_size, n_features)
  mu_t2  <- .expand_module_vector(mu_tr2,     n_modules, module_size, n_features)
  
  # Build covariance
  if (is.null(Sigma_control)) Sigma_control <- .build_block_cov(sd_within, n_modules, module_size, n_features)
  if (is.null(Sigma_tr1))     Sigma_tr1     <- .build_block_cov(sd_within, n_modules, module_size, n_features)
  if (is.null(Sigma_tr2))     Sigma_tr2     <- .build_block_cov(sd_within, n_modules, module_size, n_features)
  
  # Draw
  control    <- MASS::mvrnorm(n_samples * ratio, mu = mu_con, Sigma = Sigma_control)
  treatment1 <- MASS::mvrnorm(n_samples,         mu = mu_t1,  Sigma = Sigma_tr1)
  treatment2 <- MASS::mvrnorm(n_samples,         mu = mu_t2,  Sigma = Sigma_tr2)
  print(Sigma_control)
  # Names
  rn_c  <- paste0("Sample_", seq_len(n_samples * ratio))
  rn_t1 <- paste0("Sample_", seq_len(n_samples))
  rn_t2 <- paste0("Sample_", seq_len(n_samples))
  cn    <- .feature_names(n_modules, module_size, prefix = feature_prefix)
  
  rownames(control)    <- rn_c
  rownames(treatment1) <- rn_t1
  rownames(treatment2) <- rn_t2
  colnames(control) <- colnames(treatment1) <- colnames(treatment2) <- cn
  
  list(
    control = control,
    treatment1 = treatment1,
    treatment2 = treatment2,
    meta = list(
      n_samples = n_samples,
      n_modules = n_modules,
      module_size = module_size,
      n_features = n_features
    )
  )
}

n_samples  <- 300
n_modules  <- 7
module_sz  <- 5

mu_control_mods <- c(1, 1, 1, 1, 1, 1, 1)   # length n_modules
mu_tr1_mods <- c(2, 2, 1, 5, 1, 5, 1)   # length n_modules
mu_tr2_mods <- c(2, 2, 2, 1, 5, 1, 5)   # length n_modules

sim <- simulate_modular_data(
  n_samples   = n_samples,
  n_modules   = n_modules,
  module_size = module_sz,
  ratio       = 3,
  mu_control  = mu_control_mods,               # scalar => 0 everywhere
  mu_tr1      = mu_tr1_mods,     # vector of length n_modules
  mu_tr2      = mu_tr2_mods,     # vector of length n_modules
  sd_within   = c(1,1,1,1,1,1,1)                # homoskedastic by default
)

# ------------- Your simulation loop scaffolding -------------
n_simulations <- 1
n_folds <- 5
test_statistics <- numeric(n_simulations)
grouping_vector <- rep(1:n_modules, each = module_sz)
names(grouping_vector) <- colnames(sim$control)

for (i in 1:n_simulations) {
  set.seed(i)
  print(i)
  result <- convergence_testing(sim$control, 
                                sim$treatment1, 
                                sim$treatment2,
                                classifier_method = "group_lasso",
                                n_folds = 5,
                                group = grouping_vector,
                                lambda_type = 'lambda.1se')
  test_statistics[i] <- result$test_statistic
}

qq_df <- data.frame(
  sample_quantiles = sort(test_statistics),
  theoretical_quantiles = qnorm(ppoints(n_simulations))
)

ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
  geom_point(color = "steelblue", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Test Statistics", x = "Theoretical Quantiles (Normal)", y = "Sample Quantiles") +
  theme_minimal()

visualize_convergence_top_genes(
  fold_data = result$fold_data[1:3],
  top_n = 30
)

collect_active_features(result, 
                        group = grouping_vector,
                        group_threshold = 3)

