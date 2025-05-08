library(glue)
library(data.table)
library(grpreg)
library(HMC)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
source(glue(work_directory, 'R/group_lasso_function.R'))
source(glue(work_directory, 'R/convergence.R'))
control_subsetting_sample_size <- 500

residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_small.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

preprocess_one_setting <- function(residual_subset, treatment1_name, treatment2_name, clustering){
  
  control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
  treatment1 <- residual_subset[Guides_collapsed_by_gene == treatment1_name, -"Guides_collapsed_by_gene"]
  treatment2 <- residual_subset[Guides_collapsed_by_gene == treatment2_name, -"Guides_collapsed_by_gene"]
  
  control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]
  treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
  treatment2 <- treatment2[, match(clustering$gene_name, colnames(treatment2)), with = FALSE]
  
  control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
  return(list(control = control,
              treatment1 = treatment1,
              treatment2 = treatment2))
}

treatment1_name <- 'STAT2'
treatment2_name <- 'ADO'
process_data <- preprocess_one_setting(residual_subset, treatment1_name, treatment2_name, clustering)

test_result <- convergence_testing(process_data$control, 
                    process_data$treatment1, 
                    process_data$treatment2,
                    pca_method = "dense_pca",
                    classifier_method = "group_lasso",
                    lambda_type = 'lambda.min',
                    n_folds = 5,
                    group = clustering$cluster_index,
                    verbose = TRUE
)
test_result$p_value
active_gene_name <- collect_active_features(test_result$fold_data)
table(clustering[gene_name  %in% active_gene_name, ]$cluster_index)

