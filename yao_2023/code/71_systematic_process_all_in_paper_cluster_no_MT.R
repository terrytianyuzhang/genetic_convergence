message('this code is similar to 51_systematic_process....')
message('the MT module, module 42, is removed')
message('and I used standardized genes to perform the convergence analysis')

library(glue)
library(data.table)
library(grpreg)
library(HMC)
library(doParallel)
library(foreach)

rm(list = ls())
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
work_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'

source(glue(work_directory, 'R/convergence.R'))

residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_all_in_paper_no_MTmodule.rds'))
# hist(table(residual_subset[Guides_collapsed_by_gene!= 'non-targeting']$Guides_collapsed_by_gene), breaks = 100)
# length(unique(residual_subset$Guides_collapsed_by_gene))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- glue(work_directory, "data/module_list_df_2000_genes_no_MTmodule.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

# Subset control samples
set.seed(1)
control_subsetting_sample_size <- 500
control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]

dir.create(glue(work_directory, 'data/intermediate_data/71_pairwise/'), recursive = TRUE, showWarnings = FALSE)
dir.create(glue(work_directory, 'log/'), recursive = TRUE, showWarnings = FALSE)

preprocess_one_setting <- function(residual_subset, treatment1_name, treatment2_name, clustering) {
  treatment1 <- residual_subset[Guides_collapsed_by_gene == treatment1_name, -"Guides_collapsed_by_gene"]
  treatment2 <- residual_subset[Guides_collapsed_by_gene == treatment2_name, -"Guides_collapsed_by_gene"]
  
  treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
  treatment2 <- treatment2[, match(clustering$gene_name, colnames(treatment2)), with = FALSE]
  
  return(list(treatment1 = treatment1, treatment2 = treatment2))
}

# Parallelization setup
num_cores <- 20
cl <- makeCluster(num_cores)
registerDoParallel(cl)

all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']

log_file <- glue(work_directory, 'log/logs.txt')

# Generate all unique treatment pairs
treatment_pairs <- expand.grid(treatment1 = all_treatment_names, treatment2 = all_treatment_names)
treatment_pairs <- treatment_pairs[treatment_pairs$treatment1 != treatment_pairs$treatment2, ]

summary_results <- foreach(i = 1:nrow(treatment_pairs), .combine = rbind, .packages = c("data.table", "glue", "grpreg", "HMC")) %dopar% {
  treatment1_name <- treatment_pairs$treatment1[i]
  treatment2_name <- treatment_pairs$treatment2[i]
  
  log_message <- glue("{Sys.time()} - Processing {treatment1_name} vs {treatment2_name} in parallel\n")
  write(log_message, file = log_file, append = TRUE)
  
  process_data <- preprocess_one_setting(residual_subset, treatment1_name, treatment2_name, clustering)
  test_result <- convergence_testing(control,
                                     process_data$treatment1,
                                     process_data$treatment2,
                                     pca_method = "dense_pca",
                                     classifier_method = "group_lasso",
                                     lambda_type = 'lambda.min',
                                     n_folds = 5,
                                     group = clustering$cluster_index,
                                     standardize_feature = TRUE,
                                     verbose = TRUE)

  output_filename <- glue(work_directory, 'data/intermediate_data/71_pairwise/{treatment1_name}_{treatment2_name}.rds')
  saveRDS(test_result, file = output_filename)
  
  return(data.table(treatment1 = treatment1_name, treatment2 = treatment2_name, file_path = output_filename))
}

# Save summary file
fwrite(summary_results, glue(work_directory, 'data/intermediate_data/71_pairwise/all_results_summary.csv'))

# Stop parallel cluster
stopCluster(cl)

message("Processing completed. Summary saved.")


