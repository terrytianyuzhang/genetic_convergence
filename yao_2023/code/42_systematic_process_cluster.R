library(glue)
library(data.table)
library(grpreg)
library(HMC)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'
# work_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'  

# source(glue(work_directory, 'R/group_lasso_function.R'))
source(glue(work_directory, 'R/convergence.R'))
control_subsetting_sample_size <- 500

residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_small.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

####
clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

dir.create(glue(work_directory, 'data/intermediate_data/42_pairwise/'), recursive = TRUE, showWarnings = FALSE)
dir.create(glue(work_directory, 'log/'), recursive = TRUE, showWarnings = FALSE)

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

#####
all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']
table(residual_subset$Guides_collapsed_by_gene)

summary_results <- list()
log_file <- glue(work_directory, 'log/logs.txt')

for(treatment1_name in all_treatment_names){
  for(treatment2_name in all_treatment_names){
    if(treatment1_name == treatment2_name)
      next
    else{
      message(glue("Processing {treatment1_name} vs {treatment2_name}"))
      log_message <- glue("{Sys.time()} - Processing {treatment1_name} vs {treatment2_name}\n")
      write(log_message, file = log_file, append = TRUE)
      
      process_data <- preprocess_one_setting(residual_subset, treatment1_name, treatment2_name, clustering)
      test_result <- convergence_testing(process_data$control,
                                         process_data$treatment1,
                                         process_data$treatment2,
                                         pca_method = "dense_pca",
                                         classifier_method = "group_lasso",
                                         lambda_type = 'lambda.min',
                                         n_folds = 5,
                                         group = clustering$cluster_index,
                                         verbose = TRUE)
      
      output_filename <- glue(work_directory, 'data/intermediate_data/42_pairwise/{treatment1_name}_{treatment2_name}.rds')
      saveRDS(test_result, file = output_filename)
      
      summary_results <- append(summary_results, list(data.table(treatment1 = treatment1_name, 
                                                                 treatment2 = treatment2_name, 
                                                                 file_path = output_filename)))
    }
  }
}

# Save summary file
summary_results_dt <- rbindlist(summary_results)
fwrite(summary_results_dt, glue(work_directory, 'data/intermediate_data/42_pairwise/all_results_summary.csv'))
message("Processing completed. Summary saved.")



