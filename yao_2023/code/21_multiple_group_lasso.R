library(glue)
library(data.table)
library(grpreg)
library(HMC)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
source(glue(work_directory, 'R/group_lasso_function.R'))
control_subsetting_sample_size <- 2000

residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_small.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

####
clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)


one_pair_convergence <- function(treatment1_name, treatment2_name, clustering){
  ######
  control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
  treatment1 <- residual_subset[Guides_collapsed_by_gene == treatment1_name, -"Guides_collapsed_by_gene"]
  treatment2 <- residual_subset[Guides_collapsed_by_gene == treatment2_name, -"Guides_collapsed_by_gene"]
  
  control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]
  treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
  treatment2 <- treatment2[, match(clustering$gene_name, colnames(treatment2)), with = FALSE]
  
  ####first fit the treatment 1 group
  print("starting the first fit")
  control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
  selected_beta <- fit_grlasso(control, treatment1, clustering)
  
  ####
  selected_beta_binary <- selected_beta[-1] != 0
  retained_genes <- clustering[selected_beta_binary, ]
  group_count <- table(retained_genes$cluster_index)
  group_forward <- as.numeric(names(group_count)[group_count >= group_truncation_num])
  
  if(length(group_forward) == 0){
    warning("the first treatment does not select any groups")
    return(list(treatment1 = treatment1_name,
                treatment2 = treatment2_name,
                group_count = 0,
                group_count2 = 0,
                anchor_test_stat = 0,
                convergence_gene = NA))
  }
  
  clustering <- clustering[cluster_index %in% group_forward, ]
  filtered_gene <- clustering$gene_name
  
  control <- control[, colnames(control) %in% filtered_gene, with = FALSE]
  treatment2 <- treatment2[, colnames(treatment2) %in% filtered_gene, with = FALSE]
  
  print("starting the second fit")
  selected_beta2 <- fit_grlasso(control, treatment2, clustering)
  
  #####
  selected_beta_binary2 <- selected_beta2[-1] != 0
  retained_genes2 <- clustering[selected_beta_binary2, ]
  group_count2 <- table(retained_genes2$cluster_index)
  group_count2
  ###
  
  anchor_test <- anchored_lasso_testing(control, treatment2, pca_method = 'dense_pca')
  summarize_feature_name(anchor_test)
  
  return(list(treatment1 = treatment1_name,
              treatment2 = treatment2_name,
              group_count = group_count,
              group_count2 = group_count2,
              anchor_test_stat = anchor_test$test_statistics,
              convergence_gene = summarize_feature_name(anchor_test)))
}


#####
all_treatment_names <- unique(residual_subset$Guides_collapsed_by_gene)
all_treatment_names <- all_treatment_names[all_treatment_names != 'non-targeting']
table(residual_subset$Guides_collapsed_by_gene)

all_pairs_result <- list()
pair_count <- 1
for(treatment1_name in all_treatment_names){
  for(treatment2_name in all_treatment_names){
    if(treatment1_name == treatment2_name)
      next
    else{
      print(glue("tratment 1 is ", treatment1_name))
      print(glue("tratment 2 is ", treatment2_name))
      
      one_pair_result <- one_pair_convergence(treatment1_name, treatment2_name, clustering)
      all_pairs_result[[pair_count]] <- one_pair_result
      pair_count <- pair_count + 1
      
      saveRDS(all_pairs_result, file = glue(work_directory, 'data/intermediate_data/multiple_paris_first_pass.rds'))
    }
  }
}

