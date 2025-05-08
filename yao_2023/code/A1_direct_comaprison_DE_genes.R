library(glue)
library(data.table)
library(grpreg)
library(HMC)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
source(glue(work_directory, 'R/group_lasso_function.R'))

residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_small.rds'))
residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

#####
table(residual_subset$Guides_collapsed_by_gene)
treatment1_name <- "IRAK1"
treatment2_name <- "IRAK4"

######
control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
treatment1 <- residual_subset[Guides_collapsed_by_gene == treatment1_name, -"Guides_collapsed_by_gene"]
treatment2 <- residual_subset[Guides_collapsed_by_gene == treatment2_name, -"Guides_collapsed_by_gene"]
control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]

#####

get_DE_gene <- function(control, treatment){
  anchor_test <- anchored_lasso_testing(control, treatment, pca_method = 'dense_pca', n_folds = 10)
  DE_gene <- summarize_feature_name(anchor_test)
  return(list(anchor_test$test_statistics, DE_gene))
}

DE_gene1 <- get_DE_gene(control, treatment1)
DE_gene2 <- get_DE_gene(control, treatment2)

print(DE_gene1[[1]])
print(DE_gene2[[1]])
intersect(DE_gene1[[2]], DE_gene2[[2]])

