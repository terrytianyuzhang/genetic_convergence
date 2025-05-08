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
treatment1_name <- "IRAK4"
treatment2_name <- "ADO"

######
control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
treatment1 <- residual_subset[Guides_collapsed_by_gene == treatment1_name, -"Guides_collapsed_by_gene"]
treatment2 <- residual_subset[Guides_collapsed_by_gene == treatment2_name, -"Guides_collapsed_by_gene"]

####
clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]
treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
treatment2 <- treatment2[, match(clustering$gene_name, colnames(treatment2)), with = FALSE]

control[,1:5]
treatment1[,1:5]
treatment2[,1:5]
clustering[1:5,]

####first fit the treatment 1 group
control <- control[sample(1:nrow(control), control_subsetting_sample_size), ]
selected_beta <- fit_grlasso(control, treatment1, clustering)

####
selected_beta_binary <- selected_beta[-1] != 0
retained_genes <- clustering[selected_beta_binary, ]
group_count <- table(retained_genes$cluster_index)
group_forward <- as.numeric(names(group_count)[group_count >= group_truncation_num])

clustering <- clustering[cluster_index %in% group_forward, ]
filtered_gene <- clustering$gene_name

control <- control[, colnames(control) %in% filtered_gene, with = FALSE]
treatment2 <- treatment2[, colnames(treatment2) %in% filtered_gene, with = FALSE]

selected_beta2 <- fit_grlasso(control, treatment2, clustering)

#####
selected_beta_binary2 <- selected_beta2[-1] != 0
retained_genes2 <- clustering[selected_beta_binary2, ]
group_count2 <- table(retained_genes2$cluster_index)
group_count2
###

anchor_test1 <- anchored_lasso_testing(control, treatment2, pca_method = 'dense_pca')
anchor_test1$test_statistics
summarize_feature_name(anchor_test1)

