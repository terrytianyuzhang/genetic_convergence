rm(list = ls())
library(glue)
work_directory <- '~/Documents/genetic_convergence/yao_2023/'    
residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_all_in_paper.rds'))
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)
MT_module_gens <- clustering[cluster_index == 42, ]$gene_name
residual_subset <- residual_subset[, !..MT_module_gens]

saveRDS(residual_subset, 
        glue(work_directory, 'data/intermediate_data/residual_matrix_all_in_paper_no_MTmodule.rds'))
