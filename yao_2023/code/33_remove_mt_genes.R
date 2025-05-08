library(glue)
library(data.table)

rm(list = ls())
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'

clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering <- clustering[cluster_index != 42, ]

output_file <- glue("{work_directory}data/module_list_df_2000_genes_no_MTmodule.csv")
fwrite(clustering, output_file)
