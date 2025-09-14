library(glue)
library(data.table)
library(ggplot2)
rm(list = ls())
work_directory <- '/Users/tianyuzhang/Documents/genetic_convergence/'
source(glue(work_directory, 'R/collect_and_structure_results.R'))

high_num_pair <- fread(file = paste0(work_directory, "/yao_2023/data/intermediate_data/B1_interesting_pairs/greater_than_4_modules.csv"))

batch_name <- '51_pairwise'
processed_results <- fread(glue(work_directory, '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))
# Assuming `prossed_results` is a data.table

processed_results <- processed_results[, .(comparison, p_value, active_group)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(processed_results, processed_results, 
                        by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))

merge(x = high_num_pair,
      y = merged_results[, .(gene1, gene2, active_group.x, active_group.y)],
      all.x = TRUE)


processed_results[gene1 == 'TRAF6' & gene2 == "TRIB1", ]
processed_results[gene1 == 'TRIB1' & gene2 == "TRAF6", ]

processed_results[gene1 == 'STAT1' & gene2 == "STAT2", ]
processed_results[gene1 == 'STAT2' & gene2 == "STAT1", ]

