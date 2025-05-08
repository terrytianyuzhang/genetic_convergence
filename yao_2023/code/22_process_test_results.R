library(glue)
library(data.table)
library(grpreg)
library(HMC)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
source(glue(work_directory, 'R/group_lasso_function.R'))

all_pairs_result <- readRDS(file = glue(work_directory, 'data/intermediate_data/multiple_paris_first_pass.rds'))

structured_data <- data.frame()

for(pair_index in 1:length(all_pairs_result)){
  
  current_result <- all_pairs_result[[pair_index]]
  forward_group <- names(current_result$group_count)[current_result$group_count >= group_truncation_num]
  final_group <- names(current_result$group_count2)[current_result$group_count2 >= group_truncation_num]
  
  structured_data <- rbind(structured_data, 
                           data.frame(treatment1 = current_result$treatment1,
                                      treatment2 = current_result$treatment2,
                                      forward_group = paste(forward_group, collapse = "/"),
                                      final_group = paste(final_group, collapse = "/"),
                                      p_value = pnorm(-abs(current_result$anchor_test_stat))))
}

write.csv(structured_data, file = glue(work_directory, 'report/first_pass_result.csv'))
