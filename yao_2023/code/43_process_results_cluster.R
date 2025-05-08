library(glue)
library(data.table)

# work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'  
work_directory <- '/raid6/Tianyu/convergence_risk_gene/try_Cleary_data/'
source(glue(work_directory, 'R/collect_and_structure_results.R'))

batch_name <- '42_pairwise'
summary_file <- glue(work_directory, 'data/intermediate_data/', batch_name, '/all_results_summary.csv')

clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

# Reading and processing saved results
if (file.exists(summary_file)) {
  results_summary <- fread(summary_file)
} else {
  stop("Summary file not found! Ensure that the results have been generated.")
}

# Read all results into a list
all_results <- lapply(results_summary$file_path, load_results)
names(all_results) <- paste(results_summary$treatment1, results_summary$treatment2, sep = "_vs_")


processed_results <- lapply(all_results, process_result, group = clustering$cluster_index)
processed_results_dt <- rbindlist(processed_results, fill = TRUE, idcol = "comparison")
fwrite(processed_results_dt, glue(work_directory, 'data/intermediate_data/', batch_name, '/processed_results.csv'))
message("Processed results saved.")

