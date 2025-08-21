library(glue)
library(data.table)
library(ggplot2)
rm(list = ls())
work_directory <- '/Users/tianyuzhang/Documents/genetic_convergence/'
source(glue(work_directory, 'R/collect_and_structure_results.R'))

batch_name <- '51_pairwise'
processed_results <- fread(glue(work_directory, '/yao_2023/data/intermediate_data/', batch_name, '/processed_results.csv'))
# Assuming `prossed_results` is a data.table

processed_results <- processed_results[, .(comparison, p_value)]
processed_results[, c("gene1", "gene2") := tstrsplit(comparison, "_vs_")]

merged_results <- merge(processed_results, processed_results, 
                        by.x = c("gene1", "gene2"), by.y = c("gene2", "gene1"))
merged_results[, p_min := pmin(p_value.x, p_value.y)]
merged_results[, p_max := pmax(p_value.x, p_value.y)]

merged_results[, p_to_plot := ifelse(gene1 > gene2, p_min, p_max)]

threshold <- 0.025
merged_results[, is_significant := 0]
# merged_results[, is_significant := is_significant + (p_to_plot < threshold)]
merged_results[, is_significant := is_significant + (p_to_plot < threshold/2450)]
merged_results[, is_significant := factor(is_significant, levels = c(0, 1))]

# Plot all gene pairs with color indicating significance
p <- ggplot(merged_results, aes(x = gene1, y = gene2, color = is_significant)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(
    values = c("0" = "#EEEEEE", 
               "1" = "#013e75"
               # , "2" = "#f5b70a"
               ),
    labels = c("No", "Yes"),
    name = "Convergence"
  ) +
  theme_minimal() +
  labs(
    title = "Gene Pair Significance",
    x = "Gene 1", y = "Gene 2"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

plot(p)

# Define the file path without using glue
plot_path <- paste0("/Users/tianyuzhang/Documents/genetic_convergence/paper_plots/data/",batch_name, "_overall_significance.pdf")

# Save the plot as a transparent PDF
ggsave(filename = plot_path, plot = p, width = 9, height = 6*9/8, bg = "transparent")

