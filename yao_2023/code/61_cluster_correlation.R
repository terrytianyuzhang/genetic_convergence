# Load required libraries
library(glue)
library(data.table)
library(corrplot)
library(CSCORE)
rm(list = ls())
# Set working directory
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'   

# === Load and preprocess residual matrix ===
residual_subset <- readRDS(glue("{work_directory}data/intermediate_data/residual_matrix_small.rds"))
residual_subset <- data.table(residual_subset[, !c("ID", "Cell_cycle_phase")])
residual_subset <- residual_subset[Guides_collapsed_by_gene == 'non-targeting', ]
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

# === Load clustering file ===
clustering_file <- glue("{work_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering[, gene_name := gsub("\\.", "-", gene_name)]

# mean(clustering[, .N, by = cluster_index]$N)
# === Select gene subset based on cluster indices ===
plot_cluster_index <- c(1, 31, 40, 42)
subset_gene_name <- clustering[cluster_index %in% plot_cluster_index, gene_name]

# === Subset residual matrix based on selected genes ===
selected_column <- colnames(residual_subset) %in% subset_gene_name
residual_subset <- residual_subset[, ..selected_column]

# === Compute correlation matrix ===
cor_mat <- cor(residual_subset)

# === Load CSCORE results ===
gene_count <- 2000
CSCORE_result <- readRDS(glue("{work_directory}data/CSCORE_result_{gene_count}_gene.rds"))
CSCORE_coexp <- CSCORE_result$est
CSCORE_p <- CSCORE_result$p_value

# === Standardize gene names ===
genes_selected <- clustering$gene_name
colnames(CSCORE_coexp) <- gsub("\\.", "-", colnames(CSCORE_coexp))
rownames(CSCORE_coexp) <- gsub("\\.", "-", rownames(CSCORE_coexp))

# === Adjust p-values and threshold coexpression matrix ===
p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
# CSCORE_coexp[p_matrix_BH > 0.05] <- 0

# === Plot coexpression matrix for selected subset ===
mat <- CSCORE_coexp[subset_gene_name, subset_gene_name]
corrplot(mat)

# === Visualize sampled subset of selected clusters ===
sampled_genes <- lapply(plot_cluster_index, function(idx) {
  genes <- clustering[cluster_index == idx, gene_name]
  sample(genes, min(20, length(genes)))
})

# Combine all sampled genes into a single vector
gene_names <- unique(unlist(sampled_genes))

CSCORE_mat <- CSCORE_coexp[gene_names, gene_names]
corrplot(CSCORE_mat)

output_file <- glue("{work_directory}/report/Kathryn_presentation/CSCORE_corrplot_sampled_genes_original.pdf")

# Save the plot to PDF (size in inches)
pdf(file = output_file, width = 10, height = 10)
corrplot(CSCORE_mat)
dev.off()


original_coexp <- cor(residual_subset)
original_mat <- original_coexp[gene_names, gene_names]
corrplot(original_mat)

