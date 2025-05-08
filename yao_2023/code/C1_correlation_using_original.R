# === Load required libraries ===
library(glue)
library(data.table)
library(corrplot)
library(CSCORE)
library(Seurat)

# === Clear workspace and set working directory ===
rm(list = ls())
work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'

# === Load and preprocess residual matrix ===
residual_subset <- readRDS(glue("{work_directory}data/intermediate_data/residual_matrix_small.rds"))
residual_subset <- data.table(residual_subset[, !c("ID", "Cell_cycle_phase")])
residual_subset <- residual_subset[Guides_collapsed_by_gene == 'non-targeting', ]
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

# === Load clustering file and select gene subset ===
clustering_file <- glue("{work_directory}data/module_list_df_2000_genes.csv")
clustering <- fread(clustering_file)
clustering[, gene_name := gsub("\\.", "-", gene_name)]

plot_cluster_index <- c(1, 31, 40, 42)
subset_gene_name <- clustering[cluster_index %in% plot_cluster_index, gene_name]

# === Subset residual matrix and compute correlation ===
selected_column <- colnames(residual_subset) %in% subset_gene_name
residual_subset_small <- residual_subset[, ..selected_column]
cor_mat <- cor(residual_subset_small)

# === Load CSCORE results ===
gene_count <- 2000
CSCORE_result <- readRDS(glue("{work_directory}data/CSCORE_result_{gene_count}_gene.rds"))
CSCORE_coexp <- CSCORE_result$est
CSCORE_p <- CSCORE_result$p_value

# === Standardize gene names ===
genes_selected <- clustering$gene_name
colnames(CSCORE_coexp) <- gsub("\\.", "-", colnames(CSCORE_coexp))
rownames(CSCORE_coexp) <- gsub("\\.", "-", rownames(CSCORE_coexp))

# === Adjust p-values (BH) and optionally threshold ===
p_matrix_BH <- matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] <- p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)
# CSCORE_coexp[p_matrix_BH > 0.05] <- 0  # Optional thresholding

# === Visualize CSCORE matrix for selected subset ===
mat <- CSCORE_coexp[subset_gene_name, subset_gene_name]
corrplot(mat)

# === Sample a small subset of genes from each cluster ===
sampled_genes <- lapply(plot_cluster_index, function(idx) {
  genes <- clustering[cluster_index == idx, gene_name]
  sample(genes, min(20, length(genes)))
})
gene_names <- unique(unlist(sampled_genes))

# === Plot CSCORE and original correlation matrices for sampled genes ===
CSCORE_mat <- CSCORE_coexp[gene_names, gene_names]
corrplot(CSCORE_mat)

original_coexp <- cor(residual_subset_small)
original_mat <- original_coexp[gene_names, gene_names]
corrplot(original_mat)

# === Load Seurat object and extract non-targeting subset ===
file_name <- glue("{work_directory}data/intermediate_data/after_subsetting.rds")
Cleary <- readRDS(file_name)
Cleary_subset <- subset(Cleary, subset = Guides_collapsed_by_gene == 'non-targeting')

# === Extract and normalize expression matrix ===
expr_mat <- as.matrix(Cleary_subset@assays$RNA@counts)
expr_mat_t <- t(expr_mat)  # Cells as rows, genes as columns
colnames(expr_mat_t) <- gsub("\\.", "-", colnames(expr_mat_t))
expr_mat_t <- expr_mat_t[, subset_gene_name]
expr_mat_t_norm <- expr_mat_t / rowSums(expr_mat_t)

# === Compute and visualize correlation matrix ===
original_coexp <- cor(expr_mat_t_norm, method = "pearson")
original_mat <- original_coexp[gene_names, gene_names]
corrplot(original_mat)

# === Example gene-gene plots ===
plot(expr_mat_t[, 'SERF2'], expr_mat_t[, 'MT-ND1'])
plot(expr_mat_t_norm[, 'SERF2'], expr_mat_t_norm[, 'MT-ND1'], xlim = c(0, 0.05))
cor(expr_mat_t_norm[, 'SERF2'], expr_mat_t_norm[, 'MT-ND1'])
