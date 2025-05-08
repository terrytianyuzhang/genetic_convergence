work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'    
file_name <- glue::glue(work_directory, 'data/GSM6858447_KO_conventional.rds')

Cleary_raw <- readRDS(file_name)
# Cleary_raw@meta.data$"10X_channel" #channel informaton? in total 19 of them
treatment_cells_binary <- Cleary_raw@meta.data$Guides_collapsed_by_gene %in% c('STAT2', 'TYK2')
control_cells_binary <- Cleary_raw@meta.data$Guides_collapsed_by_gene %in% c('non-targeting')

hist(Cleary_raw@meta.data$G2M_score)
table(Cleary_raw@meta.data$Cell_cycle_phase)
Cleary_raw@assays$perturbations
perturbation_matrix <- GetAssayData(Cleary_raw, assay = "perturbations") #if a cell gets both non-targeting and a perturbation, meta data thinks it gets a perturbation only

Cleary_subset <- Cleary_raw[, Cleary_raw$Guides_collapsed_by_gene %in% c('non-targeting', 'STAT2', 'TYK2')]
saveRDS(Cleary_subset, file = glue::glue(work_directory,'/data/Cleary_subset.rds'))


mean_exp <- rowMeans(Cleary_subset@assays$RNA@counts/Cleary_subset$Total_RNA_count)
genes_selected <- names(sort.int(mean_exp, decreasing = T))[1:2000]

saveRDS(genes_selected, file = glue::glue(work_directory,'/data/genes_selected.rds'))

CSCORE_result <- CSCORE(Cleary_subset, genes = genes_selected)

CSCORE_coexp <- CSCORE_result$est

dim(CSCORE_coexp)

CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0

# Compute the adjacency matrix based on the co-expression matrix
adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 6)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
# Run hierarchical clustering as in the WGCNA workflow
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                     distM = dissTOM, 
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

names(memb) = genes_selected
memb_tab <- table(memb)
module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))
module_list

saveRDS(module_list, file = glue::glue(work_directory,'/data/gene_clustering.rds'))


####verify the distribution between two genes identified in the same cluster
RNA_matrix <- Cleary_subset@assays$RNA@counts

inds <- which(CSCORE_coexp >0. & row(CSCORE_coexp) != col(CSCORE_coexp), arr.ind = TRUE)
# Create a data frame including row names, column names, and the values
significant_pairs <- data.frame(
  row = rownames(CSCORE_coexp)[inds[, "row"]],
  col = colnames(CSCORE_coexp)[inds[, "col"]],
  value = CSCORE_coexp[inds]
)

hist(significant_pairs$value)
# Display the result
print(significant_pairs)
for(index in 1:nrow(significant_pairs)){
  plot_pair <- significant_pairs[index,]
  gene_A <- RNA_matrix[plot_pair$row,]
  gene_B <- RNA_matrix[plot_pair$col,]
  # plot(gene_A, gene_B)
  print(cor(gene_A, gene_B))
  significant_pairs$p_cor[index] <- cor(gene_A, gene_B)
}

plot(significant_pairs[1:426,]$value, significant_pairs[1:426,]$p_cor)


####check batch effect
library(ggplot2)
DefaultAssay(Cleary_subset) <- "RNA"

# (1) Data preprocessing: Normalize and scale the data
Cleary_subset <- NormalizeData(Cleary_subset)
Cleary_subset <- FindVariableFeatures(Cleary_subset)
Cleary_subset <- ScaleData(Cleary_subset)

# (2) Run PCA using the variable features
Cleary_subset <- RunPCA(Cleary_subset, features = VariableFeatures(object = Cleary_subset))

# (3) Visualize the PCA results: Plot cells in the PC1 vs. PC2 space
DimPlot(Cleary_subset, reduction = "pca", group.by = "10X_channel") +
  ggtitle("PCA Projection Colored by Batch")

DimPlot(Cleary_subset, reduction = "pca", group.by = "Cell_cycle_phase") +
  ggtitle("PCA Projection Colored by Batch") ##cell cycle shows strong batch effect

DimPlot(Cleary_subset, reduction = "pca", group.by = "Guides_collapsed_by_gene") +
  ggtitle("PCA Projection Colored by Batch")



