library(pheatmap)
library(RColorBrewer)

# Set up output directory (optional)
output_dir <- "~/Documents/genetic_convergence/FlowCharts/data/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Function to generate one heatmap with random zeros
generate_and_save_heatmap <- function(file_name, n_row = 10, n_col = 20, zero_prop = 0.2) {
  # Simulate data
  expression_matrix <- matrix(rnorm(n_row * n_col), nrow = n_row, ncol = n_col)
  
  # Randomly zero out a proportion of values
  n_total <- length(expression_matrix)
  n_zero <- floor(zero_prop * n_total)
  zero_indices <- sample(n_total, n_zero)
  expression_matrix[zero_indices] <- 0
  
  # Color palette
  heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  
  # Save heatmap as PDF (vector graphic)
  pdf(file = file.path(output_dir, file_name), width = 5, height = 2.5)
  pheatmap(expression_matrix,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           show_colnames = FALSE,
           legend = FALSE,
           color = heatmap_colors,
           border_color = NA)
  dev.off()
}

# Generate and save three heatmaps
generate_and_save_heatmap("heatmap1.pdf")
generate_and_save_heatmap("heatmap2.pdf")
generate_and_save_heatmap("heatmap3.pdf")
