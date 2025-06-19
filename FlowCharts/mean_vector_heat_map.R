library(pheatmap)
library(RColorBrewer)

# Simulate gene expression data
set.seed(123)
expression_matrix <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)

# Set a random proportion (e.g., 20%) to zero
proportion_zero <- 0.4
n_total <- length(expression_matrix)
n_zero <- floor(proportion_zero * n_total)

zero_indices <- sample(n_total, n_zero)
expression_matrix[zero_indices] <- 0

# Define color palette
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Plot minimal heatmap
pheatmap(expression_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         legend = FALSE,
         color = heatmap_colors,
         border_color = NA)





library(ggplot2)
red_length <- 4
blue_length <- 2
gray_length <- 9
# Generate a matrix with random values

color_matrix <- matrix(runif(gray_length), nrow = 1, ncol = gray_length)

# Normalize the matrix values to be between 0 and 1
color_matrix <- 1 - 0.5*(color_matrix - min(color_matrix)) / (max(color_matrix) - min(color_matrix))

# Define a function to convert a matrix of values to colors
value_to_gray <- function(x) {
    rgb(x, x, x)
}

# Convert the matrix values to colors
gray_matrix <- apply(color_matrix, c(1, 2), value_to_gray)

#########

if(red_length > 0){
  color_matrix <- matrix(runif(red_length), nrow = 1, ncol = red_length)
  # Normalize the matrix values to be between 0 and 1
  color_matrix <- 1 - 0.2*(color_matrix - min(color_matrix)) / (max(color_matrix) - min(color_matrix))
  # Define a function to convert a matrix of values to colors
  value_to_color <- function(x) {
    rgb(x, 0.5*x, 0.5*x)
  }
  
  # Convert the matrix values to colors
  red_matrix <- apply(color_matrix, c(1, 2), value_to_color)
}else{
  red_matrix <- NULL
}


if(blue_length > 0){
  color_matrix <- matrix(runif(blue_length), nrow = 1, ncol = blue_length)
  # Normalize the matrix values to be between 0 and 1
  color_matrix <- 1 - 0.2*(color_matrix - min(color_matrix)) / (max(color_matrix) - min(color_matrix))
  # Define a function to convert a matrix of values to colors
  value_to_color <- function(x) {
    rgb(0.5*x, 0.5*x, x)
  }
  
  # Convert the matrix values to colors
  blue_matrix <- apply(color_matrix, c(1, 2), value_to_color)
}
plot_matrix <- cbind(blue_matrix, red_matrix, gray_matrix)

# Create a data frame from the matrix for ggplot2
df <- expand.grid(x = 1:nrow(plot_matrix), y = 1:ncol(plot_matrix))
df$color <- as.vector(plot_matrix)

# Plot the matrix of colored cells
ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

