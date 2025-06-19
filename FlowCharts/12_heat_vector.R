library(ggplot2)
library(scales)

# Output directory
output_dir <- "~/Documents/genetic_convergence/FlowCharts/data/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Vector of group types
my_plot <- c("control", 
             "top_increase", "tail_increase",
             "top_increase", "top_increase_more",
             "top_increase", "top_decrease",
             "top_decrease", "top_mixture")

# Function to generate vector values
generate_vector_values <- function(group, vec_length = 20) {
  if (group == "control") {
    values <- rnorm(vec_length, mean = 0, sd = 0.3)
  } else if (group == "top_increase") {
    values <- c(rnorm(5, mean = 2.5, sd = 0.3), rnorm(vec_length - 5, mean = 0, sd = 0.5))
  } else if (group == "tail_increase") {
    values <- c(rnorm(vec_length - 5, mean = 0, sd = 0.5), rnorm(5, mean = 2.5, sd = 0.3))
  } else if (group == "top_decrease") {
    values <- c(rnorm(5, mean = -2.5, sd = 0.3), rnorm(vec_length - 5, mean = 0, sd = 0.5))
  } else if (group == "top_increase_more") {
    values <- c(rnorm(7, mean = 2.5, sd = 0.3), rnorm(vec_length - 7, mean = 0, sd = 0.5))
  } else if (group == "top_mixture") {
    values <- c(rnorm(2, mean = 2.5, sd = 0.3), rnorm(3, mean = -2.5, sd = 0.3), rnorm(vec_length - 5, mean = 0, sd = 0.5))
  } else {
    stop("Unrecognized group name.")
  }
  return(values)
}

# Gradient function: blue → gray → yellow
blue_gray_red <- colorRampPalette(c("#013e75", "#aaaaaa", "#dddddd", "#aaaaaa", "#f5b70a"))

# Loop through each group and save plots
for (plot_index in 1:length(my_plot)) {
  
  group <- my_plot[plot_index]
  vec_length <- 20
  values <- generate_vector_values(group, vec_length)
  
  # Normalize to [0, 1]
  norm_values <- (values + 3) / 6
  norm_values <- pmin(pmax(norm_values, 0), 1)
  
  # Map to colors
  colors <- blue_gray_red(100)[ceiling(norm_values * 99) + 1]
  
  # Create data frame
  df <- data.frame(
    x = 1:vec_length,
    y = 1,
    fill_color = colors
  )
  
  # Save to PDF
  pdf(file = file.path(output_dir, paste0("heat_vector_", plot_index, ".pdf")), 
      width = 5, height = 1)
  ggplot(df, aes(x = x, y = y, fill = fill_color)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    coord_fixed(ratio = 1)
  dev.off()
}

