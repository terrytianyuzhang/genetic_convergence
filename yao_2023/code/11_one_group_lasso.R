library(glue)
library(data.table)
library(grpreg)

work_directory <- '/Users/tianyuzhang/Documents/convergence_risk_gene/try_Cleary_data/'    
residual_subset <- readRDS(glue(work_directory, 'data/intermediate_data/residual_matrix_small.rds'))

residual_subset <- residual_subset[, -c("ID", "Cell_cycle_phase")]
residual_subset <- data.table(residual_subset)
colnames(residual_subset) <- gsub("\\.", "-", colnames(residual_subset))

table(residual_subset$Guides_collapsed_by_gene)
control <- residual_subset[Guides_collapsed_by_gene == "non-targeting", -"Guides_collapsed_by_gene"]
treatment1 <- residual_subset[Guides_collapsed_by_gene == "MEF2C", -"Guides_collapsed_by_gene"]

####
clustering_file <- glue(work_directory, "data/module_list_df_2000_genes.csv")
clustering <- data.table(read.csv(clustering_file))
clustering$gene_name <- gsub("\\.", "-", clustering$gene_name)

control <- control[, match(clustering$gene_name, colnames(control)), with = FALSE]
treatment1 <- treatment1[, match(clustering$gene_name, colnames(treatment1)), with = FALSE]
control[,1:5]
treatment1[,1:5]
clustering[1:5,]

####
control <- control[sample(1:nrow(control), 500), ]
expression <- rbind(control, treatment1)
group_label <- c(rep(0, nrow(control)), rep(1, nrow(treatment1)))
# fit <- grpreg(X = expression, 
#               y = group_label, 
#               group = clustering$gene_name, 
#               penalty = "grLasso", 
#               family = 'binomial')

cv_fit <- cv.grpreg(X = expression, 
                    y = group_label, 
                    group = clustering$gene_name, 
                    penalty = "grLasso", 
                    family = 'binomial')
selected_beta <- coef(cv_fit, s = "lambda.min")
plot(selected_beta)

selected_beta_binary <- selected_beta[-1] != 0
retained_genes <- clustering[selected_beta_binary, ]
table(retained_genes$cluster_index)
