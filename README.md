# genetic_convergence
Statistical method for genetic convergence analysis using cross-fitting

# Rep. structure 

```
.
├── R
│   ├── (to_be_remove?)group_lasso_function.R
│   ├── collect_and_structure_results.R
│   ├── convergence.R
│   └── group_lasso_function.R
├── README.md
└── yao_2023
    ├── code
    │   ├── 01_subsetting_genes.R
    │   ├── 02_remove_cell_cycle.R
    │   ├── 03_structure_residuals.R
    │   ├── 04_subsetting_perturbation.R
    │   ├── 05_remove_MT_genes.R
    │   ├── 11_one_group_lasso.R
    │   ├── 12_two_group_lasso.R
    │   ├── 13_under_the_null.R
    │   ├── 21_multiple_group_lasso.R
    │   ├── 22_process_test_results.R
    │   ├── 31_establish_clustering.R
    │   ├── 32_GO_analysis.R
    │   ├── 33_remove_mt_genes.R
    │   ├── 41_first_try_new_convergence_code.R
    │   ├── 42_systematic_process_cluster.R
    │   ├── 43_process_results_cluster.R
    │   ├── 44_visualization.R
    │   ├── 51_systematic_process_all_in_paper_cluster.R
    │   ├── 52_process_results_cluster.R
    │   ├── 53_visualization.R
    │   ├── 54_graph_visualization.R
    │   ├── 61_cluster_correlation.R
    │   ├── 71_systematic_process_all_in_paper_cluster_no_MT.R
    │   ├── A1_direct_comaprison_DE_genes.R
    │   ├── B1_testing_my_convergence.R
    │   ├── B2_test_if_normal.R
    │   ├── B3_test_if_gr_lasso_normal.R
    │   ├── B4_test_helper_function.R
    │   ├── C1_correlation_using_original.R
    │   ├── QC.zip
    │   ├── establish_clustering.R
    │   ├── first_group_lasso_classification.R
    │   └── read_in_data.R
    ├── data
    │   ├── CSCORE_result_2000_gene.rds
    │   ├── Cleary_subset.rds
    │   ├── GSM6858447_KO_conventional_FRPerturb_effect_sizes.csv.gz
    │   ├── ego_result.rds
    │   ├── ego_result_167_genes_second.rds
    │   ├── ego_result_2000_genes.rds
    │   ├── ego_result_4000_genes.rds
    │   ├── final_data
    │   ├── finer_module_list_167_genes.rds
    │   ├── finer_module_list_2000_genes.rds
    │   ├── gene_clustering.rds
    │   ├── genes_selected.rds
    │   ├── intermediate_data
    │   ├── module_list_df_167_genes.csv
    │   ├── module_list_df_2000_genes.csv
    │   ├── module_list_df_2000_genes_no_MTmodule.csv
    │   ├── top_ego_result_2000_genes.txt
    │   └── top_ego_result_4000_genes.txt
    ├── log
    └── report
        ├── 42_pairwise
        ├── 51_pairwise
        ├── 51_pairwise_igraph
        ├── Kathryn_presentation
        ├── PCA_cell_cycle_after.pdf
        ├── PCA_cell_cycle_before.pdf
        ├── PCA_chennel_before.pdf
        ├── first_pass_result.csv
        └── violin_plot_QC.pdf```