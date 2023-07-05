# zebraheart-development

This is the repository of the Master's thesis: Integrative study of zebrafish heart transcriptome through development

Each folder contains several scripts required to carry out this project:

01_data_processing: Contains the scripts that download, annotate and TMM normalize the selected data. Also containst the script that performs the exploratory analysis.

02_data_integration: Contains the script that integrates the selected samples of each into a single dataset. Also normalizes the dataset by quantile normalization.

03_batch_correction: Performs the proposed batch correction strategies and a exploratory analysis for each.

04_differential_expression: Performs the differential expression using edgeR and plot the results into a volcano plot.

05_functional_enrichment: performs an over-representation analysis for both upregulated and downregulated genes for each ontology.
