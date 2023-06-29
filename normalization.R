# Load the edgeR package
library(edgeR)

# Load the RNAseq data table
load("data/datasets/discarded_excluded_unnorm/constructed_dataset_intersect.RData")

# Create a DGEList object
dge_list <- DGEList(counts = selected_data)

# Calculate library sizes and normalize for RNA composition
dge_list <- calcNormFactors(dge_list, method = "TMM")

# Get the normalized counts
norm_counts <- cpm(dge_list, log = FALSE)

# Save the normalized counts to a file
save(norm_counts, file = "data/datasets/dataset_intersect_tmm.RData")
