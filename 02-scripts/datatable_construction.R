##-----------------------------------------------------------------------------
##
## datatable_construction.R
##
## Script that takes multiple featureCounts output files and merges them into
## a single data table. Normalizes the counts in the table
## by TMM normalization and Quantile Normalization (QN).
##
## Author: Aarón García Blázquez
##
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------
# Get featureCounts output files
files <- list.files(path = paste0(dir_data, "counts/"),
                    full.names = TRUE)

# Loop through each file and read it into a data frame
datalist <- lapply(files,
                   FUN = read.table,
                   header = TRUE)

annotations <- read.table(paste0(dir_docs, "annotations.txt"),
                          sep = "\t",
                          header = TRUE)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene Names
##-----------------------------------------------------------------------------
# Get gene names
names <- lapply(datalist, "[", , "Geneid")

# Make intersection of gene names
names <- Reduce(unique, names)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene Counts
##-----------------------------------------------------------------------------
# Create empty matrix with desired dims
datamatrix <- matrix(NA, nrow = length(names), ncol = length(datalist))

# Apply gene IDs as rownames
rownames(datamatrix) <- names

# Apply sample IDs as colnames
colnames(datamatrix) <- str_extract_all(files, "(?<=\\/)[^\\/]*(?=_)")

# Loop through files and merge all counts into a single table
for (i in seq_along(datalist)) {
  datamatrix[, i] <- datalist[[i]][, 7]
}

datamatrix <- as.data.frame(datamatrix)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## TMM normalization
##-----------------------------------------------------------------------------
# Create a DGEList object
dge_list <- DGEList(counts = datamatrix)

# Calculate library sizes and normalize for RNA composition
dge_list <- calcNormFactors(dge_list, method = "TMM")

# Get the normalized counts
tmm_counts <- cpm(dge_list, log = FALSE)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Quantile normalization
##-----------------------------------------------------------------------------
# Calculate the normalized quantiles
qn_counts <- normalize.quantiles(as.matrix(tmm_counts),
                                 copy = TRUE,
                                 keep.names = TRUE)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Save
##-----------------------------------------------------------------------------
write.table(datamatrix,
            file = paste0(dir_data, "datatables/datamatrix_raw.txt"),
            sep = "\t",
            quote = FALSE)
write.table(tmm_counts,
            file = paste0(dir_data, "datatables/datamatrix_tmm.txt"),
            sep = "\t",
            quote = FALSE)
write.table(qn_counts,
            file = paste0(dir_data, "datatables/datamatrix_qn.txt"),
            sep = "\t", 
            quote = FALSE)
##-----------------------------------------------------------------------------
