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
## Date: 25/09/2023
##
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
# Clean environment
rm(list=ls())

# Required packages
require(rstudioapi)
require(tidyverse)
require(edgeR)
require(preprocessCore)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------
# Get featureCounts output files
files <- list.files(path = getwd(), pattern = 'counts.txt', full.names=TRUE)

# Loop through each file and read it into a data frame
for (i in 1:length(files)) {
  load(files[i])
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene Names
##-----------------------------------------------------------------------------
# Save counts into list of objects
datafiles <- ls(pattern='counts')

# Get gene names
names <- lapply(datafiles, function(x) rownames(get(x)))

# Make intersection of gene names
names <- Reduce(unique, names)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene Counts
##-----------------------------------------------------------------------------
datamatrix <- matrix(NA, nrows=length(names), ncols=length(datafiles))
rownames(datamatrix) <- names
colnames(datamatrix) <- strsplit(files, split = '_counts.txt')

# Loop through files and merge all counts into a single table
for (fi in datafiles) {
    datamatrix[,fi] <- fi[,-1]
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## TMM normalization
##-----------------------------------------------------------------------------
# Create a DGEList object
dge_list <- DGEList(counts = datamatrix)

# Calculate library sizes and normalize for RNA composition
dge_list <- calcNormFactors(dge_list, method = 'TMM')

# Get the normalized counts
tmm_counts <- cpm(dge_list, log = FALSE)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Quantile normalization
##-----------------------------------------------------------------------------
# Calculate the normalized quantiles
qn_counts <- normalize.quantiles(as.matrix(tmm_counts), copy = TRUE)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Save
##-----------------------------------------------------------------------------
write.table(datamatrix, file='datamatrix.txt', sep ='\t')
write.table(tmm_counts, file='datamatrix_tmm.txt', sep='\t')
write.table(qn_counts, file='datamatrix_qn.txt', sep='\t')
##-----------------------------------------------------------------------------

# Clean environment
rm(list=ls())