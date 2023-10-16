##-----------------------------------------------------------------------------
##
## gene_filtering.R
##
## Script that takes an expression matrix and a gene list and filter the rows
## based on the gene list.
##
## Author: Aarón García Blázquez
##
##
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
# Clean environment
rm(list = ls())

# Required packages
require(rstudioapi)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
# Load expression matrix
datamatrix <- read.table("datamatrix_qn.txt", header = TRUE, sep = "\t")

# Load genes list to filter
gene_list <- read.table("heart_genes.txt")
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene filter
##-----------------------------------------------------------------------------
# Convert list to character vector
gene_list <- gene_list[, 1]

# Apply filter
datamatrix_filtered <- datamatrix[gene_list, ]

# Remove missing genes
datamatrix_filtered <- datamatrix_filtered[!is.na(datamatrix_filtered[, 1]), ]
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data save
##-----------------------------------------------------------------------------
write.table(datamatrix_filtered,
            "datamatrix_filtered.txt",
            sep = "\t",
            quote = FALSE)
##-----------------------------------------------------------------------------