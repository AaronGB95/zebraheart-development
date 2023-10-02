##-----------------------------------------------------------------------------
##
## time_series_analysis.R
##
## Script that performs a time series analysis of gene expression data
## with a time variable and optional additional factors.
##
## Author: Aarón García Blázquez
##
## Date: 25/09/2023 # nolint: commented_code_linter.
##
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
# Clean environment
rm(list = ls())

# Load required packages
require(maSigPro)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table("datamatrix.txt",
                         sep = "\t",
                         header = 1,
                         row.names = 1)

phenodata <- read.table("phenodata.txt",
                        sep = "\t",
                        header = 1,
                        row.names = 1)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Time series analysis
##-----------------------------------------------------------------------------
# Make design matrix
design <- make.design.matrix(phenodata)

# Find significant genes
significant_genes <- p.vector(datamatrix, design, counts = TRUE)

# Find significant differences
significant_differences <- T.fit(NBp)

# Get list of significant genes
get <- get.siggenes(NBt, vars = "all")
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Results plots
##-----------------------------------------------------------------------------
suma2Venn(get$summary[, c(2:4)])
##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())