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
require(plyr)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table("datamatrix_qn.txt",
                         sep = "\t",
                         header = 1,
                         row.names = 1)

phenodata <- read.table("phenodata.txt",
                        sep = "\t",
                        header = 1,
                        row.names = 1)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data processing
##-----------------------------------------------------------------------------
# Add Age as a time variable
phenodata$Time <- revalue(phenodata$Age, c("48 hpf"="0",
                                           "72 hpf"="1",
                                           "120 hpf"="2",
                                           "Adult"="3"))

# Turn variable to factor
phenodata$Time <- as.factor(as.character(phenodata$Time))

# Add a dummy variable
phenodata$Animal <- "Fish"
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Time series analysis
##-----------------------------------------------------------------------------
# Make design matrix
design <- make.design.matrix(phenodata,
                             time.col = "Time",
                             repl.col = "Set",
                             group.cols = NULL)

# Find significant genes
significant_genes <- p.vector(datamatrix, design, counts = TRUE)

# Find significant differences
significant_differences <- T.fit(significant_genes)

# Get list of significant genes
get <- get.siggenes(significant_differences, vars = "all")
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Results plots
##-----------------------------------------------------------------------------
suma2Venn(get$summary[, c(2:4)])
##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())