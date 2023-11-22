##-----------------------------------------------------------------------------
##
## time_series_analysis.R
##
## Script that performs a time series analysis of gene expression data
## with a time variable and optional additional factors.
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
rm(list = ls())

# Load required packages
require(maSigPro)
require(plyr)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table("datamatrix_tmm.txt",
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
phenodata$Time <- as.factor(as.numeric(phenodata$Time))

# Add a dummy variable
phenodata$Animal <- "Fish"

phenodata$Sample <- rownames(phenodata)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Time series analysis
##-----------------------------------------------------------------------------
# Make design matrix
design <- make.design.matrix(phenodata,
                             degree = 1,
                             time.col = 3,
                             repl.col = 1)

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

##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())