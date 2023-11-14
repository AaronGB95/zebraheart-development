##-----------------------------------------------------------------------------
##
## heatmap.R
##
## Script that takes a datamatrix and creates a heatmap with the
## expression data
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
require(ComplexHeatmap)
require(InteractiveComplexHeatmap)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------
datamatrix <- read.table("datamatrix_qn.txt", sep = "\t", header = TRUE)
phenodata <- read.table("phenodata.txt", sep = "\t", header = TRUE)

phenodata$name <- with(phenodata,
                       paste0(phenodata$Sample, sep = " - ", phenodata$Age))
colnames(datamatrix) <- phenodata$name
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Heatmap
##-----------------------------------------------------------------------------
ht <- Heatmap(datamatrix)
ht <- draw(ht)

htShiny(ht, output_ui_float = TRUE)
##-----------------------------------------------------------------------------


