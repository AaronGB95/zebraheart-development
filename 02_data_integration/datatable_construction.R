##-----------------------------------------------------------------------------
##
## Data table construction
##
## Script that takes multiple featureCounts output files and merges them into
## a single data table.
##
## Author: Aarón García Blázquez
##
## Date: 25/09/2023
##
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
rm(list=ls())

# Required packages
require(rstudioapi)
require(tidyverse)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------

# Get featureCounts output files
files <- list.files(path = getwd(), pattern = "counts.txt", full.names=TRUE)

# Loop through each file and read it into a data frame
for (i in 1:length(files)) {
  load(files[i])
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Gene Names
##-----------------------------------------------------------------------------
# Save counts into list of objects
datafiles <- ls(pattern="counts")

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
colnames(datamatrix) <- strsplit(files, split = "_counts.txt")

# Loop through files and merge all counts into a single table
for (fi in datafiles) {
    datamatrix[,fi] <- fi[,-1]
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Save
##-----------------------------------------------------------------------------
write.table(datamatrix, file="datamatrix.txt", sep ="\t")
##-----------------------------------------------------------------------------

# Clean environment
rm(list=ls())