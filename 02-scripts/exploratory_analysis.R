##-----------------------------------------------------------------------------
##
## Exploratory Analysis
##
## This script performs an exploratory analysis of an expression matrix
## through a boxplot, PCA and clustering analysis
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

# Load required packages
require(rstudioapi)
require(tidyverse)
require(data.table)
require(knitr)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(ggdendro)
require(PCAtools)
require(tidyr)
require(tibble)
require(RColorBrewer)

# Set file path to working directory
setwd("..")
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table("03-data/datamatrix.txt",
                         sep = "\t",
                         header = TRUE,
                         row.names = 1)
phenodata <- read.table("01-documentation/phenodata.txt",
                        sep = "\t",
                        header = TRUE,
                        row.names = 1)

# Order phenodata rows by datamatrix columns order
phenodata$Sample <- rownames(phenodata)
phenodata <- phenodata[match(colnames(datamatrix), rownames(phenodata)), ]

# Give file name to save plots
file <- "raw_counts"
dir_plots <- "05-results/plots"
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Boxplot
##-----------------------------------------------------------------------------
file_name <- paste0(dir_plots, "/_", file, "_age_boxplot.png")

png(filename = file_name, width = 720, height = 720)
par(mar = c(10, 4.1, 4.1, 2.1))
par(cex.axis = 1.3)
par(las = 2)

# For the dataset
group <- as.factor(phenodata$Set)

log2(datamatrix + 1) %>%
  rownames_to_column("Genes") %>%
  gather(Sample, Sample_value, -Genes) %>%
  left_join(phenodata, by = "Sample") %>%
  ggplot(aes(x = Sample, y = Sample_value, fill = Age)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Boxplot with Set
##-----------------------------------------------------------------------------
file_name <- paste(file, "_set_boxplot.png", sep = "")

png(filename = file_name, width = 720, height = 720)
par(mar = c(10, 4.1, 4.1, 2.1))
par(cex.axis = 1.3)
par(las = 2)

# For the dataset
group <- as.factor(phenodata$Set)

log2(datamatrix + 1) %>%
  rownames_to_column("Genes") %>%
  gather(Sample, Sample_value, -Genes) %>%
  left_join(phenodata, by = "Sample") %>%
  ggplot(aes(x = Sample, y = Sample_value, fill = Set)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## PCA
##-----------------------------------------------------------------------------
# Create PCA object
p <- pca(datamatrix, metadata = phenodata)

# Name for PCA plot file
file_name <- paste0(dir_plots, "/_", file, "_age_pca.png")

# PCA plot
png(filename = file_name, width = 720, height = 720)

biplot(p,
       colby = "Age",
       legendPosition = "right")

dev.off()

# Name for PCA plot file
file_name <- paste(file, "_set_pca.png", sep = "")  # Change group here

# PCA plot
png(filename = file_name, width = 720, height = 720)

biplot(p,
       colby = "Set",
       legendPosition = "right")  # Change group here

dev.off()
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Correlation clustering
##-----------------------------------------------------------------------------
file_name <- paste0(dir_plots, "/_", file, "_age_cluster.png")

png(filename = file_name, width = 720, height = 720)

par(mar=c(3.1, 0.1, 0.1, 1.1))

correlacion <- cor(datamatrix)

distancia <- as.dist((1 - correlacion) / 2)

dd <- hclust(distancia)

ddata_x <- dendro_data(dd)

ddata_x$labels <- merge(label(ddata_x),
                        phenodata,
                        by.x = "label",
                        by.y = "Sample")

# Correlation with Age
file_name <- paste(file, "_age_cluster.png", sep = "")  # Change group here

png(filename = file_name, width = 720, height = 720)

par(mar=c(3.1, 0.1, 0.1, 1.1))



dendroplot <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = label(ddata_x),
            aes(label = label, x = x, y = 0, colour = Age, hjust = 0)) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  labs(color = "Age")

dendroplot

dev.off()

# Correlation with Set
file_name <- paste(file, "_set_cluster.png", sep = "")  # Change group here

png(filename = file_name, width = 720, height = 720)

par(mar=c(3.1, 0.1, 0.1, 1.1))



dendroplot <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data = label(ddata_x),
            aes(label = label, x = x, y = 0, colour = Set, hjust = 0)) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  labs(color = "Set")

dendroplot

dev.off()
##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())