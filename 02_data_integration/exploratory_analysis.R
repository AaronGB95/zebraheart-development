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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table("datamatrix_qn.txt",
                         sep = "\t",
                         header = TRUE,
                         row.names = 1)
phenodata <- read.table("phenodata.txt",
                        sep = "\t",
                        header = TRUE,
                        row.names = 1)

# Order phenodata rows by datamatrix columns order
phenodata$Sample <- rownames(phenodata)
phenodata <- phenodata[match(colnames(datamatrix), rownames(phenodata)), ]

# Give file name to save plots
file <- "qn_counts"
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Boxplot
##-----------------------------------------------------------------------------
file_name <- paste(file, "_set_boxplot.png", sep = "") # Change group here

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
  ggplot(aes(x = Sample, y = Sample_value, fill = Set)) +   # Change group here
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
file_name <- paste(file, "_age_pca.png", sep = "")  # Change group here

# PCA plot
png(filename = file_name, width = 720, height = 720)
biplot(p,
       colby = "Age",
       legendPosition = "right")  # Change group here
dev.off()
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Clustering by correlation distance
##-----------------------------------------------------------------------------
# Correlation distance calculation
correlacion <- cor(datamatrix)
distancia <- as.dist((1 - correlacion) / 2)
hc <- hclust(distancia)
samples_order <- colnames(datamatrix)
dend_data <- dendro_data(hc, type = "rectangle")

# Name for Clustering plot
file_name <- paste(file, "_age_cluster.png", sep = "")  # Change group here

png(filename = file_name, width = 720, height = 720)

plot <- ggplot(dend_data$segments) +
  labs(title = file_name,
       x = "Sample",
       y = "Distance") +
  geom_segment(aes(x = x,
                   y = y,
                   xend = xend,
                   yend = yend),
               linewidth = 1.5) +
  theme(legend.text = element_text(size = 17),
        legend.title = element_text(size = 17),
        axis.title=element_text(size = 17))

if (!is.null(group)) {
  if (is.null(samples_order)) {
    stop("samples_order is needed")
  }
  if (length(samples_order) != length(group) ||
        length(group) != length(dend_data$labels$label)) {
    stop("lengths are not equal")
  }
  group <- group[match(as.vector(dend_data$labels$label), samples_order)]
  plot <- plot + scale_y_continuous(expand = c(0.2, 0)) +
    geom_text(data = dend_data$labels,
              aes(x,
                  y,
                  label = label,
                  group = group),
              hjust = 1,
              size = 5,
              angle = 90) +
    theme_classic()
}

plot

dev.off()
##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())