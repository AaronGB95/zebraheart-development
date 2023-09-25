##-----------------------------------------------------------------------------
##
## Exploratory Analysis
##
## This script performs an exploratory analysis of an expression matrix
## through a boxplot, PCA and clustering analysis
##
## Author: Aarón García Blázquez
##
## Date: 06-02-2023
##
## Email: aaron.garcia.blazquez@gmail.com
##
##-----------------------------------------------------------------------------

rm(list = ls())

## Required packages
##----------------------

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


## Data load
##--------------------

# Change the working directory to script path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Introduce study ID
dataset <- "dataset_combat"

dirdata <- paste0(paste0("data/datasets/",dataset))
archivos <- ".RData"

datamatrix <- paste0(dirdata, archivos)
datamatrix <- loadRData(datamatrix)
datamatrix <- as.data.frame(datamatrix)

phenodata <- "data/datasets/phenodata.RData"
phenodata <- loadRData(phenodata)

dirplots <- paste0("plots/", dataset)
dir.create(dirplots, showWarnings = FALSE)

file <- paste0(dirplots, paste0("/", dataset))


# PCA
p <- pca(datamatrix, metadata = phenodata)

# Clustering analysis by correlation distance
correlacion <- cor(datamatrix)
distancia <- as.dist((1- correlacion) / 2)
hc <- hclust(distancia)

samples_order <- colnames(datamatrix)

dend_data <- dendro_data(hc, type = "rectangle")

# Change group here
group_choice <- phenodata$Set
palette_name <- "Set3"
color_palette <-  c("#004949","#24ff24","#ff6db6","#db6d00",
                    "#490092","#006ddb","#b66dff","#6db6ff",
                    "#920000","#924900")
phenodata$colors <- color_palette[as.factor(group_choice)]



## Boxplot
##---------
  
file_name <- paste(file, "_set_boxplot.png", sep = "") # Change group here

png(filename = file_name,width = 720, height = 720)
par(mar = c(10, 4.1, 4.1, 2.1))
par(cex.axis = 1.3)
par(las=2)

# For the dataset
group <- as.factor(group_choice)
phenodata$Sample <- rownames(phenodata)

datamatrix %>%
  rownames_to_column("Genes") %>%                         
  gather(Sample, Sample_value, -Genes) %>%                 
  left_join(phenodata, by = "Sample") %>%        
  ggplot(aes(x=Sample, y=Sample_value, fill=Set)) +   # Change group here
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  scale_fill_manual(values = color_palette)

dev.off()
  



## PCA plot
##--------------------------------------

file_name <- paste(file, "_set_pca.png", sep = "")  # Change group here

png(filename = file_name ,width = 720, height = 720)

biplot(p, colby = "Set", colkey = color_palette, legendPosition = "right")  # Change group here

dev.off()




## Clustering plot
##----------------------------------------------

file_name <- paste(file, "_set_cluster.png", sep = "")  # Change group here

png(filename = file_name ,width = 720, height = 720)

plot <- ggplot(dend_data$segments) + 
  labs(title = file_name, x = "Sample", y = "Distance") + 
  scale_fill_manual(name = "Group", palette = palette_name, values = color_palette) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth=1.5) + 
  theme(legend.text = element_text(size=17), legend.title = element_text(size=17), axis.title=element_text(size=17))

if (!is.null(group)) {
  if (is.null(samples_order)) {
    stop("samples_order is needed")
  }
  if (length(samples_order) != length(group) | length(group) != length(dend_data$labels$label)) {
    stop("lengths are not equal")
  }
  group <- group[match(as.vector(dend_data$labels$label), samples_order)]
  
  plot <- plot + scale_y_continuous(expand = c(0.2, 0)) +
    geom_text(data = dend_data$labels, aes(x,y, label = label, group = group, color = group_choice), hjust = 1, size = 5, angle = 90) +
    scale_color_manual(values = color_palette) +
    theme_classic()
}

plot

dev.off()




  


