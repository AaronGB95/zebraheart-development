##-----------------------------------------------------------------------------
##
## Análisis exploratorio
##
## Este script realiza un análisis exploratorio de una matriz de expresión
## mediante un boxplot, un PCA y un análisis cluster
##
## Autor: Aarón García Blázquez
##
## Fecha de creación: 06-02-2023
##
## Email: aaron.garcia.blazquez@gmail.com
##
##-----------------------------------------------------------------------------

rm(list = ls())

## Librerías necesarias
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


## Carga de los datos
##--------------------

# Cambiamos el directorio de trabajo a la ubicacion del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Introducimos el indicador del estudio
dataset <- "dataset_full_included"

dirdata <- paste0("input/",dataset)
archivos <- ".RData"

datamatrix <- paste0(dirdata, archivos)
datamatrix <- loadRData(datamatrix)
datamatrix <- as.data.frame(datamatrix)

phenodata <- "input/dataset_qn_phenodata.RData"
phenodata <- loadRData(phenodata)

dirplots <- paste0("output/", dataset)
dir.create(dirplots, showWarnings = FALSE)

file <- paste0(dirplots, paste0("/", dataset))


# Realizamos el PCA
p <- pca(datamatrix, metadata = phenodata)

# Análisis cluster por distancia de correlación
correlacion <- cor(datamatrix)
distancia <- as.dist((1- correlacion) / 2)
hc <- hclust(distancia)

samples_order <- colnames(datamatrix)

dend_data <- dendro_data(hc, type = "rectangle")

# Cambiar aquí
group_choice <- phenodata$Set
palette_name <- "Set3"
color_palette <-  c("#004949","#24ff24","#ff6db6","#db6d00",
                    "#490092","#006ddb","#b66dff","#6db6ff",
                    "#920000","#924900")
phenodata$colors <- color_palette[as.factor(group_choice)]



## Boxplot
##---------
  
file_name <- paste(file, "_set_boxplot.png", sep = "") # Cambiar aquí

png(filename = file_name,width = 720, height = 720)
par(mar = c(10, 4.1, 4.1, 2.1))
par(cex.axis = 1.3)
par(las=2)

# Para el conjunto de los datos
group <- as.factor(group_choice)
phenodata$Sample <- rownames(phenodata)

datamatrix %>%
  rownames_to_column("Genes") %>%                         
  gather(Sample, Sample_value, -Genes) %>%                 
  left_join(phenodata, by = "Sample") %>%        
  ggplot(aes(x=Sample, y=Sample_value, fill=Set)) +   # Cambiar aquí
  geom_boxplot() +
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position = "bottom") +
  scale_fill_manual(values = color_palette)

dev.off()
  



## Gráfico del PCA
##--------------------------------------

file_name <- paste(file, "_set_pca.png",
                   sep = "")  # Cambiar aquí

png(filename = file_name,
    width = 720,
    height = 720)

biplot(p,
       colby = "Set", # Cambiar aquí
       colkey = color_palette,
       legendPosition = "bottom") 

dev.off()




## Gráfico del dendrograma
##----------------------------------------------

file_name <- paste(file, "_set_cluster.png", sep = "")  # Cambiar aquí

png(filename = file_name ,width = 720, height = 720)

plot <- ggplot(dend_data$segments) + 
  labs(title = file_name,
       x = "Sample",
       y = "Distance") + 
  scale_fill_manual(name = "Group",
                    palette = palette_name,
                    values = color_palette) +
  geom_segment(aes(x = x,
                   y = y,
                   xend = xend,
                   yend = yend),
               linewidth=1.5) + 
  theme(legend.text = element_text(size=17),
        legend.title = element_text(size=17),
        legend.position = "bottom"
        axis.title=element_text(size=17))

if (!is.null(group)) {
  if (is.null(samples_order)) {
    stop("samples_order is needed")
  }
  if (length(samples_order) != length(group) | length(group) != length(dend_data$labels$label)) {
    stop("lengths are not equal")
  }
  group <- group[match(as.vector(dend_data$labels$label), samples_order)]
  
  plot <- plot + scale_y_continuous(expand = c(0.2, 0)) +
    geom_text(data = dend_data$labels,
              aes(x,y,
                  label = label,
                  group = group,
                  color = group_choice),
              hjust = 1,
              size = 5,
              angle = 90) +
    scale_color_manual(values = color_palette) +
    theme_classic()
}

plot

dev.off()




  


