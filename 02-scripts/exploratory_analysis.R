##-----------------------------------------------------------------------------
##
## Exploratory Analysis
##
## This script performs an exploratory analysis of an expression matrix
## through a boxplot, PCA and clustering analysis
##
## Author: Aarón García Blázquez
##
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data load
##-----------------------------------------------------------------------------
datamatrix <- read.table(paste0(dir_data,"datatables/datamatrix_tmm.txt"),
                         sep = "\t",
                         header = TRUE,
                         row.names = 1)
phenodata <- read.table(paste0(dir_docs, "phenodata_5_groups.txt"),
                        sep = "\t",
                        header = TRUE,
                        row.names = 1)

# Order phenodata rows by datamatrix columns order
phenodata$Sample <- rownames(phenodata)
phenodata <- phenodata[match(colnames(datamatrix), rownames(phenodata)), ]

# Give file name to save plots
file <- "tmm_counts"
dir_plots <- "04-output/plots/"
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Boxplot
##-----------------------------------------------------------------------------
file_name <- paste0(dir_plots, "boxplots/", file, "_age_boxplot.png")

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
file_name <- paste0(dir_plots, "boxplots/", file, "_set_boxplot.png")

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
file_name <- paste0(dir_plots, "pca_plots/", file, "_age_pca.png")

# PCA plot
png(filename = file_name, width = 720, height = 720)

biplot(p,
       colby = "Age",
       legendPosition = "right")

dev.off()

# Name for PCA plot file
file_name <- paste0(dir_plots, "pca_plots/", file, "_set_pca.png")

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

correlacion <- cor(datamatrix)

distancia <- as.dist((1 - correlacion) / 2)

dd <- hclust(distancia)

ddata_x <- dendro_data(dd)

ddata_x$labels <- merge(label(ddata_x),
                        phenodata,
                        by.x = "label",
                        by.y = "Sample")

# Correlation with Age
file_name <- paste0(dir_plots, "clustering/", file, "_age_cluster.png")

png(filename = file_name, width = 720, height = 720)

par(mar=c(3.1, 0.1, 0.1, 1.1))

dendroplot <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(ddata_x),
            aes(label = label, x = x, y = 0, colour = Age, hjust = 0)) +
  coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
  labs(color = "Age")

dendroplot

dev.off()


# Correlation with Set
file_name <- paste0(dir_plots, "clustering/", file, "_set_cluster.png")

png(filename = file_name, width = 720, height = 720)

par(mar = c(3.1, 0.1, 0.1, 1.1))

dendroplot <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = label(ddata_x),
            aes(label = label, x = x, y = 0, colour = Set, hjust = 0)) +
  coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
  labs(color = "Set")

dendroplot

dev.off()
##-----------------------------------------------------------------------------