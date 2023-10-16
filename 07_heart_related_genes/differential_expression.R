##------------------------------------------------------------------------------
##
## Análisis de expresión diferencial
##
## Este script realiza el análisis de expresión diferencial dándole una
## matriz de datos y una matriz de diseño
##
## Autor: Aarón García Blázquez
##
## Fecha de creación: 17-04-2023
##
## Email: aaron.garcia.blazquez@gmail.com
##
##------------------------------------------------------------------------------

rm(list = ls())

## Librerías necesarias
##----------------------

require(rstudioapi)
require(limma)
require(edgeR)
require(DESeq2)
require(dplyr)
require(sva)
require(EnhancedVolcano)
require(ggplot2)

## Carga de los datos
##--------------------

## Cambiamos el directorio de trabajo a la ubicacion del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Cargamos los datos necesarios
datamatrix <- read.table("datamatrix_filtered.txt",
                         sep = "\t",
                         row.names = 1,
                         header = 1)
phenodata <- read.table("phenodata.txt",
                        sep = "\t",
                        header = 1)
annotations <- read.table("annotations.txt",
                          sep = "\t",
                          header = 1)
colnames(annotations) <- c("ENSEMBLID", "Gene", "ENTREZID")
annotations <- annotations[!duplicated(annotations$Gene) & !(annotations$Gene == ""), ]

## Expresión diferencial
##-----------------------

# Crea la matriz modelo y la matriz nula
mod <- model.matrix(~ 0 + phenodata$Age, data = phenodata)
mod0 <- model.matrix(~1, data = phenodata)

colnames(mod) <- c("hpf120", "hpf48", "hpf72", "Adult")

dge <- DGEList(counts = datamatrix)

dge <- estimateDisp(dge, design = mod)

fit <- glmFit(dge, design = mod)


## Contrastes
##-----------

## Creamos los contrastes
cont_mod <- makeContrasts(hpf72vshpf48 = hpf72 - hpf48,
                          hpf120vshpf72 = hpf120 - hpf72,
                          Adultvshpf120 = Adult - hpf120,
                          levels = mod)

lrt1 <- glmLRT(fit, contrast = cont_mod[, 1])
tt_72_48 <- topTags(lrt1, n = Inf, adjust.method = "fdr")
tt_72_48$table$Gene <- rownames(tt_72_48$table)
tt_72_48$table <- merge(tt_72_48$table, annotations, by = "Gene")
save(tt_72_48, file = "tt_72_48_filtered.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod[, 2])
tt_120_72 <- topTags(lrt2, n = Inf, adjust.method = "fdr")
tt_120_72$table$Gene <- rownames(tt_120_72$table)
tt_120_72$table <- merge(tt_120_72$table, annotations, by = "Gene")
save(tt_120_72, file = "tt_120_72_filtered.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod[, 3])
tt_adult_120 <- topTags(lrt3, n = Inf, adjust.method = "fdr")
tt_adult_120$table$Gene <- rownames(tt_adult_120$table)
tt_adult_120$table <- merge(tt_adult_120$table, annotations, by = "Gene")
save(tt_adult_120, file = "tt_adult_120_filtered.RData")

## Guardado de los resultados en excel

## Todos los genes

write.table(tt_72_48$table,
            file = "72_vs_48_filtered.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_120_72$table,
            file = "120_vs_72_filtered.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_adult_120$table,
            file = "Adulto_vs_120_filtered.txt",
            quote = FALSE,
            sep = "\t")

## Ups y Downs por separado

### Downs

write.table(tt_72_48$table[which(tt_72_48$table$logFC < -1.5 &
                                   tt_72_48$table$FDR < 0.05), ],
            file = "72_vs_48_filtered_downs.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_120_72$table[which(tt_120_72$table$logFC < -1.5 &
                                    tt_120_72$table$FDR < 0.05), ],
            file = "120_vs_72_filtered_downs.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_adult_120$table[which(tt_adult_120$table$logFC < -1.5 &
                                       tt_adult_120$table$FDR < 0.05), ],
            file = "Adulto_vs_120_filtered_downs.txt",
            quote = FALSE,
            sep = "\t")

### Ups

write.table(tt_72_48$table[which(tt_72_48$table$logFC > 1.5 &
                                   tt_72_48$table$FDR < 0.05), ],
            file = "72_vs_48_filtered_ups.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_120_72$table[which(tt_120_72$table$logFC > 1.5 &
                                    tt_120_72$table$FDR < 0.05), ],
            file = "120_vs_72_filtered_ups.txt",
            quote = FALSE,
            sep = "\t")
write.table(tt_adult_120$table[which(tt_adult_120$table$logFC > 1.5 &
                                       tt_adult_120$table$FDR < 0.05), ],
            file = "Adulto_vs_120_filtered_ups.txt",
            quote = FALSE,
            sep = "\t")


##------------------------------------------------------------------------------
##
##  Volcano Plots
##
##------------------------------------------------------------------------------

# Función que realiza los volcano plots
volcanoplot <- function(tt, contrname) {
  res <- as.data.frame(tt$table)

  keyvals <- ifelse(
    res$logFC < -1.5 & res$PValue < 0.05, "royalblue",
    ifelse(res$logFC > 1.5 & res$PValue < 0.05, "firebrick2", "black")
  )
  keyvals[is.na(keyvals)] <- "black"
  names(keyvals)[keyvals == "firebrick2"] <- "up"
  names(keyvals)[keyvals == "black"] <- "n.s."
  names(keyvals)[keyvals == "royalblue"] <- "down"

  file_name <- paste0("volcano_filtered_", paste0(contrname, ".png"))

  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = "logFC",
                       y = "FDR",
                       subtitle = NULL,
                       pCutoff = 0.05,
                       selectLab = rownames(res)[which(names(keyvals) %in%
                                                         c("up", "down"))],
                       labSize = 3,
                       pointSize = 1.0,
                       title = contrname,
                       titleLabSize = 12,
                       colCustom = keyvals,
                       colAlpha = 1,
                       legendPosition = "right",
                       legendLabSize = 8,
                       legendIconSize = 2.0,
                       axisLabSize = 8)

  ggsave(filename = file_name, plot = p, dpi = 300)

}

contrname <- "72 hpf vs 48 hpf"
volcanoplot(tt_72_48, contrname = contrname)

contrname <- "120 hpf vs 72 hpf"
volcanoplot(tt_120_72, contrname = contrname)

contrname <- "Adult vs 120 hpf"
volcanoplot(tt_adult_120, contrname = contrname)