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
require(xlsx)
require(ggplot2)

## Carga de los datos
##--------------------

## Cambiamos el directorio de trabajo a la ubicacion del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Cargamos los datos necesarios
load("data/datasets/dataset_full_qn.RData")
load("data/datasets/phenodata.RData")
load("data/studies_per_gene.RData")

rownames(counts) <- counts$ind


## Expresión diferencial
##-----------------------

# Crea la matriz modelo y la matriz nula
mod <- model.matrix(~ 0 + new_pData$Age, data = new_pData)
mod0 <- model.matrix(~1, data = new_pData)

colnames(mod) <- c("hpf120", "hpf48", "hpf72", "Adult")

dge <- DGEList(counts = selected_data)

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
tt_72_48$table <- left_join(tt_72_48$table %>%
                              mutate(ind = rownames(tt_72_48$table)),
                            counts,
                            by = "ind")
save(tt_72_48, file = "data/differential_expression/tt_72_48.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod[, 2])
tt_120_72 <- topTags(lrt2, n = Inf, adjust.method = "fdr")
tt_120_72$table <- left_join(tt_120_72$table %>%
                               mutate(ind = rownames(tt_120_72$table)),
                             counts,
                             by = "ind")
save(tt_120_72, file = "data/differential_expression/tt_120_72.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod[, 3])
tt_adult_120 <- topTags(lrt3, n = Inf, adjust.method = "fdr")
tt_adult_120$table <- left_join(tt_adult_120$table %>%
                                  mutate(ind = rownames(tt_adult_120$table)),
                                counts,
                                by = "ind")
save(tt_adult_120, file = "data/differential_expression/tt_adult_120.RData")


## Guardado de los resultados en excel

rownames(tt_72_48$table) <- tt_72_48$table$ind
rownames(tt_120_72$table) <- tt_120_72$table$ind
rownames(tt_adult_120$table) <- tt_adult_120$table$ind

## Todos los genes

write.xlsx2(tt_72_48$table[, c("logFC", "FDR")],
            file = "data/results/72_vs_48_todos.xlsx")
write.xlsx2(tt_120_72$table[, c("logFC", "FDR")],
            file = "data/results/120_vs_72_todos.xlsx")
write.xlsx2(tt_adult_120$table[, c("logFC", "FDR")],
            file = "data/results/Adulto_vs_120_todos.xlsx")

## Ups y Downs por separado

### Downs

write.xlsx2(tt_72_48$table[which(tt_72_48$table$logFC < -1.5 &
                                   tt_72_48$table$FDR < 0.05),
                           c("logFC", "FDR")],
            file = "data/results/72_vs_48_downs.xlsx")
write.xlsx2(tt_120_72$table[which(tt_120_72$table$logFC < -1.5 &
                                    tt_120_72$table$FDR < 0.05),
                            c("logFC", "FDR")],
            file = "data/results/120_vs_72_downs.xlsx")
write.xlsx2(tt_adult_120$table[which(tt_adult_120$table$logFC < -1.5 &
                                       tt_adult_120$table$FDR < 0.05),
                               c("logFC", "FDR")],
            file = "data/results/Adulto_vs_120_downs.xlsx")

### Ups

write.xlsx2(tt_72_48$table[which(tt_72_48$table$logFC > 1.5 &
                                   tt_72_48$table$FDR < 0.05),
                           c("logFC", "FDR")],
            file = "data/results/72_vs_48_ups.xlsx")
write.xlsx2(tt_120_72$table[which(tt_120_72$table$logFC > 1.5 &
                                    tt_120_72$table$FDR < 0.05),
                            c("logFC", "FDR")],
            file = "data/results/120_vs_72_ups.xlsx")
write.xlsx2(tt_adult_120$table[which(tt_adult_120$table$logFC > 1.5 &
                                       tt_adult_120$table$FDR < 0.05),
                               c("logFC", "FDR")],
            file = "data/results/Adulto_vs_120_ups.xlsx")


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

  file_name <- paste0("plots/volcano_", paste0(contrname, ".png"))

  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = "logFC",
                       y = "PValue",
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