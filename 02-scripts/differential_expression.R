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

## Cargamos los datos necesarios
phenodata <- read.table(paste0(dir_docs, "phenodata_4_groups.txt"),
                        sep = "\t",
                        header = 1)
annotations <- read.table(paste0(dir_docs, "annotations.txt"),
                          sep = "\t",
                          header = 1)

differential_expression <- function(datamatrix,
                                    phenodata,
                                    annotations,
                                    normalization) {

  ## Expresión diferencial
  ##-----------------------

  # Crea la matriz modelo y la matriz nula
  mod <- model.matrix(~ 0 + phenodata$Age, data = phenodata)
  mod0 <- model.matrix(~1, data = phenodata)

  colnames(mod)[1:4] <- c("hpf120", "hpf48", "hpf72", "Adult")

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
  save(tt_72_48, file = paste0(dir_output,
                               "differential_expression/",
                               normalization,
                               "_control_tt_72_48.RData"))

  lrt2 <- glmLRT(fit, contrast = cont_mod[, 2])
  tt_120_72 <- topTags(lrt2, n = Inf, adjust.method = "fdr")
  tt_120_72$table$Gene <- rownames(tt_120_72$table)
  tt_120_72$table <- merge(tt_120_72$table, annotations, by = "Gene")
  save(tt_120_72, file = paste0(dir_output,
                                "differential_expression/",
                                normalization,
                                "_control_tt_120_72.RData"))

  lrt3 <- glmLRT(fit, contrast = cont_mod[, 3])
  tt_adult_120 <- topTags(lrt3, n = Inf, adjust.method = "fdr")
  tt_adult_120$table$Gene <- rownames(tt_adult_120$table)
  tt_adult_120$table <- merge(tt_adult_120$table, annotations, by = "Gene")
  save(tt_adult_120, file = paste0(dir_output,
                                   "differential_expression/",
                                   normalization,
                                   "_control_tt_adult_120.RData"))

  # Función que realiza los volcano plots
  volcanoplot <- function(tt, contrname, normalization) {
    res <- as.data.frame(tt$table)

    keyvals <- ifelse(
      res$logFC < -1.5 & res$PValue < 0.05, "royalblue",
      ifelse(res$logFC > 1.5 & res$PValue < 0.05, "firebrick2", "black")
    )
    keyvals[is.na(keyvals)] <- "black"
    names(keyvals)[keyvals == "firebrick2"] <- "up"
    names(keyvals)[keyvals == "black"] <- "n.s."
    names(keyvals)[keyvals == "royalblue"] <- "down"

    file_name <- paste0(dir_output,
                        "plots/volcano_plots/",
                        normalization,
                        "_control_volcano_",
                        contrname,
                        ".png")

    p <- EnhancedVolcano(res,
                         lab = "Gene",
                         x = "logFC",
                         y = "FDR",
                         subtitle = NULL,
                         pCutoff = 0.05,
                         selectLab = rownames(res)[which(names(keyvals) %in%
                                                           c("up", "down"))],
                         labSize = 3,
                         pointSize = 1.0,
                         title = paste0(contrname, " ", normalization),
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
  volcanoplot(tt_72_48,
              contrname = contrname,
              normalization = normalization)

  contrname <- "120 hpf vs 72 hpf"
  volcanoplot(tt_120_72,
              contrname = contrname,
              normalization = normalization)

  contrname <- "Adult vs 120 hpf"
  volcanoplot(tt_adult_120,
              contrname = contrname,
              normalization = normalization)

}

datamatrix <- read.table(paste0(dir_data,
                                "datatables/",
                                "datamatrix_tmm.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

differential_expression(datamatrix = datamatrix,
                       phenodata = phenodata,
                       annotations = annotations,
                       normalization = "TMM")

datamatrix <- read.table(paste0(dir_data,
                                "datatables/",
                                "datamatrix_qn.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

differential_expression(datamatrix = datamatrix,
                       phenodata = phenodata,
                       annotations = annotations,
                       normalization = "QN")

datamatrix <- read.table(paste0(dir_data,
                                "datatables/",
                                "datamatrix_tmm_ncrna_removed.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

differential_expression(datamatrix = datamatrix,
                        phenodata = phenodata,
                        annotations = annotations,
                        normalization = "TMM_ncRNA")

datamatrix <- read.table(paste0(dir_data,
                                "datatables/",
                                "datamatrix_qn_ncrna_removed.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

differential_expression(datamatrix = datamatrix,
                        phenodata = phenodata,
                        annotations = annotations,
                        normalization = "QN_ncRNA")

