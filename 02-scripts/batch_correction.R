##------------------------------------------------------------------------------
##
## batch_correction.R
##
## Este script corrige el posible efecto lote en un conjunto de datos
## mediante ComBat o SVA
##
## Author: Aarón García Blázquez
##
## Email: aaron.garcia.blazquez@gmail.com
##
##------------------------------------------------------------------------------

rm(list = ls())

## Required libraries
require(rstudioapi)
require(sva)
require(limma)
require(edgeR)
require(DESeq2)
require(dplyr)
require(ggplot2)
require(EnhancedVolcano)

## Data load
##--------------------

## Change working directory to script path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Data load
datamatrix <- read.table("datamatrix_qn.txt", sep = "\t")
phenodata <- read.table("phenodata.txt",
                        sep = "\t",
                        header = 1,
                        row.names = 1)

# Variable of interest as factor
phenodata$Age <- as.factor(phenodata$Age)

# Create model matrix and null matrix
mod <- as.matrix(model.matrix(~ 0 + phenodata$Age, data = phenodata))
mod0 <- as.matrix(model.matrix(~1, data = phenodata))

colnames(mod) <- c("hpf120", "hpf48", "hpf72", "Adult")

datamatrix <- as.matrix(log2(datamatrix + 1))

##------------------------------------------------------------------------------
##
## Surrogate Variable Analysis (SVA)
##
##------------------------------------------------------------------------------

# Estima el número de variables sustitutas para el modelo
n_sv <- num.sv(datamatrix,
               mod,
               method = "leek")

# Estimamos los coeficientes de 2 variables sustitutas
svobj <- svaseq(datamatrix,
                as.matrix(mod),
                as.matrix(mod0),
                n.sv = n_sv)

# Eliminamos la variabilidad de los datos para el análisis exploratorio
dataset_adjusted <- removeBatchEffect(datamatrix,
                                      batch = svobj$sv[, 1],
                                      batch2 = svobj$sv[, 2],
                                      design = mod)

# Guardamos la tabla con los datos transformados
write.table(dataset_adjusted,
            file = "datamatrix_sva.txt",
            sep = "\t",
            quote = FALSE)

# Introduce las variables sustitutas en las matrices modelo y nula
mod_sv <- cbind(mod, svobj$sv)
mod0_sv <- cbind(mod0, svobj$sv)

colnames(mod_sv)[5:ncol(mod_sv)] <- paste0("sv", 1:n_sv)

# Definimos los contrastes
cont_mod_sv <- makeContrasts(hpf72vshpf48 = hpf72 - hpf48,
                             hpf120vshpf72 = hpf120 - hpf72,
                             Adultvshpf120 = Adult - hpf120,
                             levels = mod_sv)

# Convertimos la matriz de conteos en una lista de expresión digital
dge <- DGEList(counts = datamatrix)

# Estimamos la dispersión de los datos
dge <- estimateDisp(dge, design = mod_sv)

plotBCV(dge)

# Ajustamos el modelo con las variables sustitutas
fit <- glmFit(dge, design = mod_sv)


# Realizamos el análisis de expresión diferencial
lrt1 <- glmLRT(fit, contrast = cont_mod_sv[, 1])
tt_72_48_sva <- topTags(lrt1, n = Inf, adjust.method = "fdr")

save(tt_72_48_sva, file = "tt_72_48_sva.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod_sv[, 2])
tt_120_72_sva <- topTags(lrt2, n = Inf, adjust.method = "fdr")

save(tt_120_72_sva, file = "tt_120_72_sva.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod_sv[, 3])
tt_adult_120_sva <- topTags(lrt3, n = Inf, adjust.method = "fdr")

save(tt_adult_120_sva, file = "tt_Adult_120_sva.RData")


##------------------------------------------------------------------------------
##
## Corrección del efecto lote con ComBat
##
##------------------------------------------------------------------------------


modcombat <- model.matrix(~1, data = phenodata)

combat_data <- ComBat_seq(datamatrix,
                          batch = phenodata$Set,
                          group = NULL,
                          full_mod = FALSE)

mod <- model.matrix(~0 + phenodata$Age, data = phenodata)
colnames(mod) <- c("hpf120", "hpf48", "hpf72", "Adult")

dge <- DGEList(counts = combat_data)

dge <- estimateDisp(dge, design = mod)

plotBCV(dge)

fit <- glmFit(dge, design = mod)

cont_mod <- makeContrasts(hpf72vshpf48 = hpf72 - hpf48,
                          hpf120vshpf72 = hpf120 - hpf72,
                          Adultvshpf120 = Adult - hpf120,
                          levels = mod)

lrt1 <- glmLRT(fit, contrast = cont_mod[, 1])
tt_72_48_combat <- topTags(lrt1, n = Inf, adjust.method = "fdr")

save(tt_72_48_combat,
     file = "tt_72_48_combat.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod[, 2])
tt_120_72_combat <- topTags(lrt2, n = Inf, adjust.method = "fdr")

save(tt_120_72_combat,
     file = "tt_120_72_combat.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod[, 3])
tt_adult_120_combat <- topTags(lrt3, n = Inf, adjust.method = "fdr")

save(tt_adult_120_combat,
     file = "tt_Adult_120_combat.RData")

write.table(combat_data,
            file = "datamatrix_combat.txt",
            sep = "\t",
            quote = FALSE)


##------------------------------------------------------------------------------
##
## Volcano plots
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
  file_name <- paste0("volcano_", paste0(contrname, ".png"))
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


## Plots con ComBat

contrname <- "72 hpf vs 48 hpf combat"
volcanoplot(tt_72_48_combat, contrname = contrname)

contrname <- "120 hpf vs 72 hpf combat"
volcanoplot(tt_120_72_combat, contrname = contrname)

contrname <- "Adult vs 120 hpf combat"
volcanoplot(tt_adult_120_combat, contrname = contrname)


## Plots con SVA

contrname <- "72 hpf vs 48 hpf sva"
volcanoplot(tt_72_48_sva, contrname = contrname)

contrname <- "120 hpf vs 72 hpf sva"
volcanoplot(tt_120_72_sva, contrname = contrname)

contrname <- "Adult vs 120 hpf sva"
volcanoplot(tt_adult_120_sva, contrname = contrname)
