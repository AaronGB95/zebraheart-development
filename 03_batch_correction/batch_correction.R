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
datamatrix <- as.matrix(read.table("datamatrix_qn.txt", sep = "\t"))
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
dataset_adjusted <- removeBatchEffect(edata,
                                      batch = svobj$sv[, 1],
                                      batch2 = svobj$sv[, 2],
                                      design = mod)

# Guardamos la tabla con los datos transformados
save(dataset_adjusted, file = "dataset_sva.RData")

# Introduce las variables sustitutas en las matrices modelo y nula
mod_sv <- cbind(mod, svobj$sv)
mod0_sv <- cbind(mod0, svobj$sv)

colnames(mod_sv)[5:5] <- paste0("sv", 1:n_sv)

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
tt_72_48_sva <- left_join(tt_72_48_sva$table %>%
                            mutate(ind = rownames(tt_72_48_sva$table)),
                          counts,
                          by = "ind")
save(tt_72_48_sva, file = "data/differential_expression/tt_72_48_sva.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod_sv[, 2])
tt_120_72_sva <- topTags(lrt2, n = Inf, adjust.method = "fdr")
tt_120_72_sva <- left_join(tt_120_72_sva$table %>%
                             mutate(ind = rownames(tt_120_72_sva$table)),
                           counts,
                           by = "ind")
save(tt_120_72_sva, file = "data/differential_expression/tt_120_72_sva.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod_sv[, 3])
tt_adult_120_sva <- topTags(lrt3, n = Inf, adjust.method = "fdr")
tt_adult_120_sva <- left_join(tt_adult_120_sva$table %>%
                                mutate(ind = rownames(tt_adult_120_sva$table)),
                              counts,
                              by = "ind")
save(tt_adult_120_sva,
     file = "data/differential_expression/tt_Adult_120_sva.RData")


##------------------------------------------------------------------------------
##
## Corrección del efecto lote con ComBat
##
##------------------------------------------------------------------------------


modcombat <- model.matrix(~1, data = phenodata)

combat_data <- ComBat_seq(edata,
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
tt_72_48_combat <- left_join(tt_72_48_combat$table %>%
                               mutate(ind = rownames(tt_72_48_combat$table)),
                             counts,
                             by = "ind")
save(tt_72_48_combat,
     file = "data/differential_expression/tt_72_48_combat.RData")

lrt2 <- glmLRT(fit, contrast = cont_mod[, 2])
tt_120_72_combat <- topTags(lrt2, n = Inf, adjust.method = "fdr")
tt_120_72_combat <- left_join(tt_120_72_combat$table %>%
                                mutate(ind = rownames(tt_120_72_combat$table)),
                              counts,
                              by = "ind")
save(tt_120_72_combat,
     file = "data/differential_expression/tt_120_72_combat.RData")

lrt3 <- glmLRT(fit, contrast = cont_mod[, 3])
tt_adult_120_combat <- topTags(lrt3, n = Inf, adjust.method = "fdr")
tt_adult_120_combat <- left_join(tt_adult_120_combat$table %>%
  mutate(
         ind = rownames(tt_adult_120_combat$table)),
counts,
by = "ind")
save(tt_adult_120_combat,
     file = "data/differential_expression/tt_Adult_120_combat.RData")

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

  res <- as.data.frame(tt)
  res$logFC <- as.numeric(res$logFC)
  res$PValue <- as.numeric(res$PValue)
  rownames(res) <- res$ind

  res$color <- ifelse(
                      res$logFC < -2.5, "royalblue",
                      ifelse(res$logFC > 2.5, "firebrick2", "black"))
  res$color[is.na(res$color)] <- "black"
  res$group[res$color == "firebrick2"] <- "up"
  res$group[res$color == "black"] <- "n.s."
  res$group[res$color == "royalblue"] <- "down"

  file_name <- paste0("plots/volcano_", paste0(contrname, ".png"))

  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = res$logFC,
                       y = res$FDR,
                       subtitle = NULL,
                       pCutoff = 0.05,
                       pointSize = 1.0,
                       title = contrname,
                       titleLabSize = 12,
                       colCustom = color,
                       colAlpha = 1,
                       legendPosition = "right",
                       legendLabSize = 8,
                       legendIconSize = 2.0,
                       axisLabSize = 8)

  ggsave(filename = file_name, plot = p, dpi = 300)

}


## Plots con ComBat

contrname <- "72 hpf vs 48 hpf"
volcanoplot(tt_72_48_combat, contrname = contrname)

contrname <- "120 hpf vs 72 hpf"
volcanoplot(tt_120_72_combat, contrname = contrname)

contrname <- "Adult vs 120 hpf"
volcanoplot(tt_Adult_120_combat, contrname = contrname)


## Plots con SVA

contrname <- "72 hpf vs 48 hpf"
volcanoplot(tt_72_48_sva, contrname = contrname)

contrname <- "120 hpf vs 72 hpf"
volcanoplot(tt_120_72_sva, contrname = contrname)

contrname <- "Adult vs 120 hpf"
volcanoplot(tt_Adult_120_sva, contrname = contrname)
