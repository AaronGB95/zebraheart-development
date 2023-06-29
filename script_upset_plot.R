##------------------------------------------------------------------------------
##
## Upset plot
##
## Realiza el upset plot de un conjunto de archivos de expresión
##
## Autor: Aarón García Blázquez
##
## Fecha 09/02/2023
##
##------------------------------------------------------------------------------

rm(list = ls())


## Librerías necesarias
##---------------------

library(rstudioapi)
library(ComplexHeatmap)

## Carga de los datos
##-------------------

# Cambiamos el directorio de trabajo a la ubicacion del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

files <- dir(path = "data/studies", pattern = "annot")

# Con los archivos anotados
list_files <- lapply(files, function(x) paste0("data/studies/", x))
lapply(list_files,load,.GlobalEnv)


## Guardamos los datos en una tabla
##---------------------------------

EMTAB10860_mat_annot <- expression_matrix_annot
rm(expression_matrix_annot)

studies <- ls(pattern = "_mat_annot")
studies <- unlist(strsplit(unlist(studies), "_mat_annot"))
studies <- unique(studies)

list_files <- mget(ls(pattern = 'annot$'))

dims <- NULL

for (fi in list_files) {
  dims <- c(dims, nrow(fi))
}

df <- matrix(NA, nrow = max(dims), ncol = length(studies))

for (i in 1:length(list_files)) {
  values <- rownames(list_files[[i]])
  length(values) <- nrow(df)
  df[,i] <- values
}

colnames(df) <- studies

df <- as.data.frame(df)


## Hacemos la matriz de combinación
##---------------------------------

# Creamos una lista a partir de los nombres de cada estudio
lt = list(EMTAB10860 = df$EMTAB10860,
          GSE79585 = df$GSE79585,
          GSE103169 = df$GSE103169,
          GSE107228 = df$GSE107228,
          GSE133130 = df$GSE133130,
          GSE143346 = df$GSE143346,
          GSE152389 = df$GSE152389,
          GSE160107 = df$GSE160107,
          GSE160398 = df$GSE160398,
          GSE189934 = df$GSE189934,
          GSE206948 = df$GSE206948)

# Eliminamos los NAs
lt = lapply(lt, function(x) x[!is.na(x)])

comb_matrix = make_comb_mat(lt)
dim(comb_matrix)


# Filtramos intersecciones
cond1 <- colSums(comb_matrix) == 1
table(cond1)

cond2 <- colSums(comb_matrix) >= (nrow(comb_matrix) - 1)
table(cond2)

mask <- (cond1 | cond2)

m2 <- comb_matrix[,mask]


png(filename = "plots/upset_annotated.png", width = 720, height = 720)
UpSet(m2, pt_size = unit(5, "mm"),
      comb_order = order(comb_size(m2), decreasing = TRUE),
      comb_col = c("firebrick2", "royalblue", "black")[as.factor(comb_degree(m2))],
      top_annotation = upset_top_annotation(m2, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m2, add_numbers = TRUE))

dev.off()

write.csv(df, file = "data/datasets/gene_names.csv", row.names = FALSE)




