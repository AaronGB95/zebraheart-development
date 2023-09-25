

rm(list = ls())

require(rstudioapi)
require(dplyr)
require(edgeR)
require(stringr)
require(preprocessCore)

# Establish document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Data load


file_list <- list.files(path = "input/", pattern = "annot.txt$", full.names = TRUE)

# Loop through each file and read it into a data frame
for (i in 1:length(file_list)) {
  load(file_list[i])
}

file_list <- list.files("data/phenodata", pattern = "\\.RData$", full.names = TRUE)

# Loop through each file and read it into a data frame
for (i in 1:length(file_list)) {
  load(file_list[i])
}

# Common genes intersection


# we select expression matrix
matrices <- ls(pattern = "mat")

# Get gene names
nombres <- lapply(matrices, function(x) rownames(get(x)))

# Intersection of all names
genes <- Reduce(unique, nombres)


# Data normalization

tmm <- function(datamatrix) {
  
  # Change NAs to 0s
  datamatrix[is.na(datamatrix)] <- 0
  
  # Create a DGEList object
  dge_list <- DGEList(counts = datamatrix)

  # Calculate library sizes and normalize for RNA composition
  dge_list <- calcNormFactors(dge_list, method = "TMM")

  # Get the normalized counts
  norm_counts <- cpm(dge_list, log = FALSE)
  
  # Return the normalized counts
  return(norm_counts)
}


# Relevant Sample Selection {.tabset}

## E-MTAB-10860


# Check study variables
head(pData)

# Interest variable is group
# We choose "none" samples by group variable
selection <- pData$group == "none"
selected_samples <- rownames(pData[selection,])

# Table filter using sample selection and gene names list
expression_matrix_annot <- as.data.frame(expression_matrix_annot)

# TMM normalization
expression_matrix_annot <- tmm(expression_matrix_annot)

selected_data <- expression_matrix_annot[genes,selected_samples]



# Build new phenotypic data table
selected_pdata <- as.data.frame(matrix(NA, nrow = length(selected_samples), ncol = 3))
rownames(selected_pdata) <- selected_samples
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "E-MTAB-10860"
selected_pdata$Age <- "120 hpf"
selected_pdata$Norm <- "FPKM"

new_pData <- selected_pdata


## GSE79585


# Data load
load("data/annot_norm/GSE79585_tmm.RData")

# Check study variables
head(GSE79585_pData)

# All samples are relevant, we only filter by gene names list
GSE79585_sel <- GSE79585_mat_annot[genes,]

# Estos datos ya están normalizados mediante TMM, así que no los modificamos

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE79585_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE79585_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE79585_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE79585"
selected_pdata$Age <- "48-56 hpf"
selected_pdata$Norm <- "TMM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE103169


# Vemos las variables de este estudio
head(GSE103169_pData)

# En este caso nos interesan los wild-type, así que seleccionamos ese grupo
selection <- GSE103169_pData$group == "wt"
selected_samples <- rownames(GSE103169_pData[selection,])

# Filtramos las muestras y los genes
GSE103169_sel <- GSE103169_mat_raw[genes,selected_samples]

# Normalizamos mediante TMM
GSE103169_sel <- tmm(GSE103169_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE103169_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE103169_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE103169_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE103169"
selected_pdata$Age <- "48-56 hpf"
selected_pdata$Norm <- "TC"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE107228


'# Vemos las variables de este estudio
head(GSE107228_pData)

# En este caso nos interesan los wild-type, así que seleccionamos ese grupo
selection <- GSE107228_pData$group == "Control"
selected_samples <- rownames(GSE107228_pData[selection,])

# Filtramos las muestras y los genes
GSE107228_sel <- GSE107228_mat_annot[genes,selected_samples]

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE107228_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]'



'# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE107228_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE107228_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE107228"
selected_pdata$Age <- "Adult"
selected_pdata$Norm <- "FPKM"'



'new_pData <- bind_rows(new_pData, selected_pdata)'


## GSE133130


# Vemos las variables de este estudio
head(GSE133130_pData)

# En este caso nos interesan los wild-type, así que seleccionamos ese grupo
selection <- GSE133130_pData$group == "wt"
selected_samples <- rownames(GSE133130_pData[selection,])

# Filtramos las muestras y los genes
GSE133130_mat_annot <- as.data.frame(GSE133130_mat_annot)
GSE133130_sel <- GSE133130_mat_annot[genes,selected_samples]

# Normalizamos mediante TMM
GSE133130_sel <- tmm(GSE133130_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE133130_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE133130_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE133130_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE133130"
selected_pdata$Age <- "Adult"
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE143346


# Vemos las variables de este estudio
head(GSE143346_pData)

# En este caso nos interesan todas las muestras, así que solo filtramos por genes
GSE143346_sel <- GSE143346_mat_raw[genes,]

# Normalizamos mediante TMM
GSE143346_sel <- tmm(GSE143346_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE143346_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE143346_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE143346_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE143346"
selected_pdata$Age <- "Adult"
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE152389


# Vemos las variables de este estudio
head(GSE152389_pData)

# En este caso nos interesan los controles, así que seleccionamos ese grupo
selection <- GSE152389_pData$group == "Control"
selected_samples <- rownames(GSE152389_pData[selection,])

# Filtramos las muestras y los genes
GSE152389_sel <- GSE152389_mat_annot[genes,selected_samples]

# Normalizamos mediante TMM
GSE152389_sel <- tmm(GSE152389_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE152389_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE152389_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE152389_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE152389"
selected_pdata$Age <- "120 hpf"
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE160107


# Vemos las variables de este estudio
head(GSE160107_pData)

# En este caso nos interesan las muestras GFP negativas, así que seleccionamos ese grupo
selection <- str_detect(GSE160107_pData$group, "GFP-.")
selected_samples <- rownames(GSE160107_pData[selection,])

# Filtramos las muestras y los genes
GSE160107_mat_annot <- as.data.frame(GSE160107_mat_annot)
GSE160107_sel <- GSE160107_mat_annot[genes,selected_samples]

# Normalizamos mediante TMM
GSE160107_sel <- tmm(GSE160107_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE160107_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE160107_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE160107_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE160107"
selected_pdata$Age <- c("48-56 hpf", "48-56 hpf", "48-56 hpf", "72 hpf", "72 hpf", "72 hpf")
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE160398


# Vemos las variables de este estudio
head(GSE160398_pData)

# En este caso nos interesan las muestras negativas, así que seleccionamos ese grupo
selection <- GSE160398_pData$group == "negative"
selected_samples <- rownames(GSE160398_pData[selection,])

# Filtramos las muestras y los genes
GSE160398_sel <- GSE160398_mat_raw[genes,selected_samples]

# Normalizamos mediante TMM
GSE160398_sel <- tmm(GSE160398_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE160398_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE160398_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE160398_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE160398"
selected_pdata$Age <- "72 hpf"
selected_pdata$Norm <- "TPM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE189934


# Vemos las variables de este estudio
head(GSE189934_pData)

# En este caso nos interesan las muestras control, así que seleccionamos ese grupo
selection <- GSE189934_pData$group == "Control"
selected_samples <- rownames(GSE189934_pData[selection,])

# Filtramos las muestras y los genes
GSE189934_sel <- GSE189934_mat_raw[genes,selected_samples]

# Normalizamos mediante TMM
GSE189934_sel <- tmm(GSE189934_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE189934_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE189934_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE189934_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE189934"
selected_pdata$Age <- "48-56 hpf"
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


## GSE206948


# Vemos las variables de este estudio
head(GSE206948_pData)

# En este caso nos interesan las muestras control, así que seleccionamos ese grupo
selection <- GSE206948_pData$group == "Control"
selected_samples <- rownames(GSE206948_pData[selection,])

# Filtramos las muestras y los genes
GSE206948_sel <- GSE206948_mat_raw[genes,selected_samples]

# Normalizamos mediante TMM
GSE206948_sel <- tmm(GSE206948_sel)

# Lo añadimos al dataset nuevo
selected_data <- merge(selected_data, GSE206948_sel, by = "row.names", all = TRUE)
rownames(selected_data)=selected_data$Row.names
selected_data=selected_data[,2:ncol(selected_data)]



# Construimos la tabla de datos fenotípicos nueva
selected_pdata <- as.data.frame(matrix(NA, nrow = ncol(GSE206948_sel), ncol = 3))
rownames(selected_pdata) <- colnames(GSE206948_sel)
colnames(selected_pdata) <- c("Set", "Age", "Norm")
selected_pdata$Set <- "GSE206948"
selected_pdata$Age <- "Adult"
selected_pdata$Norm <- "FPKM"



new_pData <- bind_rows(new_pData, selected_pdata)


# {-}

# Normalización por cuantiles


normalized_data <- normalize.quantiles(as.matrix(selected_data), copy = TRUE)
normalized_data <- as.data.frame(normalized_data)

rownames(normalized_data) <- rownames(selected_data)
colnames(normalized_data) <- colnames(selected_data)

normalized_data[is.na(normalized_data)] <- 0

selected_data <- normalized_data

new_names <- unique(unlist(lapply(nombres, unique)))

selected_data <- selected_data[new_names,]
selected_data[is.na(selected_data)] = 0



# Guardado de los datos


save(selected_data, file = "data/datasets/dataset_full_qn.RData")
save(new_pData, file = "data/datasets/dataset_qn_phenodata.RData")


# Conteo de genes por estudio


counts <- stack(table(unlist(nombres)))

save(counts, file = "data/studies_per_gene.RData")











