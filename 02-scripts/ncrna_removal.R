##-----------------------------------------------------------------------------
## ncRNA genes removal
##-----------------------------------------------------------------------------

rm(list = ls())

require("rjson")
require("rstudioapi")
require("stringr")
require("dplyr")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

genes <- fromJSON(file = "zebrafish_rrna_genes.json")

ids <- NULL

for (gene in seq_along(genes)) {
  ids <- c(ids, genes[[gene]]["description"])
}

ids <- unlist(ids)

ids <- str_extract(ids, pattern = "ENSDARG[0-9]+")

#ENSEMBL to gene symbols
annotations <- read.table("annotations.txt",
                          sep = "\t",
                          header = 1)
colnames(annotations) <- c("ENSEMBLID", "Gene", "ENTREZID")
annotations <- annotations[!duplicated(annotations$Gene) &
                             !(annotations$Gene == ""), ]

ncrna_genes <- subset(annotations, ENSEMBLID %in% ids)

##-----------------------------------------------------------------------------
## Apply to expression matrix
##-----------------------------------------------------------------------------

datamatrix <- read.table("datamatrix_tmm.txt",
                         sep = "\t",
                         row.names = 1,
                         header = 1)

datamatrix_filtered <- datamatrix[!ids, ]
