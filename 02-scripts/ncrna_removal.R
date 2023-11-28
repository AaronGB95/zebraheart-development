##-----------------------------------------------------------------------------
## ncRNA genes removal
##-----------------------------------------------------------------------------

require("rjson")
require("rstudioapi")
require("stringr")
require("dplyr")

genes <- fromJSON(file = paste0(dir_docs, "zebrafish_rrna_genes.json"))

ids <- NULL

for (gene in seq_along(genes)) {
  ids <- c(ids, genes[[gene]]["description"])
}

ids <- unlist(ids)

ids <- str_extract(ids, pattern = "ENSDARG[0-9]+")

#ENSEMBL to gene symbols
annotations <- read.table(paste0(dir_docs, "annotations.txt"),
                          sep = "\t",
                          header = 1)
colnames(annotations) <- c("ENSEMBLID", "Gene", "ENTREZID")
annotations <- annotations[!duplicated(annotations$Gene) &
                             !(annotations$Gene == ""), ]

ncrna_genes <- subset(annotations, ENSEMBLID %in% ids)

##-----------------------------------------------------------------------------
## Apply to expression matrix
##-----------------------------------------------------------------------------

datamatrix <- read.table(paste0(dir_data, "datatables/","datamatrix_tmm.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

datamatrix_filtered <- datamatrix[!(rownames(datamatrix) %in% ncrna_genes$Gene), ]

write.table(datamatrix_filtered,
            file = paste0(dir_data,
                          "datatables/",
                          "datamatrix_tmm_ncrna_removed.txt"),
            sep = "\t",
            quote = FALSE)

datamatrix <- read.table(paste0(dir_data, "datatables/","datamatrix_qn.txt"),
                         sep = "\t",
                         row.names = 1,
                         header = 1)

datamatrix_filtered <- datamatrix[!(rownames(datamatrix) %in% ncrna_genes$Gene), ]

write.table(datamatrix_filtered,
            file = paste0(dir_data,
                          "datatables/",
                          "datamatrix_qn_ncrna_removed.txt"),
            sep = "\t",
            quote = FALSE)
