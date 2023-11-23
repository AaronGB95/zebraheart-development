mart <- useEnsembl(biomart = "ensembl",
                   dataset = "drerio_gene_ensembl",
                   version = 102)

drerio_attributes <- c("ensembl_gene_id",
                       "zfin_id_symbol",
                       "entrezgene_id")

annotations <- getBM(attributes = drerio_attributes, 
                     values = TRUE,
                     mart = mart)

colnames(annotations) <- c("ENSEMBLID", "Gene", "ENTREZID")

# Eliminar genes duplicados y sin anotación
annotations <- annotations[!duplicated(annotations$Gene) &
                             !(annotations$Gene == "") &
                             !is.na(annotations$Gene), ]

# Guardar los resultados en un archivo
write.table(annotations,
            file = paste0(dir_docs, "annotations.txt"),
            sep = "\t",
            quote = FALSE)

