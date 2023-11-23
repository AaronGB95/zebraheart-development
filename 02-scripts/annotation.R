##-----------------------------------------------------------------------------
## Annotation
##-----------------------------------------------------------------------------
datamatrix <- read.table(file = "", header = TRUE)

datamatrix$ENSEMBLID <- rownames(datamatrix)

datamatrix <- merge(datamatrix, annotations, by = "ENSEMBLID")

datamatrix <- datamatrix[!is.na(datamatrix$Gene), ]

rownames(datamatrix) <- datamatrix$Gene

columns_to_remove <- c("ENSEMBLID", "Gene", "ENTREZID")

datamatrix <- datamatrix[, !names(datamatrix) %in% columns_to_remove]

write.table(datamatrix,
            file = paste0("annotated_", args[1]),
            sep = "\t",
            quote = FALSE)
##-----------------------------------------------------------------------------