
# Required packages
library(biomaRt)
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

keys <- keys(org.Dr.eg.db, keytype = "SYMBOL")

# GO terms
mart <- useMart(biomart = "ensembl",
                dataset = "drerio_gene_ensembl")

drerio_go <- getBM(attributes = c("external_gene_name",
                                  "go_id",
                                  "namespace_1003"),
                   mart = mart)

# All GO Terms
drerio_go <- drerio_go[!(drerio_go$external_gene_name == "" |
                           drerio_go$go_id == ""), ]

write.table(drerio_go[, c("go_id", "external_gene_name")],
            file = paste0(dir_docs, "zebrafish_go_terms.txt"),
            sep = "\t",
            quote = FALSE)

# Biological Process
drerio_bp <- drerio_go[drerio_go$namespace_1003 == "biological_process", ]
drerio_bp <- drerio_bp[, c("go_id", "external_gene_name")]

write.table(drerio_bp,
            file = paste0(dir_docs, "zebrafish_go_bp.txt"),
            sep = "\t",
            quote = FALSE)

# Cellular Component
drerio_cc <- drerio_go[drerio_go$namespace_1003 == "cellular_component", ]
drerio_cc <- drerio_cc[, c("go_id", "external_gene_name")]

write.table(drerio_cc,
            file = paste0(dir_docs, "zebrafish_go_cc.txt"),
            sep = "\t",
            quote = FALSE)

# Molecular Function
drerio_mf <- drerio_go[drerio_go$namespace_1003 == "molecular_function", ]
drerio_mf <- drerio_mf[, c("go_id", "external_gene_name")]

write.table(drerio_mf,
            file = paste0(dir_docs, "zebrafish_go_mf.txt"),
            sep = "\t",
            quote = FALSE)

# KEGG pathways
drerio_kegg <- AnnotationDbi::select(org.Dr.eg.db,
                                     keys = keys,
                                     columns = c("PATH", "SYMBOL"),
                                     keytype = "SYMBOL")

drerio_kegg <- na.exclude(drerio_kegg)

write.table(drerio_kegg,
            file = paste0(dir_docs, "zebrafish_kegg_terms.txt"),
            sep = "\t",
            quote = FALSE)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Human annotations

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")

# GO terms
mart <- useMart(biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl")

hsapiens_go <- getBM(attributes = c("external_gene_name",
                                    "go_id",
                                    "namespace_1003"),
                     mart = mart)

# All GO Terms
hsapiens_go <- hsapiens_go[!(hsapiens_go$external_gene_name == "" |
                             hsapiens_go$go_id == ""), ]

write.table(hsapiens_go[, c("go_id", "external_gene_name")],
            file = paste0(dir_docs, "human_go_terms.txt"),
            sep = "\t",
            quote = FALSE)

# Biological Process
hsapiens_bp <- hsapiens_go[hsapiens_go$namespace_1003 == "biological_process", ]
hsapiens_bp <- hsapiens_bp[, c("go_id", "external_gene_name")]

write.table(hsapiens_go,
            file = paste0(dir_docs, "human_go_bp.txt"),
            sep = "\t",
            quote = FALSE)

# Cellular Component
hsapiens_cc <- hsapiens_go[hsapiens_go$namespace_1003 == "cellular_component", ]
hsapiens_cc <- hsapiens_cc[, c("go_id", "external_gene_name")]

write.table(hsapiens_cc,
            file = paste0(dir_docs, "human_go_cc.txt"),
            sep = "\t",
            quote = FALSE)

# Molecular Function
hsapiens_mf <- hsapiens_go[hsapiens_go$namespace_1003 == "molecular_function", ]
hsapiens_mf <- hsapiens_mf[, c("go_id", "external_gene_name")]

write.table(hsapiens_mf,
            file = paste0(dir_docs, "human_go_mf.txt"),
            sep = "\t",
            quote = FALSE)

# KEGG pathways
hsapiens_kegg <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = keys,
                                       columns = c("PATH", "SYMBOL"),
                                       keytype = "SYMBOL")

hsapiens_kegg <- na.exclude(hsapiens_kegg)

write.table(hsapiens_kegg,
            file = paste0(dir_docs, "human_kegg_terms.txt"),
            sep = "\t",
            quote = FALSE)

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Mouse annotations

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

keys <- keys(org.Mm.eg.db, keytype = "SYMBOL")

# GO terms
mart <- useMart(biomart = "ensembl",
                dataset = "mmusculus_gene_ensembl")

mmusculus_go <- getBM(attributes = c("external_gene_name",
                                     "go_id",
                                     "namespace_1003"),
                      mart = mart)

# All GO Terms
mmusculus_go <- mmusculus_go[!(mmusculus_go$external_gene_name == "" |
                               mmusculus_go$go_id == ""), ]

write.table(mmusculus_go[, c("go_id", "external_gene_name")],
            file = paste0(dir_docs, "mouse_go_terms.txt"),
            sep = "\t",
            quote = FALSE)

# Biological Process
mmusculus_bp <- mmusculus_go[mmusculus_go$namespace_1003 == "biological_process", ]
mmusculus_bp <- mmusculus_bp[, c("go_id", "external_gene_name")]

write.table(mmusculus_go,
            file = paste0(dir_docs, "mouse_go_bp.txt"),
            sep = "\t",
            quote = FALSE)

# Cellular Component
mmusculus_cc <- mmusculus_go[mmusculus_go$namespace_1003 == "cellular_component", ]
mmusculus_cc <- mmusculus_cc[, c("go_id", "external_gene_name")]

write.table(mmusculus_cc,
            file = paste0(dir_docs, "mouse_go_cc.txt"),
            sep = "\t",
            quote = FALSE)

# Molecular Function
mmusculus_mf <- mmusculus_go[mmusculus_go$namespace_1003 == "molecular_function", ]
mmusculus_mf <- mmusculus_mf[, c("go_id", "external_gene_name")]

write.table(mmusculus_mf,
            file = paste0(dir_docs, "mouse_go_mf.txt"),
            sep = "\t",
            quote = FALSE)

# KEGG pathways
mmusculus_kegg <- AnnotationDbi::select(org.Mm.eg.db,
                                        keys = keys,
                                        columns = c("PATH", "SYMBOL"),
                                        keytype = "SYMBOL")

mmusculus_kegg <- na.exclude(mmusculus_kegg)

write.table(mmusculus_kegg,
            file = paste0(dir_docs, "mouse_kegg_terms.txt"),
            sep = "\t",
            quote = FALSE)