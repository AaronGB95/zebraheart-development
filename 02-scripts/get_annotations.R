
# Required packages
library(biomaRt)
library(org.Dr.eg.db)
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
