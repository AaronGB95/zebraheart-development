load(paste0(dir_output,
            "gene_set_enrichment_analysis/",
            "human_TMM_GSEA_DB_GOBP_adult_120.RData"))

sig_terms <- gsea_GO_DB@result[gsea_GO_DB@result$p.adjust < 0.05, ]

dotplot(gsea_GO_DB,
        showCategory = 20,
        split = ".sign",
        font.size = 10,
        label_format = function(x)
          stringr::str_wrap(x, width = 100)) +
  facet_grid(.~.sign)

heart_development <- read.table(file = paste0(dir_docs,
                                              "heart_development.txt"))

cardiac_contraction <- read.table(file = paste0(dir_docs,
                                                "cardiac_contraction.txt"))

calcium_transport <- read.table(file = paste0(dir_docs,
                                              "calcium_transport.txt"))

dev_terms <- sig_terms$ID %in% heart_development$V1

table(dev_terms)

dev_terms <- sig_terms$ID %in% cardiac_contraction$V1

table(dev_terms)

dev_terms <- sig_terms$ID %in% calcium_transport$V1

table(dev_terms)




