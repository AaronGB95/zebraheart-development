
library(mdgsa)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(fgsea)
library(dplyr)
library(org.Dr.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(simplifyEnrichment)


# Contrastes

drerio_go <- read.delim(file = paste0(dir_docs, "zebrafish_go_terms.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_bp <- read.delim(file = paste0(dir_docs, "zebrafish_go_bp.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_cc <- read.delim(file = paste0(dir_docs, "zebrafish_go_cc.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_mf <- read.delim(file = paste0(dir_docs, "zebrafish_go_mf.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_kegg <- read.delim(file = paste0(dir_docs, "zebrafish_kegg_terms.txt"),
                          sep = "\t",
                          header = TRUE)

## 72 hpf vs 48 hpf
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_72_48.RData"))

diffexp <- tt_72_48$table

### Gene Ontology
#------------------------------------------------------------------------------

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

### GSEA with own gene universe
gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = drerio_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

### Add GO term description
matching_terms <- match(gsea_result@result$ID, drerio_go$go_id)
gsea_result@result$Description <- drerio_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "zebrafish_TMM_GSEA_OG_GOBP_72_48.RData"))

### Dotplot for all significant terms
gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with own genes. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Upregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_72_48_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Downregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_72_48_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Dr.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_TMM_GSEA_OG_GOBP_72_48.RData"))

### GSEA with annotation package
#------------------------------------------------------------------------------
gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Dr.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

### Plot for all terms
gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with database package. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_DB_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "zebrafish_TMM_GSEA_DB_GOBP_72_48.RData"))

### KEGG
#------------------------------------------------------------------------------
## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$ENTREZID
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'dre',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for KEGG. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_gsea_result_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "zebrafish_TMM_gsea_result_72_48.RData"))

## 120 hpf vs 72 hpf

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_120_72.RData"))

diffexp <- tt_120_72$table

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = drerio_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, drerio_go$go_id)

gsea_result@result$Description <- drerio_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "zebrafish_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with own genes. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Upregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_120_72_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Downregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_120_72_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Dr.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Dr.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with database package. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_DB_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "zebrafish_TMM_GSEA_DB_GOBP_120_72.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$ENTREZID
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'dre',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for KEGG. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_gsea_result_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "zebrafish_TMM_gsea_result_120_72.RData"))

## Adult vs 120 hpf

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_adult_120.RData"))

diffexp <- tt_adult_120$table

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = drerio_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, drerio_go$go_id)

gsea_result@result$Description <- drerio_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "zebrafish_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with own genes. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Upregulated terms. Adults vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_adult_120_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with own genes, Downregulated terms. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_TMM_GSEA_OG_GOBP_adult_120_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Dr.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Dr.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with database package. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_GSEA_DB_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "zebrafish_TMM_GSEA_DB_GOBP_adult_120.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$ENTREZID
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'dre',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for KEGG. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_TMM_gsea_result_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "zebrafish_TMM_gsea_result_adult_120.RData"))

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Human Orthologs

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Contrastes

hsapiens_go <- read.delim(file = paste0(dir_docs, "human_go_terms.txt"),
                          sep = "\t",
                          header = TRUE)

hsapiens_bp <- read.delim(file = paste0(dir_docs, "human_go_bp.txt"),
                          sep = "\t",
                          header = TRUE)

hsapiens_cc <- read.delim(file = paste0(dir_docs, "human_go_cc.txt"),
                          sep = "\t",
                          header = TRUE)

hsapiens_mf <- read.delim(file = paste0(dir_docs, "human_go_mf.txt"),
                          sep = "\t",
                          header = TRUE)

hsapiens_kegg <- read.delim(file = paste0(dir_docs, "human_kegg_terms.txt"),
                            sep = "\t",
                            header = TRUE)

## 72 hpf vs 48 hpf

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_72_48.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$HGNC_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = hsapiens_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, hsapiens_go$go_id)

gsea_result@result$Description <- hsapiens_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_OG_GOBP_72_48.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human own genes. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Upregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_72_48_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Downregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_72_48_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Hs.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_human_TMM_GSEA_OG_GOBP_72_48.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human database package. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_DB_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_DB_GOBP_72_48.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Human
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'hsa',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for human KEGG pathways. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_gsea_result_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "human_TMM_gsea_result_72_48.RData"))

## 120 hpf vs 72 hpf

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_120_72.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$HGNC_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = hsapiens_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, hsapiens_go$go_id)

gsea_result@result$Description <- hsapiens_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human own genes. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Upregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_120_72_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Downregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_120_72_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Hs.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_human_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human database package. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_DB_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_DB_GOBP_120_72.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Human
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'hsa',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for human KEGG pathways. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_gsea_result_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "human_TMM_gsea_result_120_72.RData"))

## Adult vs 120 hpf

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_adult_120.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$HGNC_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = hsapiens_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, hsapiens_go$go_id)

gsea_result@result$Description <- hsapiens_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human own genes. Adults vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Upregulated terms. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_adult_120_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with human own genes, Downregulated terms. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_human_TMM_GSEA_OG_GOBP_adult_120_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Hs.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_human_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Hs.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       font.size = 10,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with human database package. Adults vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_GSEA_DB_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "human_TMM_GSEA_DB_GOBP_adult_120.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Human
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)

gsea_result <- gseKEGG(geneList     = geneList_kegg,
                     organism     = 'hsa',
                     keyType      = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose      = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for human KEGG pathways. Adults vs 120 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_human_TMM_gsea_result_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "human_TMM_gsea_result_adult_120.RData"))

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Mouse Orthologs

##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

# Contrastes

mmusculus_go <- read.delim(file = paste0(dir_docs, "mouse_go_terms.txt"),
                           sep = "\t",
                           header = TRUE)

mmusculus_bp <- read.delim(file = paste0(dir_docs, "mouse_go_bp.txt"),
                           sep = "\t",
                           header = TRUE)

mmusculus_cc <- read.delim(file = paste0(dir_docs, "mouse_go_cc.txt"),
                           sep = "\t",
                           header = TRUE)

mmusculus_mf <- read.delim(file = paste0(dir_docs, "mouse_go_mf.txt"),
                           sep = "\t",
                           header = TRUE)

mmusculus_kegg <- read.delim(file = paste0(dir_docs, "mouse_kegg_terms.txt"),
                             sep = "\t",
                             header = TRUE)

## 72 hpf vs 48 hpf

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_72_48.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$MGI_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = mmusculus_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, mmusculus_go$go_id)

gsea_result@result$Description <- mmusculus_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_OG_GOBP_72_48.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse own genes. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Upregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_72_48_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Downregulated terms. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_72_48_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Mm.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_mouse_TMM_GSEA_OG_GOBP_72_48.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse database package. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_DB_GOBP_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_DB_GOBP_72_48.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Mouse
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'mmu',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for mouse KEGG pathways. 72 hpf vs 48 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_gsea_result_72_48.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "mouse_TMM_gsea_result_72_48.RData"))

## 120 hpf vs 72 hpf

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_120_72.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$MGI_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = mmusculus_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, mmusculus_go$go_id)

gsea_result@result$Description <- mmusculus_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse own genes. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Upregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_120_72_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Downregulated terms. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_120_72_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Mm.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_mouse_TMM_GSEA_OG_GOBP_120_72.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse database package. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_DB_GOBP_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_DB_GOBP_120_72.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Mouse
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'mmu',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      font.size = 10,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for mouse KEGG pathways. 120 hpf vs 72 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_gsea_result_120_72.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "mouse_TMM_gsea_result_120_72.RData"))

## Adult vs 120 hpf

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_adult_120.RData"))

diffexp <- annotated_data

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$MGI_Symbol
geneList <- sort(geneList, decreasing = TRUE)

gsea_result <- GSEA(geneList = geneList,
                   TERM2GENE = mmusculus_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 1,
                   eps = 0)

matching_terms <- match(gsea_result@result$ID, mmusculus_go$go_id)

gsea_result@result$Description <- mmusculus_go$name_1006[matching_terms]

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse own genes. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter overexpressed terms
gsea_result_up <- gsea_result
gsea_result_up@result <- gsea_result_up@result[gsea_result_up@result$NES > 0, ]

### Overexpressed terms plot
gsea_result_up_plot <- enrichplot::dotplot(gsea_result_up,
                                          x = "GeneRatio",
                                          color = "p.adjust",
                                          showCategory = 20,
                                          font.size = 10,
                                          orderBy = "x",
                                          label_format = function(x)
                                            stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Upregulated terms. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_up_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_adult_120_up.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Filter underexpressed terms
gsea_result_down <- gsea_result
gsea_result_down@result <- gsea_result_down@result[gsea_result_down@result$NES < 0, ]

### Underexpressed terms plot
gsea_result_down_plot <- enrichplot::dotplot(gsea_result_down,
                                            x = "GeneRatio",
                                            color = "p.adjust",
                                            showCategory = 20,
                                            font.size = 10,
                                            orderBy = "x",
                                            label_format = function(x)
                                              stringr::str_wrap(x, width = 60)) +
  ggtitle("GSEA with mouse own genes, Downregulated terms. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = gsea_result_down_plot,
       filename = "dotplot_mouse_TMM_GSEA_OG_GOBP_adult_120_down.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

### Simplify GSEA result
gsea_result_matrix <- GO_similarity(gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"],
                                   ont = "BP",
                                   db = "org.Mm.eg.db")
gsea_result_sum <- simplifyGO(gsea_result_matrix)
save(gsea_result_sum,
     file = paste0(dir_output,
                   "gene_set_enrichment_analysis/",
                   "sum_mouse_TMM_GSEA_OG_GOBP_adult_120.RData"))

gsea_result <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Mm.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH",
                    eps = 0)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                       x = "GeneRatio",
                                       color = "p.adjust",
                                       showCategory = 20,
                                       orderBy = "x",
                                       label_format = function(x)
                                         stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for GSEA with mouse database package. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 50))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_GSEA_DB_GOBP_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                               "gene_set_enrichment_analysis/",
                               "mouse_TMM_GSEA_DB_GOBP_adult_120.RData"))

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg <- diffexp$logFC
names(geneList_kegg) <- diffexp$EntrezGeneID_Mouse
geneList_kegg <- sort(geneList_kegg, decreasing = TRUE)


gsea_result <- gseKEGG(geneList = geneList_kegg,
                     organism = 'mmu',
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     verbose = TRUE)

gsea_result_plot <- enrichplot::dotplot(gsea_result,
                                      x = "GeneRatio",
                                      color = "p.adjust",
                                      showCategory = 20,
                                      orderBy = "x",
                                      label_format = function(x)
                                        stringr::str_wrap(x, width = 60)) +
  ggtitle("Dotplot for mouse KEGG pathways. Adult vs 120 hpf") +
  theme(plot.title = element_text(hjust = 1),
        plot.margin = margin(l = 40))

ggsave(plot = gsea_result_plot,
       filename = "dotplot_mouse_TMM_gsea_result_adult_120.jpg",
       path = paste0(dir_output, "plots/gsea/"),
       dpi = 300)

save(gsea_result, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "mouse_TMM_gsea_result_adult_120.RData"))

