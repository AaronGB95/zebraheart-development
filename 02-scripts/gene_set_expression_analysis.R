
library(mdgsa)
library(clusterProfiler)
library(enrichplot)
library(fgsea)
library(org.Dr.eg.db)
library(ggplot2)
library(pathview)
library(KEGGREST)
library(DT)


# Contrastes

drerio_go <- read.table(file = paste0(dir_docs, "zebrafish_go_terms.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_bp <- read.table(file = paste0(dir_docs, "zebrafish_go_bp.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_cc <- read.table(file = paste0(dir_docs, "zebrafish_go_cc.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_mf <- read.table(file = paste0(dir_docs, "zebrafish_go_mf.txt"),
                        sep = "\t",
                        header = TRUE)

drerio_kegg <- read.table(file = paste0(dir_docs, "zebrafish_kegg_terms.txt"),
                          sep = "\t",
                          header = TRUE)

## 72 hpf vs 48 hpf

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_72_48.RData"))

diffexp <- tt_72_48$table

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_GO_OG <- GSEA(geneList = geneList,
                   TERM2GENE = drerio_bp,
                   verbose = T,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pAdjustMethod = "BH",
                   eps = 0)

gsea_GO_DB <- gseGO(geneList = geneList,
                    ont = "BP",
                    OrgDb = org.Dr.eg.db,
                    keyType = 'SYMBOL',
                    verbose = T,
                    minGSSize = 10,
                    maxGSSize = 500,
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH")

save(gsea_BH, file = paste0(dir_output,
                            "gene_set_enrichment_analysis/",
                            "TMM_GSEA_GOBP_72_48.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_BH@result$pvalue < 0.05))))

dotplot(gsea_BH, showCategory = 20, font.size = 8) +
  ggtitle("GSEA GO:BP 72 hpf vs 48 hpf")

### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg = diffexp$logFC
names(geneList_kegg) = diffexp$ENTREZID
geneList_kegg = sort(geneList_kegg, decreasing = TRUE)


gsea_kegg = gseKEGG(geneList     = geneList_kegg,
                    organism     = 'dre',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.05,
                    verbose      = TRUE)

save(gsea_kegg, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "TMM_GSEA_KEGG_72_48.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_kegg@result$pvalue < 0.05))))

dotplot(gsea_kegg, showCategory = 20, font.size = 8) +
  ggtitle("GSEA KEGG 72 hpf vs 48 hpf")


## 120 hpf vs 72 hpf

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_120_72.RData"))

diffexp <- tt_120_72$table

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_BH <- gseGO(geneList  = geneList,
                 ont       = "BP",
                 OrgDb     = org.Dr.eg.db,
                 keyType   = 'SYMBOL',
                 verbose   = T,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH")

save(gsea_BH, file = paste0(dir_output,
                            "gene_set_enrichment_analysis/",
                            "TMM_GSEA_GOBP_120_72.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_BH@result$pvalue < 0.05))))

dotplot(gsea_BH, showCategory = 20, font.size = 8) +
  ggtitle("GSEA GO:BP 120 hpf vs 72 hpf")


### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg = diffexp$logFC
names(geneList_kegg) = diffexp$ENTREZID
geneList_kegg = sort(geneList_kegg, decreasing = TRUE)

gsea_kegg = gseKEGG(geneList     = geneList_kegg,
                    organism     = 'dre',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.05,
                    verbose      = TRUE)

save(gsea_kegg, file = paste0(dir_output,
                              "gene_set_enrichment_analysis/",
                              "TMM_GSEA_KEGG_120_72.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_kegg@result$pvalue < 0.05))))

dotplot(gsea_kegg, showCategory = 20, font.size = 8) +
  ggtitle("GSEA KEGG 120 hpf vs 72 hpf")


## Adult vs 120 hpf

load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_adult_120.RData"))

diffexp <- tt_adult_120$table

### Gene Ontology

geneList <- diffexp$logFC
names(geneList) <- diffexp$Gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_BH <- gseGO(geneList  = geneList,
                 ont       = "BP",
                 OrgDb     = org.Dr.eg.db,
                 keyType   = 'SYMBOL',
                 verbose   = T,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH")

save(gsea_BH, file = paste0(dir_output,
                            "gene_set_enrichment_analysis",
                            "TMM_GSEA_GOBP_adult_120.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_BH@result$pvalue < 0.05))))

dotplot(gsea_BH, showCategory = 20, font.size = 8) +
  ggtitle("GSEA GO:BP Adult vs 120 hpf")


### KEGG

## Ranked list. Id must be Entrez id
geneList_kegg = diffexp$logFC
names(geneList_kegg) = diffexp$ENTREZID
geneList_kegg = sort(geneList_kegg, decreasing = TRUE)

gsea_kegg = gseKEGG(geneList     = geneList_kegg,
                    organism     = 'dre',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.05,
                    verbose      = TRUE)

save(gsea_kegg, file = paste0(dir_output,
                              "gene_set_enrichment_analysis",
                              "TMM_GSEA_KEGG_adult_120.RData"))

print(paste("Número de términos significativos:",
            length(which(gsea_kegg@result$pvalue < 0.05))))

dotplot(gsea_kegg, showCategory = 20, font.size = 8) +
  ggtitle("GSEA KEGG Adult vs 120 hpf")

