
overrepresentation <- function(sig_genes, ontology, contrast, group) {
  if (ontology == "KEGG") {
    kegg_item <- enricher(gene = sig_genes,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          TERM2GENE = r_kegg,
                          TERM2NAME = dre_pathways)
    save(kegg_item,
         file = paste0(dir_outputs,
                       "over_representation_analysis/",
                       ontology,
                       "_",
                       group,
                       "_",
                       contrast,
                       ".RData"))
    
    if (dim(kegg_item)[1] > 0) {
      kegg_table <- kegg_item@result[kegg_item@result$p.adjust < 0.05,
                                     c("ID","p.adjust")]
      write.table(kegg_table,
                  file = paste0(dir_outputs,
                                "over_representation_analysis/",
                                ontology,
                                "_",
                                group,
                                "_",
                                contrast,
                                ".txt"),
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
      
      png(filename = paste0(dir_outputs,
                            "plots/ora/",
                            ontology,
                            "_",
                            group,
                            "_",
                            contrast,
                            ".png"))
      dotplot(keggTop,
              showCategory = 20,
              x = "geneRatio",
              font.size = 6)
      dev.off()
    }
  } else {
    ora_item <- enrichGO(gene = bottom_genes,
                         OrgDb = org.Dr.eg.db,
                         keyType = 'SYMBOL',
                         ont = ontology,
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)
    save(ora_item,
         file = paste0(dir_outputs,
                       "over_representation_analysis/",
                       ontology,
                       "_",
                       group,
                       "_",
                       contrast,
                       ".RData"))
    
    if (dim(ora_item)[1] > 0) {
      ora_table <- ora_item@result[ora_item@result$p.adjust < 0.05,
                                     c("ID","p.adjust")]
      write.table(ora_table,
                  file = paste0(dir_outputs,
                                "over_representation_analysis/",
                                ontology,
                                "_",
                                group,
                                "_",
                                contrast,
                                ".txt"),
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
      
      png(filename = paste0(dir_outputs,
                            "plots/ora/",
                            ontology,
                            "_",
                            group,
                            "_",
                            contrast,
                            ".png"))
      dotplot(ora_item,
              showCategory = 20,
              x = "geneRatio",
              font.size = 6)
      dev.off()
    }    
  }
}

threshold <- 1.5

ontologies <- c("BP", "MF", "CC", "KEGG")

if (length(ls(pattern = "tt_72_48")) == 0) {
  load(paste0(dir_outputs, "differential_expression/TMM_", "tt_72_48.RData"))
}

diffexp <- tt_72_48$table

sig.genes <- diffexp[diffexp$FDR < 0.05, ]

# Top genes
top <- sig.genes[sig.genes$logFC > threshold, ]
top_genes <- row.names(top)

# Bottom genes
bottom <- sig.genes[sig.genes$logFC < -threshold, ]
bottom_genes <- row.names(bottom)

sig_genes <- list(top_genes, bottom_genes)
names(sig_genes) <- c("top", "bottom")

for (ontology in ontologies) {
  for (group in names(sig_genes)) {
    overrepresentation(sig_genes = sig_genes[[group]],
                       ontology = ontology,
                       contrast = "72vs48")
  }
  
}

