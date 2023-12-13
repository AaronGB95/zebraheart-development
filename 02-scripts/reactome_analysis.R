require(ReactomePA)

reactomeAnalysis <- function(diffexp, contrast, organism) {
  if (organism == "zebrafish") {
    geneList <- diffexp$logFC
    names(geneList) <- diffexp$ENTREZID
  } else if (organism == "human") {
    geneList <- diffexp[!is.na(diffexp$EntrezGeneID_Human), "logFC"]
    names(geneList) <- diffexp$EntrezGeneID_Human[!is.na(diffexp$EntrezGeneID_Human)]
  } else if (organism == "mouse") {
    geneList <- diffexp[!is.na(diffexp$EntrezGeneID_Mouse), "logFC"]
    names(geneList) <- diffexp$EntrezGeneID_Mouse[!is.na(diffexp$EntrezGeneID_Mouse)]
    }
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  gsea_reac <- gsePathway(geneList = geneList,
                          pvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          verbose = FALSE,
                          organism = organism,
                          eps = 0)
  
  save(gsea_reac, file = paste0(dir_output,
                                "gene_set_enrichment_analysis/",
                                organism,
                                "_TMM_GSEA_Reactome_",
                                contrast,
                                ".RData"))
}

##------------------------------------------------------------------------------
## 72 hpf vs 48 hpf
##------------------------------------------------------------------------------
load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_72_48.RData"))
reactomeAnalysis(tt_72_48$table, "72_48", "zebrafish")

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_72_48.RData"))
reactomeAnalysis(annotated_data, "72_48", "human")

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_72_48.RData"))
reactomeAnalysis(annotated_data, "72_48", "mouse")

##------------------------------------------------------------------------------
## 120 hpf vs 72 hpf
##------------------------------------------------------------------------------
load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_120_72.RData"))
reactomeAnalysis(tt_120_72$table, "120_72", "zebrafish")

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_120_72.RData"))
reactomeAnalysis(annotated_data, "72_48", "human")

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_120_72.RData"))
reactomeAnalysis(annotated_data, "120_72", "mouse")

##------------------------------------------------------------------------------
## Adult vs 120 hpf
##------------------------------------------------------------------------------
load(paste0(dir_output,
            "differential_expression/",
            "TMM_control_tt_adult_120.RData"))
reactomeAnalysis(tt_adult_120$table, "adult_120", "zebrafish")

load(paste0(dir_output,
            "differential_expression/",
            "human_annotated_tt_adult_120.RData"))
reactomeAnalysis(annotated_data, "adult_120", "human")

load(paste0(dir_output,
            "differential_expression/",
            "mouse_annotated_tt_adult_120.RData"))
reactomeAnalysis(annotated_data, "adult_120", "mouse")