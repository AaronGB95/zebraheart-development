##-----------------------------------------------------------------------------
##
## runall.R
##
## script that executes the whole project
##
##-----------------------------------------------------------------------------

dir_docs <- "01-documentation/"
dir_src <- "02-scripts/"
dir_data <- "03-data/"
dir_outputs <- "04-output/"
dir_reports <- "05-reports/"

##-----------------------------------------------------------------------------
## Required packages
##-----------------------------------------------------------------------------
require(rstudioapi)
require(stringr)
require(tidyverse)
require(edgeR)
require(preprocessCore)
require(biomaRt)
require(data.table)
require(knitr)
require(dplyr)
require(ggplot2)
require(ggpubr)
require(ggdendro)
require(PCAtools)
require(tidyr)
require(tibble)
require(RColorBrewer)
require(org.Dr.eg.db)
require(org.Hs.eg.db)
require(org.Mm.eg.db)
require(clusterProfiler)
require(enrichplot)
require(fgsea)
require(mdgsa)
require(pathview)
require(KEGGREST)
require(limma)
require(DT)
require(sva)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Get gene annotations
##-----------------------------------------------------------------------------
source(paste0(dir_src, "get_annotations.R"))
##-----------------------------------------------------------------------------

##-----------------------------------------------------------------------------
## Create Expression Matrix
##-----------------------------------------------------------------------------
source(paste0(dir_src, "datatable_construction.R"))

##-----------------------------------------------------------------------------
## Batch correction
##-----------------------------------------------------------------------------
source(paste0(dir_src, "batch_correction.R"))

##-----------------------------------------------------------------------------
## Exploratory Analysis
##-----------------------------------------------------------------------------
source(paste0(dir_src, "exploratory_analysis.R"))

##-----------------------------------------------------------------------------
## Differential Expression
##-----------------------------------------------------------------------------
source(paste0(dir_src, "differential_expression.R"))

##-----------------------------------------------------------------------------
## Over Representation Analysis
##-----------------------------------------------------------------------------

# Get KEGG pathways descriptions
keys <- keys(org.Dr.eg.db, keytype="SYMBOL")
r_kegg <- AnnotationDbi::select(org.Dr.eg.db, keys=keys,
                                columns=c("SYMBOL", "PATH"),
                                keytype="SYMBOL")
r_kegg <- na.exclude(r_kegg)
r_kegg <- r_kegg[,c(2,1)]

dre_pathways <- keggList("pathway", "dre") %>% 
  tibble(pathway = names(.), description = .)

path_ids <- str_split_i(dre_pathways$pathway, "dre", 2)
path_descriptions <- str_split_i(dre_pathways$description, " -", 1)
dre_pathways$pathway <- path_ids
dre_pathways$description <- path_descriptions

source(paste0(dir_src, "ora.R"))

##-----------------------------------------------------------------------------
## Gene Set Enrichment Analysis
##-----------------------------------------------------------------------------
source(paste0(dir_src, "gsea.R"))

##-----------------------------------------------------------------------------
## Get orthologs
##-----------------------------------------------------------------------------

# Mouse
source(paste0(dir_src, "mouse_annotation.R"))

# Human
source(paste0(dir_src, "human_annotation.R"))







