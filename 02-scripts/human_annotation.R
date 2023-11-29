##-----------------------------------------------------------------------------
##
## human_annotation.R
##
## Script that takes a differential expression result object and annotates
## the genes with the human orthologs.
##
## Author: Aaron Garcia Blazquez
##
## Date: 08/11/2023
##
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
require(biomaRt)
require(dplyr)
require(curl)
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------

contrasts <- list.files(path = paste0(dir_output, "differential_expression"),
                        all.files = TRUE,
                        pattern = "^TMM_control",
                        recursive = TRUE,
                        full.names = TRUE)

for (file in contrasts) {
  load(file)
}

##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Annotation Function
##-----------------------------------------------------------------------------
convert_danio_to_human <- function(m) {

  dre <- useEnsembl(biomart = "ensembl",
                    dataset = "drerio_gene_ensembl",
                    version = 102)

  hsap <- useEnsembl(biomart = "ensembl",
                     dataset = "hsapiens_gene_ensembl",
                     version = 102)

  danio_ensembl_id_attributes <- c("ensembl_gene_id",
                                   "hsapiens_homolog_perc_id")
  #"mmusculus_homolog_goc_score") #orthology score to mouse from danio side

  danio_zfin_id_attributes <- c("ensembl_gene_id",
                                "entrezgene_id",
                                "zfin_id_symbol",
                                "description")

  hsap_attributes <- c("hgnc_symbol",
                       "ensembl_gene_id",
                       "entrezgene_id")

  # make a converted dataframe with ensemble ids and mouse orthology scores
  genes_conv_ensembl <- getLDS(values = m,
                               filters = "ensembl_gene_id",
                               attributes = danio_ensembl_id_attributes,
                               mart = dre, # use danio mart
                               attributesL = hsap_attributes,
                               martL = hsap, # use human mart
                               uniqueRows = TRUE)


  # make dataframe for dre ensembl id and zfin
  genes_conv_zfin <- getBM(attributes = danio_zfin_id_attributes,
                           filters = "ensembl_gene_id",
                           values = m,
                           mart = dre)


  # merge conversion df with zfin id df
  gene_merge <- merge(genes_conv_ensembl,
                      genes_conv_zfin,
                      by.x = "Gene.stable.ID",
                      by.y = "ensembl_gene_id",
                      all.y = TRUE)

  # rename column names to readable format
  colnames(gene_merge)[colnames(gene_merge) == "Gene.stable.ID"] <- "Ensembl_Danio"
  colnames(gene_merge)[colnames(gene_merge) == "X.id..target.Human.gene.identical.to.query.gene"] <- "Human_homology_percentage"
  colnames(gene_merge)[colnames(gene_merge) == "NCBI.gene..formerly.Entrezgene..ID"] <- "EntrezGeneID_Human"
  colnames(gene_merge)[colnames(gene_merge) == "entrezgene_id"] <- "EntrezGeneID_Danio"
  colnames(gene_merge)[colnames(gene_merge) == "Gene.stable.ID.1"] <- "Ensembl_Human"
  colnames(gene_merge)[colnames(gene_merge) == "HGNC.symbol"] <- "HGNC_Symbol"
  colnames(gene_merge)[colnames(gene_merge) == "zfin_id_symbol"] <- "Zfin_Symbol"
  colnames(gene_merge)[colnames(gene_merge) == "description"] <- "Description"


  # select the top ortholog with max homology percentage
  gene_merge <- gene_merge %>% group_by(Ensembl_Danio) %>%
    arrange(desc(Human_homology_percentage), .by_group = TRUE) %>%
    mutate(num_orthologues = n()) %>%
    distinct(Ensembl_Danio, .keep_all = TRUE) %>%
    mutate(unique_ortho_check = n()) %>%
    arrange(desc(num_orthologues))
  cat("If worked try again switching old host=\
  \"https://www.ensembl.org\" and mirror = \"useast\" ")
  return(gene_merge)
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Apply Annotation Function
##-----------------------------------------------------------------------------

human_annotation <- function(datatable) {

  name <- deparse(substitute(datatable))

  annotations <- convert_danio_to_human(datatable$table$ENSEMBLID)

  annotated_data <- merge(x = annotations,
                          y = datatable$table |>
                            distinct(ENSEMBLID,
                                   .keep_all = TRUE),
                          by.x = "Ensembl_Danio",
                          by.y = "ENSEMBLID",
                          all.x = TRUE)

  save(annotated_data,
       file = paste0(dir_output,
                     "differential_expression/",
                     "human_annotated_",
                     name,
                     ".RData"))
}

##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data save
##-----------------------------------------------------------------------------
human_annotation(tt_72_48)
human_annotation(tt_120_72)
human_annotation(tt_adult_120)
##-----------------------------------------------------------------------------