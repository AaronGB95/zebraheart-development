##-----------------------------------------------------------------------------
##
## mouse_annotation.R
##
## Script that takes a differential expression result object and annotates
## the genes with the mouse orthologs.
##
## Author: Aaron Garcia Blazquez
##
## Date: 08/11/2023
##
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Setup
##-----------------------------------------------------------------------------
# Clean environment
rm(list = ls())

require(rstudioapi)
require(biomaRt)
require(dplyr)
require(curl)

# Set document path as working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data Load
##-----------------------------------------------------------------------------

contrasts <- list.files(path = ".",
                        all.files = TRUE,
                        pattern = "^tt_",
                        recursive = TRUE,
                        full.names = TRUE)

for (file in contrasts) {
  load(file)
}

##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Annotation Function
##-----------------------------------------------------------------------------
convertDanio2Mouse <- function(m){
  
  dre <- useEnsembl(biomart = 'ensembl', 
                    dataset = 'drerio_gene_ensembl',
                    version = 102)
  
  mmus <- useEnsembl(biomart = 'ensembl', 
                     dataset = 'mmusculus_gene_ensembl',
                     version = 102)
  
  danio_ensembl_id_attributes <-c("ensembl_gene_id",
                                  "mmusculus_homolog_perc_id")#,# attributes to query from danio database (query organism) 
  #"mmusculus_homolog_goc_score") #orthology score to mouse from danio side
  
  danio_zfin_id_attributes <-c("ensembl_gene_id",
                               "entrezgene_id",
                               "zfin_id_symbol",
                               "description") # attributes to convert ensembl id to zfin id 
  
  mmus_attributes <- c("mgi_symbol",
                       "ensembl_gene_id",
                       "entrezgene_id") # attributes that you want from mouse database
  
  # make a converted dataframe with ensemble ids and mouse orthology scores 
  genes_conv_ensembl <- getLDS(values = m ,
                               filters = "ensembl_gene_id",  # query using ensemble ids
                               attributes = danio_ensembl_id_attributes, # query these attributes from danio mart
                               mart = dre, # use danio mart
                               attributesL = mmus_attributes, # convert to attributes of mouse mart 
                               martL = mmus,# use mouse mart
                               uniqueRows=T) 
  
  
  # make dataframe for dre ensembl id and zfin
  genes_conv_zfin <- getBM(attributes = danio_zfin_id_attributes,
                           filters = "ensembl_gene_id",
                           values = m,
                           mart = dre)
  
  
  # merge conversion df with zfin id df
  gene_merge <- merge(genes_conv_ensembl,
                      genes_conv_zfin,
                      by.x= "Gene.stable.ID",
                      by.y="ensembl_gene_id",
                      all.y=TRUE)
  
  # rename column names to readable format
  colnames(gene_merge)[colnames(gene_merge)== "Gene.stable.ID"] <- "Ensembl_Danio"
  colnames(gene_merge)[colnames(gene_merge)== "X.id..target.Mouse.gene.identical.to.query.gene"] <- "Mouse_homology_percentage"
  colnames(gene_merge)[colnames(gene_merge)== "NCBI.gene..formerly.Entrezgene..ID"] <- "EntrezGeneID_Mouse"
  colnames(gene_merge)[colnames(gene_merge)== "entrezgene_id"] <- "EntrezGeneID_Danio"
  colnames(gene_merge)[colnames(gene_merge)== "Gene.stable.ID.1"] <- "Ensembl_Mouse"
  colnames(gene_merge)[colnames(gene_merge)== "MGI.symbol"] <- "MGI_Symbol"
  colnames(gene_merge)[colnames(gene_merge)== "zfin_id_symbol"] <- "Zfin_Symbol"
  colnames(gene_merge)[colnames(gene_merge)== "description"] <- "Description"
  
  
  # select the top ortholog with max homology percentage 
  
  gene_merge <- gene_merge %>% group_by(Ensembl_Danio) %>% # group by ensembl id (since ensemble ids are the one used for query)
    arrange(desc(Mouse_homology_percentage),.by_group = TRUE) %>% #arrange each of them by homology score then by GOC score
    mutate(num_orthologues=n()) %>% #write in a column how many orthologues are possible
    distinct(Ensembl_Danio,.keep_all = TRUE) %>%   # if two orthologues have same homology percentage and goc then choose the top one
    mutate(unique_ortho_check=n()) %>%
    arrange(desc(num_orthologues)) # arrange by the genes which have maximum orthogues
  cat('If worked try again switching old host=\"https://www.ensembl.org\" and mirror = \"useast\" ')
  return(gene_merge)
}
##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Apply Annotation Function
##-----------------------------------------------------------------------------

mouseAnnotation <- function(datatable) {
  
  name <- deparse(substitute(datatable))
  
  annotations <- convertDanio2Mouse(datatable$table$ENSEMBLID)

  annotated_data <- merge(x = annotations,
                        y = datatable$table |>
                          distinct(ENSEMBLID,
                                   .keep_all = TRUE),
                        by.x = "Ensembl_Danio",
                        by.y = "ENSEMBLID",
                        all.x = TRUE)
  
  save(annotated_data, file = paste0("mouse_annotated_", name, ".RData"))
}

##-----------------------------------------------------------------------------


##-----------------------------------------------------------------------------
## Data save
##-----------------------------------------------------------------------------
mouseAnnotation(tt_72_48)
mouseAnnotation(tt_120_72)
mouseAnnotation(tt_adult_120)
##-----------------------------------------------------------------------------

# Clean environment
rm(list = ls())





