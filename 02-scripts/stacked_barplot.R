require(ggplot2)
require(dplyr)
require(tidyverse)
require(hrbrthemes)
require(reshape)

##------------------------------------------------------------------------------
## General Terms
##------------------------------------------------------------------------------
# Load child terms in every general term
ont1 <- "Heart development"
heart_development <- read.table(file = paste0(dir_docs,
                                              "heart_development.txt"))

heart_development <- as.character(heart_development$V1)

ont2 <- "Cardiac contraction"
cardiac_contraction <- read.table(file = paste0(dir_docs,
                                                "cardiac_contraction.txt"))

cardiac_contraction <- as.character(cardiac_contraction$V1)

ont3 <- "Calcium transport"
calcium_transport <- read.table(file = paste0(dir_docs,
                                              "calcium_transport.txt"))

calcium_transport <- as.character(calcium_transport$V1)

onts <- list(heart_development, cardiac_contraction, calcium_transport)

names(onts) <- c(ont1, ont2, ont3)

contrasts_names <- c("72 hpf vs 48 hpf",
                     "120 hpf vs 72 hpf",
                     "Adult vs 120 hpf")
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------
## Data transformation function
##------------------------------------------------------------------------------
df_creation <- function(organism, dataset, ontology) {

  # Initialize an empty data frame
  result_df <- matrix(nrow = length(onts) + 1,
                      ncol = length(contrasts_names))
  
  rownames(result_df) <- c(names(onts), "Other")
  
  colnames(result_df) <- contrasts_names
  
  result_df <- replace(result_df, is.na(result_df), 0)
  
  # Get dataframe of significative terms in 72 vs 48 hpf
  load(paste0(dir_output,
              "gene_set_enrichment_analysis/",
              organism,
              "_TMM_GSEA_",
              dataset,
              "_",
              ontology,
              "_72_48.RData"))
  
  sig_terms_1 <- gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"]
  
  # Get dataframe of significative terms in 120 vs 72 hpf
  load(paste0(dir_output,
              "gene_set_enrichment_analysis/",
              organism,
              "_TMM_GSEA_",
              dataset,
              "_",
              ontology,
              "_120_72.RData"))
  
  sig_terms_2 <- gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"]
  
  # Get dataframe of significative terms in Adult vs 120 hpf
  load(paste0(dir_output,
              "gene_set_enrichment_analysis/",
              organism,
              "_TMM_GSEA_",
              dataset,
              "_",
              ontology,
              "_adult_120.RData"))
  
  sig_terms_3 <- gsea_result@result[gsea_result@result$p.adjust < 0.05, "ID"]
  
  terms <- list(sig_terms_1, sig_terms_2, sig_terms_3)
  
  names(terms) <- contrasts_names
  
  for (contrast in contrasts_names) {
    for (id in terms[[contrast]]) {
      if (id %in% heart_development) {
        result_df["Heart development", contrast] = result_df["Heart development", contrast] + 1
      } else if (id %in% cardiac_contraction) {
        result_df["Cardiac contraction", contrast] = result_df["Cardiac contraction", contrast] + 1
      } else if (id %in% calcium_transport) {
        result_df["Calcium transport", contrast] = result_df["Calcium transport", contrast] + 1
      } else {
        result_df["Other", contrast] = result_df["Other", contrast] + 1
      }
      
    }
    
  }
  
  return(result_df)
  }
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
## Apply Data Transformation
##------------------------------------------------------------------------------

zebrafish_og_bp <- df_creation("zebrafish", "OG", "GOBP")
zebrafish_db_bp <- df_creation("zebrafish", "DB", "GOBP")

human_og_bp <- df_creation("human", "OG", "GOBP")
human_db_bp <- df_creation("human", "DB", "GOBP")

mouse_og_bp <- df_creation("mouse", "OG", "GOBP")
mouse_db_bp <- df_creation("mouse", "DB", "GOBP")

##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
## Stacked Barplots
##------------------------------------------------------------------------------

stacked_barplot <- function(dataframe, name) {
  dataframe <- melt(dataframe)
  colnames(dataframe) <- c("term", "contrast", "count")
  
  plot <- ggplot(dataframe, aes(fill = term, y = count, x = contrast)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_viridis_d(name = "GO term") +
    xlab("Comparision") +
    ylab("Number of terms") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(plot = plot,
         filename = paste0("04-output/plots/functional_classification/",
                           name,
                           ".png"))
}

##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
## Create plots
##------------------------------------------------------------------------------

stacked_barplot(zebrafish_og_bp, "zebrafish_OG_GOBP")

stacked_barplot(zebrafish_db_bp, "zebrafish_DB_GOBP")

stacked_barplot(human_og_bp, "human_OG_GOBP")

stacked_barplot(human_db_bp, "human_DB_GOBP")

stacked_barplot(mouse_og_bp, "mouse_OG_GOBP")

stacked_barplot(mouse_db_bp, "mouse_DB_GOBP")

##------------------------------------------------------------------------------
