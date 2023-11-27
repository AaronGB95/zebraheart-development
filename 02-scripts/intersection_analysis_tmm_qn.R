require(dplyr)
require(ggplot2)
require(gridExtra)
require(VennDiagram)
require(MetBrewer)

contrast_names <- c("72 hpf vs 48 hpf",
                    "120 hpf vs 72 hpf",
                    "Adult vs 120 hpf")

normalizations <- c("TMM", "QN")

# Data Load

## TMM normalization

### 72 hpf vs 48 hpf
load(paste0(dir_output, "differential_expression/", "TMM_tt_72_48.RData"))

#### Upregulated genes
tmm_72_48_ups <- tt_72_48$table[tt_72_48$table$logFC > 1.5 &
                                  tt_72_48$table$FDR < 0.05,
                                "Gene"]

#### Downregulated genes
tmm_72_48_downs <- tt_72_48$table[tt_72_48$table$logFC < -1.5 &
                                    tt_72_48$table$FDR < 0.05,
                                  "Gene"]

### 120 hpf vs 72 hpf
load(paste0(dir_output, "differential_expression/", "TMM_tt_120_72.RData"))

#### Upregulated genes
tmm_120_72_ups <- tt_120_72$table[tt_120_72$table$logFC > 1.5 &
                                    tt_120_72$table$FDR < 0.05,
                                  "Gene"]

#### Downregulated genes
tmm_120_72_downs <- tt_120_72$table[tt_120_72$table$logFC < -1.5 &
                                      tt_120_72$table$FDR < 0.05,
                                    "Gene"]

### Adults vs 120 hpf
load(paste0(dir_output, "differential_expression/", "TMM_tt_adult_120.RData"))

#### Upregulated genes
tmm_adult_120_ups <- tt_adult_120$table[tt_adult_120$table$logFC > 1.5 &
                                          tt_adult_120$table$FDR < 0.05,
                                        "Gene"]

#### Downregulated genes
tmm_adult_120_downs <- tt_adult_120$table[tt_adult_120$table$logFC < -1.5 &
                                            tt_adult_120$table$FDR < 0.05,
                                          "Gene"]

## TMM + QN normalization

### 72 hpf vs 48 hpf
load(paste0(dir_output, "differential_expression/", "QN_tt_72_48.RData"))

#### Upregulated genes
qn_72_48_ups <- tt_72_48$table[tt_72_48$table$logFC > 1.5 &
                                  tt_72_48$table$FDR < 0.05,
                                "Gene"]

#### Downregulated genes
qn_72_48_downs <- tt_72_48$table[tt_72_48$table$logFC < -1.5 &
                                    tt_72_48$table$FDR < 0.05,
                                  "Gene"]

### 120 hpf vs 72 hpf
load(paste0(dir_output, "differential_expression/", "QN_tt_120_72.RData"))

#### Upregulated genes
qn_120_72_ups <- tt_120_72$table[tt_120_72$table$logFC > 1.5 &
                                    tt_120_72$table$FDR < 0.05,
                                  "Gene"]

#### Downregulated genes
qn_120_72_downs <- tt_120_72$table[tt_120_72$table$logFC < -1.5 &
                                      tt_120_72$table$FDR < 0.05,
                                    "Gene"]

### Adults vs 120 hpf
load(paste0(dir_output, "differential_expression/", "QN_tt_adult_120.RData"))

#### Upregulated genes
qn_adult_120_ups <- tt_adult_120$table[tt_adult_120$table$logFC > 1.5 &
                                          tt_adult_120$table$FDR < 0.05,
                                        "Gene"]

#### Downregulated genes
qn_adult_120_downs <- tt_adult_120$table[tt_adult_120$table$logFC < -1.5 &
                                            tt_adult_120$table$FDR < 0.05,
                                          "Gene"]


# Intersection analysis

dir_plot <- paste0(dir_output, "plots/venn_diagrams/")

draw_venn_diagram <- function(genes_list, groups, file_name) {
  VennDiagram::venn.diagram(x = genes_list,
                            category.names = normalizations,
                            filename = file_name,
                            disable.logging = TRUE,
                            output = TRUE,
                            imagetype = "png",
                            height = 480,
                            width = 720,
                            resolution = 300,
                            compression = "lzw",
                            lwd = 1,
                            col = c("#440154ff", '#fde725ff'),
                            fill = c(alpha("#440154ff",0.3),
                                     alpha('#fde725ff',0.3)),
                            cex = 0.7,
                            fontfamily = "sans",
                            cat.pos = c(325, 35),
                            cat.cex = 0.7,
                            cat.default.pos = "outer",
                            cat.dist = c(0.05, 0.05),
                            cat.fontfamily = "sans",
                            cat.col = "black")
}

## 72 hpf vs 48 hpf

### Upregulated genes

D_list <- list(tmm_72_48_ups, qn_72_48_ups)

file_name <- paste0(dir_plot, "venn_tmm_qn_72_48_up.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

### Dowregulated genes

D_list <- list(tmm_72_48_downs, qn_72_48_downs)

file_name <- paste0(dir_plot, "venn_tmm_qn_72_48_down.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

## 120 hpf vs 72 hpf

### Upregulated genes

D_list <- list(tmm_120_72_ups, qn_120_72_ups)

file_name <- paste0(dir_plot, "venn_tmm_qn_120_72_up.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

### Downregulated genes

D_list <- list(tmm_120_72_downs, qn_120_72_downs)

file_name <- paste0(dir_plot, "venn_tmm_qn_120_72_down.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

## Adults vs 120 hpf

### Upregulated genes

D_list <- list(tmm_adult_120_ups, qn_adult_120_ups)

file_name <- paste0(dir_plot, "venn_tmm_qn_adult_120_up.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

### Downregulated genes

D_list <- list(tmm_adult_120_downs, qn_adult_120_downs)

file_name <- paste0(dir_plot, "venn_tmm_qn_adult_120_down.png")

draw_venn_diagram(genes_list = D_list,
                  groups = normalizations,
                  file_name = file_name)

