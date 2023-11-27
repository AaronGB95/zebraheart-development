
require(dplyr)
require(ggplot2)
require(gridExtra)
require(VennDiagram)
require(MetBrewer)


load("tt_72_48.RData")
tt_72_48$tablel$ind <- rownames(tt_72_48$table)

load("tt_120_72.RData")
tt_120_72$tablel$ind <- rownames(tt_120_72$table)

load("tt_adult_120.RData")
tt_adult_120$tablel$ind <- rownames(tt_adult_120$table)

# Resultados empleando Combat

load("tt_72_48_combat.RData")
tt_72_48_combat$tablel$ind <- rownames(tt_72_48_combat$table)

load("tt_120_72_combat.RData")
tt_120_72_combat$tablel$ind <- rownames(tt_120_72_combat$table)

load("tt_adult_120_combat.RData")
tt_adult_120_combat$tablel$ind <- rownames(tt_adult_120_combat$table)

# Resultados empleando SVA

load("tt_72_48_sva.RData")
tt_72_48_sva$tablel$ind <- rownames(tt_72_48_sva$table)

load("tt_120_72_sva.RData")
tt_120_72_sva$tablel$ind <- rownames(tt_120_72_sva$table)

load("tt_adult_120_sva.RData")
tt_adult_120_sva$tablel$ind <- rownames(tt_adult_120_sva$table)

n_ups <- NULL
n_downs <- NULL
contrasts_names <- c("72hpf vs 48hpf", "120hpf vs 72hpf", "Adult vs 120hpf")
batch_names <- c("Control", "ComBat", "SVA")

diff_df <- tt_72_48$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

control_1_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
control_1_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_120_72$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

control_2_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
control_2_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_adult_120$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

control_3_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
control_3_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_72_48_combat$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

combat_1_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
combat_1_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_120_72_combat$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

combat_2_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
combat_2_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_adult_120_combat$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

combat_3_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
combat_3_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_72_48_sva$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

sva_1_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
sva_1_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_120_72_sva$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

sva_2_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
sva_2_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

diff_df <- tt_adult_120_sva$table[,c("logFC", "FDR")]
diff_df$gene_name <- rownames(diff_df)

diff_df$group <- "N.S."

diff_df[which(diff_df['FDR'] < 0.05 & diff_df['logFC'] > 1.5 ),"group"] <- "Up"
diff_df[which(diff_df['FDR'] < 0.05 &
diff_df['logFC'] <  -1.5 ),"group"] <- "Down"

n_ups <- c(n_ups,length(which(diff_df$group == "Up")))
n_downs <- c(n_downs,length(which(diff_df$group == "Down")))

sva_3_ups <- rownames(diff_df[which(diff_df["group"] == "Up"),])
sva_3_downs <- rownames(diff_df[which(diff_df["group"] == "Down"),])

# Intersección de genes significativos {.tabset}

## Control {.tabset}

### Upregulated
D_ups <- list(control_1_ups, control_2_ups, control_3_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = contrasts_names,
                          filename = "venn_control_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_control_up.png")


### Downregulated


D_downs <- list(control_1_downs, control_2_downs, control_3_downs)

VennDiagram::venn.diagram(x = D_downs,
                          category.names = contrasts_names,
                          filename = "venn_control_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_control_down.png")


## ComBat {.tabset}

### Upregulated
D_ups <- list(combat_1_ups, combat_2_ups, combat_3_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = contrasts_names,
                          filename = "venn_combat_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_combat_up.png")


### Downregulated
D_downs <- list(combat_1_downs, combat_2_downs, combat_3_downs)

VennDiagram::venn.diagram(x = D_downs,
                          category.names = contrasts_names,
                          filename = "venn_combat_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_combat_down.png")


## SVA {.tabset}

### Upregulated
D_ups <- list(sva_1_ups, sva_2_ups, sva_3_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = contrasts_names,
                          filename = "venn_sva_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_sva_up.png")
```

### Downregulated

D_downs <- list(sva_1_downs, sva_2_downs, sva_3_downs)

VennDiagram::venn.diagram(x = D_downs,
                          category.names = contrasts_names,
                          filename = "venn_sva_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_sva_down.png")

# Intersección por método de corrección del efecto batch {.tabset}

## 72 hpf vs 48 hpf {.tabset}

### Upregulated
D_ups <- list(control_1_ups, combat_1_ups, sva_1_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_72_48_batch_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_72_48_batch_up.png")

### Downregulated

D_ups <- list(control_1_downs, combat_1_downs, sva_1_downs)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_72_48_batch_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_72_48_batch_down.png")

## 120 hpf vs 72 hpf {.tabset}

### Upregulated
D_ups <- list(control_2_ups, combat_2_ups, sva_2_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_120_72_batch_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_120_72_batch_up.png")


### Downregulated
D_ups <- list(control_2_downs, combat_2_downs, sva_2_downs)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_120_72_batch_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_120_72_batch_down.png")

## Adult vs 120 hpf {.tabset}

### Upregulated
D_ups <- list(control_3_ups, combat_3_ups, sva_3_ups)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_adult_120_batch_up.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_adult_120_batch_up.png")


### Downregulated
D_ups <- list(control_3_downs, combat_3_downs, sva_3_downs)

VennDiagram::venn.diagram(x = D_ups,
                          category.names = batch_names,
                          filename = "venn_adult_120_batch_down.png",
                          output = TRUE,
                          imagetype = "png",
                          height = 480,
                          width = 720,
                          resolution = 300,
                          compression = "lzw",
                          lwd = 1,
                          col = c("#440154ff", '#21908dff', '#fde725ff'),
                          fill = c(alpha("#440154ff",0.3),
                                   alpha('#21908dff',0.3),
                                   alpha('#fde725ff',0.3)),
                          cex = 0.7,
                          fontfamily = "sans",
                          cat.cex = 0.7,
                          cat.default.pos = "outer",
                          cat.pos = c(-27, 27, 135),
                          cat.dist = c(0.055, 0.055, 0.055),
                          cat.fontfamily = "sans",
                          cat.col = "black",
                          rotation = 1)

knitr::include_graphics("venn_adult_120_batch_down.png")
