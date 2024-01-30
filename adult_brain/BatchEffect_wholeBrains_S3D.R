library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)

# loading subset of the whole dataset, subset to samples coming from whole, (non-dissected, non-sorted) brains only
Zebrafish_brain_wholeBrain <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/brain_pool_corr_reg_wholeBrainOnly.rds")

Zebrafish_brain_wholeBrain@meta.data

## whole brain, non-sorted libraries, all cell types ####

colscale <- brewer.pal(name = "Set1", n = 7)

library_cols = c("b8_mult" = colscale[1],
                 "b11_mult" = colscale[2],
                 "b12_mult" = colscale[3],
                 "b13_mult" = colscale[4],
                 "b15_mult" = colscale[5],
                 "b17_tdmr" = colscale[6],
                 "b19_tdmr" = colscale[7])


DimPlot(Zebrafish_brain_wholeBrain, group.by = "ident", cols = library_cols, pt.size = 0.6) + ggtitle("")

ggsave(filename = "20240117_radialgliaRev_wholeBrainOnly_libraries.pdf",
       plot = last_plot(), units = "in", height = 7, width = 8)



mct_colors <- c("Neurons" = "#80b1d3", #blue
                "Radial glia" = "#bebada", #purple
                "Oligodendroglia" = "#fccde5", #pink
                "Immune cells" = "#8dd3c7", #green
                "Ependymal cells" = "#ffed6f", #yellow
                "Erythrocytes" = "#fb8072", #red
                "Epithelial cells" = "#d9d9d9", #grey
                "Endothelial cells" = "#fdb462" #orange
)

DimPlot(Zebrafish_brain_wholeBrain, group.by = "major_celltypes_pc28_res0.6", label = TRUE, cols = mct_colors, 
        label.size = 5.5, label.box = TRUE, pt.size = 0.6, repel = TRUE, raster = FALSE) + NoLegend() + ggtitle("")

ggsave(filename = "20240117_radialgliaRev_wholeBrainOnly_MainCellTypes.pdf",
       plot = last_plot(), units = "in", height = 7, width = 7)

## whole brain, non-sorted libraries, RG only ####

Zebrafish_brain_wholeBrain_RG <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/rg_pool_reg_sub.rds")
Zebrafish_brain_wholeBrain_RG@meta.data

wholeBrainLibs <- as.data.frame(table(Zebrafish_brain_wholeBrain@meta.data$orig.ident))$Var1

Zebrafish_brain_wholeBrain_RG <- subset(Zebrafish_brain_wholeBrain_RG, orig.ident %in% wholeBrainLibs)

DimPlot(Zebrafish_brain_wholeBrain_RG)

DimPlot(Zebrafish_brain_wholeBrain_RG, group.by = "orig.ident", cols = library_cols, pt.size = 0.6) + ggtitle("") +
  scale_colour_manual(values = library_cols, labels = c("b11 (5prime_v1)", "b12 (5prime_v2)", "b13 (5prime_v2)", "b15 (5prime_v2)", "b17 (3prime_v3.1)", "b19 (3prime_v3.1)", "b8 (3prime_v2)"))

ggsave(filename = "20240117_radialgliaRev_wholeBrain_RG-Only_libraries.pdf",
       plot = last_plot(), units = "in", height = 7, width = 9)

df_RG_origins <- as.data.frame(table(Zebrafish_brain_wholeBrain_RG@meta.data$orig.ident, Zebrafish_brain_wholeBrain_RG@meta.data$pc28_res0.8_named_detailed))
# these cell types have no cells assigned to them
df_RG_origins <- subset(df_RG_origins, Var2 != "Proliferating cells 5'")
df_RG_origins <- subset(df_RG_origins, Var2 != "Radial glia atp1b1b+ (tel/dien)")
df_RG_origins <- subset(df_RG_origins, Var2 != "Radial glia cd74a+")


# one cell type was renamed in the course of the analysis, but not in this specific source file
df_RG_origins <- df_RG_origins %>% 
  mutate(Var2 = str_replace(Var2, "Radial glia gfap", "Radial glia her4"))

ggplot(df_RG_origins, mapping = aes(y = Freq, x = Var2, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = library_cols, name = "", labels = c("b11 (5prime_v1)", "b12 (5prime_v2)", "b13 (5prime_v2)", "b15 (5prime_v2)", "b17 (3prime_v3.1)", "b19 (3prime_v3.1)", "b8 (3prime_v2)")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("number of cells\ncontributing to cluster") +
  xlab("")

ggsave(filename = "20240117_radialgliaRev_wholeBrain_RG-Only_clusters_origins.pdf",
       plot = last_plot(), units = "in", height = 4, width = 9)

col <- brewer.pal(4, "Purples") 
rg_col <- colorRampPalette(col)(18)

umap_num <- DimPlot(Zebrafish_brain_wholeBrain_RG, label = TRUE, cols = rg_col, 
                    label.size = 5.5, label.box = TRUE, pt.size = 0.6, raster = FALSE) + NoLegend()
umap_num

ggsave(filename = "20240117_radialgliaRev_wholeBrain_RG-Only_numberedClusters.pdf",
       plot = last_plot(), units = "in", height = 7, width = 7)

