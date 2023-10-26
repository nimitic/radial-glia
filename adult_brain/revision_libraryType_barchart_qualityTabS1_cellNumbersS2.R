library(Seurat)
library(ggplot2)
library(reshape2)

all_cells <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/brain_pool_corr_reg.rds")
metadata_all_cells <- all_cells@meta.data

radialglia <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/rg_pool_reg_sub.rds")
metadata_RG <- radialglia@meta.data

## cell numbers for table S2
as.data.frame(table(metadata_RG$pc28_res0.8_named_detailed))

## quality numbers for table S1
sample_IDs <- unique(metadata_all_cells$orig.ident)

results <- t(data.frame(row.names = c("sample", "percent.mt", "nUMIs", "nGenes"), c("test",0,0,0)))
results <- results[-1,]

i <- 1
for (i in 1:length(sample_IDs)) {
  percent_mt <- round(mean(subset(metadata_all_cells, orig.ident == sample_IDs[i])$percent.mt), 2)
  nUMIs <- round(mean(subset(metadata_all_cells, orig.ident == sample_IDs[i])$nCount_RNA), 2)
  nGenes <- round(mean(subset(metadata_all_cells, orig.ident == sample_IDs[i])$nFeature_RNA), 2)
  results <- rbind(results, c(as.character(sample_IDs[i]), percent_mt, nUMIs, nGenes))
  i <- i+1
}

results

write.csv(results, file = "20230830_radialgliaRev_libraryQuality.csv")

#### notch inhibition datasets from separate object:

notchInhib <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/perturb_pool_soupx_dbrm.rds")
notchMetadata <- notchInhib@meta.data

sample_IDs_n <- unique(notchMetadata$orig.ident)

results_n <- t(data.frame(row.names = c("sample", "percent.mt", "nUMIs", "nGenes"), c("test",0,0,0)))
results_n <- results_n[-1,]

i <- 1
for (i in 1:length(sample_IDs_n)) {
  percent_mt <- round(mean(subset(notchMetadata, orig.ident == sample_IDs_n[i])$percent.mt), 2)
  nUMIs <- round(mean(subset(notchMetadata, orig.ident == sample_IDs_n[i])$nCount_RNA), 2)
  nGenes <- round(mean(subset(notchMetadata, orig.ident == sample_IDs_n[i])$nFeature_RNA), 2)
  results_n <- rbind(results_n, c(as.character(sample_IDs_n[i]), percent_mt, nUMIs, nGenes))
  i <- i+1
}

results_n

write.csv(results_n, file = "20230830_radialgliaRev_libraryQuality_notchSamples.csv")

## library types for radialglia subtypes
metadata_all_cells <- subset(metadata_all_cells, major_celltypes_pc28_res0.6 == "Radial glia")
metadata_all_cells <- subset(metadata_all_cells, rownames(metadata_all_cells) %in% rownames(metadata_RG))




RD.by.celltype <- table(metadata_all_cells$celltype_detailed, metadata_all_cells$library_type)
RD.by.celltype <- round(prop.table(RD.by.celltype, margin = 1), digits = 10) * 100
RD.by.celltype <- melt(RD.by.celltype, ID = 0)

libtype_colors <- c("dissected" = "#66c2a5", 
                    "sorted" = "#e78ac3", 
                    "whole" = "#8da0cb"
)

ggplot(data = RD.by.celltype, mapping = aes(x = Var1, y = value, fill = Var2))+
  geom_bar(stat = "identity", position = "fill")+
  coord_cartesian(ylim = c(0,1)) + 
  labs(x="", y="Fraction", fill = "Library type") +  
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        plot.title = element_text(size = 18, hjust = 0.5), legend.title = element_text(size = 18), 
        legend.text = element_text(size = 16), legend.justification = "top",
        axis.text.x = element_text(angle = 55, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"))+
  scale_fill_manual(values = libtype_colors) 

ggsave(filename = "RG_barplot_libtype.pdf",
       plot = last_plot(), units = "in", height = 8, width = 20)


## snap25 and ggctb, as well as the generic radial glia markers fabp7a and gfap, across radial glia

#### there is a column called 'ident' in the metadata, which prevents the 'label' argument in FeaturePlot from working correctly.
#### renaming that column fixes the issue - so we are working on an edited copy of the original object here
radialglia_edit <- radialglia

radialglia_edit@meta.data$ident_1 <- radialglia_edit@meta.data$ident
radialglia_edit@meta.data$ident <- NULL


Vln1 <- VlnPlot(radialglia, features = c("fabp7a", "gfap"))
Vln2 <- VlnPlot(radialglia, features = c("snap25a", "ggctb"))

Fea1 <- FeaturePlot(radialglia_edit, features = c("fabp7a", "gfap", "snap25a", "ggctb"), label = TRUE)

cluster_names <- as.data.frame(levels(radialglia_edit@meta.data$pc28_res0.8_named_detailed))
cluster_names$number <- rownames(cluster_names)
cluster_names$cell_type <- cluster_names$`levels(radialglia_edit@meta.data$pc28_res0.8_named_detailed)`
cluster_names$`levels(radialglia_edit@meta.data$pc28_res0.8_named_detailed)` <- NULL
cluster_names <- cluster_names[c(1:18),]

cluster_names[10,2] <- "Radial glia her4++"
cluster_names[6,2] <- "Bergmann glia (prss35+; rhom)"

#install.packages("ggpmisc")
library("ggpmisc")

annotation_1 <- ggplot()+
  annotate(geom = "table",
         x = 0,
         y = 0,
         label = list(cluster_names[c(1:6),])) +
  theme_void()

annotation_2 <- ggplot()+
  annotate(geom = "table",
           x = 0,
           y = 0,
           label = list(cluster_names[c(7:12),])) +
  theme_void()

annotation_3 <- ggplot()+
  annotate(geom = "table",
           x = 0,
           y = 0,
           label = list(cluster_names[c(13:18),])) +
  theme_void()

#radialglia_edit@meta.data$pc28_res0.8_clusters_numeric_edited, 

library(ggpubr)
plots1 <- ggarrange(Vln1, Vln2, nrow = 2, ncol = 1)
annotations <- ggarrange(annotation_1, annotation_2, annotation_3, nrow = 1)
R1 <- ggarrange(plots1, Fea1, annotations, nrow = 3, ncol = 1, heights = c(4.5,4.5,1))


pdf(file = "20230831_radialgliaRev_markersInRGsubtypes.pdf", width = 15, height = 25)
R1
dev.off()

table(radialglia_edit@meta.data$pc28_res0.8_clusters_numeric, radialglia_edit@meta.data$pc28_res0.8_named_detailed)
