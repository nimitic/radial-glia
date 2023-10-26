library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

brain_pool <- readRDS(file = "/home/aneusch/data_analysis/20230721_radialgliaRev/brain_pool_corr_reg.rds")

brain_pool@meta.data

libtype_colors <- c("dissected" = "#66c2a5", 
                    "sorted" = "#e78ac3", 
                    "whole" = "#8da0cb"
)

S2b_new <- ggplot(brain_pool@meta.data, mapping = aes(x = major_celltypes_pc28_res0.6, fill = library_type)) +
  geom_bar(position="fill") +
  scale_fill_manual(values = libtype_colors) +
  coord_cartesian(ylim = c(0,1)) + 
  labs(x="", y="Fraction", fill = "Library type") +  
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 18), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 20), legend.text = element_text(size = 18), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"))


ggsave(filename = "20230927_radialgliaRev_S2B_non-normalised.pdf",
       plot = S2b_new, units = "in", height = 8, width = 12)
