library(Seurat)
library(ggplot2)

radial_glia <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/rg_pool_reg_sub.rds")


DotPlot(radial_glia, features = c("snap25a", "ggctb", "AL954697.1", "crabp1b", "prss35", "apof"), group.by = "pc28_res0.8_named_detailed", cols = "RdBu") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Radial glia subtypes")

ggsave(last_plot(), filename = "20231025_RadialgliaRev_RNAscopeProbes_RNAseq-expression_radialglia.pdf", width = 13, height = 5)


neurons <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/neu_pool_reg.rds")
neurons@meta.data

DotPlot(neurons, features = c("snap25a", "ggctb"), group.by = "pc30_res0.8_biological_names", cols = "RdBu") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("Neuronal subtypes")

ggsave(last_plot(), filename = "20231025_RadialgliaRev_RNAscopeProbes_RNAseq-expression_neurons.pdf", width = 13, height = 5)

rm(neurons)
rm(radial_glia)
gc()

all <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/brain_pool_corr_reg.rds")

all@meta.data

# done on the server
all_dissected <- subset(x = all, subset = library_type == "dissected")

all_dissected <- readRDS("/home/aneusch/data_analysis/20230721_radialgliaRev/20231025_brain_pool_corr_reg_dissectedOnly.rds")

table(all_dissected@meta.data$region_detailed)

all_dissected_sub <- subset(all_dissected, subset = region_detailed %in% c("Dien", "Mes", "Rhom", "Tel"))

# for peace of mind regarding the failed cebpd probe (showed expression in diencephalon and mesencephalon)
VlnPlot(all_dissected_sub, features = "cebpd")
