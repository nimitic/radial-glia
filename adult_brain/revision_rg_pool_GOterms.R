library(reshape2)
library(ggplot2)
library(DESeq2)
library(gprofiler2)
library(dplyr)

RG_marker_genes <- read.csv("/home/aneusch/journals/analysisjournal/20230727/radialglia_revision/20230727_rg_pool_sub_markers.csv")


RG_marker_genes$X <- NULL

RG_marker_genes <- subset(RG_marker_genes, p_val_adj <= 0.01 )

RG_marker_genes




# Adding column based on other column:
RG_marker_genes <- RG_marker_genes %>%
  mutate(celltype = case_when(
    cluster == 1 ~ "Radial glia id2b+",
    cluster == 2 ~ "Radial glia snap25a+",
    cluster == 3 ~ "Proliferating cells mki67+",
    cluster ==4 ~ "Radial glia nppc+",
    cluster ==5 ~ "Radial glia apof+",
    cluster ==6 ~ "Radial glia prss35+",
    cluster ==7 ~ "Radial glia id3+",
    cluster ==8 ~ "Proliferating cells mcm6+",
    cluster ==9 ~ "Radial glia enkur+",
    cluster ==10 ~ "Radial glia her4++",
    cluster ==11 ~ "Radial glia foxn4+",
    cluster ==12 ~ "Radial glia crabp1b+",
    cluster ==13 ~ "Radial glia ntn1b+",
    cluster ==14 ~ "Radial glia nrg1+",
    cluster ==15 ~ "Neuroepithelial cells",
    cluster ==16 ~ "Radial glia RP high",
    cluster ==17 ~ "Radial glia stat2+",
    cluster ==18 ~ "Radial glia stra6+"
  ))


write.csv(RG_marker_genes, file = "20230726_RG_pool_markers.csv")




celltype = "Radial glia id2b+"

RunGost <- function(RG_marker_genes, celltype){
  gostres <- gost(query = subset(RG_marker_genes, RG_marker_genes$celltype == celltype)$gene,
                  organism = "drerio",
                  significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.01, correction_method = "fdr",
                  highlight = TRUE)
  gostres$result$parents <- NULL
  write.csv(gostres$result, file = paste0("20230727_RG_pool_", celltype, "_GOterms_DiffExprGenes_p001.csv"))
}


celltypes <- unique(RG_marker_genes$celltype)

i <- 1
for (i in 1:length(celltypes)) {
  print(celltypes[i])
  RunGost(RG_marker_genes, celltype = celltypes[i])
  i <- i + 1
}





          