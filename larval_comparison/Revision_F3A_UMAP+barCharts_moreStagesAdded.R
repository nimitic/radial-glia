library(Seurat)
library(ggplot2)

# use Set2 as palette for stage
stage_colors <-  c("adult" = "#66c2a5",
                   "larvae_15dpf" = "#fc8d62",
                   "larvae_2dpf" = "#8da0cb",
                   "larvae_3dpf" = "#e78ac3",
                   "larvae_5dpf" = "#a6d854",
                   "larvae_8dpf" = "#ffd92f",
                   "larvae_36hpf" = "#E5C494",
                   "larvae_24hpf" = "#B3B3B3",
                   "larvae_20hpf" = "#296d8b")
library(scales)
show_col(stage_colors)

# load dataset, mase using seu_raj_2020_all_v2_wEarlierStages in /local/users/aneusch/jupyterlab/radialglia_revision/
prog.combined <- readRDS(file = "/home/aneusch/data_analysis/20230721_radialgliaRev/int_prog_2dpf-ad_3000vf_wEarlierStages_AN.rds")

# plotting F3A UMAP inclusing all stages

int_umap <- DimPlot(prog.combined, reduction = "umap", label = FALSE, pt.size = 0.6, group.by = "stage", cols = stage_colors, shuffle = TRUE)# +
  #scale_color_manual(values= stage_colors, labels = c("adult", "larvae_15dpf", "larvae_2dpf", "")
int_umap

pdf(file = "20230721_radialgliaRev_3A_UMAP.pdf", width = 8, height = 7)
int_umap
dev.off()

DimPlot(prog.combined, label = TRUE)
DimPlot(prog.combined, group.by = "seurat_clusters", label = TRUE)

# make a plot highlighting proliferation marker genes
# taken from table S2

# cluster 3:	Proliferating cells : mki67+	hmgn2, hmgb2a, tuba8l, mki67, top2a, nusap1, cenpf, pcna, stmn1a, hmgb2b
# cluster 8:	Proliferating cells mcm6+	mcm2, hells, mcm6, cdca7a, mcm5, npm1a, pcna, stmn1a, hmgb2b

prolif_markers <- unique(c("mki67",	"hmgn2", "hmgb2a", "tuba8l", "mki67", "top2a", "nusap1", "cenpf", "pcna", "stmn1a", "hmgb2b",
                    "mcm6",	"mcm2", "hells", "mcm6", "cdca7a", "mcm5", "npm1a", "pcna", "stmn1a", "hmgb2b"))

## percentageFeatureSet uses the raw counts stored in RNA slot
prog.combined[["proliferation_markers"]] <- PercentageFeatureSet(object = prog.combined, features = prolif_markers, assay = "RNA")

umap_prolife <- FeaturePlot(prog.combined, features = "proliferation_markers", cols = c("lightgrey", "#30005c")) +
  ggtitle("proliferation marker genes")
umap_prolife

pdf(file = "20230727_radialgliaRev_3A_UMAP_proliferation.pdf", width = 4, height = 3.5)
umap_prolife
dev.off()


## cell type fraction bar plots start here

table(prog.combined$stage)

# with normalisation
celltype_ad.by.cluster <- table(prog.combined$pc28_res0.8_named_detailed, prog.combined$seurat_clusters)
stage.by.cluster <- table(prog.combined$stage, prog.combined$seurat_clusters)

# how many cells of which identity are in which cluster?
# adult annotation
# prop.table(m, 1) : The value of each cell divided by the sum of the row cells (i.e. normalising by total number of neuroepithelial cells)
celltype_ad.by.cluster <- round(prop.table(celltype_ad.by.cluster, margin = 1), digits = 10) * 100
# how many cells of which stage are in which cluster?
stage.by.cluster <- round(prop.table(stage.by.cluster, margin = 1), digits = 10) * 100 #0

### plan B : what's the dominating adult cell type in each cluster?
### make sure each cluster has a significant fraction of adult cells comprising the above
### otherwise annotate as larval and find closest match to annotate as progenitors of X
### add feature plot prolif. markers

celltype_ad.in.each.cluster <- round(prop.table(celltype_ad.by.cluster, margin = 2), digits = 10) * 100
write.csv(celltype_ad.in.each.cluster, file = "20230724_radialgliaRev_adultCelltypesInEachCluster_larvaeInt.csv")

stage.In.each.cluster <- round(prop.table(stage.by.cluster, margin = 2), digits = 10) * 100
write.csv(stage.In.each.cluster, file = "20230724_radialgliaRev_StagesInEachCluster_larvaeInt.csv")

## manual test case before setting up a function

cumsum_vector <- cumsum(sort(celltype_ad.by.cluster[c("Radial glia snap25a+"),]))
cumsum_vector

cutoff_0.1 <- 0.1*sum(celltype_ad.by.cluster[c("Radial glia snap25a+"),])
cutoff_0.1

# here we figure out which clusters have a siginificant amount of the cell type of interest
clusters_snap25 <- names(cumsum_vector[which(cumsum_vector > cutoff_0.1)])
clusters_snap25

fraction_clusters_snap25 <- celltype_ad.by.cluster[c("Radial glia snap25a+"),clusters_snap25]

# and here we check in the other table, what these clusters are comprised of
sub_snap25 <- stage.by.cluster[, clusters_snap25]
sub_snap25

# this is weighting the clusters by adult contribution (not perfect, but maybe slightly better than straight up using the whole cluster)
i <- 1
for (i in 1:length(colnames(sub_snap25))) {
  y <- clusters_snap25[i]
  sub_snap25[,y] <- sub_snap25[,y] * (fraction_clusters_snap25[[y]]/100)
  i <- i+1
}

sum_snap25 <- rowSums(sub_snap25)

sum_snap25

df_snap25 <- data.frame("Stage" = names(sum_snap25),
                        "Nr_cells" = sum_snap25,
                        "Celltype_adult" = "Radial glia snap25a+")
df_snap25$Stage <- factor(df_snap25$Stage, levels = c("larvae_20hpf", "larvae_24hpf", "larvae_36hpf",
                                                      "larvae_2dpf", "larvae_3dpf", "larvae_5dpf",
                                                      "larvae_8dpf", "larvae_15dpf", "adult"))

head(df_snap25)

ggplot(data = df_snap25, aes(x = Celltype_adult, y = Nr_cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill")

celltype <- "Radial glia snap25a+"

## and here is the function - this one weights the clusters by adult contribution (e.g. 25% of cells are added if 25% of adult cell type are in that cluster)
makeBarplotByStage <- function(celltype_ad.by.cluster, stage.by.cluster, celltype){
  cumsum_vector <- cumsum(sort(celltype_ad.by.cluster[c(celltype),]))
  cutoff_0.1 <- 0.1*sum(celltype_ad.by.cluster[c(celltype),])
  subset_clusters <- names(cumsum_vector[which(cumsum_vector > cutoff_0.1)])
  fraction_clusters <- celltype_ad.by.cluster[celltype,subset_clusters]
  subset_stage <- stage.by.cluster[, subset_clusters]
  i <- 1
  for (i in 1:length(subset_clusters)) {
    y <- subset_clusters[i]
    subset_stage[,y] <- subset_stage[,y] * (fraction_clusters[[y]]/100)
    i <- i+1
  }
  sum_celltype <- rowSums(subset_stage)
  df_celltype <- data.frame("Stage" = names(sum_celltype),
                            "Nr_cells" = sum_celltype,
                            "Celltype_adult" = celltype)
  df_celltype$Stage <- factor(df_celltype$Stage, levels = c("larvae_20hpf", "larvae_24hpf", "larvae_36hpf",
                                                        "larvae_2dpf", "larvae_3dpf", "larvae_5dpf",
                                                        "larvae_8dpf", "larvae_15dpf", "adult"))
  ggsave(filename = paste0("/home/aneusch/journals/analysisjournal/20230727/radialglia_revision/", celltype, "_weightedByAdultContribution.pdf"),
         plot = ggplot(data = df_celltype, aes(x = Celltype_adult, y = Nr_cells, fill = Stage)) +
           geom_bar(stat = "identity", position = "fill") +
           scale_fill_manual(values = stage_colors) +
           ggtitle(celltype) +
           theme_minimal() +
           xlab("") +
           ylab("fraction of cells") +
           theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 text = element_text(size = 20),
                 plot.margin = margin(1,1,1.5,1.2, "cm")) +
           NoLegend(),
         units = "cm", height = 12, width = 12)
}

# run function manually
makeBarplotByStage(celltype_ad.by.cluster, stage.by.cluster, "Radial glia snap25a+")


# run function in a loop for all cell types
all_adult_celltypes <- rownames(celltype_ad.by.cluster)

i <- 1
for (i in c(1:length(all_adult_celltypes))) {
  print(all_adult_celltypes[i])
  makeBarplotByStage(celltype_ad.by.cluster, stage.by.cluster, all_adult_celltypes[i])
  i <- i+1
}


## here we are trying to put all in one plot to save space


df_celltype_all <- t(data.frame(c(1,2,3), row.names = c("Stage", "Nr_cells", "Celltype_adult")))
df_celltype_all <- df_celltype_all[-1,]

i <- 1
for (i in c(1:length(all_adult_celltypes))) {
  
  celltype <- all_adult_celltypes[i]
  
  cumsum_vector <- cumsum(sort(celltype_ad.by.cluster[c(celltype),]))
  cutoff_0.1 <- 0.1*sum(celltype_ad.by.cluster[c(celltype),])
  subset_clusters <- names(cumsum_vector[which(cumsum_vector > cutoff_0.1)])
  fraction_clusters <- celltype_ad.by.cluster[celltype,subset_clusters]
  subset_stage <- stage.by.cluster[, subset_clusters]
  i <- 1
  for (i in 1:length(subset_clusters)) {
    y <- subset_clusters[i]
    subset_stage[,y] <- subset_stage[,y] * (fraction_clusters[[y]]/100)
    i <- i+1
  }
  sum_celltype <- rowSums(subset_stage)
  df_celltype <- data.frame("Stage" = names(sum_celltype),
                            "Nr_cells" = sum_celltype,
                            "Celltype_adult" = celltype)
  df_celltype$Stage <- factor(df_celltype$Stage, levels = c("larvae_20hpf", "larvae_24hpf", "larvae_36hpf",
                                                            "larvae_2dpf", "larvae_3dpf", "larvae_5dpf",
                                                            "larvae_8dpf", "larvae_15dpf", "adult"))
  df_celltype_all <- rbind(df_celltype_all, df_celltype)
  i <- i +1
}

## new cluster name for gfap++
df_celltype_all$Celltype_adult <- gsub(df_celltype_all$Celltype_adult, pattern = "Radial glia gfap++", replacement = "Radial glia her4++", fixed = TRUE)



ggplot(data = df_celltype_all, aes(x = Celltype_adult, y = Nr_cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = stage_colors) +
  theme_minimal() +
  xlab("") +
  ylab("fraction of cells") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size = 20),
        plot.margin = margin(1,1,1.5,1.2, "cm"),
        axis.text.x=element_text(angle=90, hjust=1)) +
  NoLegend()

pdf(file = "20230727_RG_larvaeIntegration_WeightedByAdult.pdf", width = 10, height = 6)
last_plot()
dev.off()

### the same without adult cells (no weighting!)


# with normalisation
celltype_ad.by.cluster <- table(prog.combined$pc28_res0.8_named_detailed, prog.combined$seurat_clusters)
stage.by.cluster <- table(prog.combined$stage, prog.combined$seurat_clusters)

# how many cells of which identity are in which cluster?
# adult annotation
# prop.table(m, 1) : The value of each cell divided by the sum of the row cells (i.e. normalising by total number of neuroepithelial cells)
celltype_ad.by.cluster <- round(prop.table(celltype_ad.by.cluster, margin = 1), digits = 10) * 100
# how many cells of which stage are in which cluster?
stage.by.cluster <- round(prop.table(stage.by.cluster, margin = 1), digits = 10) * 1000

# this is what removes the adult cells
stage.by.cluster <- stage.by.cluster[-1,]

makeBarplotByStageNoAdult <- function(celltype_ad.by.cluster, stage.by.cluster, celltype){
  cumsum_vector <- cumsum(sort(celltype_ad.by.cluster[c(celltype),]))
  cutoff_0.1 <- 0.1*sum(celltype_ad.by.cluster[c(celltype),])
  subset_clusters <- names(cumsum_vector[which(cumsum_vector > cutoff_0.1)])
  subset_stage <- stage.by.cluster[, subset_clusters]
  sum_celltype <- rowSums(subset_stage)
  df_celltype <- data.frame("Stage" = names(sum_celltype),
                            "Nr_cells" = sum_celltype,
                            "Celltype_adult" = celltype)
  df_celltype$Stage <- factor(df_celltype$Stage, levels = c("larvae_20hpf", "larvae_24hpf", "larvae_36hpf",
                                                            "larvae_2dpf", "larvae_3dpf", "larvae_5dpf",
                                                            "larvae_8dpf", "larvae_15dpf"))
  ggsave(filename = paste0("/home/aneusch/journals/analysisjournal/20230721/radialglia_revision/", celltype, ".pdf"),
         plot = ggplot(data = df_celltype, aes(x = Celltype_adult, y = Nr_cells, fill = Stage)) +
           geom_bar(stat = "identity", position = "fill") +
           scale_fill_manual(values = stage_colors) +
           ggtitle(celltype) +
           theme_minimal() +
           xlab("") +
           ylab("fraction of cells") +
           theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank()) +
           NoLegend(),
         units = "in", height = 4, width = 3)
}

# run function in a loop for all cell types
all_adult_celltypes <- rownames(celltype_ad.by.cluster)

i <- 1
for (i in c(1:length(all_adult_celltypes))) {
  print(all_adult_celltypes[i])
  makeBarplotByStageNoAdult(celltype_ad.by.cluster, stage.by.cluster, all_adult_celltypes[i])
  i <- i+1
}


