library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

options(future.globals.maxSize=(110*1024^3)) # set max size as 60G
setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/disease comparison/")


# 1. integrating datasets 
pbmc_list = list()
sepsis_raw = LoadH5Seurat("./data/sepsis_PBMC_raw.h5seurat")
ms_raw = LoadH5Seurat("./data/ms_PBMC_raw_merged.h5Seurat")

# a loop for reading 5 COVID-19 PBMC datasets
Dataset_names = c("COVID-19_Aru","COVID-19_Guo","COVID-19_Lee","COVID-19_Sch","COVID-19_Wil")
raw_files = list.files("../raw data/",pattern = ".h5seurat")
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  pbmc_data = LoadH5Seurat(paste0("../raw data/",raw_file))
  pbmc_data$Dataset = Dataset_names[i]
  
  pbmc_list[[i]] = pbmc_data
}

sepsis_raw$Dataset = "Sepsis"
pbmc_list[[6]] = sepsis_raw
ms_raw$Dataset = "Multiple Sclerosis"
pbmc_list[[7]] = ms_raw

# 2.  normalization and variable genes
pbmc_list = lapply(X = pbmc_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = pbmc_list)
pbmc_list = lapply(X = pbmc_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = pbmc_list,reduction = "rpca", dim=1:30)
pbmc_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. downstream analysis
pbmc_integrated = ScaleData(pbmc_integrated)
pbmc_integrated = RunPCA(pbmc_integrated, npcs=30)
pbmc_integrated = RunUMAP(pbmc_integrated, reduction = "pca",dims=1:30)
pbmc_integrated = FindNeighbors(pbmc_integrated, dims = 1:30)
pbmc_integrated = FindClusters(pbmc_integrated, resolution = 1.0)
pbmc_integrated = FindClusters(pbmc_integrated, resolution = 0.5)
pbmc_integrated = FindClusters(pbmc_integrated, resolution = 2.0)
saveRDS(pbmc_integrated,"PBMC_multiDiseases.RDS")

# 6. read cell annotations
df_cells = read.table("ctype_integrated_concise_v3_adjusted.txt",sep="\t",header=1,row.names = 1)
pbmc_integrated = pbmc_integrated[,rownames(df_cells)]
pbmc_integrated@meta.data = df_cells
DimPlot(pbmc_integrated, group.by = "Lineage")

pbmc_integrated$Disease_concise = mapvalues(pbmc_integrated$Disease,
                                            from = c("Leuk-UTI","URO","Int-URO","ICU-NoSEP","ICU-SEP","Bac-SEP","COVID-19 Moderate/Mild/Early stage/NonVent","COVID-19 Severe/Late stage/Vent","COVID-19 Remission"),
                                            to = c(rep("Sepsis",6),"COVID-19 Mild","COVID-19 Severe","COVID-19 Mild"))
pbmc_integrated$Disease_concise = factor(pbmc_integrated$Disease_concise, levels = c("Control","IIF","Multiple Sclerosis","Sepsis","Influenza Severe","COVID-19 Mild","COVID-19 Severe"))

DimPlot(pbmc_integrated,group.by = "Disease_concise",cols = c("#F2FAF8","#DBFFFF","#336600","#FFCC99","#CC0000","#F9F99B","#FFABD5")) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=24,family = "Helvetica",face="bold")) + ggtitle("")

pbmc_integrated$Cell.class_v2_factor = factor(pbmc_integrated$Cell.class_v2,
                                              levels = c("B naive","B memory","B intermediate","Plasmablast",
                                                         "CD4+ CTL","CD4+ T naive","CD4+ Tcm","CD4+ Tem","CD4+ T activated","Treg",
                                                         "CD8+ T naive","CD8+ Tcm","CD8+ Tem","MAIT","CD8+ T activated",
                                                         "dn T","gd T","T/NK proliferative",
                                                         "NK","NK CD56bright","NK activated","ILC",
                                                         "CD14+ Monocyte","CD16+ Monocyte","Neutrophil","immature Neutrophil","cDC","cDC1","cDC2","tDC","pDC",
                                                         "HSPC","RBC","Platelet"))

colors_class = c("#FFD6C7","#FFFF66","#FF6666","#FFA500",
                 "#B0E0E6","#4682B4","#CCE5FF","#00008B","#FF69B4","#CC0066",
                 "#FF3399","#D8BFD8","#CDB5E5","#9933FF","#0066CC",
                 "#F9842B","#666600","#336600",
                 "#FFFACD","#F4A460","#CD853F","#330066",
                 "#B0C4DE","#778899","#CD5C5C","#DEA4A4","#BDB76B","#FFCCE5","#BDB76B","#663300","#999900",
                 "#9ACD32","#0BAF0B","#FF0000")

DimPlot(pbmc_integrated,group.by = "Cell.class_v2_factor",cols = colors_class) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=24,family = "Helvetica",face="bold")) + ggtitle("")



