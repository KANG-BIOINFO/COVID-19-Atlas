library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)

options(future.globals.maxSize=(90*1024^3)) # set max size as 40G
setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/Plasmablast/PBMC_BAL//")

# load cell metadata
df_all = read.table("/data/aronow/Kang/data for toppcell/covid19 Atlas/All/ctype_integrated_v2_subcluster_added_indexAdjusted_v2.txt",
                    sep="\t",header=1, row.names = 1)
df_pbmc_bal = df_all[df_all$Location %in% c("PBMC","BAL"),]
df_PB_pbmc_bal = df_pbmc_bal[df_pbmc_bal$Cell.group == "Plasmablast",]

### 1. read data
# 1.1 read PBMC raw data
raw_files = list.files("../../raw data/",pattern = ".h5seurat")
PB_list = list()
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  
  PB_data = LoadH5Seurat(paste0("../../raw data/",raw_file))
  common_PBs = intersect(colnames(PB_data), rownames(df_PB_pbmc_bal))
  PB_data = PB_data[,common_PBs]
  
  PB_list[[i]] = PB_data
}

# 1.2 read BAL raw data
raw_files_bal = list.files("../../../BAL Atlas/raw data/", pattern = ".h5seurat")
for (i in 6:7)
{
  raw_file = raw_files_bal[i-5]
  
  PB_data = LoadH5Seurat(paste0("../../../BAL Atlas/raw data/",raw_file))
  common_PBs = intersect(colnames(PB_data), rownames(df_PB_pbmc_bal))
  PB_data = PB_data[,common_PBs]
  
  PB_list[[i]] = PB_data
}


### 2. integrate
PB_list = lapply(X = PB_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = PB_list)
PB_list = lapply(X = PB_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = PB_list,reduction = "rpca", dim=1:30)
PB_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. clustering
#PB_integrated = PB_integrated[,!PB_integrated$integrated_snn_res.0.5 %in% c(5,7,9,10)]
PB_integrated = ScaleData(PB_integrated)
PB_integrated = RunPCA(PB_integrated, npcs=30)
PB_integrated = RunUMAP(PB_integrated, reduction = "pca",dims=1:30)

PB_integrated = FindNeighbors(PB_integrated, dims = 1:30)
PB_integrated = FindClusters(PB_integrated,resolution = 2)
PB_integrated = FindClusters(PB_integrated,resolution = 1)
PB_integrated = FindClusters(PB_integrated,resolution = 0.5)
PB_integrated = FindClusters(PB_integrated,resolution = 0.3)

df_PB_pbmc_bal = df_PB_pbmc_bal[colnames(PB_integrated),]
PB_integrated@meta.data = PB_integrated@meta.data[,c("integrated_snn_res.0.5","integrated_snn_res.1","integrated_snn_res.2")]
for (col in colnames(df_PB_pbmc_bal)){PB_integrated[[col]] = df_PB_pbmc_bal[[col]]}
saveRDS(PB_integrated, file = "PB_integrated.RDS")
saveRDS(PB_integrated, file = "PB_integrated_doubletRemoved.RDS")

# 6. visualization
PB_integrated$subtype = mapvalues(PB_integrated$integrated_snn_res.0.3,
                                         from = c(0,2,4,1,3),
                                         to = c(rep("Plasmablast",3),rep("Developping Plasmablast",2)))
DimPlot(PB_integrated, group.by = "integrated_snn_res.0.3",label=T,pt.size = 1,label.size = 8.5) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

PB_integrated$Location = factor(PB_integrated$Location, levels = c("PBMC","BAL"))
DimPlot(PB_integrated, group.by = "Location",pt.size = 1, cols = c("#FFE5CC","#3399FF")) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(PB_integrated, group.by = "sub_cluster_v2",pt.size = 1,label=T, label.size = 10) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(PB_integrated, group.by = "disease.group_v2",pt.size = 1, cols = c("#5062AA","#F2D7A6","#EC3284")) +
        theme(legend.position = c(0.02,0.1),
              legend.text = element_text(family = "Helvetica", size = 20),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DefaultAssay(PB_integrated) = "RNA"
FeaturePlot(PB_integrated, features = c("IGHG1","IGHG2","IGHA1","IGHA2","IGHM","IGHD"),
            cols = c("#F5F5F5","#901026"))
FeaturePlot(PB_integrated, features = c("BCL2","CD79A","PLK1","MKI67"),
            cols = c("#F5F5F5","#901026"))

# draw DE analysis
library(EnhancedVolcano)

setwd("./covid19 Atlas/PBMC Atlas/Plasmablast/PBMC_BAL/")

PB_vs_DevPB = read.table("./DEGs_DevPB_vs_PB.txt", sep="\t", header=1)
PB_vs_DevPB[PB_vs_DevPB$pvals_adj < 1e-100,"pvals_adj"] = 1e-100
PB_vs_DevPB$names = as.character(PB_vs_DevPB$names)

EnhancedVolcano(PB_vs_DevPB, lab = PB_vs_DevPB$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("MKI67","CDC7","MCM2","PLK1","CDC6","CDCA2","CENPA","CENPE","CENPF",
                              "SRSF1","SRSF3","HLA-DOB","IGHA1","IGHA2","IGKC","IGLC2",
                              "IL2RG","TNFRSF13C","BCL2","LAX1","PRDM1","CD79A","CD40","CD27","ITM2A","IGF1"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "PB vs. Developping PB", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-4,8))
# Dev barplot
df_DevPB = read.table("./Dev_PB_enrichment.txt",sep="\t",header=1)
barplot(df_DevPB$Gene.Enrichment.Score,horiz=T, col = "#00BFC4")

# PB barplot
df_PB = read.table("PB_enrichment.txt",sep="\t",header=1)
barplot(df_PB$Gene.Enrichment.Score,horiz=T, col = "#F8766D")

