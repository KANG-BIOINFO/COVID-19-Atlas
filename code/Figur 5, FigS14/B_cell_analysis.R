library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

setwd("./covid19 Atlas/PBMC Atlas/B cells/B cell PBMC_BAL/")

# load cell metadata
df_all = read.table("/data/aronow/Kang/data for toppcell/covid19 Atlas/All/ctype_integrated_v2_subcluster_added_indexAdjusted_v2.txt",
                    sep="\t",header=1, row.names = 1)
df_pbmc_bal = df_all[df_all$Location %in% c("PBMC","BAL"),]
df_B_pbmc_bal = df_pbmc_bal[df_pbmc_bal$Cell.group == "B cell",]

### 1. read data
# 1.1 read PBMC raw data
raw_files = list.files("../../raw data/",pattern = ".h5seurat")
B_list = list()
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  
  B_data = LoadH5Seurat(paste0("../../raw data/",raw_file))
  common_PBs = intersect(colnames(B_data), rownames(df_B_pbmc_bal))
  B_data = B_data[,common_PBs]
  
  B_list[[i]] = B_data
}

# 1.2 read BAL raw data
raw_files_bal = list.files("../../../BAL Atlas/raw data/", pattern = ".h5seurat")
for (i in 6:7)
{
  raw_file = raw_files_bal[i-5]
  
  B_data = LoadH5Seurat(paste0("../../../BAL Atlas/raw data/",raw_file))
  common_PBs = intersect(colnames(B_data), rownames(df_B_pbmc_bal))
  B_data = B_data[,common_PBs]
  
  B_list[[i]] = B_data
}


### 2. integrate
B_list = lapply(X = B_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = B_list)
B_list = lapply(X = B_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = B_list,reduction = "rpca", dim=1:30)
B_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. clustering
B_integrated = ScaleData(B_integrated)
B_integrated = RunPCA(B_integrated, npcs=30)
B_integrated = RunUMAP(B_integrated, reduction = "pca",dims=1:30)

B_integrated = FindNeighbors(B_integrated, dims = 1:30)
B_integrated = FindClusters(B_integrated,resolution = 2)
B_integrated = FindClusters(B_integrated,resolution = 1)
B_integrated = FindClusters(B_integrated,resolution = 0.5)
B_integrated = FindClusters(B_integrated,resolution = 0.4)

df_B_pbmc_bal = df_B_pbmc_bal[colnames(B_integrated),]
B_integrated@meta.data = B_integrated@meta.data[,c("integrated_snn_res.0.5","integrated_snn_res.1","integrated_snn_res.2")]
for (col in colnames(df_B_pbmc_bal)){B_integrated[[col]] = df_B_pbmc_bal[[col]]}
saveRDS(B_integrated, file = "B_integrated.RDS")
write.table(B_integrated@meta.data, file = "ctype_B.txt",sep="\t")
write.table(B_integrated@reductions)

# 6. visualization
DimPlot(B_integrated, group.by = "integrated_snn_res.0.4",label=T,pt.size = 1,label.size = 8.5) +
      theme(legend.position = "none",
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

B_integrated$Location = factor(B_integrated$Location, levels = c("PBMC","BAL"))
DimPlot(B_integrated, group.by = "Location",pt.size = 1, cols = c("#FFE5CC","#3399FF")) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(B_integrated, group.by = "Cell.class_v2",pt.size = 1) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(B_integrated, group.by = "disease.group_v2",pt.size = 1,
        cols = c("#5062AA","#F2D7A6","#EC3284")) +
          theme(legend.position = c(0.25,0.1),
                legend.text = element_text(family = "Helvetica", size = 20),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DefaultAssay(B_integrated) = "RNA"
FeaturePlot(B_integrated, features = c("IGHM","IGHD","IGHG1","IGHG2","IGHA1","IGHA2"),
            cols = c("#F5F5F5","#901026"))

FeaturePlot(B_integrated, features = c("CD79A","MS4A1","BANK1","BLK"),
            cols = c("#F5F5F5","#901026"))

FeaturePlot(B_integrated, features = c("OAS1","OAS3","IFIT1","IFIT3"),
            cols = c("#F5F5F5","#901026"))

df_B_format = read.table("./ctype_B_formatted.txt", sep ="\t", header=1, row.names = 1)
B_integrated = B_integrated[,rownames(df_B_format)]
B_integrated@meta.data = df_B_format
DefaultAssay(B_integrated) = "integrated"
VlnPlot(B_integrated, group.by = "disease.group_v2", features = c("CXCR4","CXCR5"), cols = c("#5062AA", "#F2D7A6","#EC3284"), 
        pt.size = 0.01) + 
        xlab("Disease Condition") + ylab("Scaled Expression")

VlnPlot(B_integrated, group.by = "sub_cluster", 
        features = c("CXCR5"), pt.size = 0.01) + 
        ylab("Scaled Expression") +
        theme(legend.position = "none")

