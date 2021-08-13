library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

DC_list = list()
df_DCs = read.table("./covid19 Atlas/BAL Atlas/cDC/ctype_cDC_v1.txt",
                    sep="\t",header = 1,row.names = 1)
all_DCs = rownames(df_DCs)

count_datasets = 1
setwd("./covid19 Atlas/PBMC Atlas/raw data/")
raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Arunachalam","Guo","Lee","Schulte-Schrepping","Wilk")

# 1. read PBMC data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  DC_data = LoadH5Seurat(raw_file)
  DC_data$Dataset = Dataset_names[i]
  
  common_DCs = intersect(colnames(DC_data), all_DCs)
  DC_data = DC_data[,common_DCs]
  print(sample_name); print(length(common_DCs))
  
  if (length(common_DCs) > 50)
  {
    DC_list[[count_datasets]] = DC_data
    count_datasets = count_datasets + 1
  }
}

# read BAL data
setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/BAL Atlas/raw data/")
raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Grant","Liao")
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  DC_data = LoadH5Seurat(raw_file)
  DC_data$Dataset = Dataset_names[i]
  
  common_DCs = intersect(colnames(DC_data), all_DCs)
  DC_data = DC_data[,common_DCs]
  print(sample_name); print(length(common_DCs))
  
  if (length(common_DCs) > 50)
  {
    DC_list[[count_datasets]] = DC_data
    count_datasets = count_datasets + 1
  }
}


# 2. normalization and variable genes
DC_list = lapply(X = DC_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = DC_list)
DC_list = lapply(X = DC_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = DC_list,reduction = "rpca", dim=1:30)
DC_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. Visualization
DC_integrated = ScaleData(DC_integrated)
DC_integrated = RunPCA(DC_integrated, npcs=30)
DC_integrated = RunUMAP(DC_integrated, reduction = "pca",dims=1:30)

DC_integrated = FindNeighbors(DC_integrated, dims = 1:30)
DC_integrated = FindClusters(DC_integrated, resolution = 1)
DC_integrated = FindClusters(DC_integrated, resolution = 2)
DC_integrated = FindClusters(DC_integrated, resolution = 0.5)

setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/BAL Atlas/cDC/")
saveRDS(DC_integrated, file = "DC_integrated_Rep_PCA.RDS")

DC_integrated$cell.class_v2 = mapvalues(DC_integrated$Dataset, from = c("Arunachalam","Lee","Schulte-Schrepping","Wilk","Grant","Liao"),
                                        to = c(rep("cDC_blood",4), rep("cDC_lung",2)))
df_DCs = df_DCs[colnames(DC_integrated),]
DC_integrated$disease.group_v2 = df_DCs$disease_v3

write.table(DC_integrated_v2@meta.data, file = "ctype_cDC_v2_clustered.txt" , sep="\t")
write.table(DC_integrated_v2@reductions$umap@cell.embeddings, file = "umap.txt",sep="\t")
DimPlot(DC_integrated, group.by = "cell.class_v2")+
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.position = c(0.7,0.93), legend.text = element_text(size=24))

DimPlot(DC_integrated, group.by = "integrated_snn_res.0.5",label=T, label.size = 10)+
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.position = "none")

DC_integrated_v2 = DC_integrated[,DC_integrated$disease.group_v2!="Influenza Severe"]
DC_integrated_v2$disease.group_v2 = as.factor(as.character(DC_integrated_v2$disease.group_v2))
DC_integrated_v2$disease.group_v2 = factor(DC_integrated_v2$disease.group_v2, 
                                           levels = c("Healthy Control_blood","COVID-19_Mild/Remission_blood","COVID-19_Severe_blood",
                                                      "Healthy Control_lung","COVID-19_lung"))
colors_condition = c("#87CEEB","#FFE4C4","#99004C","#6495ED","#FF1493")
DimPlot(DC_integrated_v2, group.by = "disease.group_v2",cols = colors_condition)+
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=20))

DefaultAssay(DC_integrated) = "RNA"
FeaturePlot(DC_integrated, features = c("CCR7","FCER1A","IL12B","OAS1","IFIT1"), cols = c("#F5F5F5","#901026"),ncol = 4)
FeaturePlot(DC_integrated, features = c("HLA-DQA1","HLA-DRB1"), cols = c("#F5F5F5","#901026"),ncol = 4)
FeaturePlot(DC_integrated, features = c("CCR7","CSF2RA","FCER1A","CD1C",
                                        "HLA-DQA2","MKI67","OAS1","IFIT1"), cols = c("#F5F5F5","#901026"),ncol = 4)

# box
Idents(DC_integrated) = paste0("DC_",Idents(DC_integrated))
DimPlot(DC_integrated, label = T, label.size = 7,pt.size = 1) + 
      theme(legend.position = "none",
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.border = element_rect(colour = "#FCB7B3",fill=NA,size=5)) + ggtitle("")

DC_integrated = DC_integrated[,DC_integrated$disease.group_v2 != "Influenza Severe"]
DC_integrated$disease.group_v2 = as.factor(as.character(DC_integrated$disease.group_v2))
DC_integrated$disease.group_v2 = factor(DC_integrated$disease.group_v2, 
                                           levels = c("Healthy Control_blood","COVID-19_Mild/Remission_blood","COVID-19_Severe_blood",
                                                      "Healthy Control_lung","COVID-19_lung"))
colors_condition = c("#87CEEB","#FFE4C4","#99004C","#6495ED","#FF1493")
DimPlot(DC_integrated,group.by = "disease.group_v2",pt.size = 1,cols = colors_condition) + 
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "#FCB7B3",fill=NA,size=5)) + ggtitle("")

