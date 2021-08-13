library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("./covid19 Atlas/PBMC Atlas/Platelet/")

# 1. load data 
df_all = read.table("../ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2_diseaseFixed_myeloidReformatted_Bcell-InfluenzaExcluded_concise_1129_lineageAdded_subclusterAdded.txt",
                    sep="\t",header=1,row.names = 1)
platelet_cells = rownames(df_all[df_all$Cell.class_concise == "Platelet",])

raw_files = list.files("../raw data/",pattern = ".h5seurat")
Dataset_names = c("Aru","Guo","Lee","Sch","Wil")
platelet_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  platelet_data = LoadH5Seurat(paste0("../raw data/",raw_file))
  platelet_data$Dataset = Dataset_names[i]
  
  platelet_shared = intersect(platelet_cells, colnames(platelet_data))
  platelet_data = platelet_data[,platelet_shared]
  
  platelet_list[[i]] = platelet_data
}

# 2. normalization and variable genes
platelet_list = lapply(X = platelet_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = platelet_list)
platelet_list = lapply(X = platelet_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = platelet_list,reduction = "rpca", dim=1:30)
platelet_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. pca and umap
platelet_integrated = ScaleData(platelet_integrated)
platelet_integrated = RunPCA(platelet_integrated, npcs=30)
platelet_integrated = RunUMAP(platelet_integrated, reduction = "pca",dims=1:30)

# 6. clustering
df_platelet = df_all[colnames(platelet_integrated),]
platelet_integrated@meta.data = df_platelet

platelet_integrated = FindNeighbors(platelet_integrated, dims = 1:30)
platelet_integrated = FindClusters(platelet_integrated, resolution = 1)
platelet_integrated = FindClusters(platelet_integrated, resolution = 2)
platelet_integrated = FindClusters(platelet_integrated, resolution = 0.5)
platelet_integrated = FindClusters(platelet_integrated, resolution = 0.1)

write.table(platelet_integrated@meta.data, "ctype_platelet_v1.txt",sep="\t")
write.table(platelet_integrated@reductions$umap@cell.embeddings, "umap_v1.txt",sep="\t")

# remove multiplets
platelet_integrated_v2 = platelet_integrated[,platelet_integrated$integrated_snn_res.0.1 %in% c(0,1,2)]
platelet_integrated_v2 = RunPCA(platelet_integrated_v2, npcs=30)
platelet_integrated_v2 = RunUMAP(platelet_integrated_v2, reduction = "pca",dims=1:30)
platelet_integrated_v2 = FindNeighbors(platelet_integrated_v2, dims = 1:30)
platelet_integrated_v2 = FindClusters(platelet_integrated_v2, resolution = 1)
platelet_integrated_v2 = FindClusters(platelet_integrated_v2, resolution = 0.5)
platelet_integrated_v2 = FindClusters(platelet_integrated_v2, resolution = 0.3)

saveRDS(platelet_integrated_v2, "platelet.RDS")
write.table(platelet_integrated_v2@meta.data, "ctype_platelet_v2.txt",sep="\t")
write.table(platelet_integrated_v2@reductions$umap@cell.embeddings, "umap_v2.txt",sep="\t")

# 7. visualization
platelet_integrated_v2$Cluster = paste0("PLT_",platelet_integrated_v2$integrated_snn_res.0.3)

DimPlot(platelet_integrated_v2, group.by = "disease.group_v4",pt.size=1,cols = c("#4169E1","#F5DEB3", "#C68C8D", "#FF1493"))+
  theme(legend.position = c(0.50,0.86),
        legend.text = element_text(family = "Helvetica",size=25),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")
DimPlot(platelet_integrated_v2, group.by = "Cluster",label=T,pt.size = 1,label.size = 8.5) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DefaultAssay(platelet_integrated_v2) = "RNA"
# inteacting genes
FeaturePlot(platelet_integrated_v2, features = c("ITGA2B","ITGB3","ITGB5","CCL5", "GP1BA", "GP9",
                                                 "JAM3", "PDGF4", "PF4", "PPBP", "SELP", "TGFB1"),
            cols = c("#F5F5F5","#901026"),ncol=4)

# integrins
integrins = c("ITGA1","ITGA2","ITGA3","ITGA4","ITGA5","ITGA6","ITGA7","ITGA8","ITGA9","ITGA10","ITGA11",
              "ITGAD","ITGAE","ITGAL","ITGAM","ITGAV","ITGA2B","ITGAX","ITGB1","ITGB2","ITGB3","ITGB4",
              "ITGB2","ITGB3","ITGB4","ITGB5","ITGB6","ITGB7","ITGB8")
FeaturePlot(platelet_integrated_v2, features = integrins, cols = c("#F5F5F5","#901026"), ncol = 6)

VlnPlot(platelet, group.by = "Cluster", features = c("ITGA2B"))
VlnPlot(platelet, group.by = "Cluster", features = c("ITGA2B","ITGB1","ITGB3","ITGB5"))
DefaultAssay(platelet) = "RNA"
DotPlot(platelet, group.by = "Cluster",features = c("ITGA2B","ITGB1","ITGB3","ITGB5")) + 
          theme(panel.background = element_rect(fill=NA,colour = "black"),
                axis.title = element_blank(),
                axis.text.x = element_text(family = "Helvetica",angle = 45,vjust=0.5, size = 16),
                axis.text.y = element_text(family = "Helvetica", size = 16)) +
          scale_colour_gradient2(low = "#004C99", mid = "#FFFFFF", high = "#750B22")

DotPlot(platelet, group.by = "Cluster",features = c("ITGA2B","ITGB1","ITGB3","ITGB5","PF4","PPBP","GP9",
                                                    "TGFB1","THBS1","JAM3","PDGFA","SELP",
                                                    "CCL5","GP1BA","PDGFB")) + 
        theme(panel.background = element_rect(fill=NA,colour = "black"),
              axis.title = element_blank(),
              axis.text.x = element_text(family = "Helvetica",angle = 45,vjust=1, hjust = 1, size = 16),
              axis.text.y = element_text(family = "Helvetica", size = 16)) +
        scale_colour_gradient2(low = "#004C99", mid = "#FFFFFF", high = "#750B22")

DotPlot(platelet, group.by = "disease.group_v3",features = c("ITGA2B","ITGB1","ITGB3","ITGB5","PF4","PPBP","GP9",
                                                    "TGFB1","THBS1","JAM3","PDGFA","SELP",
                                                    "CCL5","GP1BA","PDGFB")) + 
        theme(panel.background = element_rect(fill=NA,colour = "black"),
              axis.title = element_blank(),
              axis.text.x = element_text(family = "Helvetica",angle = 45,vjust=1, hjust = 1, size = 16),
              axis.text.y = element_text(family = "Helvetica", size = 16)) +
        scale_colour_gradient2(low = "#004C99", mid = "#FFFFFF", high = "#750B22")

VlnPlot(platelet, group.by = "Cluster",features = c("ITGA2B","ITGB1","ITGB3","ITGB5","PF4","PPBP","GP9",
                                                    "TGFB1","THBS1","JAM3","PDGFA","SELP",
                                                    "CCL5","GP1BA","PDGFB"))

