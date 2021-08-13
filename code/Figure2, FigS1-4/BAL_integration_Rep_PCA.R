# BAL integration
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("./covid19 Atlas/BAL Atlas/integration//")

raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Grant","Liao")
BAL_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  BAL_data = LoadH5Seurat(raw_file)
  BAL_data$Dataset = Dataset_names[i]
  
  BAL_list[[i]] = BAL_data
}

# 2. normalization and variable genes
BAL_list = lapply(X = BAL_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = BAL_list)
BAL_list = lapply(X = BAL_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = BAL_list,reduction = "rpca", dim=1:30)
BAL_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. Visualization
BAL_integrated = ScaleData(BAL_integrated)
BAL_integrated = RunPCA(BAL_integrated, npcs=30)
BAL_integrated = RunUMAP(BAL_integrated, reduction = "pca",dims=1:30)

BAL_integrated = FindNeighbors(object = BAL_integrated, dim = 1:30)
BAL_integrated = FindClusters(BAL_integrated, resolution = 0.5)
BAL_integrated = FindClusters(BAL_integrated, resolution = 1)
BAL_integrated = FindClusters(BAL_integrated, resolution = 2)

BAL_integrated$disease.group.g = mapvalues(BAL_integrated$disease.group,from = c("COVID-19","Control","Healthy Control","Severe","Mild"),
                                           to = c("Severe","Control","Control","Severe","Mild"))

DimPlot(BAL_integrated,group.by = "Dataset", cols = c("#E0E0E0","#FFB266"))+
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(BAL_integrated,group.by = "Cell.class",label=T) + theme(legend.position = "")
DimPlot(BAL_integrated,group.by = "Cell.class",label=T)
DimPlot(BAL_integrated,group.by = "disease.group.g",cols = c("#66B2FF","#FF9999"))+
      theme(axis.text = element_text(size=35,family="Helvetica"),
            axis.title = element_text(size=35,family="Helvetica"))

# A new round 
BAL_integrated = BAL_integrated[,! BAL_integrated$Cell.class_v2 %in% c("undefined", "undefined Epi")]
BAL_integrated$Cell.class_v2 = as.factor(as.character(BAL_integrated$Cell.class_v2))
BAL_integrated$Cell.class_v3 = factor(BAL_integrated$Cell.class_v2, 
                                      levels = c("AT1/AT2","Basal/Club","Ciliated","transitional Epi",
                                                 "B cell","Plasmablast",
                                                 "CD4+ T_1","CD4+ T_2","CD8+ T","NK","proliferating T/NK",
                                                 "Neutrophil", "MoAM1","MoAM2","MoAM3","MoAM4","MoAM5","transitional Macro","TRAM1","TRAM2","TRAM3","cDC","pDC","proliferating Myeloid cells"
                                                 ))
colors_class = c("#FFD6C7","#FFFF66","#FF6666","#FFA500",
                 "#B0E0E6","#4682B4","#1E90FF","#00008B",
                 "#FF3399","#D8BFD8","#800080","#9933ff",
                 "#FFFACD","#FFFF33","#336600","#FF69B4",
                 "#D2691E","#B266FF","#CD853F",
                 "#B0C4DE","#CD5C5C","#DEA4A4","#BDB76B","#808000")#,
                 #"#9ACD32","#0BAF0B")
DimPlot(BAL_integrated, group.by = "Cell.class_v3",cols = colors_class) + 
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        axis.line = element_line(colour = "black",size=2))+ggtitle("")
