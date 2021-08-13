# BAL integration
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/BAL Atlas/T cell///")
options(future.globals.maxSize=(90*1024^3)) # set max size as 40G

df_cells = read.table("../integration/ctype_BAL_integrated_v3.txt", sep="\t", header=1, row.names = 1)
df_T_cells = df_cells[df_cells$integrated_snn_res.1 %in% c(5,6,13,14,20),]
T_NK_cells = rownames(df_T_cells)

raw_files = list.files("../raw data//",pattern = ".h5seurat")
Dataset_names = c("Grant","Liao")
BAL_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  BAL_data = LoadH5Seurat(paste0("../raw data/",raw_file))
  BAL_data$Dataset = Dataset_names[i]
  
  overlap = intersect(T_NK_cells, colnames(BAL_data))
  BAL_data = BAL_data[, overlap]
  
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

saveRDS(BAL_integrated, "T_integrated.RDS")

DefaultAssay(T_integrated) = "RNA"
FeaturePlot(T_integrated, features = c("MKI67","CD3D","CD3E","CD4","CD8A","CD8B","FOXP3","NKG7","PDCD1"))
write.table(T_integrated@meta.data, file = "T_subclusters_v2.txt",sep="\t")
write.table(T_integrated@reductions$umap@cell.embeddings, file = "umap.txt",sep="\t")

# 6. visualization
df_cellxgene = read.table("./ctype_T_cellxgene_0210.csv", sep=",", header=1, row.names = 1)
df_cellxgene = df_cellxgene[colnames(T_integrated),]
T_integrated$Cell.class_v2 = df_cellxgene$Cell_class_v2

DimPlot(T_integrated, group.by = "Cell.class_v2",pt.size = 0.5) +
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(T_integrated, group.by = "Dataset",pt.size = 0.5) +
        theme(legend.position = c(0.75,0.2),
              legend.text = element_text(family = "Helvetica", size = 28),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

T_integrated$disease.group_v2 = mapvalues(T_integrated$disease.group,
                                          from = c("Healthy Control","COVID-19"),
                                          to = c("Control","Severe"))
DimPlot(T_integrated, group.by = "disease.group_v2",pt.size = 0.5, cols = c("#5062AA","#F2D7A6","#EC3284")) +
        theme(legend.position = c(0.68,0.2),
              legend.text = element_text(family = "Helvetica", size = 28),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DefaultAssay(T_integrated) = "RNA"
exhausted_t_markers = c("PDCD1","LAG3","HAVCR2","CD160",
                        "TPM2", "CTLA4", "TIGIT", "SNX9", "RBPJ")
exhausted_t_markers = c("PDCD1","HAVCR2","CTLA4","LAG3")
isg_markers = c("OAS1","IFIT1","IFIT3","OAS3")

FeaturePlot(T_integrated, features = exhausted_t_markers, cols = c("#F5F5F5","#901026"))
FeaturePlot(T_integrated, features = isg_markers, cols = c("#F5F5F5","#901026"))

VlnPlot(T_integrated, group.by = "Cell.class_v2", features = c("CXCL13")) + 
  ylab("Scaled Expression") + xlab("Cell class") +
  theme(legend.position = "none",
        plot.title = element_text(family = "Arial", size = 30),
        axis.text = element_text(family = "Arial", size = 25),
        axis.title = element_text(family = "Arial", size = 25))

VlnPlot(T_integrated, group.by = "disease.group_v2", features = c("CXCL13")) + 
  ylab("Scaled Expression") + xlab("Disease Condition") +
  theme(legend.position = "none",
        plot.title = element_text(family = "Arial", size = 30),
        axis.text = element_text(family = "Arial", size = 25),
        axis.title = element_text(family = "Arial", size = 25))

