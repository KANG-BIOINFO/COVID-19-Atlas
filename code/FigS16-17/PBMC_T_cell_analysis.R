library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

setwd("./covid19 Atlas/PBMC Atlas/T cells//")

# remove sepsis cells and influenza cells
# T_cell_integrated = readRDS("./T_cell_integrated_Rep_PCA_all_Tcells.RDS")
SS_meta = T_cell_integrated[,T_cell_integrated$Dataset=="Schulte-Schrepping"]@meta.data
Lee_meta = T_cell_integrated[,T_cell_integrated$Dataset=="Lee"]@meta.data
sepsis_cells = rownames(SS_meta[SS_meta$orig.ident == "sepsis",])
influenza_cells = rownames(Lee_meta[Lee_meta$disease == "influenza",])
removal_cells = c(sepsis_cells, influenza_cells)

# integrate datasets
setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/raw data/")
df_T_cells = read.table("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/T cells/ctype_T_cell_integrated_clean_sampleAdded.txt",
                        sep="\t", header=1, row.names = 1)
df_T_cells = df_T_cells[setdiff(rownames(df_T_cells), removal_cells),]

raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Arunachalam","Guo","Lee","Schulte-Schrepping","Wilk")
T_cell_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  T_cell_data = LoadH5Seurat(raw_file)
  T_cell_data$Dataset = Dataset_names[i]
  
  common_T_cells = intersect(colnames(T_cell_data), rownames(df_T_cells))
  T_cell_data = T_cell_data[,common_T_cells]
  print(sample_name); print(length(common_T_cells))
  
  T_cell_list[[i]] = T_cell_data
}

# 2. normalization and variable genes
T_cell_list = lapply(X = T_cell_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = T_cell_list)
T_cell_list = lapply(X = T_cell_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = T_cell_list,reduction = "rpca", dim=1:30)
T_cell_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. Visualization
T_cell_integrated = ScaleData(T_cell_integrated)
T_cell_integrated = RunPCA(T_cell_integrated, npcs=30)
T_cell_integrated = RunUMAP(T_cell_integrated, reduction = "pca",dims=1:30)

T_cell_integrated = FindNeighbors(T_cell_integrated, dims = 1:30)
T_cell_integrated = FindClusters(T_cell_integrated,resolution = 2)
T_cell_integrated = FindClusters(T_cell_integrated,resolution = 1)
T_cell_integrated = FindClusters(T_cell_integrated,resolution = 0.5)


df_cells = read.table("ctype_T_cell_integrated_clean_sampleAdded_Reannotated.txt",sep="\t",header=1,row.names = 1)
df_cells = df_cells[colnames(T_cell_integrated),]
for (i in colnames(df_cells))
{
  T_cell_integrated@meta.data[[i]] = df_cells[[i]]
}
saveRDS(T_cell_integrated, file = "T_cell_integrated_Rep_PCA_Tcells_COVID-19_only_multiplets_removed.RDS")

# 7. visualization
T_cell_integrated$Cell.class_v2 = as.factor(as.character(T_cell_integrated$Cell.class_v2))

colors_class = c("#000000","#A63F48","#CEE8F2","#FF7F50","#EA979A","#F5DEB3","#4174D9",
           "#FFD700","#DC143C","#20B2AA","#D8BFD8","#800080","#BA55D3","#FFC0CB",
           "#90EE90")
DimPlot(T_cell_integrated, group.by = "Cell.class_v2",cols = colors_class) +
          theme(legend.position = "none",
                legend.text = element_text(family = "Helvetica",size=25),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

colors_condition = c("#FFE4C4","#FFD700","#FA8072","#87CEEB")
DimPlot(T_cell_integrated, group.by = "disease.group_v2",cols=colors_condition) +
          theme(legend.position = "none",
                legend.text = element_text(family = "Helvetica",size=25),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

DimPlot(T_cell_integrated, group.by = "Dataset") +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.position = c(0.02,0.90),
        legend.text = element_text(family = "Helvetica",size=20))

DefaultAssay(T_cell_integrated) = "RNA"
FeaturePlot(T_cell_integrated, features = c("OAS1","PDCD1","IFNG","MKI67"),
            cols = c("#F5F5F5","#901026"))

FeaturePlot(T_cell_integrated, features = c("XAF1","IRF1","TP53","BCL2L11","CASP3"),
            cols = c("#F5F5F5","#901026"))


# dot plot
# organize metadata of T cells
T_cell_integrated_v2 = T_cell_integrated[,T_cell_integrated$disease.group_v2 != "COVID-19 Remission"]
T_cell_integrated_v2$disease.group_v2 = as.factor(as.character(T_cell_integrated_v2$disease.group_v2))
T_cell_integrated_v2$disease.group_v2 = mapvalues(T_cell_integrated_v2$disease.group_v2, from = c("Healthy/Control","COVID-19 Severe/Late stage/Vent","COVID-19 Moderate/Mild/Early stage/NonVent"),
                                                  to = c("H","S","M"))
T_cell_integrated_v2$Cell.class_v3 = paste0(T_cell_integrated_v2$Cell.class_v2,"_",T_cell_integrated_v2$disease.group_v2)
T_cell_integrated_v2 = T_cell_integrated_v2[, !T_cell_integrated_v2$Cell.class_v3 %in% c("CD4+ T activated_H","CD8+ T activated_H","NK activated_H")]
T_cell_integrated_v2$Cell.class_v3 = as.factor(as.character(T_cell_integrated_v2$Cell.class_v3))
cell_class_v3_order = read.table("./cell_class_v3.txt",sep="\t",header=1,row.names = 1)
cell_class_v3_order$color = paste0("#",cell_class_v3_order$color)

T_cell_integrated_v2$Cell.class_v3_ordered = factor(T_cell_integrated_v2$Cell.class_v3, levels = rownames(cell_class_v3_order))

shown_features = c("CD3D","CD3E","CD3G",
                   "CD4","CCR7","IL7R","ISG15","IFI44L","XAF1","IL2RA","CD27","FOXP3","CCR10",
                   "CD8A","CD8B","CCL5","NKG7","SLC4A10",
                   "GNLY","GZMA","GZMB","FCGR3A","SELL","TNFRSF18","MKI67","TYMS",
                   "TRDC","TRGC1","TRGC2","TRAC","TRBC1","TRBC2")
DotPlot(T_cell_integrated_v2,features = shown_features,
        group.by = c("Cell.class_v3_ordered"))+
        theme(panel.background = element_rect(fill=NA,colour = "black"),
              axis.title = element_blank(),
              axis.text.x = element_text(family = "Helvetica",angle = 45,vjust=0.5, size = 16),
              axis.text.y = element_text(family = "Helvetica", size = 16, colour = cell_class_v3_order$color)) +
        scale_colour_gradient2(low = "#004C99", mid = "#FFFFFF", high = "#750B22")





