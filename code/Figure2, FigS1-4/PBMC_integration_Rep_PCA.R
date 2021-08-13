library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("./covid19 Atlas/PBMC Atlas/raw data/")

raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Aru","Guo","Lee","Sch","Wil")
pbmc_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  pbmc_data = LoadH5Seurat(raw_file)
  pbmc_data$Dataset = Dataset_names[i]
  
  pbmc_list[[i]] = pbmc_data
}

# 2. normalization and variable genes
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

# 5. Visualization
pbmc_integrated = ScaleData(pbmc_integrated)
pbmc_integrated = RunPCA(pbmc_integrated, npcs=30)
pbmc_integrated = RunUMAP(pbmc_integrated, reduction = "pca",dims=1:30)
DimPlot(pbmc_integrated,group.by = "Dataset")
DimPlot(pbmc_integrated,group.by = "Cell.class")
DimPlot(pbmc_integrated,group.by = "Cell.class",label=T)

# 5.1 B cells
highlight_cells = colnames(pbmc_integrated[,pbmc_integrated$Cell.class %in% c("B","B cell","B cells_1","B cells_2","B cells_3")])
DimPlot(pbmc_integrated,group.by = "Cell.class", cells.highlight = highlight_cells)

# 5.2 cd4+ T cells
highlight_cells = colnames(pbmc_integrated[,pbmc_integrated$Cell.class %in% c("CD4 Memory T cell","CD4m T","CD4 Naive T cell","CD4n T",
                                                                              "CD4 T","CD4 T cell","CD4+ T cells_1","CD4+ T cells_2",
                                                                              "CD4+ T cells_3")])
DimPlot(pbmc_integrated,group.by = "Cell.class", cells.highlight = highlight_cells)

# 5.3 cd8+ T cells
highlight_cells = colnames(pbmc_integrated[,pbmc_integrated$Cell.class %in% c("CD8 Effective T cell","CD8 Effector T cell","CD8eff T",
                                                                              "CD8 Memory T cell","CD8m T","CD8+ T cells_1","CD8+ T cells_2",
                                                                              "CD8+ T cells_3")])
DimPlot(pbmc_integrated,group.by = "Cell.class", cells.highlight = highlight_cells)

# 5.4 plasmablast
highlight_cells = colnames(pbmc_integrated[,pbmc_integrated$Cell.class %in% c("Plasmablast","Plasmablasts","IgA PB","IgG PB")])
DimPlot(pbmc_integrated,group.by = "Cell.class", cells.highlight = highlight_cells)

# 5.5 Covid condition
DimPlot(pbmc_integrated, group.by = "disease.group")
pbmc_integrated$disease.group_concise = mapvalues(pbmc_integrated$disease.group, from = unique())


# 5.6 Others
df_all = read.table("../ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2_diseaseFixed_myeloidReformatted.txt",sep="\t",header=1,row.names = 1)
overlap_cells = intersect(rownames(df_all),colnames(pbmc_integrated))
df_all = df_all[overlap_cells,]
pbmc_integrated = pbmc_integrated[,overlap_cells]
pbmc_integrated@meta.data = df_all

pbmc_integrated = pbmc_integrated[,pbmc_integrated$Cell.class != "B cell"]
pbmc_integrated$Cell.class = mapvalues(pbmc_integrated$Cell.class, from = c("IgA PB","IgG PB"), to = c("Plasmablast","Plasmablast"))
pbmc_integrated$Cell.class = as.factor(as.character(pbmc_integrated$Cell.class))
pbmc_integrated$Cell.class = factor(pbmc_integrated$Cell.class,
                                    levels = c("B naive","B intermediate","B memory","Plasmablast",
                                               "CD4+ T naive","CD4+ CTL","CD4+ Tcm","CD4+ T activated",
                                               "CD8+ T naive","CD8+ Tem","CD8+ T activated","MAIT",
                                               "Treg","dn T","gd T","T/NK proliferative",
                                               "NK","NK CD56bright","NK activated",
                                               "Classical Monocyte","Non-classical Monocyte","Neutrophil","immature Neutrophil","cDC","pDC",
                                               "HSPC","Platelet","RBC"))
colors_class = c("#FFD6C7","#FFFF66","#FF6666","#FFA500",
                 "#B0E0E6","#4682B4","#1E90FF","#00008B",
                 "#FF3399","#D8BFD8","#800080","#9933ff",
                 "#D2691E","#FFFF33","#336600","#FF69B4",
                 "#FFFACD","#F4A460","#CD853F",
                 "#B0C4DE","#778899","#CD5C5C","#DEA4A4","#BDB76B","#808000",
                 "#9ACD32","#0BAF0B","#FF0000")
DimPlot(pbmc_integrated,group.by = "Cell.class",cols = colors_class) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=20))

DimPlot(pbmc_integrated,group.by = "Cell.class",cols = colors_class) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.position = "none")

DimPlot(pbmc_integrated,group.by = "Dataset",cols = colors_class) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=24,family = "Helvetica",face="bold"))

colors_disease = c("#4169E1","#F5DEB3","#FFA500","#FF1493")
pbmc_integrated = pbmc_integrated[,pbmc_integrated$disease.group_v2!="Influenza Severe"]
pbmc_integrated$disease.group_v2 = as.factor(as.character(pbmc_integrated$disease.group_v2)) 
pbmc_integrated$disease.group_v2 = factor(pbmc_integrated$disease.group_v2, levels = c("Healthy/Control","COVID-19 Moderate/Mild/Early stage/NonVent",
                                                                                       "COVID-19 Remission","COVID-19 Severe/Late stage/Vent"))
DimPlot(pbmc_integrated,group.by = "disease.group_v2",cols = colors_disease) +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"),
        legend.text = element_text(size=12,family = "Helvetica",face="bold"))

DefaultAssay(pbmc_integrated) = "RNA"
FeaturePlot(pbmc_integrated, features = c("IFNG","TNF","CSF2","IL6", "IL2", "IL10", "IL1B","CSF3","PDCD1"),
            cols = c("#F5F5F5","#901026"),ncol = 3)
saveRDS(pbmc_integrated,file="/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/integration/Integrated_RepPCA_v3.RDS")
