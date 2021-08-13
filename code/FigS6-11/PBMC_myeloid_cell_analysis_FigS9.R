library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

### integrate datasets
setwd("./covid19 Atlas/PBMC Atlas/raw data/")

df_myeloid_cells = read.table("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded.txt",
                        sep="\t", header=1, row.names = 1)
all_myeloid_cells = rownames(df_myeloid_cells[df_myeloid_cells$Cell.group %in% c("cDC","Monocyte","pDC","Neutrophil"),])

raw_files = list.files(".",pattern = ".h5seurat")
Dataset_names = c("Arunachalam","Guo","Lee","Schulte-Schrepping","Wilk")
myeloid_cell_list = list()
# 1. read data
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  myeloid_cell_data = LoadH5Seurat(raw_file)
  myeloid_cell_data$Dataset = Dataset_names[i]
  
  common_myeloid_cells = intersect(colnames(myeloid_cell_data), all_myeloid_cells)
  myeloid_cell_data = myeloid_cell_data[,common_myeloid_cells]
  print(sample_name); print(length(common_myeloid_cells))
  
  myeloid_cell_list[[i]] = myeloid_cell_data
}

# 2. normalization and variable genes
myeloid_cell_list = lapply(X = myeloid_cell_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = myeloid_cell_list)
myeloid_cell_list = lapply(X = myeloid_cell_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = myeloid_cell_list,reduction = "rpca", dim=1:30)
myeloid_cell_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. Visualization
myeloid_cell_integrated = ScaleData(myeloid_cell_integrated)
myeloid_cell_integrated = RunPCA(myeloid_cell_integrated, npcs=30)
myeloid_cell_integrated = RunUMAP(myeloid_cell_integrated, reduction = "pca",dims=1:30)

# add annotation
df_myeloid_cells = df_myeloid_cells[colnames(myeloid_cell_integrated),]
myeloid_cell_integrated$Cell.class = df_myeloid_cells$Cell.class
myeloid_cell_integrated$Cell.group = df_myeloid_cells$Cell.group
myeloid_cell_integrated$sample = df_myeloid_cells$sample
myeloid_cell_integrated$disease.group = df_myeloid_cells$disease.group

# add annotations again
df_new_myeloid = read.table("/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2.txt",
                            sep = "\t", header = 1, row.names = 1)
overlapping = intersect(rownames(df_myeloid_cells),colnames(myeloid_cell_integrated))
df_new_myeloid = df_new_myeloid[overlapping,]; myeloid_cell_integrated = myeloid_cell_integrated[, overlapping]
myeloid_cell_integrated$disease.group_v2 = df_new_myeloid$disease.group_v2

saveRDS(myeloid_cell_integrated,"/data/aronow/Kang/data for toppcell/covid19 Atlas/PBMC Atlas/Myeloid cells/myeloid_cell_integrated.RDS")

# draw UMAPs
DimPlot(myeloid_cell_integrated,group.by = "Cell.class") +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(myeloid_cell_integrated,group.by = "disease.group_v2") +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(myeloid_cell_integrated,group.by = "Dataset") +
  theme(axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

# dot plots
saveRDS(myeloid_subset,"myeloid_cell_integrated_subset.RDS")
Idents(myeloid_subset) = "Cell.class_v3"
myeloid_subset$Cell.class_v3 = factor(myeloid_subset$Cell.class_v3, levels = rev(c("immature Neu","Neu","pDC","cDC",
                                                                                   "cMono1","cMono2","cMono3","cMono4","cMono_influenza","ncMono")))
DotPlot(myeloid_subset,assay = "RNA",
        #col.min = -1.2,col.max = 2.4,
        features = c("HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5",
                     "IFITM1","IFITM3","IFIT1","ISG15","OAS1",
                     "CD14","FCN1","S100A8","S100A9","IL6","IL10","IL17A","CCL2","CXCL1","IFNG","IFNA1",
                     "MPO","MKI67"))+
  theme(panel.background = element_rect(fill=NA,colour = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(family = "Helvetica",angle = 45,vjust=0.5, size = 16),
        axis.text.y = element_text(family = "Helvetica", size = 16)) +
  scale_colour_gradient2(low = "#004C99", mid = "#FFFFFF", high = "#750B22")