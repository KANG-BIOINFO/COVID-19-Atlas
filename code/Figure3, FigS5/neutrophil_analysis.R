# neutrophil integration & visualization
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)
library(EnhancedVolcano)

setwd("./covid19 Atlas/PBMC Atlas/neutrophil integration//")

# 1. data loading and extraction
# get neutrophils from PBMC
df_cells = read.table("../ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2_diseaseFixed_myeloidReformatted.txt",
                      sep="\t",header=1,row.names = 1)
neutrophil_cells = rownames(df_cells[df_cells$Cell.group=="Neutrophil",])

raw_files = c("Covid-19_PBMC_Schulte-Schrepping et al._Cell_10X_raw.h5seurat",
              "Covid-19_PBMC_Wilk et al._Nature Medicine_raw.h5seurat")
Dataset_names = c("Schulte-Schrepping","Wilk")
Neu_list = list()
for (i in 1:length(raw_files))
{
  raw_file = raw_files[i]
  sample_name = strsplit(raw_file,".h5")[[1]][1]
  print(sample_name)
  
  Neu_data = LoadH5Seurat(paste0("../raw data/",raw_file))
  overlap_cells = intersect(colnames(Neu_data), neutrophil_cells)
  Neu_data = Neu_data[,overlap_cells]
  df_specific = df_cells[overlap_cells,]
  
  Neu_data@meta.data = df_specific
  Neu_list[[i]] = Neu_data
}
# get neutrophils from BAL
df_bal_cells = read.table("../../BAL Atlas/integration/ctype_BAL_integrated_v1.txt",
                          sep="\t",header=1,row.names = 1)
bal_neutrophil_cells = rownames(df_bal_cells[df_bal_cells$integrated_snn_res.1=="21",]) # 806 neutrophils from Liao et al.
Neu_data = LoadH5Seurat("../../BAL Atlas/raw data/Covid-19_BAL_Liao et al._Nature Medicine_raw.h5seurat")
overlap_cells = intersect(colnames(Neu_data), bal_neutrophil_cells)
Neu_data = Neu_data[,overlap_cells]
df_bal_cells_neu = df_bal_cells[overlap_cells,]
Neu_data@meta.data = df_bal_cells_neu

Neu_list[[3]] = Neu_data


# 2. normalization and variable genes
Neu_list = lapply(X = Neu_list, FUN = function(x){
  x = NormalizeData(x); x = FindVariableFeatures(x)
})

# 3. select features and run pca
features = SelectIntegrationFeatures(object.list = Neu_list)
Neu_list = lapply(X = Neu_list, FUN = function(x){
  x = ScaleData(x, features = features); x = RunPCA(x, features = features)
})

# 4. anchors and integration
anchors = FindIntegrationAnchors(object.list = Neu_list,reduction = "rpca", dim=1:30)
Neu_integrated = IntegrateData(anchorset = anchors, dims = 1:30)

# 5. Visualization
Neu_integrated = ScaleData(Neu_integrated)
Neu_integrated = RunPCA(Neu_integrated, npcs=30)
Neu_integrated = RunUMAP(Neu_integrated, reduction = "pca",dims=1:30)

# 6. Mapping annotations
Neu_integrated$Cell.class = mapvalues(Neu_integrated$Cell.class, from=c("B","Macrophages_FABP4+","Macrophages_FCN1 high","Macrophages_FCN1+SPP1+","mDC"),
                                                             to = rep("Neutrophil_Lung",5))

Neu_integrated = FindNeighbors(Neu_integrated, dims = 1:30)
Neu_integrated = FindClusters(Neu_integrated,resolution = 1)
Neu_integrated = FindClusters(Neu_integrated,resolution = 0.5)
Neutro_0.5 = FindAllMarkers(Neu_integrated, assay = "RNA",slot = "data", "test.use" = "t")

saveRDS(Neu_integrated,"Neutrophil_integrated.RDS")

# 7. Visualization
DimPlot(Neu_integrated,group.by = "Cell.class_v2",label=T,label.size = 10) + 
  theme(legend.position = "none",
        axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(Neu_integrated,group.by = "integrated_snn_res.0.5",label=T,label.size = 12) + 
  theme(legend.position = "none",
        axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(Neu_integrated,group.by = "Dataset_Class") + 
  theme(legend.position = c(0.03,0.15),
        legend.text = element_text(family = "Helvetica", size=12),
        axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DimPlot(Neu_integrated,group.by = "disease.group_v2") + 
  theme(legend.position = c(0.03,0.15),
        legend.text = element_text(family = "Helvetica", size=12),
        axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

Neu_integrated$disease.group_v2 = mapvalues(Neu_integrated$disease.group_v2, from = NA, to = "COVID-19_Lung")
colors_disease = c("#FF6666","#FFCC99","#FFCCE5","#3399FF")
DimPlot(Neu_integrated,group.by = "disease.group_v2", cols = colors_disease) +
  theme(legend.position = c(0.03,0.15),
        legend.text = element_text(family = "Helvetica", size=12),
        axis.text = element_text(size=35,family="Helvetica"),
        axis.title = element_text(size=35,family="Helvetica"))

DefaultAssay(Neu_integrated) = "RNA"
FeaturePlot(Neu_integrated,features = c("MKI67","MPO","CXCR1","CXCR2","IFIT1","IFITM1","S100A8","S100A9","CXCL1"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("FCN1","CD14"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("S100A8","S100A9"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("CEACAM1","CEACAM6","CEACAM8"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("ELANE","DEFA4","PRTN3"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("CCL2","CCL3","CCL4"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("HLA-DRA","HLA-DMA","HLA-DRB1"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("CD14","CD68"),cols = c("#F5F5F5","#901026"))

FeaturePlot(Neu_integrated,features = c("ARG1"),cols = c("#F5F5F5","#901026"))
FeaturePlot(Neu_integrated,features = c("CD14","FCN1","S100A8",
                                        "ELANE","DEFA4","PRTN3",
                                        "IFIT1","IFITM1","OAS1",
                                        "MKI67","MPO","CXCR1"),cols = c("#F5F5F5","#901026"),ncol=4)

Idents(neutrophil) = "integrated_snn_res.0.5"
FeaturePlot(neutrophil,features = c("ELANE","DEFA4","MKI67","MPO",
                                    "CXCL1","CXCR1","IFIT1","IFITM1"),cols = c("#F5F5F5","#901026"),ncol=4)

# markers
Cluster2_markers = FindMarkers(Neu_integrated, slot = "data", ident.1 = 2, ident.2 = 4, only.pos = T)
Cluster4_markers = FindMarkers(Neu_integrated, slot = "data", ident.1 = 4, ident.2 = 2, only.pos = T)

# remove cluster 5,6
neutrophil = neutrophil[,! neutrophil$integrated_snn_res.0.5 %in% c(5,6)]
neutrophil$integrated_snn_res.0.5 = as.factor(as.character(neutrophil$integrated_snn_res.0.5))
neutrophil$integrated_snn_res.0.5 = mapvalues(neutrophil$integrated_snn_res.0.5, from = c(0,1,2,3,4),
                                              to = c("Neu_0","Neu_1","Neu_2","Neu_3","Neu_4"))
DimPlot(neutrophil, group.by = "integrated_snn_res.0.5", label = T, label.size = 10,pt.size = 2) + 
        theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

neutrophil$Cell.class = mapvalues(neutrophil$Cell.class, from = c("immature Neutrophil"),to=c("immature Neu"))
DimPlot(neutrophil, group.by = "Cell.class", label = T, label.size = 10,pt.size = 2,
        cols =c("#FF6666","#FF8800","#999900")) + 
        theme(legend.position = "none",
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "#0270C0",fill=NA,size=5)) + ggtitle("")

c("#FF6666","#FFCC99","#FFCCE5","#3399FF")
neutrophil$disease.group_v2 = mapvalues(neutrophil$disease.group_v2,
                                        from = c("Healthy/Control", "COVID-19 Moderate/Mild/Early stage/NonVent", "COVID-19 Severe/Late stage/Vent",NA),
                                        to = c("Healthy", "Mild/Remission","Severe","Severe"))
DimPlot(neutrophil, group.by = "disease.group_v3",pt.size = 2, cols = c("#3399FF","#FFCC99","#FF6666")) + 
        theme(legend.position = c(0.05, 0.15),
              legend.text = element_text(family = "Helvetica", size = 35),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# main figure
DefaultAssay(neutrophil) = "RNA"
FeaturePlot(neutrophil, features = c("FCGR3B","CEACAM1","MMP8",
                                     "MKI67","ELANE","DEFA4",
                                     "IFITM1","IFIT1"), cols = c("#F5F5F5","#901026"),ncol=4)
# supplementary figure
FeaturePlot(neutrophil, features = c("FCGR3B","CEACAM1","MMP8",
                                     "MKI67","ELANE","DEFA4", "MPO",
                                     "IFITM1","IFIT1", "S100A8","S100A9","MMP9","CEACAM8","PADI4","AZU1","RNASE2",
                                     "PRTN3","TLR5","CSF3R","IL4R","ARG1","CXCL1","CXCL8","CXCR2",
                                     "TNF","CCL2","CCL3"), cols = c("#F5F5F5","#901026"),ncol=9)

# new supplementary figure
FeaturePlot(neutrophil, features = c("S100A8","S100A9","MMP9","CEACAM8","PADI4","AZU1","RNASE2",
                                     "MPO", "PRTN3","TLR5","TLR6","TLR10",
                                     "CSF3R","IL4R","ARG1","CXCL1","CXCL8","CXCR2",
                                     "TNF","CCL2","CCL3"), cols = c("#F5F5F5","#901026"),ncol=7)

# volcano plot
neu2_vs_neu4 = read.table("./DEGs/DEGs_Neu2_vs_Neu4.txt", sep="\t", header=1)
neu2_vs_neu4[neu2_vs_neu4$pvals_adj < 1e-75,"pvals_adj"] = 1e-75
neu2_vs_neu4$names = as.character(neu2_vs_neu4$names)
EnhancedVolcano(neu2_vs_neu4, lab = neu2_vs_neu4$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("CEACAM1","FCN1","IGITM2","IL1R2","MMP8","MMP9","S100A8","S100A9","NCF1","S100A8","S100A9",
                              "CEACAM8","DEFA4","AZU1","ELANE","MPO","MKI67","PRTN3","CEACAM6","CDK6","CCNA2","RNASE2",
                              "HP", "PGLYRP1","CYBB","CAMP","LTF","LCN2"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "Neu4 vs. Neu2", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-15,15))

# volcano plot for Neu1,2 vs. Neu 3
neu3_vs_neu1_2 = read.table("./DEGs/DEGs_Neu3_vs_Neu0,1.txt", sep="\t", header=1)
neu3_vs_neu1_2[neu3_vs_neu1_2$pvals_adj < 1e-100,"pvals_adj"] = 1e-100
neu3_vs_neu1_2$names = as.character(neu3_vs_neu1_2$names)
EnhancedVolcano(neu3_vs_neu1_2, lab = neu3_vs_neu1_2$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("MMP9","S100A8","S100A9","S100A8","S100A9","CCL2","CCL3","CCL4","CCL8","CCR1",
                              "CXCL10","CXCL2","CXCL8","CXCL16","CXCR4","FTL","HLA-DRA","IFI27","IFI30","IFIH1",
                              "IFIT1","IFIT2","IFIT5","ISG15","ISG20","OAS1","SLC3A2",
                              "CR1","MME","MMP9","PTPRC","LYZ"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "Neu1,2 vs. Neu3", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-10,11))

# Neu2 barplot
df_neu2 = read.table("./DEGs/Neu2_enrichment.txt",sep="\t",header=1)
barplot(df_neu2$Gene.Enrichment.Score,horiz=T, col = "#48BA7B")

# Neu4 barplot
df_neu4 = read.table("./DEGs/neu4_enrichment.txt",sep="\t",header=1)
barplot(df_neu4$Gene.Enrichment.Score,horiz=T, col = "#B587A7")


# new genes from Bruce (3.31)
neutrophil = readRDS("./Neutrophil_integrated_v2.RDS")
gene_list = c("HP","PGLYRP1", "ATP5F1E", "CYBB", "CAMP", "HMGN2", "LTF", "LCN2", 
              "RETN", "MMP8", "ANXA1", "CD24", "NKG7")
DefaultAssay(neutrophil) = "RNA"
FeaturePlot(neutrophil, features = gene_list, cols = c("#F5F5F5","#901026"))

DimPlot(neutrophil, group.by = "Dataset",pt.size = 2, cols = c("#3399FF","#FFCC99","#CCCC00")) + 
        theme(legend.position = c(0.05, 0.15),
              legend.text = element_text(family = "Helvetica", size = 30),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# Neu0,1 barplot
df_neu0_1 = read.table("./DEGs/Neu0,1_enrichment.txt",sep="\t",header=1)
barplot(df_neu0_1$Gene.Enrichment.Score,horiz=T, col = "#F78A22")

# Neu4 barplot
df_neu3 = read.table("./DEGs/Neu3_enrichment.txt",sep="\t",header=1)
barplot(df_neu3$Gene.Enrichment.Score,horiz=T, col = "#00C4F1")
