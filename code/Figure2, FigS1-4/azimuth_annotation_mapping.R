library(Seurat)
library(SeuratDisk)

# Azimuth reference was used for cell label projection in each indivial PBMC and BAL dataset. 
# This is an example of Azimuth application in the data of Schulte-Schrepping et al. (Cell). Same code was used for other datasets.
setwd("./covid19 Atlas/Covid-19_PBMC_Schulte-Schrepping et al._Cell/")
# load reference data
reference = LoadH5Seurat("pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# load query data
pbmc = LoadH5Seurat("./10x_PBMC/raw_10x_PBMC.h5seurat")
# SCTransform
pbmc = SCTransform(pbmc)

# Find transfer anchors
anchors = FindTransferAnchors(
    reference = reference,
    query = pbmc,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
)

pbmc <- TransferData(
  anchorset = anchors, 
  reference = reference,
  query = pbmc,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT")
)
write.table(pbmc@meta.data,"pbmc_azimuth.txt",sep="\t")
pbmc = AddMetaData(object = pbmc,metadata = MappingScore(anchors = anchors),col.name="mapping.score")
write.table(pbmc@meta.data,"pbmc_azimuth_v2.txt",sep="\t")
