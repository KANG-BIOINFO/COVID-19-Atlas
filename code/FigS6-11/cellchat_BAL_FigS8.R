library(CellChat)
library(ggalluvial)
library(svglite)
library(ggplot2)
library(patchwork)
library(plyr)

# Generate cellchat objects
bal.integrated = readRDS("./integration/BAL_integrated.RDS")
bal.liao = bal.integrated[,bal.integrated$Dataset == "Liao"]

for (condition in c("Healthy Control","Mild","Severe"))
{
  specific_bal = bal.liao[,bal.liao$disease.group == condition]
  for (col in colnames(specific_bal@meta.data))
  {
    if (! col %in% c("nCount_RNA","nFeature_RNA"))
    {
      specific_bal@meta.data[[col]] = as.factor(as.character(specific_bal@meta.data[[col]]))
    }
  }
  
  # 1. load cell chat data
  data.input <- GetAssayData(specific_bal, assay = "RNA", slot = "data") # normalized data matrix
  labels <- specific_bal$Cell.class
  identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat <- createCellChat(data = data.input)
  cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  
  # 2. set interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
  # showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  # dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  # 3. preprocessing the expression table
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  saveRDS(cellchat, paste0("cellchat_", condition, "_bal.RDS"))
}

# load cellchat objects
bal_healthy = readRDS("./cellchat_Healthy Control_bal.RDS")
bal_mild = readRDS("./cellchat_Mild_bal.RDS")
bal_severe = readRDS("./cellchat_Severe_bal.RDS")

cellchat = bal_healthy

# visualization
groupSize = as.numeric(table(cellchat@idents))
pathways.show <- c("CCL") 
vertex.receiver = seq(1,9) # a numeric vector
# Hierarchy plot
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- c("CCL") 
netVisual_aggregate(bal_healthy, signaling = pathways.show, layout = "circle", vertex.size = groupSize_healthy)
netVisual_aggregate(bal_severe, signaling = pathways.show, layout = "circle", vertex.size = groupSize_severe)


cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(cellchat, signaling = pathways.show)

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
# incoming
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

# 4.3 Manifold and classification learning of signaling networks
cellchat = computeNetSimilarity(cellchat, type = "functional")
cellchat = netEmbedding(cellchat, type = "functional")
cellchat = netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional")
netVisual_embeddingZoomIn(cellchat, type = "functional")

cellchat = computeNetSimilarity(cellchat, type = "structural")
cellchat = netEmbedding(cellchat, type = "structural")
cellchat = netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural")
netVisual_embeddingZoomIn(cellchat, type = "structural")

### integrate
cellchat.obj1 <- readRDS("cellchat_Healthy Control_bal.RDS")
cellchat.obj2 <- readRDS("cellchat_Severe_bal.RDS")
cellchat <- mergeCellChat(list(cellchat.obj1, cellchat.obj2), add.names = c("Healthy","Severe"))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization
netVisual_embeddingPairwise(cellchat, type = "structural")
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural")

rankSimilarity(cellchat, type = "structural")
rankNet(cellchat, mode = "comparison")

