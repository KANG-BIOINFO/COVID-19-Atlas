library(plyr)
library(ggplot2)
require(gridExtra)
library(ggpubr)
library(tibble)
library(stats)

df_all = read.table("GSE158127_04all-lineages-metadata_NU_Parenchyma_only_macrophageFixed_Adjusted.txt",sep="\t",header=1,row.names = 1)
df_all$Dataset = "Ankit"

pval_generator <- function(df = df_dataset, celltype = "cDC", 
                           celltype_col = "Cell.type", 
                           condition_col = "Diagnosis", 
                           donor_col = "orig.ident")
{
  # compute the numerate for each Donor and mapping relations
  # for (col in colnames(df)){df[[col]] = as.factor(as.character(df[[col]]))}
  donor_totalCells = table(df[[donor_col]])
  donor2condition = df[,condition_col]
  names(donor2condition) = df[,donor_col]
  
  
  # get the specific cell type
  df_celltype = df[df[[celltype_col]] == celltype,]
  
  # compute cell distribution across donors in that cell type
  donor_specificCells = table(df_celltype[[donor_col]])
  
  # get the dataframe
  donor_pecentage = c()
  condition_vals = c()
  for (d in names(donor_specificCells))
  {
    donor_pecentage[d] = ((donor_specificCells[d] + 0.01) / donor_totalCells[d]) * 100
    condition_vals = append(condition_vals, as.character(donor2condition[d]))
  }
  donor_pecentage = as.data.frame(donor_pecentage)
  
  donor_pecentage$Condition = condition_vals
  donor_pecentage$Donor = rownames(donor_pecentage)
  names(donor_pecentage)[1] <- "Proportion"
  
  # iterate for each comparison
  Comparisons <- combn(unique(donor_pecentage$Condition), 2, simplify=FALSE)
  compRecorder <- c()
  pvalRecoder <- c()
  k = 0
  for (i in 1:length(Comparisons)) {
    pair <- Comparisons[[i]]
    tmp <- rbind(donor_pecentage[donor_pecentage$Condition==pair[1], ], donor_pecentage[donor_pecentage$Condition==pair[2], ])
    tmp <- as.data.frame(tmp)
    tmp$Condition = factor(tmp$Condition, levels=unique(tmp$Condition))
    res <- wilcox.test(Proportion ~ Condition, data=tmp)
    p.val = res$p.value
    # add pvals into list
    k = k + 1
    # print(pair)
    # print(p.val)
    pair = paste(pair, collapse = " vs. ")
    print(pair)
    compRecorder = append(compRecorder, pair)
    pvalRecoder = append(pvalRecoder, p.val)
  }
  
  return(list(compRecorder, pvalRecoder))
}

compRecoder_all = c()
pvalRecoder_all = c()
datasetRecoder_all = c()
celltypeRecoder_all = c()

for (dataset in unique(df_all$Dataset)) 
{
  # get data per dataset
  df_dataset = df_all[df_all$Dataset == dataset,]
  print(dataset)
  for (col in colnames(df_dataset)){df_dataset[[col]] = as.factor(as.character(df_dataset[[col]]))}
  # get significance for each cell type
  for (celltype in unique(df_dataset$Cell.type))
  {
    print(celltype)
    output = pval_generator(df_dataset, celltype = celltype)
    compRecoder = paste0(dataset, ", ", celltype, ", ", output[[1]])
    pvalRecoder = output[[2]]
    compRecoder_all = append(compRecoder_all, compRecoder)
    pvalRecoder_all = append(pvalRecoder_all, pvalRecoder)
    datasetRecoder_all = append(datasetRecoder_all, rep(dataset, length(pvalRecoder)))
    celltypeRecoder_all = append(celltypeRecoder_all, rep(celltype, length(pvalRecoder)))
  }
}

# create a data frame for all information
df_output = data.frame(row.names = compRecoder_all,
                       "dataset" = datasetRecoder_all,
                       "cell class" = celltypeRecoder_all,
                       "pval" = pvalRecoder_all)

# multi-correction (global)
pval_adjusted_BH = p.adjust(pvalRecoder_all, method = "BH")
pval_adjusted_bonferroni = p.adjust(pvalRecoder_all, method = "bonferroni")
df_output["pval_BH_global"] = pval_adjusted_BH; df_output["pval_bonferroni_global"] = pval_adjusted_bonferroni

# multi-correction (for each dataset)
for (dataset in unique(df_output$dataset))
{
  df_sub = df_output[df_output$dataset == dataset, ]
  comps = rownames(df_sub)
  pvals_sub = df_sub[["pval"]]
  pval_adjusted_BH_local = p.adjust(pvals_sub, method = "BH")
  pval_adjusted_bonferroni_local = p.adjust(pvals_sub, method = "bonferroni")
  df_output[comps, "pval_BH_dataset"] = pval_adjusted_BH_local
  df_output[comps, "pval_bonferroni_dataset"] = pval_adjusted_bonferroni_local
}

# create output table
write.table(df_output, "LungBiopsy_DA_pval_adjusted.txt", sep = "\t")

