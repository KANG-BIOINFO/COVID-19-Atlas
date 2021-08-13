library(dplyr)
library(plyr)
library(ggplot2)
require(gridExtra)
library(ggpubr)
library(tibble)

# 1. generate the boxplot
box_generator <- function(df = df_bal, celltype = "Developing Neutrophil", celltype_col = "cell_type_fines", condition_col = "Ventilated", 
                          donor_col = "Donor_full", donor_color, comparison_pair)
{
  # compute the numerate for each Donor and mapping relations
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
  
  # order the position of 3 or more conditions in a reasonable way
  d = as.character(unique(df[["Dataset"]]))
  if (d == "Liao") {donor_pecentage$Condition = factor(donor_pecentage$Condition, levels = c("Healthy Control","Mild","Severe"))}
  else {donor_pecentage$Condition = factor(donor_pecentage$Condition, levels = c("Control","COVID-19"))}
  
  # draw box plots
  colour <- donor_color[donor_pecentage$Donor]
  p = ggplot(donor_pecentage, aes(x = Condition, y = Proportion)) + geom_boxplot(width = 0.5, lwd = 0.8, outlier.shape = NA) +
    geom_jitter(colour = colour, size = 5, width = 0.15) + 
    #geom_point(colour = colour, size = 5) + 
    ylab("") + xlab("") +ggtitle(celltype) + ylim(c(-0.0001,NA)) + 
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 2.5),
          axis.title = element_text(family = "Helvetica", size = 20),
          axis.text.y = element_text(family = "Helvetica", size = 15, hjust = -0.1),
          axis.text.x = element_text(family = "Helvetica", size = 12, vjust = 0.6, angle = 15),
          plot.title = element_text(family = "Helvetica", size = 20, hjust = 0.5, vjust = 0.2),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey", size = 0.25)) +
    stat_compare_means(comparisons = comparison_pair, method = "wilcox.test")
  
  return (p)
}

# 2. generate comparison list of different diseage stages
comparison_generator = function(df, disease_col = "disease.group")
{
  disea_pairs = list()
  disea_vals = unique(as.character(df[[disease_col]]))
  counts = 1
  for (i in 1:(length(disea_vals)-1))
  {
    for (j in (i+1):length(disea_vals))
    {
      disea_pairs[[counts]] = c(disea_vals[i], disea_vals[j])
      counts = counts + 1
    }
  }
  return (disea_pairs)
}

# 3. color code generator
color_generator = function(df, donor_col = "sample",condition_col = "disease.group")
{
  # get 1-1 relation between sample and disease condition
  donor2condition = df[,condition_col]
  names(donor2condition) = df[,donor_col]
  color_map = c("Healthy"= "#4169E1", "control"= "#4169E1","Healthy donor" = "#4169E1","Healthy Control" = "#4169E1","Control"="#4169E1",
                "Mild"="#F5DEB3","COVID-19 Moderate"="#F5DEB3", "mild COVID-19" = "#F5DEB3", "mild COVID-19\n(asymptomatic)" = "#F5DEB3","mild_early"="#F5DEB3","mild_late"="#F5DEB3","NonVent" ="#F5DEB3", 
                "COVID-19"="#FF2F93","Severe"="#FF2F93","COVID-19 Severe" = "#FF2F93","severe" = "#FF2F93", "severe COVID-19" = "#FF2F93", "severe influenza" = "#FF2F93","severe_early"="#FF2F93","severe_late"="#FF2F93","Vent"="#FF2F93",
                "remission" = "#DC9DDC")
  colors = c()
  donors = c()
  for (i in unique(df[[donor_col]]))
  {
    donors = append(donors, i)
    colors = append(colors, as.character(color_map[as.character(donor2condition[i])]))
  }
  names(colors) = donors
  return (colors)
}

# 4. run the box_generator function and generate box plots for BAL
df_all = read.table("./covid19 Atlas/BAL Atlas/integration/ctype_BAL_integrated_v2.txt",sep="\t",header=1, row.names = 1)
setwd("./covid19 Atlas/BAL Atlas/cell proportion/")
df_all$cell.group = df_all$Cell.class_v2
df_all$cell.group = mapvalues(df_all$Cell.class_v2, from = c("MoAM1","MoAM2","MoAM3","MoAM4","MoAM5",
                                                             "TRAM1","TRAM2","TRAM3",
                                                            "CD8+ T","CD4+ T_1","CD4+ T_2",
                                                             "Ciliated","Basal/Club","AT1/AT2","undefined Epi","transitional Epi"),
                              to = c(rep("MoAM",5), rep("TRAM",3),rep("T cell",3), rep("Epithelial cell",5)))

# run the code for each dataset
for (dataset in c("Grant","Liao"))
{
  df_dataset = df_all[df_all$Dataset == dataset,]
  for (col in colnames(df_dataset)){df_dataset[[col]] = as.factor(as.character(df_dataset[[col]]))}
  
  # run box_generator and save figure per cell type per dataset
  for (celltype in unique(df_dataset$cell.group))
  {
    p = box_generator(df_dataset, celltype = celltype, celltype_col = "cell.group",condition_col = "disease.group",
                      donor_col = "sample", donor_color = color_generator(df_dataset), comparison_pair = comparison_generator(df_dataset))
    outputFile_name = sub("/","_",paste0(dataset,"_",celltype,".png"))
    ggsave(p, file = outputFile_name, width = 11, height = 11, units = "cm")
  }
}

# 5. run the box_generator function and generate box plots for PBMC
df_all = read.table("ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2_diseaseFixed.txt",sep="\t",header=1,row.names = 1)
for (dataset in c("Lee","Wilk","Arunachalam","Guo","Schulte-Schrepping")) 
{
  df_dataset = df_all[df_all$Dataset == dataset,]
  for (col in colnames(df_dataset)){df_dataset[[col]] = as.factor(as.character(df_dataset[[col]]))}
  
  # run box_generator and save figure per cell type per dataset
  for (celltype in unique(df_dataset$Cell.group))
  {
    p = box_generator(df_dataset, celltype = celltype, celltype_col = "Cell.group",condition_col = "disease.group",
                      donor_col = "sample", donor_color = color_generator(df_dataset))
    outputFile_name = sub("/","_",paste0(dataset,"_",celltype,".png"))
    ggsave(p, file = paste0("Cell proportion_0127/cell group/",outputFile_name), width = 6, height = 6, units = "cm")
  }
}

# 6. run the box_generator function and generate box plots for Lung Biopsy
df_all = read.table("GSE158127_04all-lineages-metadata_NU_Parenchyma_only_macrophageFixed_Adjusted.txt",sep="\t",header=1, row.names = 1)
df_dataset = df_all
for (col in colnames(df_dataset)){df_dataset[[col]] = as.factor(as.character(df_dataset[[col]]))}
for (celltype in unique(df_dataset$Cell.type))
{
  p = box_generator(df_dataset, celltype = celltype, celltype_col = "Cell.type",condition_col = "Diagnosis",
                    donor_col = "orig.ident", donor_color = color_generator(df_dataset, "orig.ident", "Diagnosis"))
  outputFile_name = sub("/","_",paste0(celltype,".png"))
  ggsave(p, file = paste0("Cell proportion_0207_Parenchyma/",outputFile_name), width = 6, height = 6, units = "cm")
}


