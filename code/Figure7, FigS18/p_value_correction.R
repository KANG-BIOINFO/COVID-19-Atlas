library(dplyr)
library(plyr)
library(ggplot2)

df_all = read.table("../Covid Atlas/disease comparison/ctype_integrated_concise_v2_adjusted.txt",sep="\t",header=1)
colnames(df_all)

#######################################################################
# function to create a single ratio vector per cell class per dataset #
#######################################################################
single_box = function(df, cell_class, disease_control, disease_target)
{
  df = df[df$Disease %in% c(disease_control, disease_target),]
  df = df[df$cell.class_concise == cell_class, ]
  df$sample = as.factor(as.character(df$sample))
  df$cell.class_concise = as.factor(as.character(df$cell.class_concise))
  df$Disease = as.factor(as.character(df$Disease))
  ratio_vector_HH = c()
  ratio_vector_DH = c()
  proportion_control = df[df$Disease == disease_control,"proportion.of.cells"]
  proportion_target = df[df$Disease == disease_target,"proportion.of.cells"]
  for (i in 1:length(proportion_control))
  {
    for (j in 1:length(proportion_target))
    {
      ratio_vector_DH = append(ratio_vector_DH,log2(proportion_target[j] / proportion_control[i]))
    }
    for (k in 1:length(proportion_control))
    {
      if (i != k) {ratio_vector_HH = append(ratio_vector_HH, log2(proportion_control[i] / proportion_control[k]))}
    }
  }
  significance = ks.test(ratio_vector_DH, ratio_vector_HH)$p.value
  return (list(ratio_vector_DH, significance))
}


#############################################################################################
# create multiple ratio vectors using all datasets and user-customized cell class selection #
#############################################################################################
multiple_box = function(df = df_all, cell_class_selection, control_map, target_map)
{
  # loop for each dataset
  for (dataset in unique(df$Dataset))
  {
    df_dataset = df[df$Dataset == dataset,]
    # get control and target level according to each dataset
    control_level = as.character(control_map[dataset])
    target_levels = as.character(target_maps[names(target_maps) == dataset])
    for (col in colnames(df_dataset)){df_dataset[[col]] = as.factor(as.character(df_dataset[[col]]))}
    
    ##################################################
    #### create abundance table for each dataset  ####
    ##################################################
    df_cellabundance_dataset = as.data.frame(table(df_dataset$sample, df_dataset$Cell.class_concise))
    colnames(df_cellabundance_dataset) = c("sample","cell.class_concise","num.of.cells")
    # find sample abundance as denominator and calculate the proportion
    df_sample_abundance = as.numeric(table(df_dataset$sample))
    names(df_sample_abundance) = names(table(df_dataset$sample))
    df_cellabundance_dataset$proportion.of.cells = 0
    for (i in rownames(df_cellabundance_dataset))
    {
      sample_name = as.character(df_cellabundance_dataset[i,1])
      df_cellabundance_dataset[i,4] = df_cellabundance_dataset[i,3] / as.numeric(df_sample_abundance[sample_name])
    }
    # find the sample to disease relation and add disease column
    temp = as.data.frame(table(df_dataset$sample,df_dataset$Disease))
    temp = temp[temp$Freq != 0, ]
    sample_to_disease = temp$Var2
    names(sample_to_disease) = temp$Var1
    df_cellabundance_dataset$Disease = as.character(sample_to_disease[as.character(df_cellabundance_dataset$sample)])

    ##################################
    #### loop for each cell class ####
    ##################################
    pvalRecorder = c()
    datasetRecorder = c()
    celltypeRecorder = c()
    compRecorder = c()
    print(target_levels)
    for (target_level in target_levels)
    {
      ratio_vector_all = c() # column ratio
      cell_class_vector_all = c() # columns cell class
      print(dataset)
      print(target_level)
      print(control_level)
      for (cell_class in cell_class_selection)
      {
        # run single_box function
        #list(ratio_vector, p_val) = single_box(df_cellabundance_dataset, cell_class = cell_class,disease_control = control_level,disease_target = target_level)
        out_temp = single_box(df_cellabundance_dataset, cell_class = cell_class,disease_control = control_level,disease_target = target_level)
        ratio_vector = out_temp[[1]]
        #print(ratio_vector)
        pval = out_temp[[2]]
        # add data to recorder
        pvalRecorder = append(pvalRecorder, pval)
        datasetRecorder = append(datasetRecorder, dataset)
        celltypeRecorder = append(celltypeRecorder, cell_class)
        comp = paste0(target_level, " vs. ", control_level, ", ", cell_class)
        compRecorder = append(compRecorder, comp)

      }
    }
  }
  return (list(compRecorder, datasetRecorder, celltypeRecorder, pvalRecorder))
}

################################################
###############   run our dataset ##############
################################################

selected_cell_classes = levels = rev(c("CD14+ Monocyte","CD16+ Monocyte","cDC","pDC","CD4+ T naive","CD4+ Tcm","CD8+ T naive","CD8+ Tem","NK","B naive","B memory","Plasmablast","Platelet"))
control_maps = c("Healthy/Control","Healthy/Control","IIF","Control")
names(control_maps) = c("Influenza","integrated COVID-19","multiple sclerosis","sepsis")

target_maps = c("Influenza Severe","COVID-19 Moderate/Mild/Early stage/NonVent", 
                "COVID-19 Severe/Late stage/Vent","Multiple Sclerosis","Leuk-UTI","ICU-SEP")
names(target_maps) = c("Influenza","integrated COVID-19","integrated COVID-19",
                       "multiple sclerosis","sepsis","sepsis")

# ms
df_ms = df_all[df_all$Dataset == "multiple sclerosis",]
for (col in colnames(df_ms)){ df_ms[[col]] = as.factor(as.character(df_ms[[col]]))}
df_output_ms = multiple_box(df_ms, cell_class_selection = selected_cell_classes, control_map = control_maps, target_map = target_maps)

# influenza
df_flu = df_all[df_all$Dataset == "Influenza",]
for (col in colnames(df_flu)){ df_flu[[col]] = as.factor(as.character(df_flu[[col]]))}
df_output_flu = multiple_box(df_flu, cell_class_selection = selected_cell_classes, control_map = control_maps, target_map = target_maps)

# sepsis
df_sepsis = df_all[df_all$Dataset == "sepsis",]
for (col in colnames(df_sepsis)){ df_sepsis[[col]] = as.factor(as.character(df_sepsis[[col]]))}
cell_class_selection_sepsis = rev(c("CD14+ Monocyte","CD16+ Monocyte","cDC","pDC","CD4+ T naive","CD4+ Tcm","CD8+ T naive","CD8+ Tem","NK","B naive","B memory","Plasmablast"))
df_output_sepsis = multiple_box(df_sepsis, cell_class_selection = cell_class_selection_sepsis, control_map = control_maps, target_map = target_maps)

# integrated covid-19 (select Lee et al. only)
df_covid = df_all[df_all$Dataset == "integrated COVID-19",]
for (col in colnames(df_covid)){ df_covid[[col]] = as.factor(as.character(df_covid[[col]]))}
df_covid_Lee = df_covid[df_covid$sample %in% c("Normal 1","Normal 2","Normal 3","Normal 4",
                                               "nCoV 1","nCoV 2","nCoV 3","nCoV 4","nCoV 5","nCoV 6","nCoV 7","nCoV 8","nCoV 9","nCoV 10","nCoV 11"),]
for (col in colnames(df_covid_Lee)){ df_covid_Lee[[col]] = as.factor(as.character(df_covid_Lee[[col]]))}
df_output_covid = multiple_box(df_covid_Lee, cell_class_selection = selected_cell_classes, control_map = control_maps, target_map = target_maps)


##########################
## organize the output ###
##########################
organize_output = function(output)
{
  compRecoder_all = output[[1]]
  datasetRecoder_all = output[[2]]
  celltypeRecoder_all = output[[3]]
  pvalRecoder_all = output[[4]]
  # create a data frame for all information
  df_output = data.frame(row.names = compRecoder_all,
                         "dataset" = datasetRecoder_all,
                         "cell group" = celltypeRecoder_all,
                         "pval" = pvalRecoder_all)
  
  # multi-correction (for each dataset)
  pval_adjusted_BH_local = p.adjust(pvalRecoder_all, method = "BH")
  pval_adjusted_bonferroni_local = p.adjust(pvalRecoder_all, method = "bonferroni")
  df_output[, "pval_BH_dataset"] = pval_adjusted_BH_local
  df_output[, "pval_bonferroni_dataset"] = pval_adjusted_bonferroni_local
  
  return (df_output)
}

df_output_all = rbind(organize_output(df_output_ms),
                      organize_output(df_output_flu),
                      organize_output(df_output_sepsis),
                      organize_output(df_output_covid))
df_output_all[, "pval_BH_global"] = p.adjust(df_output_all$pval, method = "BH")
df_output_all[, "pval_bonferroni_global"] = p.adjust(df_output_all$pval, method = "bonferroni")
write.table(df_output_all, "DiseaseComparison_DA_pval_adjusted.txt", sep = "\t")
save.image("pval_corr_fig7b.rda")
