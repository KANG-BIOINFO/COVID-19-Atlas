library(dplyr)
library(plyr)
library(ggplot2)

### 1. load data
df_cells = read.table("./covid19 Atlas/PBMC Atlas/ctype_cellClass_manual_B_T_NK_cell_predicted_fixed_RemoveConfused_sampleAdded_v2_diseaseFixed_myeloidReformatted.txt",sep="\t",header=1,row.names = 1)

### 2. draw an overall barplot
# add abbreviation to cell class
df_cells$Dataset_abbr = mapvalues(df_cells$Dataset, from = as.character(unique(df_cells$Dataset)), to = c("Aru","Guo","Lee","Sch","Wilk"))
df_cells$Cell.class_v2 = paste0(df_cells$Dataset_abbr,"_", df_cells$Cell.class)
write.table(df_cells, file = "ctype_merged_1129.txt",sep="\t")

# add cell group level
df_anno = read.table("./ctype_merged_1219_anno.txt",sep="\t",header=1, row.names = 1) # cell.class_v2 metadata
df_cells$Cell.group = df_anno[as.character(df_cells$Cell.class_v2),"Cell.group"] # add cell group level

# order cell.class_v2 in a reasonable order
df_anno = read.table("./ctype_merged_1219_anno.txt",sep="\t",header=1)
df_cells = df_cells[df_cells$Cell.class_v2 %in% unique(df_anno$Cell.class_v2),]
df_cells$Cell.class_v2 = factor(df_cells$Cell.class_v2, levels = as.character(df_anno$Cell.class_v2))

# map COVID19 conditions to consistent criteria
unique(df_cells$disease.group)
df_cells$disease.group_v3 = mapvalues(df_cells$disease.group_v2,from = as.character(unique(df_cells$disease.group_v2)), 
                                      to = c("Healthy","Mild/Remission","Severe","Mild/Remission"))
for (col in colnames(df_cells)){df_cells[[col]] = as.factor(as.character(df_cells[[col]]))}

dist_cells = as.data.frame(table(df_cells$Cell.class_v2, df_cells$disease.group_v3))
colnames(dist_cells) <- c("cell class","disease group","num of cells")
dist_cells$`cell class` = factor(dist_cells$`cell class`,levels = as.character(df_anno$Cell.class_v2))
brks <- c(0,0.25,0.5,0.75,1)
cellgroup_cols = c(rep("#E9967A",5), rep("#FA8072",5),
                   rep("#BDB76B",5), rep("#ADFF2F",11), rep("#228B22",15),rep("#808000",18), 
                   rep("#5F9EA0",10), rep("#6495ED",5),rep("#1E90FF",15),rep("#6A5ACD",4),rep("#9370DB",6),
                   rep("#DDA0DD",5), rep("#FF69B4",3),rep("#CC0066",10))
ggplot(dist_cells,
      aes(fill=`disease group`,y=`num of cells`,x=`cell class`)) +
      geom_bar(stat = "identity",position="fill",width = 0.90,binwidth=0) + 
      ylab("Fraction of Cells") + labs(fill="Disease group")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1.0, family = "Helvetica",size=10,face="bold",color = cellgroup_cols),
            panel.background = element_blank(),
            axis.text.y = element_text(family = "Helvetica",size=20),
            axis.line = element_line(),
            axis.title = element_blank(),
            legend.text = element_text(size=12),
            legend.title = element_text(size=12),
            legend.position = "none",
            plot.margin = unit(c(1,1,1,1),"cm")) + 
      scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks)) + 
      scale_fill_manual("Disease group", values = c("Healthy"="#4169E1","Mild/Remission"="#F5DEB3","Severe"="#FF1493"))

write.table(df_cells, file = "ctype_merged_1129_v2.txt",sep="\t")

### 3. draw distribution of different disease stages in each datasets
dist_cells = as.data.frame(table(df_cells$Dataset, df_cells$disease.group_v3))
colnames(dist_cells) <- c("dataset","disease group","num of cells")
brks <- c(0,0.25,0.5,0.75,1)
ggplot(dist_cells,
       aes(fill=`disease group`,y=`num of cells`,x=dataset)) +
  geom_bar(stat = "identity",position="fill",width = 0.8) + 
  ylab("Fraction of Cells") + labs(fill="Disease group")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1, family = "Helvetica",size=20,face="bold"),
        panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),
        axis.line = element_line(),
        axis.title = element_blank(),
        legend.text = element_text(size=15,family = "Helvetica"),
        legend.title = element_text(size=28,family = "Helvetica",face="bold"),
        legend.position = "top",
        plot.margin = unit(c(1,1,1,1),"cm")) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks)) + 
  scale_fill_manual("Disease group", values = c("Healthy"="#4169E1","Mild/Remission"="#F5DEB3","Severe"="#FF1493"))
