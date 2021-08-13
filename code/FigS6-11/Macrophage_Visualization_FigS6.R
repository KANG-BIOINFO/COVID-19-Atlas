library(EnhancedVolcano)
library(ggplot2)

setwd("/data/aronow/Kang/data for toppcell/covid19 Atlas/BAL Atlas/DEGs/")
# 1. MoAM1,2,5 vs. MoAM3,4
# volcano plot
moam = read.table("./DEGs_MoAM1_MoAM2_MoAM5_vs_MoAM3_MoAM4.txt", sep="\t", header=1)
moam[moam$pvals_adj < 1e-300,"pvals_adj"] = 1e-300
moam$names = as.character(moam$names)
EnhancedVolcano(moam, lab = moam$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("ISG20","ISG15","IFIT1","IFITM1","IFIT3","CCL8","CXCL10","STAT1","CSF2RB","CCL3","CCL4","IL1RN","IL1R2","IL27","FCER1A",
                              "IL7R","APOE","PPARG","CCL13","CCL18","PPARG","FLNA","ENO1","MARCO","THBD","ITGAM","DAB2"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "MoAM3,4 vs. MoAM1,2,5", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-6,10))

# 2. TRAM1,2 vs. TRAM3
# volcano plot
tram = read.table("./DEGs_TRAM3_vs_TRAM1_TRAM2.txt", sep="\t", header=1)
tram[tram$pvals_adj < 1e-300,"pvals_adj"] = 1e-300
tram$names = as.character(tram$names)
EnhancedVolcano(tram, lab = tram$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("C1QC","CCR1","IFI6","IFI16","IFIT1","IFITM1","ISG15","IRF1","STAT1","TNFSF1O","C3AR1","IL10RA","IL11RA","IL1RN",
                              "S100A8","S100A9","FCN1"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "TRAM1,2 vs. TRAM3", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-6,15))

# 2. TRAM3 vs. MoAM1,2,5
# volcano plot
moam_tram = read.table("./DEGs_MoAM1_2_5_vs_TRAM3.txt", sep="\t", header=1)
moam_tram[moam_tram$pvals_adj < 1e-300,"pvals_adj"] = 1e-300
moam_tram$names = as.character(moam_tram$names)
EnhancedVolcano(moam_tram, lab = moam_tram$names, x = "logfoldchanges", y = "pvals_adj",
                selectLab = c("CCL2","CCL7","FCN1","CCL8","IFITM1","CCL3","CCL3L1","CCR5","IFITM2","ISG20",
                              "C1QA","C1QB","HLA-DQA2","TNFSF13","CXCL16","HLA-DRA"),
                boxedlabels = TRUE, drawconnectors = TRUE, colConnectors = "black",
                title = "TRAM3 vs. MoAM1,2,5", subtitle = NULL, caption = NULL,transcriptLabSize=4.5,
                xlim = c(-15,15))
