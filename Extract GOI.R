##########################
#Extract gene of interest
library(tidyverse)
library(tidyr)
#Wnt Target
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/")
Wnt_tar <- read.csv("WNT target.txt", h=F)

dds_EpiSC_WNT <- dds_EpiSC[c(rownames(dds_EpiSC) %in% Wnt_tar$V1), ]
rownames(dds_EpiSC_WNT)

WNT_tar_dat <- t(assay(dds_EpiSC_WNT))

CellLine <- substring(rownames(WNT_tar_dat), 8, 15)

TABLE <- cbind.data.frame(WNT_tar_dat, CellLine)

MEAN_TPM_run <- tibble(x=1:length(unique(CellLine)))

for (i in 1:(length(colnames(TABLE))-1)){
  Genes = colnames(TABLE)[i]
  temp_mean <- TABLE %>% dplyr::select (GOI = Genes, CellLine)
  Mean_Dat <- aggregate(GOI ~ CellLine, temp_mean, mean)
  #colnames(Mean_Dat) <- c("CellLine", Genes)
  MEAN_TPM_run <- add_column(MEAN_TPM_run, Mean_Dat[,2] )
  names(MEAN_TPM_run)[ncol(MEAN_TPM_run)] <- paste0("mean_", Genes)
}
MEAN_TPM_run <- cbind.data.frame(unique(CellLine), MEAN_TPM_run)

######################
#Plot it
MeanToPlot <- MEAN_TPM_run[c(3,7,8),-2]

MeanToPlot.reshape <- reshape(MeanToPlot, direction="long",  idvar = "unique(CellLine)", 
                        varying=c(colnames(MeanToPlot)[2:84]), sep="_")

HeatMapDat <- cbind.data.frame(MeanToPlot.reshape, 
                               c(rep(substring(colnames(MeanToPlot[2:84]),6, nchar(colnames(MeanToPlot[2:84])))
                                     , each=nrow(MeanToPlot))))
colnames(HeatMapDat) <- c("CellLine", "time", "Mean", "Gene")



########################
#heatmap. 3 -- load the function prior to run that code https://gist.github.com/tleja/7783150
# or here: C:\Users\Pierre\Desktop\Bioinformatics
require(gtools)
require(RColorBrewer)

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

heatmap.3(myMat, trace="none", scale="column", zlim=c(-3,3), reorder=FALSE,na.rm=T,
          distfun=function(x) dist(x,method = 'manhattan'), hclustfun=hclustAvg,col=rev(cols), symbreak=FALSE, cexRow = 0.75)
