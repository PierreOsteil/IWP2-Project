##########################
#Extract gene of interest
library(tidyverse)
library(tidyr)
#Wnt Target
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/")
Wnt_tar <- read.csv("WNT target.txt", h=F)

dds_EpiSC_CWT1 <- dds_EpiSC[c(rownames(dds_EpiSC) %in% colnames(CL1_WN1_TG1_clean) #put "#" before %in% to have all genes from RNAseq
                              ), ]
rownames(dds_EpiSC_CWT1)

CWT1_dat <- t(assay(dds_EpiSC_CWT1))

CellLine <- substring(rownames(CWT1_dat), 8, 15)

TABLE <- cbind.data.frame(CWT1_dat, CellLine)

MEAN_TPM_run <- tibble(x=1:length(unique(CellLine)))

for (i in 1:(length(colnames(TABLE))-1)){
  Genes = colnames(TABLE)[i]
  temp_mean <- TABLE %>% dplyr::select (GOI = Genes, CellLine)
  Mean_Dat <- aggregate(GOI ~ CellLine, temp_mean, mean)
  #colnames(Mean_Dat) <- c("CellLine", Genes)
  MEAN_TPM_run <- add_column(MEAN_TPM_run, Mean_Dat[,2] )
  names(MEAN_TPM_run)[ncol(MEAN_TPM_run)] <- paste0("mean_", Genes)
}
MEAN_TPM_run2 <- cbind.data.frame(unique(CellLine), MEAN_TPM_run)

######################
#Select data
MeanToPlot <- MEAN_TPM_run[c(3,7,8),-2]



########################
#heatmap. 3 -- load the function prior to run that code https://gist.github.com/tleja/7783150
# or here: C:\Users\Pierre\Desktop\Bioinformatics
require(gtools)
require(RColorBrewer)

myMat <- as.matrix(MeanToPlot[,-1])
rownames(myMat) <- MeanToPlot$`unique(CellLine)`
colnames(myMat) <- substring(colnames(myMat), 6, nchar(colnames(myMat)))
myMat  <- myMat [, colSums(myMat) > 1 ] #remove genes with no expression

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

distCor <- function(x) as.dist(1-cor(t(x)))
hclustAvg <- function(x) hclust(x, method="average")

heatmap.3(myMat, trace="none", scale="column", zlim=c(-3,3), reorder=FALSE,na.rm=T,
          distfun=function(x) dist(x,method = 'manhattan'), hclustfun=hclustAvg,col=rev(cols), symbreak=FALSE, cexRow = 0.75)


#######################
#validation with Biomark Data
TABLE <- CL1_WN1_TG1_clean[grep("D0", WN1$Name),]
TABLE <- TABLE[-grep("FIw", TABLE$Name),]
TABLE <- TABLE[-grep("IwF", TABLE$Name),]
#TABLE <- TABLE[-grep("Wnt3a", TABLE$Name),]

CellLineCWT1 <- substring(TABLE$Name, 1,nchar(TABLE$Name)-4)
TABLE <- cbind.data.frame(TABLE, CellLineCWT1)

MEAN_CWT1_run <- tibble(x=1:length(unique(CellLineCWT1)))

for (i in 2:(length(colnames(TABLE))-1)){
  Genes = colnames(TABLE)[i]
  temp_mean <- TABLE %>% dplyr::select (GOI = Genes, CellLineCWT1)
  Mean_Dat <- aggregate(GOI ~ CellLineCWT1, temp_mean, mean)
  #colnames(Mean_Dat) <- c("CellLine", Genes)
  MEAN_CWT1_run <- add_column(MEAN_CWT1_run, Mean_Dat[,2] )
  names(MEAN_CWT1_run)[ncol(MEAN_CWT1_run)] <- paste0("mean_", Genes)
}
MEAN_CWT1_run <- cbind.data.frame(unique(CellLineCWT1), MEAN_CWT1_run)

myMat2 <- as.matrix(MEAN_CWT1_run[,-(1:2)])
rownames(myMat2) <- MEAN_CWT1_run$`unique(CellLineCWT1)`
colnames(myMat2) <- substring(colnames(myMat2), 6, nchar(colnames(myMat2)))




################
#heatmap.3
heatmap.3(as.matrix(myMat2), trace="none", scale="column", zlim=c(-3,3), reorder=FALSE,na.rm=T,
          distfun=function(x) dist(x,method = 'minkowski'), hclustfun=hclustAvg,col=rev(cols), symbreak=FALSE, cexRow = 0.75)


################
myMat <- as.data.frame(myMat)
myMat2 <- as.data.frame(myMat2)
myMat2.1 <- myMat2[ , c(which(colnames(myMat2) %in% colnames(myMat)))]

RNAseq_ratio <- myMat[3,]   - myMat[1,]    > 0  # for reminder myMat is RNAseq data
qPCR_ratio   <- myMat2.1[3,]- myMat2.1[1,] > 0  # for reminder myMat2.1 is qPCR data

sum(RNAseq_ratio == qPCR_ratio)
sum(RNAseq_ratio == qPCR_ratio)/ncol(myMat)*100

###############
#only select the genes showing similar results
myMat.2 <- myMat[, c(which(RNAseq_ratio == qPCR_ratio))]
myMat.2 <- as.matrix(myMat.2)
