library(tidyverse)

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
TABLE <- data.frame(files)
rownames(TABLE) <- files
TABLE <- t(TABLE)

################### Change all sample names ##################################################
#clean the names
colnames(TABLE) <- c(gsub("_.R1.clean","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".R1.clean","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq_quant.sf","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq.gz_quant.sf","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq.bz2_quant.sf","", colnames(TABLE)))
###Embryos
#Chang samples
colnames(TABLE) <- c(gsub("E7.5AP","Embryo_Cha_E7.5_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5AD","Embryo_Cha_E7.0_AD", colnames(TABLE)))#Warning
colnames(TABLE) <- c(gsub("E7.5P","Embryo_Cha_E7.5_Po", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E7.0AP","Embryo_Cha_E7.0_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.0AD","Embryo_Cha_E7.5_AD", colnames(TABLE)))#Warning
colnames(TABLE) <- c(gsub("E7.0P","Embryo_Cha_E7.0_Po", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E6.5A2","Embryo_Cha_E6.5_An", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E6.5P2","Embryo_Cha_E6.5_Po", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E6Epi4","Embryo_Cha_E6.0_Al", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E5.5Epi","Embryo_Cha_E5.5_Al", colnames(TABLE)))

colnames(TABLE) <- c(gsub("Chang_EpiSC.sf","EpiSCs_Cha_F-F-", colnames(TABLE)))
colnames(TABLE) <- c(gsub("Chang_EpiSC-SF1.sf","EpiSCs_Cha_F-SB_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("Chang_EpiSC-SF2.sf","EpiSCs_Cha_F-SB_2", colnames(TABLE)))

#Wu samples
colnames(TABLE) <- c(gsub("SRR1660273","Embryo_Wu._E6.5_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660274","Embryo_Wu._E6.5_AD", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660275","Embryo_Wu._E6.5_PP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660276","Embryo_Wu._E6.5_PD", colnames(TABLE)))

colnames(TABLE) <- c(gsub("SRR1660269_1","EpiSCs_Wu._AFAF_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660270_1","EpiSCs_Wu._AFAF_2", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660271_1","EpiSCs_Wu._IwIw_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660272_1","EpiSCs_Wu._IwIw_2", colnames(TABLE)))

#Kurek
colnames(TABLE) <- c(gsub("SRR1609334","EpiSCs_Kur_AFIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609335","EpiSCs_Kur_Wnt3", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609336","EpiSCs_Kur_BMP4", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609337","EpiSCs_Kur_WnBM", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609338","EpiSCs_Kur_BmIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609339","EpiSCs_Kur_AFAF", colnames(TABLE)))

#Osteil
colnames(TABLE) <- c(gsub("EpiSC-DKKnul.","EpiSCs_Ost_DKKn", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I_AF.","EpiSCs_Ost_IwAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I.","EpiSCs_Ost_IwIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-AF.","EpiSCs_Ost_AFAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FASI.","EpiSCs_Ost_FASI_", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FSI.","EpiSCs_Ost_F-SI_", colnames(TABLE)))

colnames(TABLE)
#Change the name
Lab_stage <- substring(colnames(TABLE), 8, 15)

##Allocate the new names to the samples
SampleTable <- cbind.data.frame(substring(colnames(TABLE),1,2), Lab_stage)
rownames(SampleTable) <- colnames(TABLE)
colnames(SampleTable) <- c("TissueType", "LabStage")
colnames(txi$counts) <- rownames(SampleTable)


#Dds from txi files
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi, SampleTable, ~LabStage) #extract txi data
save(dds, file ="dds.RData")

###
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
load("dds.RData")


dds <- dds[ rowSums(counts(dds)) > 1, ]


#transform the value
#log
rld <- rlog(dds, blind=FALSE)
save(rld, file="rld.Rdata")

#Select gene of interest
rld_Chang_Frag <- rld[,grep("Embryo_Cha", colnames(rld))]
colnames(rld_Chang_Frag)


#Negative binomiale normalisation
library("edgeR")
z_rld <- apply(assay(rld_Chang_Frag) , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))
head(z_rld)



pca <- prcomp(assay(rld_Chang_Frag), center = T)
PCAdata <- as.data.frame(pca$rotation)


library("ggplot2")
library("ggrepel")
ggplot(PCAdata, aes(PCAdata[,1], PCAdata[,2]))+
  geom_point( 
             size=3) +
  geom_text_repel(aes(label=colnames(rld_Chang_Frag)))+
  #xlab(paste0("PC1: ",round((pca_StageGene$sdev[1])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  #ylab(paste0("PC2: ",round((pca_StageGene$sdev[2])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  theme(aspect.ratio=1)


###################################
GenesOnPCs = as.data.frame(predict(pca)) #extract genes contributing to each PC``

axis2 <- GenesOnPCs[order(GenesOnPCs$PC2),] #PC2 is the best at predicting timing
EarlyGenes <- rownames(axis2[c(1:100),]) #top100 genes left PC2
LateGenes <- rownames(axis2[c((length(rownames(axis2)) - 100):length(rownames(axis2))),])
rld_StageGene <- rld[c(which(rownames(rld) %in% c(paste(EarlyGenes), paste(LateGenes)))),]

Paper <- substring(colnames(assay(rld)), 8,10)

#PCA
pca_StageGene <- prcomp(assay(rld_StageGene), center = T)
PCAdata <- as.data.frame(pca_StageGene$rotation)

PCAdata <- cbind.data.frame(PCAdata, Paper)
ggplot(PCAdata, aes(PCAdata[,1], PCAdata[,2]))+
  geom_point(aes(colour = Paper), 
             size=3) +
  geom_text_repel(aes(label=colnames(rld_StageGene)), cex=3)+
  xlab(paste0("PC1: ",round((pca_StageGene$sdev[1])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  ylab(paste0("PC2: ",round((pca_StageGene$sdev[2])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  theme(aspect.ratio=1)
