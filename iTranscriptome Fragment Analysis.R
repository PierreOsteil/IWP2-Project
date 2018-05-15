# Embryonic fragment analysis
# This code requires the use of "Open_quant_sf_files.R" code to generate "txi" and "files"

#DESeq2

SampleTable <- data.frame(substring(files, 1, 1))
colnames(SampleTable) <- "TissueType"
rownames(SampleTable) <- colnames(txi$counts)


library("DESeq2")
dds <- DESeqDataSetFromTximport(txi, SampleTable, ~TissueType)

dds <- dds[ rowSums(counts(dds)) > 1, ]

TABLE <- assay(dds)
head(TABLE)

################### Change all sample names ##################################################
#clean the names
colnames(TABLE) <- c(gsub("_.R1.clean","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".R1.clean","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq_quant.sf","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq.gz_quant.sf","", colnames(TABLE)))
colnames(TABLE) <- c(gsub(".fastq.bz2_quant.sf","", colnames(TABLE)))

###Embryos
#Chang samples
colnames(TABLE) <- c(gsub("AVE","Embryo_Cha_E7.5_AV", colnames(TABLE)))
colnames(TABLE) <- c(gsub("PVE","Embryo_Cha_E7.5_PV", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5ExE","Embryo_Cha_E7.5_EX", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5AP","Embryo_Cha_E7.5_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5AD","Embryo_Cha_E7.5_AD", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5P","Embryo_Cha_E7.5_Po", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E7.0ExE","Embryo_Cha_E7.0_EX", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.0AP","Embryo_Cha_E7.0_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.0AD","Embryo_Cha_E7.0_AD", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.0P","Embryo_Cha_E7.0_Po", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E6.5A2","Embryo_Cha_E6.5_An", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E6.5P2","Embryo_Cha_E6.5_Po", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E6.5ExE2","Embryo_Cha_E6.5_EX", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E6Epi4","Embryo_Cha_E6.0_Al", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E6ExE4","Embryo_Cha_E6.0_ExE", colnames(TABLE)))

colnames(TABLE) <- c(gsub("E5.5Epi","Embryo_Cha_E5.5_Al", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E5.5ExE","Embryo_Cha_E5.5_ExE", colnames(TABLE)))

#Peng Samples
colnames(TABLE) <- c(gsub("E5.5_2.","Embryo_Pen_E5.5_", colnames(TABLE)))# 1 and 2 are similar dataset
colnames(TABLE) <- c(gsub("E6.0_2.","Embryo_Pen_E6.0_", colnames(TABLE)))# 2 contains ExE tissues
colnames(TABLE) <- c(gsub("E6.5_5.","Embryo_Pen_E6.5_", colnames(TABLE)))# The ref embryo Guangdun is using
colnames(TABLE) <- c(gsub("E7.0_1.","Embryo_Pen_E7.0_N", colnames(TABLE)))# N for New data set (ref as well)
colnames(TABLE) <- c(gsub("E7.5_5.","Embryo_Pen_E7.5_5_", colnames(TABLE))) 
colnames(TABLE) <- c(gsub("E7.5_2","Embryo_Pen_E7.5_2", colnames(TABLE)))
colnames(TABLE) <- c(gsub("E7.5_3","Embryo_Pen_E7.5_3", colnames(TABLE)))
#Wu samples
colnames(TABLE) <- c(gsub("SRR1660273","Embryo_Wu._E6.5_AP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660274","Embryo_Wu._E6.5_AD", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660275","Embryo_Wu._E6.5_PP", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660276","Embryo_Wu._E6.5_PD", colnames(TABLE)))

###EpiSCs
#Osteil
colnames(TABLE) <- c(gsub("DKKminusIW.","EpiSCs_Ost_IwAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("DKKplusIW.","EpiSCs_Ost_IwIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("YC8plusIW.","EpiSCs_Ost_AFIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("YC8minusIW.","EpiSCs_Ost_AFAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-DKKnul.","EpiSCs_Ost_DKKn", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I_AF.","EpiSCs_Ost_IwAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I.","EpiSCs_Ost_IwIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-AF.","EpiSCs_Ost_AFAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FASI.","EpiSCs_Ost_FASI", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FSI.","EpiSCs_Ost_F-SI", colnames(TABLE)))

#Kurek
colnames(TABLE) <- c(gsub("SRR1609334","EpiSCs_Kur_AFIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609335","EpiSCs_Kur_Wnt3", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609336","EpiSCs_Kur_BMP4", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609337","EpiSCs_Kur_WnBM", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609338","EpiSCs_Kur_BmIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1609339","EpiSCs_Kur_AFAF", colnames(TABLE)))

#Wu
colnames(TABLE) <- c(gsub("SRR1660269_1","EpiSCs_Wu._AFAF_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660270_1","EpiSCs_Wu._AFAF_2", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660271_1","EpiSCs_Wu._IwIw_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("SRR1660272_1","EpiSCs_Wu._IwIw_2", colnames(TABLE)))

#Peng DevCell Data
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/Salmon")
PengDevCell <- read.table("sample DevCell iTra.txt", sep= "\t")$V1
PengDevCell_dat <- data.frame(TABLE[, c(grep("SRR", colnames(TABLE)))])
colnames(PengDevCell_dat) <- PengDevCell
woPengDevCell_dat <- data.frame(TABLE[, -c(grep("SRR", colnames(TABLE)))])
TABLE_clean <- cbind(woPengDevCell_dat,PengDevCell_dat)
#Check
colnames(TABLE_clean)


################### TABLE REFERENCE EMBRYO ###################
E5.5_dat <- TABLE_clean[,c(grep("Pen_E5.5", colnames(TABLE_clean)))]
E6.0_dat <- TABLE_clean[,c(grep("Pen_E6.0", colnames(TABLE_clean)))]
E6.5_dat <- TABLE_clean[,c(grep("Pen_E6.5", colnames(TABLE_clean)))]
E7.0_dat <- TABLE_clean[,c(grep("Pen_E7.0_N", colnames(TABLE_clean)))]
E7.5_dat <- TABLE_clean[,c(grep("Pen_E7.5_5_", colnames(TABLE_clean)))]

TABLEPeng <- cbind(E5.5_dat, E6.0_dat, E6.5_dat, E7.0_dat, E7.5_dat)

#only Peng_embryos
colnames(TABLEPeng) <- c(gsub("__","_",colnames(TABLEPeng)))
TABLEPeng <- TABLEPeng[,-c(grep("M", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("EnA", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("EnP", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("V", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("EA", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("EP", colnames(TABLEPeng)))]
TABLEPeng <- TABLEPeng[,-c(grep("E7.0_N9R", colnames(TABLEPeng)))]


##### Sample annotation
Stage <- substring(colnames(TABLEPeng), 12, 15)
Lab <- substring(colnames(TABLEPeng), 8, 10)



##################################
#PCA analysis
library("FactoMineR")
library("rafalib")
library("ggplot2")
library("ggrepel")


#plot PCA
mypar(1,2)

PCA_TOT=PCA(t(TABLEPeng) , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)

PCAcoord12 <- cbind.data.frame(PCAcoord[,3], PCAcoord[,2])

PCA_data <- cbind.data.frame(PCAcoord12, Stage)
colnames(PCA_data) <- c("PC1", "PC2",  "Stage")

PCA <- ggplot(PCA_data, aes(PC1, PC2, color=Stage)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  #geom_text_repel(aes(label=colnames(TABLEPeng)))+
  theme_bw()
PCA

## 3D PCA
library("rgl")
plot3d(PCAcoord[,1],PCAcoord[,2],PCAcoord[,4], col=as.integer(PCA_data$Stage))


######################################################################################
#extract genes contributing to Axis 3
PCAvar <- as.data.frame(PCA_TOT$var)
GeneOnPCs <- data.frame(PCAvar[,3], PCAvar[,2])#PC3 is the best at predicting timing
rownames(GeneOnPCs) <- rownames(TABLEPeng)

axis3 <- GeneOnPCs[order(GeneOnPCs$PCAvar...3.),] 
EarlyGenes <- rownames(axis3[c(1:200),]) #top100 genes left PC3
LateGenes <- rownames(axis3[c((length(rownames(axis3))-199):length(rownames(axis3))),]) #top100 genes right PC2
Gene400 <- c(EarlyGenes, LateGenes)
write.csv(x = Gene400, "Gene400.csv")


##############
#Select samples
TABLEAll <- cbind(TABLEPeng, TABLE_clean[,c(grep("Ost", colnames(TABLE_clean)),
                                            grep("Wu", colnames(TABLE_clean)),
                                            grep("Cha", colnames(TABLE_clean)),
                                            grep("Kur", colnames(TABLE_clean))
)])


#With the 400 genes
TABLEStageGene <- TABLEAll[c(which(rownames(TABLEAll) %in% c(Gene400))),]

#Negative binomiale normalisation
library("edgeR")
TABLEAll_znorm <- apply(TABLEStageGene , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))
TABLEAll_znorm <- data.frame(TABLEAll_znorm)

#sample annotation
Stage <- substring(colnames(TABLEAll_znorm),12, 15)
Lab <- substring(colnames(TABLEStageGene), 8, 10)
SizeText <- c(rep("A", ncol(TABLEStageGene)))

####PCA
mypar(1,2)
PCA_TOT=PCA(t(TABLEAll_znorm) , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)
PCAcoord12 <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])


PCA_data <- cbind.data.frame(PCAcoord12, Stage, Lab, SizeText)

colnames(PCA_data) <- c("PC1", "PC2",  "Parameter", "Lab", "SizeText")

PCA <- ggplot(PCA_data, aes(PC1, PC2, color=Parameter, shape=Lab)) +
  geom_point(size=5) +
  scale_shape_manual(values=c(15,16,17,18,25))+
  scale_size_discrete()+
  scale_colour_manual(values=c("black", "blue", "darkorange3", "grey", "purple",
                               "cyan","lightseagreen","green", "green4", "darkolivegreen", 
                               "red", "pink", "magenta", "darkorange3", "grey", "grey"))+
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  geom_text_repel(aes(label=colnames(TABLEAll_znorm), size=SizeText))+
  theme_bw()
PCA

plot3d(PCAcoord[,1],PCAcoord[,2],PCAcoord[,3], col=as.integer(PCA_data$Parameter))




### TO REVISE





###################### Double Heatmap ############################

sampleDistMatrix <- as.matrix(dist(t(TABLE)))

DKKn <- sampleDistMatrix[c(grep("6.5", rownames(sampleDistMatrix))),
                         c(grep("DKK", colnames(sampleDistMatrix)), 
                           grep("Ost_AFAF", colnames(sampleDistMatrix)))]
DKKn <- DKKn[-c(grep("EP", rownames(DKKn))),]
DKKn <- DKKn[-c(grep("EA", rownames(DKKn))),]
DKKn <- DKKn[-c(grep("Cha", rownames(DKKn))),]
max(DKKn)
