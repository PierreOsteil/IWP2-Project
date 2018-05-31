#########
# This code is to extract the DEGs of Peng's data at each stage for A and P
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
colnames(TABLE) <- c(gsub("DKKplusIW.","EpiSCs_Ost_IwIw_", colnames(TABLE)))
colnames(TABLE) <- c(gsub("YC8plusIW.","EpiSCs_Ost_AFIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("YC8minusIW.","EpiSCs_Ost_AFAF_", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-DKKnul.","EpiSCs_Ost_DKKn", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I_AF.","EpiSCs_Ost_IwAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-I.","EpiSCs_Ost_IwIw", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-AF.","EpiSCs_Ost_AFAF", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FASI.","EpiSCs_Ost_FASI_", colnames(TABLE)))
colnames(TABLE) <- c(gsub("EpiSC-FSI.","EpiSCs_Ost_F-SI_", colnames(TABLE)))

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
PengDevCell_dat <- t(data.frame(TABLE[, c(grep("SRR", colnames(TABLE)))]))
colnames(PengDevCell_dat) <- PengDevCell
woPengDevCell_dat <- data.frame(TABLE[, -c(grep("SRR", colnames(TABLE)))])
TABLE_name <- c(colnames(t(woPengDevCell_dat)),colnames(PengDevCell_dat))

#Change the name
Lab_stage <- substring(TABLE_name, 12, 17)
AP <-  substring(TABLE_name, nchar(TABLE_name), nchar(TABLE_name))
  #AP <-  substring(TABLE_name, 18, 18) #for embryo E6.5
LabStageAP <- paste(Lab_stage, AP)

##Allocate the new names to the samples
SampleTable <- cbind.data.frame(substring(TABLE_name,1,2), LabStageAP)
rownames(SampleTable) <- TABLE_name
colnames(SampleTable) <- c("TissueType", "LabStageAP")
colnames(txi$counts) <- rownames(SampleTable)


#Dds from txi files
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi, SampleTable, ~LabStageAP) #extract txi data

dds_temp <- dds[,grep("Embryo_Pen_E7.5_3", colnames(dds))] #choose the embryo to analyse in the cohort. 

dds_temp$LabStageAP <- droplevels(dds_temp$LabStageAP) #Remove sample annotation for unwanted samples

dds2 <- DESeq(dds_temp) # shouldn't take more than a few minute

res5.5 <- results(dds2, contrast = c("LabStageAP", "E5.5_ A", "E5.5_ P"))
res6.0 <- results(dds2, contrast = c("LabStageAP", "E6.0_ A", "E6.0_ B"))
res6.5 <- results(dds2, contrast = c("LabStageAP", "E6.5_ A", "E6.5_ P"))
res7.0 <- results(dds2, contrast = c("LabStageAP", "E7.0_ A", "E7.0_ R"))
res7.5 <- results(dds2, contrast = c("LabStageAP", "E7.5_3 A", "E7.5_3 P"))


#extracting DEGs

DEG5.5 <- cbind.data.frame(rownames(res5.5), res5.5$log2FoldChange, res5.5$pvalue)
colnames(DEG5.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG5.5 <- DEG5.5 %>% filter ( pvalue < 0.1 , abs(log2FoldChange) > 1 )
nrow(DEG5.5)


DEG6.0 <- cbind.data.frame(rownames(res6.0), res6.0$log2FoldChange, res6.0$pvalue)
colnames(DEG6.0) <- c("Gene", "log2FoldChange", "pvalue")
DEG6.0 <- DEG6.0 %>% filter ( pvalue < 0.05, abs(log2FoldChange) > 1 )
nrow(DEG6.0)


DEG6.5 <- cbind.data.frame(rownames(res6.5), res6.5$log2FoldChange, res6.5$pvalue)
colnames(DEG6.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG6.5 <- DEG6.5 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG6.5)


DEG7.0 <- cbind.data.frame(rownames(res7.0), res7.0$log2FoldChange, res7.0$pvalue)
colnames(DEG7.0) <- c("Gene", "log2FoldChange", "pvalue")
DEG7.0 <- DEG7.0 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG7.0)


DEG7.5 <- cbind.data.frame(rownames(res7.5), res7.5$log2FoldChange, res7.5$pvalue)
colnames(DEG7.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG7.5 <- DEG7.5 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG7.5)




############################################################################################
#extracting top genes
DEG7.5 <- DEG7.5[order(DEG7.5[,2]), ]
DEG7.5 <- DEG7.5[-c(grep("Gm", DEG7.5[,1])), ]
nrow(DEG7.5)
antetop_7.5 <- DEG7.5[(nrow(DEG7.5)-99):nrow(DEG7.5), ]
posttop_7.5 <- DEG7.5[1:100, ]
topDEG7.5 <- rbind(antetop_7.5, posttop_7.5)

DEG7.0 <- DEG7.0[order(DEG7.0[,2]), ]
DEG7.0 <- DEG7.0[-c(grep("Gm", DEG7.0[,1])), ]
nrow(DEG7.0)
ante100_7.0 <- DEG7.0[(nrow(DEG7.0)-99):nrow(DEG7.0), ]
post100_7.0 <- DEG7.0[1:100, ]

DEG6.5 <- DEG6.5[order(DEG6.5[,2]), ]
DEG6.5 <- DEG6.5[-c(grep("Gm", DEG6.5[,1])), ]
nrow(DEG6.5)
ante100_6.5 <- DEG6.5[(nrow(DEG6.5)-79):nrow(DEG6.5), ] #only 80 are DEGs in the anterior part. 
post100_6.5 <- DEG6.5[1:100, ]

DEG6.0 <- DEG6.0[order(DEG6.0[,2]), ]
DEG6.0 <- DEG6.0[-c(grep("Gm", DEG6.0[,1])), ]
nrow(DEG6.0)
ante100_6.0 <- DEG6.0[(nrow(DEG6.0)-199):nrow(DEG6.0), ]
post100_6.0 <- DEG6.0[1:200, ]

DEG5.5 <- DEG5.5[order(DEG5.5[,2]), ]
DEG5.5 <- DEG5.5[-c(grep("Gm", DEG5.5[,1])), ]
nrow(DEG5.5)
ante100_5.5 <- DEG5.5[(nrow(DEG5.5)-199):nrow(DEG5.5), ]
post100_5.5 <- DEG5.5[1:200, ]

#################################

###PCA for mapping
dds_E7.5 <- as.data.frame(assay(dds))
rownames(dds_E7.5)
dds_E7.5 <- dds_E7.5 %>% filter (rownames(dds_E7.5) %in% topDEG7.5$Gene)
dds_E7.5 <- dds_E7.5[,c(c(grep("E7.5_3", colnames(dds_E7.5))),
                        c(grep("Ost_IwIw", colnames(dds_E7.5))),
                        c(grep("Ost_AFAF_", colnames(dds_E7.5))
                        ))]
dds_E7.5 <- dds_E7.5[,-c(c(grep("R", colnames(dds_E7.5))),
                        c(grep("L", colnames(dds_E7.5))))]

#Negative binomiale normalisation
library("edgeR")
z_dds_E7.5 <- apply(dds_E7.5 , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))

mypar(1,2)
PCA_TOT=PCA(t(z_dds_E7.5) , scale.unit=T,ncp=5, axes = c(1,2))
PCAcoord <- as.data.frame(PCA_TOT$ind)
PCA_data <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])

colnames(PCA_data) <- c("PC1", "PC2")

PCA <- ggplot(PCA_data, aes(PC1, PC2)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  geom_text_repel(aes(label=colnames(dds_E7.5)))+
  theme_bw()
PCA



#heatmap

hmcols <- colorRampPalette(c("blue","white","red"), space = "Lab")(256)
heatmap(as.matrix(z_dds_E7.5), col =hmcols, cexCol = 0.5
        #,Rowv = NA
        #, Colv= NA
)


####################################################################
##DEGs for different EpiSCs 
Lab_Cond <- substring(TABLE_name, 8, 16)
##Allocate the new names to the samples
SampleTable <- cbind.data.frame(substring(TABLE_name,1,2), Lab_Cond)
rownames(SampleTable) <- TABLE_name
colnames(SampleTable) <- c("TissueType", "Lab_Cond")
colnames(txi$counts) <- rownames(SampleTable)

dds <- DESeqDataSetFromTximport(txi, SampleTable, ~Lab_Cond) #extract txi data

dds_EpiSC <- dds[,grep("EpiSC", colnames(dds))] #choose the samples to analyse in the cohort. 

#DESeq2
dds_EpiSC$Lab_Cond <- droplevels(dds_EpiSC$Lab_Cond) #Remove sample annotation for unwanted samples

dds2 <- DESeq(dds_EpiSC) # shouldn't take more than a few minute

resIwIwvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_IwIw_", "Ost_AFAF_"))
resDKKvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_DKKn", "Ost_AFAF"))
resFSIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_F-SI_", "Ost_AFAF_"))
resFASIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_FASI_", "Ost_AFAF_"))

#extracting DEGs
#IwIw versus AFAF
DEGIwIwvsAFAF <- cbind.data.frame(rownames(resIwIwvsAFAF), resIwIwvsAFAF$log2FoldChange, resIwIwvsAFAF$padj)
colnames(DEGIwIwvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGIwIwvsAFAF <- DEGIwIwvsAFAF %>% filter ( padj < 0.005 , abs(log2FoldChange) > 2 )
nrow(DEGIwIwvsAFAF)

DEGIwIwvsAFAF <- DEGIwIwvsAFAF[order(DEGIwIwvsAFAF[,2]), ]
DEGIwIwvsAFAF <- DEGIwIwvsAFAF[-c(grep("Gm", DEGIwIwvsAFAF[,1])), ]
nrow(DEGIwIwvsAFAF)

#DKKnull versus AFAF
DEGDKKvsAFAF <- cbind.data.frame(rownames(resDKKvsAFAF), resDKKvsAFAF$log2FoldChange, resDKKvsAFAF$padj)
colnames(DEGDKKvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGDKKvsAFAF <- DEGDKKvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGDKKvsAFAF)

DEGDKKvsAFAF <- DEGDKKvsAFAF[order(DEGDKKvsAFAF[,2]), ]
DEGDKKvsAFAF <- DEGDKKvsAFAF[-c(grep("Gm", DEGDKKvsAFAF[,1]),
                                grep("-ps", DEGDKKvsAFAF[,1])), ]
nrow(DEGDKKvsAFAF)



#FSI versus AFAF
DEGFSIvsAFAF <- cbind.data.frame(rownames(resFSIvsAFAF), resFSIvsAFAF$log2FoldChange, resFSIvsAFAF$padj)
colnames(DEGFSIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFSIvsAFAF <- DEGFSIvsAFAF %>% filter ( padj < 0.005 , abs(log2FoldChange) > 2 )
nrow(DEGFSIvsAFAF)

DEGFSIvsAFAF <- DEGFSIvsAFAF[order(DEGFSIvsAFAF[,2]), ]
DEGFSIvsAFAF <- DEGFSIvsAFAF[-c(grep("Gm", DEGFSIvsAFAF[,1]),
                                grep("-ps", DEGFSIvsAFAF[,1])), ]
nrow(DEGFSIvsAFAF)

#FASI versus AFAF
DEGFASIvsAFAF <- cbind.data.frame(rownames(resFASIvsAFAF), resFASIvsAFAF$log2FoldChange, resFASIvsAFAF$padj)
colnames(DEGFASIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFASIvsAFAF <- DEGFASIvsAFAF %>% filter ( padj < 0.005 , abs(log2FoldChange) > 2 )
nrow(DEGFASIvsAFAF)

DEGFASIvsAFAF <- DEGFASIvsAFAF[order(DEGFASIvsAFAF[,2]), ]
DEGFASIvsAFAF <- DEGFASIvsAFAF[-c(grep("Gm", DEGFASIvsAFAF[,1]),
                                grep("-ps", DEGFASIvsAFAF[,1])), ]
nrow(DEGFASIvsAFAF)




###############################################################################################
#I think Guangdun swapped R and P regiond for E7.0
#to prove this I checked  few genes
Mixl1 <- assay(dds2)[grep("T", rownames(dds2)),]
write.csv(Mixl1, "T_data.csv")

#to prove that I will try to do a correlation between his FPKM data and my TPM

FPKM_E7.0 <- read.table("GSE65924_E1.gene.expression.txt", h=T)

# check gene name
FPKM_E7.0 <- FPKM_E7.0 %>% filter (FPKM_E7.0$Genes %in% rownames(dds))
dds_FPKMGenes <- as.data.frame(assay(dds)) %>% filter (rownames(dds) %in% FPKM_E7.0$Genes)
rownames(dds_FPKMGenes) <- FPKM_E7.0$Genes

toplot <- cbind.data.frame(FPKM_E7.0$Genes, dds_FPKMGenes$Embryo_Pen_E7.0_1_4P, FPKM_E7.0$E1.4P)
toplot <- cbind.data.frame(FPKM_E7.0$Genes, dds_FPKMGenes$Embryo_Pen_E7.0_1_4R, FPKM_E7.0$E1.4P)

cor(dds_FPKMGenes$Embryo_Pen_E7.0_1_5R, FPKM_E7.0$E1.5P)

ggplot(toplot, aes(toplot[,2], toplot[,3]))+
  geom_point()

# seems that correlation is higher despite being close 


