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

colnames(TABLE) <- c(gsub("Chang_EpiSC.sf","EpiSCs_Cha_F-F-", colnames(TABLE)))
colnames(TABLE) <- c(gsub("Chang_EpiSC-SF1.sf","EpiSCs_Cha_F-SB_1", colnames(TABLE)))
colnames(TABLE) <- c(gsub("Chang_EpiSC-SF2.sf","EpiSCs_Cha_F-SB_2", colnames(TABLE)))

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
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
PengDevCell <- read.table("DevCell_Embryo1_issue2.csv", sep= ",")$V12
PengDevCell_dat <- t(data.frame(TABLE[, c(grep("SRR", colnames(TABLE)))]))
colnames(PengDevCell_dat) <- PengDevCell[-1]
woPengDevCell_dat <- data.frame(TABLE[, -c(grep("SRR", colnames(TABLE)))])
TABLE_name <- c(colnames(t(woPengDevCell_dat)),colnames(PengDevCell_dat))

#Change the name
Lab_stage <- substring(TABLE_name, 12, 17)

##Allocate the new names to the samples
SampleTable <- cbind.data.frame(substring(TABLE_name,1,2), Lab_stage)
rownames(SampleTable) <- TABLE_name
colnames(SampleTable) <- c("TissueType", "LabStageAP")
colnames(txi$counts) <- rownames(SampleTable)


#Dds from txi files
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi, SampleTable, ~LabStageAP) #extract txi data



####################################################################
##DEGs for different EpiSCs 
Lab_Cond <- substring(TABLE_name, 8, 15)
##Allocate the new names to the samples
SampleTable <- cbind.data.frame(substring(TABLE_name,1,2), Lab_Cond)
rownames(SampleTable) <- TABLE_name
colnames(SampleTable) <- c("TissueType", "Lab_Cond")
colnames(txi$counts) <- rownames(SampleTable)

dds <- DESeqDataSetFromTximport(txi, SampleTable, ~Lab_Cond) #extract txi data
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/DEGs vs AFAF")
save(dds,file="dds.RData")

dds_EpiSC <- dds[,grep("EpiSCs", colnames(dds))] #choose the samples to analyse in the cohort. 
dds_EpiSC <- dds_EpiSC[,-c(grep("Ost_AFAF_", colnames(dds_EpiSC)), #to remove samples from the first sequencing. 
                     grep("Ost_IwIw_", colnames(dds_EpiSC)),
                     grep("Ost_IwAF1", colnames(dds_EpiSC)),
                     grep("Ost_IwAF2", colnames(dds_EpiSC)),
                     grep("Ost_IwAF3", colnames(dds_EpiSC)),
                     grep("Ost_AFIw", colnames(dds_EpiSC))
                     )]

#DESeq2
dds_EpiSC$Lab_Cond <- droplevels(dds_EpiSC$Lab_Cond) #Remove sample annotation for unwanted samples


dds2 <- DESeq(dds_EpiSC) # shouldn't take more than a few minute

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/DEGs vs AFAF")
save(dds2,file="dds2.RData")
load("dds2.RData")


resIwIwvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_IwIw", "Ost_AFAF"))
resDKKvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_DKKn", "Ost_AFAF"))
resFSIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_F-SI", "Ost_AFAF"))
resFASIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_FASI", "Ost_AFAF"))
resIwAFvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_IwAF", "Ost_AFAF"))



#extracting DEGs
#IwIw versus AFAF
DEGIwIwvsAFAF <- cbind.data.frame(rownames(resIwIwvsAFAF), resIwIwvsAFAF$log2FoldChange, resIwIwvsAFAF$padj)
colnames(DEGIwIwvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGIwIwvsAFAF <- DEGIwIwvsAFAF %>% filter ( padj < 0.05, abs(log2FoldChange) > 1 )
nrow(DEGIwIwvsAFAF)

DEGIwIwvsAFAF <- DEGIwIwvsAFAF[order(DEGIwIwvsAFAF[,2]), ]
DEGIwIwvsAFAF <- DEGIwIwvsAFAF[-c(grep("Gm", DEGIwIwvsAFAF[,1]),
                                  grep("-ps", DEGIwIwvsAFAF[,1]),
                                  grep("Olfr", DEGIwIwvsAFAF[,1])), ]

nrow(DEGIwIwvsAFAF)

IwIwgene <- DEGIwIwvsAFAF %>% filter( log2FoldChange > 0)
write.csv(IwIwgene, "IwIw_up_p05_newSeq.csv")
AFAFgene <- DEGIwIwvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "IwIw_down_p05_newSeq.csv")


#DKKnull versus AFAF
DEGDKKvsAFAF <- cbind.data.frame(rownames(resDKKvsAFAF), resDKKvsAFAF$log2FoldChange, resDKKvsAFAF$padj)
colnames(DEGDKKvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGDKKvsAFAF <- DEGDKKvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGDKKvsAFAF)

DEGDKKvsAFAF <- DEGDKKvsAFAF[order(DEGDKKvsAFAF[,2]), ]
DEGDKKvsAFAF <- DEGDKKvsAFAF[-c(grep("Gm", DEGDKKvsAFAF[,1]),
                                grep("-ps", DEGDKKvsAFAF[,1])), ]
nrow(DEGDKKvsAFAF)

DKKgene <- DEGDKKvsAFAF %>% filter( log2FoldChange > 0)
write.csv(DKKgene, "DKK_up_p05_newSeq.csv")
AFAFgene <- DEGDKKvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "DKK_down_p05_newSeq.csv")


#FSI versus AFAF
DEGFSIvsAFAF <- cbind.data.frame(rownames(resFSIvsAFAF), resFSIvsAFAF$log2FoldChange, resFSIvsAFAF$padj)
colnames(DEGFSIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFSIvsAFAF <- DEGFSIvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGFSIvsAFAF)

DEGFSIvsAFAF <- DEGFSIvsAFAF[order(DEGFSIvsAFAF[,2]), ]
DEGFSIvsAFAF <- DEGFSIvsAFAF[-c(grep("Gm", DEGFSIvsAFAF[,1]),
                                grep("-ps", DEGFSIvsAFAF[,1])), ]
nrow(DEGFSIvsAFAF)

FSIgene <- DEGFSIvsAFAF %>% filter( log2FoldChange > 0)
write.csv(FSIgene, "FSI_up_p05_newSeq.csv")
AFAFgene <- DEGFSIvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "FSI_down_p05_newSeq.csv")


#FASI versus AFAF
DEGFASIvsAFAF <- cbind.data.frame(rownames(resFASIvsAFAF), resFASIvsAFAF$log2FoldChange, resFASIvsAFAF$padj)
colnames(DEGFASIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFASIvsAFAF <- DEGFASIvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGFASIvsAFAF)

DEGFASIvsAFAF <- DEGFASIvsAFAF[order(DEGFASIvsAFAF[,2]), ]
DEGFASIvsAFAF <- DEGFASIvsAFAF[-c(grep("Gm", DEGFASIvsAFAF[,1]),
                                  grep("-ps", DEGFASIvsAFAF[,1])), ]
nrow(DEGFASIvsAFAF)

FASIgene <- DEGFASIvsAFAF %>% filter( log2FoldChange > 0)
write.csv(FASIgene, "FASI_up_p05_newSeq.csv")
AFAFgene <- DEGFASIvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "FASI_down_p05_newSeq.csv")


#IwAF versus AFAF
DEGIwAFvsAFAF <- cbind.data.frame(rownames(resIwAFvsAFAF), resIwAFvsAFAF$log2FoldChange, resIwAFvsAFAF$padj)
colnames(DEGIwAFvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGIwAFvsAFAF <- DEGIwAFvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGIwAFvsAFAF)

DEGIwAFvsAFAF <- DEGIwAFvsAFAF[order(DEGIwAFvsAFAF[,2]), ]
#DEGIwAFvsAFAF <- DEGIwAFvsAFAF[-c(grep("Gm", DEGIwAFvsAFAF[,1]),
 #                                grep("-ps", DEGIwAFvsAFAF[,1]))
  #                            ]
nrow(DEGIwAFvsAFAF)

IwAFgene <- DEGIwAFvsAFAF %>% filter( log2FoldChange > 0)
write.csv(IwAFgene, "IwAF_up_p05_newSeq.csv")
AFAFgene <- DEGIwAFvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "IwAF_down_p05_newSeq.csv")


##########################################################
# before need to reannotate the gene lists
# This function can take any of the columns(org.Hs.eg.db) as type and keys as long as the row names are in the format of the keys argument
getMatrixWithSelectedIds <- function(df, type, keys){
  require("AnnotationDbi")
  require("org.Mm.eg.db")
  
  geneSymbols <- mapIds(org.Mm.eg.db, keys=rownames(df), column=type, keytype=keys, multiVals="first")
  
  # get the entrez ids with gene symbols i.e. remove those with NA's for gene symbols
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  
  # subset your data frame based on the found_genes
  df2 <- df[names(found_genes), ]
  rownames(df2) <- found_genes
  return(df2)
}

# for example, going from SYMBOL to ENTREZID
ddsentrenz <- getMatrixWithSelectedIds(assay(dds2), type="ENTREZID", keys="SYMBOL")
ddsentrenz <- as.matrix(ddsentrenz)
head(assay(dds))
head(ddsentrenz)
ddse<-as.data.frame(ddsentrenz)

#GO analysis
library("qusage")

setwd("C:/Users/Pierre/Desktop/IWP2 Project")
load("mouse_c2_v5p2.rdata") #To download here: http://bioinf.wehi.edu.au/software/MSigDB/index.html

#or the curated list of genes http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2
c5<-read.gmt("c5.bp.v6.1.symbols.gmt")
c2<-read.gmt("c2.all.v6.1.symbols.gmt")

# Have a look at the first few gene sets
names(Mm.c2)[1:5]

# Number of gene sets in C2
length(Mm.c2)


c2.ind <- ids2indices(c2, identifiers = rownames(ddse))
head(c2.ind)
c5.ind <- ids2indices(c5, rownames(ddse))
head(c5.ind)


########################################################
dds2_test <- dds_EpiSC[,c(grep("Ost_AFAF", colnames(dds_EpiSC)), #to remove samples from the first sequencing. 
                           grep("Ost_IwIw", colnames(dds_EpiSC))
                          )]
Lab_Cond_dds2_test <- substring(colnames(assay(dds2_test)), 8, 15)
Lab_Cond_dds2_test


design <- model.matrix(~Lab_Cond_dds2_test)
design

gst.camera.c2 <- camera(assay(dds2_test),index=c2.ind,  design=design,  contrast=2,  inter.gene.cor=0.05)
gst.camera.c2[1:5,]
gst.camera.c5 <- camera(assay(dds2_test),index=c5.ind,  design=design,  contrast=2,  inter.gene.cor=0.05)
gst.camera[1:5,]


