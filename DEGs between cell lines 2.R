library(DESeq2)
library(tidyverse)


setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
load("dds.RData")

dds_EpiSC <- dds[,grep("EpiSCs", colnames(dds))]
colnames(dds_EpiSC)
nrow(dds_EpiSC)
dds_EpiSC <- dds_EpiSC[-c(grep("Gm", rownames(dds_EpiSC))),] #remove fucking Gm genes
nrow(dds_EpiSC)


#DESeq2

dds_EpiSC$LabStage <- droplevels(dds_EpiSC$LabStage)

DES <- DESeq(dds_EpiSC) # shouldn't take more than a few minute

#extract results
resIwIwvsAFAF <- results(DES, contrast = c("LabStage", "Ost_IwIw", "Ost_AFAF"))
resDKKvsAFAF <- results(DES, contrast = c("LabStage", "Ost_DKKn", "Ost_AFAF"))
resFSIvsAFAF <- results(DES, contrast = c("LabStage", "Ost_F-SI", "Ost_AFAF"))
resFASIvsAFAF <- results(DES, contrast = c("LabStage", "Ost_FASI", "Ost_AFAF"))
resIwAFvsAFAF <- results(DES, contrast = c("LabStage", "Ost_IwAF", "Ost_AFAF"))
resIwAFvsIwIw <- results(DES, contrast = c("LabStage", "Ost_IwAF", "Ost_IwIw"))
resWu_IwIwvsAFAF <- results(DES, contrast = c("LabStage", "Wu._IwIw", "Wu._AFAF"))
resKu_AFIwvsAFAF <- results(DES, contrast = c("LabStage", "Kur_AFIw", "Kur_AFAF"))
resCha_FSIvsAFAF <- results(DES, contrast = c("LabStage", "Cha_F-SB", "Cha_F-F-"))

resFSIvsFASI <- results(DES, contrast = c("LabStage", "Ost_F-SI", "Ost_FASI"))

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/DEGs vs AFAF")

#extracting DEGs
#function
DEGs_lists <- function(res, pval, lFC, line) {
  require(tidyverse)
  DEG <- cbind.data.frame(rownames(res), res$log2FoldChange, res$padj)
  colnames(DEG) <- c("Gene", "log2FoldChange", "padj")
  DEG <- DEG %>% filter ( padj < pval, abs(log2FoldChange) > lFC )
  
  DEG <- DEG[order(DEG[,2]), ]
  DEG <- DEG %>% filter(! str_detect (Gene, c("Gm"))) %>% filter(! str_detect (Gene, c("-ps"))) %>% filter(! str_detect (Gene, c("Olfr")))
  
  UPgene <- DEG %>% filter( log2FoldChange > 0)
  write.csv(UPgene, file = paste(line, "UP.csv"))
  DOWNgene <- DEG %>% filter( log2FoldChange < 0)
  write.csv(DOWNgene, file = paste(line,"DOWN.csv"))
  return(nrow(DEG))
}

#lists of DEGs
DEGs_lists(res = resIwIwvsAFAF, pval = 0.05, lFC= 1, line = "IwIwvsAFAF")
DEGs_lists(res = resIwAFvsAFAF, pval = 0.05, lFC= 1, line = "IwAFvsAFAF")
DEGs_lists(res = resDKKvsAFAF, pval = 0.05, lFC= 1, line = "DKKnvsAFAF")
DEGs_lists(res = resFSIvsAFAF, pval = 0.05, lFC= 1, line = "F-SIvsAFAF")
DEGs_lists(res = resFASIvsAFAF, pval = 0.05, lFC= 1, line = "FASIvsAFAF")
DEGs_lists(res = resWu_IwIwvsAFAF , pval = 0.05, lFC= 1, line = "Wu_IwIwvsAFAF")
DEGs_lists(res = resKu_AFIwvsAFAF , pval = 0.05, lFC= 1, line = "Ku_AFIwvsAFAF")
DEGs_lists(res = resCha_FSIvsAFAF , pval = 0.05, lFC= 1, line = "Cha_F-SIvsAFAF")
DEGs_lists(res = resFSIvsFASI , pval = 0.05, lFC= 1, line = "FSIvsFASI")
DEGs_lists(res = resIwIwvsIwAF , pval = 0.05, lFC= 1, line = "IwIwvsIwAF")


##########################################################
#GO analysis
library("qusage")

setwd("C:/Users/Pierre/Desktop/IWP2 Project")
load("mouse_c2_v5p2.rdata") #To download here: http://bioinf.wehi.edu.au/software/MSigDB/index.html

#or the curated list of genes http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2
c5<-read.gmt("c5.bp.v6.1.symbols.gmt")

# Have a look at the first few gene sets
names(Mm.c2)[1:5]

# Number of gene sets in C2
length(Mm.c2)

rownames(dds_EpiSC) <- toupper(rownames(dds_EpiSC)) #need to capitalise the gene name for GO to work

c5.ind <- ids2indices(c5, rownames(dds_EpiSC)) #allocate genes to each ontology
head(c5.ind)



########################################################
#GO results
dds2_test <- dds_EpiSC[,c(grep("Cha_F-F-", colnames(dds_EpiSC)), 
                           grep("Cha_F-SB", colnames(dds_EpiSC))
                          )]

Lab_Cond_dds2_test <- substring(colnames(assay(dds2_test)), 8, 15)
#Lab_Cond_dds2_test

design <- model.matrix(~Lab_Cond_dds2_test)
#design

IwIwvsAFAF_GO <- camera(assay(dds2_test), index=c5.ind,  design=design,  contrast=2,  inter.gene.cor=0.05)

write.csv(IwIwvsAFAF_GO, file = "IwIwvsIwAF_GO.csv")

dev <- grep(c("DEVELOPMENT|DIFFERENTIATION|GENESIS"),names(c5.ind)) #only developmental one
#names(c5.ind)[dev]
GO.dev <- camera(assay(dds2_test),index=c5.ind[dev],  design=design,  contrast=2,  inter.gene.cor=0.05)
write.csv(GO.dev, file = "ChaFFvsVSB_GO-DEV.csv")



###############################
#Barcode plot
GO_data <- resCha_FSIvsAFAF
GO_toplot <- c5.ind$GO_GASTRULATION
title <- "GO_GASTRULATION"

barcodeplot(GO_data$log2FoldChange, GO_toplot, main= title , quantiles= c(-1,1))




NEURO <- c5.ind$GO_NEUROGENESIS
barcodeplot(resWu_IwIwvsAFAF$log2FoldChange, NEURO, main="NEUROGENESIS", quantiles= c(-1,1))
HIND <- c5.ind$GO_HINDBRAIN_DEVELOPMENT
barcodeplot(resWu_IwIwvsAFAF$log2FoldChange, HIND, main="HINDBRAIN DEVELOPMENT", quantiles= c(-1,1))
ECTO <- c5.ind$GO_ECTODERM_DEVELOPMENT
barcodeplot(resDKKvsAFAF$log2FoldChange, ECTO, main="ECTODERM_DEVELOPMENT", quantiles= c(-1,1))
ENDO <- c5.ind$GO_ENDODERM_DEVELOPMENT
barcodeplot(resWu_IwIwvsAFAF$log2FoldChange, ENDO, main="ENDODERM_DEVELOPMENT", quantiles= c(-1,1))
MESO <- c5.ind$GO_MESODERM_DEVELOPMENT
barcodeplot(resWu_IwIwvsAFAF$log2FoldChange, MESO, main="MESODERM_DEVELOPMENT", quantiles= c(-1,1))
WNT <- c5.ind$GO_WNT_SIGNALING_PATHWAY
barcodeplot(resIwIwvsAFAF$log2FoldChange, WNT, main="WNT_SIGNALING_PATHWAY", quantiles= c(-1,1))
GAS <- c5.ind$GO_GASTRULATION
barcodeplot(resIwIwvsAFAF$log2FoldChange, WNT, main="WNT_SIGNALING_PATHWAY", quantiles= c(-1,1))
