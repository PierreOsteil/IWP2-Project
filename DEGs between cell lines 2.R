library(DESeq2)
library(tidyverse)


setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
load("dds.RData")

dds_EpiSC <- dds[,grep("EpiSCs", colnames(dds))]
colnames(dds_EpiSC)

#DESeq2

dds_EpiSC$LabStage <- droplevels(dds_EpiSC$LabStage)

DES <- DESeq(dds_EpiSC) # shouldn't take more than a few minute



#extract results
resIwIwvsAFAF <- results(DES, contrast = c("LabStage", "Ost_IwIw", "Ost_AFAF"))
resDKKvsAFAF <- results(DES, contrast = c("LabStage", "Ost_DKKn", "Ost_AFAF"))
resFSIvsAFAF <- results(DES, contrast = c("LabStage", "Ost_F-SI", "Ost_AFAF"))
resFASIvsAFAF <- results(DES, contrast = c("LabStage", "Ost_FASI", "Ost_AFAF"))
resIwAFvsAFAF <- results(DES, contrast = c("LabStage", "Ost_IwAF", "Ost_AFAF"))

resWu_IwIwvsAFAF <- results(DES, contrast = c("LabStage", "Wu._IwIw", "Wu._AFAF"))


setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/DEGs vs AFAF")
#extracting DEGs
#function
DEGs_lists <- function(res, pval, lFC, line) {
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


