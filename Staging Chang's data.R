setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
load("dds.RData")


dds <- dds[ rowSums(counts(dds)) > 1, ]

dds_Chang_Frag <- dds[,grep("Embryo_Cha", colnames(dds))]
dds_Chang_Frag <- dds_Chang_Frag[,-c(grep("ExE", colnames(dds_Chang_Frag))
                                    ,grep("EX", colnames(dds_Chang_Frag))
                                    ,grep("V", colnames(dds_Chang_Frag)))]
colnames(dds_Chang_Frag)


#Negative binomiale normalisation
library("edgeR")
z_dds <- apply(assay(dds_Chang_Frag) , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))
head(z_dds)



GenesOnPCs = as.data.frame(predict(pca)) #extract genes contributing to each PC``
ggplot(GenesOnPCs, aes(GenesOnPCs[,2], GenesOnPCs[,3]))+
  geom_point(size=3) +
  geom_text_repel(aes(label=rownames(GenesOnPCs)))


axis2 <- GenesOnPCs[order(GenesOnPCs$PC2),] #PC2 is the best at predicting timing
EarlyGenes <- rownames(axis2[c(1:100),]) #top100 genes left PC2
LateGenes <- rownames(axis2[c(401:500),]) #top100 genes right PC2



#Control
all(rownames(coldata_allnoExE) == colnames(cts_allnoExE)) #has to be TRUE


dds <- DESeqDataSetFromMatrix(countData = cts_allnoExE,
                              colData = coldata_allnoExE,
                              design = ~ Sequencing)
#Prefiltering
dds <- dds[ rowSums(counts(dds)) > 1, ]

#log transform
rld <- rlog(dds, blind=FALSE)
rld_dat <- assay(rld)


############################################################################################
#Select gene of interest
rld_StageGene <- rld_dat[c(which(rownames(rld_dat) %in% c(paste(EarlyGenes), paste(LateGenes)))),]

#23 common genes
Common23 <- read.table("23 common genes.txt")
rld_Common23 <- rld_dat[c(which(rownames(rld_dat) %in% c(paste(Common23$V1)))),]

#PCA
pca_StageGene <- prcomp(rld_StageGene, center = T)
PCAdata <- as.data.frame(pca_StageGene$rotation)
Paper <- coldata_allnoExE$Paper
PCAdata <- cbind.data.frame(PCAdata, Paper)
ggplot(PCAdata, aes(PCAdata[,1], PCAdata[,2]))+
  geom_point(aes(colour = Paper), 
             size=3) +
  geom_text_repel(aes(label=colnames(rld_StageGene)))+
  xlab(paste0("PC1: ",round((pca_StageGene$sdev[1])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  ylab(paste0("PC2: ",round((pca_StageGene$sdev[2])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  theme(aspect.ratio=1)


#correlation heatmap
sampleDists <- dist(t(rld_Common23))

library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Spearman correlation
SpearCor <- cor(rld_StageGene, method = "spearman")
library("corrplot")
corrplot(SpearCor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = SpearCor, col = col, symm = TRUE)

corEpiSC <- SpearCor[-c(1:10),1:10] # just cell lines versus EpISC fragments
write.csv(corEpiSC, "corEpiSC2.csv")


#TSNE
library("Rtsne")

TSNE1 <- Rtsne(rld_StageGene, perplexity = 1) #cluster genes

ggplot(as.data.frame(TSNE1$Y), aes(TSNE1$Y[,1], TSNE1$Y[,2]) )+
  geom_point(size=3) +
  geom_text_repel(aes(label=rownames(rld_StageGene)))

#Negative binomiale normalisation
library("edgeR")
z_rld <- apply(rld_dat , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))

#averaging replicates
z_rld_ave <- cbind.data.frame(apply(z_rld[,c(2:4)] , 1, function(x) mean(x)),
                              apply(z_rld[,c(5:7)] , 1, function(x) mean(x)),
                              apply(z_rld[,c(33:35)] , 1, function(x) mean(x)), 
                              z_rld[,c(8:17)])

colnames(z_rld_ave)[c(1:3)] <- c("EpiSC-I/AF", "EpiSC-I", "EpiSC-AF")

#heatmap
rld_Common23 <- rld_Common23[,c(6,7,25,26,27,28,34,35)]


path <- "C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/Gene list"
tarGeneFile <- paste(path, "/23 common genes.txt",  sep = "" )
tar <- read.table(tarGeneFile)
tar_dat <- z_rld_ave[which(rownames(z_rld_ave) %in% c(paste(tar$V1))),]
tar_dat <- tar_dat[match(paste(tar$V1), rownames(tar_dat)),match(c( "EpiSC-AF", "EpiSC-I/AF", "EpiSC-I"), colnames(tar_dat))]
tar_dat <- tar_dat[,-c(19:32)]
tar_dat[is.na(tar_dat)] <- 0

hmcols <- colorRampPalette(c("blue","white","red"), space = "Lab")(256)
heatmap(as.matrix(tar_dat), col =hmcols, cexCol = 0.5
        #,Rowv = NA, Colv= NA
)

##############################################################################
##Adding I transcriptome data

iTra_dat <- read.table( "all.samples.genes.fpkm.txt", h=T , row.names = 1, sep="\t")

rld_dat_complete <- cbind.data.frame(iTra_dat[which(rownames(iTra_dat) %in% rownames(rld_dat)),], 
                                     rld_dat[which(rownames(rld_dat) %in% rownames(iTra_dat)),])

#Select gene of interest
rld_StageGene <- rld_dat_complete[c(which(rownames(rld_dat_complete) %in% c(paste(EarlyGenes), paste(LateGenes)))),]

#PCA
pca_StageGene <- prcomp(rld_StageGene, center = T)
PCAdata <- as.data.frame(pca_StageGene$rotation)
Paper <- c(rep("iTrans", 101), coldata_allnoExE$Paper)
PCAdata <- cbind.data.frame(PCAdata, Paper)
ggplot(PCAdata, aes(PCAdata[,4], PCAdata[,5]))+
  geom_point(aes(colour = Paper), 
             size=3) +
  geom_text_repel(aes(label=colnames(rld_StageGene)))+
  xlab(paste0("PC4: ",round((pca_StageGene$sdev[1])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  ylab(paste0("PC5: ",round((pca_StageGene$sdev[2])^2 / sum(pca_StageGene$sdev^2)*100, digits = 2) )) +
  theme(aspect.ratio=1)
