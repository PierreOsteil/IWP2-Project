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

dds_temp <- dds[,grep("Embryo_Pen_E5.5_", colnames(dds))] #choose the embryo to analyse in the cohort.
dds_temp <- dds[,grep("Embryo_Pen_E6.0_", colnames(dds))]
dds_temp <- dds[,grep("Embryo_Pen_E6.5_", colnames(dds))]
dds_temp <- dds[,grep("Embryo_Pen_E7.0_1", colnames(dds))]  
dds_temp <- dds[,grep("Embryo_Pen_E7.5_3", colnames(dds))]

dds_temp$LabStageAP <- droplevels(dds_temp$LabStageAP) #Remove sample annotation for unwanted samples

dds2 <- DESeq(dds_temp) # shouldn't take more than a few minute
dds2 <- dds2[ rowSums(counts(dds2)) > 1, ]


res5.5 <- results(dds2, contrast = c("LabStageAP", "E5.5__ A", "E5.5__ P"))
res6.0 <- results(dds2, contrast = c("LabStageAP", "E6.0__ A", "E6.0__ B"))
res6.5 <- results(dds2, contrast = c("LabStageAP", "E6.5_ A", "E6.5_ P"))
res7.0 <- results(dds2, contrast = c("LabStageAP", "E7.0_1 A", "E7.0_1 P"))
res7.5 <- results(dds2, contrast = c("LabStageAP", "E7.5_3 A", "E7.5_3 P"))


#extracting DEGs

DEG5.5 <- cbind.data.frame(rownames(res5.5), res5.5$log2FoldChange, res5.5$pvalue)
colnames(DEG5.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG5.5 <- DEG5.5 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG5.5)


DEG6.0 <- cbind.data.frame(rownames(res6.0), res6.0$log2FoldChange, res6.0$pvalue)
colnames(DEG6.0) <- c("Gene", "log2FoldChange", "pvalue")
DEG6.0 <- DEG6.0 %>% filter ( pvalue < 0.05, abs(log2FoldChange) > 1 )
nrow(DEG6.0)


DEG6.5 <- cbind.data.frame(rownames(res6.5), res6.5$log2FoldChange, res6.5$pvalue)
colnames(DEG6.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG6.5 <- DEG6.5 %>% filter ( pvalue < 0.1 , abs(log2FoldChange) > 1 )
nrow(DEG6.5)


DEG7.0 <- cbind.data.frame(rownames(res7.0), res7.0$log2FoldChange, res7.0$pvalue)
colnames(DEG7.0) <- c("Gene", "log2FoldChange", "pvalue")
DEG7.0 <- DEG7.0 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG7.0)


DEG7.5 <- cbind.data.frame(rownames(res7.5), res7.5$log2FoldChange, res7.5$pvalue)
colnames(DEG7.5) <- c("Gene", "log2FoldChange", "pvalue")
DEG7.5 <- DEG7.5 %>% filter ( pvalue < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEG7.5)

##########################
#PCA staging

DEG_5stage <-rbind.data.frame(DEG5.5, DEG6.0, DEG6.5, DEG7.0, DEG7.5)


#then do select Peng's data
Mapping_dat_Peng <- dds[,c(grep("Pen", colnames(dds))
                        , grep("EpiSC", colnames(dds))
                        , grep("Embryo_Wu", colnames(dds))
                        , grep("Embryo_Cha", colnames(dds)))]
Mapping_dat_Peng <- Mapping_dat_Peng[,-c(grep("R", colnames(Mapping_dat_Peng)),
                                        grep("L", colnames(Mapping_dat_Peng)),
                                        grep("VE", colnames(Mapping_dat_Peng)),
                                        grep("EP", colnames(Mapping_dat_Peng)),
                                        grep("EA", colnames(Mapping_dat_Peng)),
                                        grep("M", substring(colnames(Mapping_dat_Peng),20,20))
                                        )]
Mapping_dat_Peng <- Mapping_dat_Peng[,c(grep("E5.5", colnames(Mapping_dat_Peng)),
                                        grep("E6.0", colnames(Mapping_dat_Peng)),
                                        grep("E6.5", colnames(Mapping_dat_Peng)),
                                        grep("E7.0_N", colnames(Mapping_dat_Peng)),
                                        grep("E7.5_3", colnames(Mapping_dat_Peng)),
                                        grep("EpiSC", colnames(Mapping_dat_Peng))
)]
dds_Peng <- as.data.frame(assay(Mapping_dat_Peng))
colnames(dds_Peng)
dds_Peng <- dds_Peng %>% filter (rownames(dds_Peng) %in% DEG_5stage$Gene)
stage <- substring(colnames(dds_Peng), 12, 15)

#binomiale transformation
library("edgeR")
z_dds_Peng <- apply(dds_Peng , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))

library("rafalib")
library("FactoMineR")
library("ggrepel")
library("ggplot2")

mypar(1,2)
#sample annotation
Stage <- substring(colnames(z_dds_Peng),12, 15)
Lab <- substring(colnames(z_dds_Peng), 8, 10)
SizeText <- c(rep("A", ncol(z_dds_Peng)))

####PCA

PCA_TOT=PCA(t(z_dds_Peng) , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)
PCAcoord12 <- cbind.data.frame(PCAcoord[,1], PCAcoord[,3])


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
  ylab(paste("PC3", "(",round(PCA_TOT$eig[3,2], 2), "% )"))+
  geom_text_repel(aes(label=colnames(z_dds_Peng), size=SizeText))+
  theme_bw()
PCA




############################################################################################
#extracting top genes
#7.5
DEG7.5 <- DEG7.5[order(DEG7.5[,2]), ]
DEG7.5 <- DEG7.5[-c(grep("Gm", DEG7.5[,1])), ]
nrow(DEG7.5)
ante_7.5 <- DEG7.5 %>% filter(log2FoldChange > 0)
post_7.5 <- DEG7.5 %>% filter(log2FoldChange < 0)

write.csv(ante_7.5, "ante_7.5.csv")
write.csv(post_7.5, "post_7.5.csv")

#7.0
DEG7.0 <- DEG7.0[order(DEG7.0[,2]), ]
DEG7.0 <- DEG7.0[-c(grep("Gm", DEG7.0[,1])), ]
nrow(DEG7.0)
ante_7.0 <- DEG7.0 %>% filter(log2FoldChange > 0)
post_7.0 <- DEG7.0 %>% filter(log2FoldChange < 0)

write.csv(ante_7.0, "ante_7.0.csv")
write.csv(post_7.0, "post_7.0.csv")

#6.5
DEG6.5 <- DEG6.5[order(DEG6.5[,2]), ]
DEG6.5 <- DEG6.5[-c(grep("Gm", DEG6.5[,1])), ]
nrow(DEG6.5)
ante_6.5  <- DEG6.5 %>% filter(log2FoldChange > 0)
post_6.5 <- DEG6.5 %>% filter(log2FoldChange < 0)

write.csv(ante_6.5, "ante_6.5.csv")
write.csv(post_6.5, "post_6.5.csv")

#E6.0
DEG6.0 <- DEG6.0[order(DEG6.0[,2]), ]
DEG6.0 <- DEG6.0[-c(grep("Gm", DEG6.0[,1])), ]
nrow(DEG6.0)
ante_6.0  <- DEG6.0 %>% filter(log2FoldChange > 0)
post_6.0 <- DEG6.0 %>% filter(log2FoldChange < 0)

write.csv(ante_6.0, "ante_6.0.csv")
write.csv(post_6.0, "post_6.0.csv")

#E5.5
DEG5.5 <- DEG5.5[order(DEG5.5[,2]), ]
DEG5.5 <- DEG5.5[-c(grep("Gm", DEG5.5[,1])), ]
nrow(DEG5.5)
ante_5.5  <- DEG5.5 %>% filter(log2FoldChange > 0)
post_5.5 <- DEG5.5 %>% filter(log2FoldChange < 0)

write.csv(ante_5.5, "ante_5.5.csv")
write.csv(post_5.5, "post_5.5.csv")







###PCA for mapping
dds_map <- as.data.frame(assay(dds))
dds_map <- dds_map %>% rownames_to_column("Gene") %>% filter (rownames(dds_map) %in% E7.5_zipcode$V1)%>% column_to_rownames("Gene")
                     
dds_map <- dds_map[,c(c(grep("E7.5_3", colnames(dds_map)))
                       ,c(grep("Ost", colnames(dds_map)))
                        )]
dds_map <- dds_map[,-c(
                          grep("C", substring(colnames(dds_map), nchar(colnames(dds_map)), 
                                                                   nchar(colnames(dds_map))))
                        , grep("D", colnames(dds_map))
                        , grep("EA", colnames(dds_map))
                        , grep("M", colnames(dds_map))
                        , grep("EP", colnames(dds_map))
                        , grep("AFAFB", colnames(dds_map))
                        , grep("AFAFC", colnames(dds_map))
                        , grep("En", colnames(dds_map))
                        )]
dds_map <- dds_map[ rowSums(dds_map) > 1, ]


#Negative binomiale normalisationfor PCA
library("edgeR")
z_dds_map <- apply(dds_map , 2, function(x) zscoreNBinom(x, size = 10, mu =mean(x)))

mypar(1,2)
PCA_TOT=PCA(t(dds_map) , scale.unit=T,ncp=5, axes = c(1,3))
PCAcoord <- as.data.frame(PCA_TOT$ind)
PCA_data <- cbind.data.frame(PCAcoord[,1], -PCAcoord[,2], substring(colnames(z_dds_map), 12, 15),
                                                         substring(colnames(z_dds_map), nchar(colnames(z_dds_map)), 
                                                                                        nchar(colnames(z_dds_map))))

colnames(PCA_data) <- c("PC1", "PC2", "Stage","Position")

PCA <- ggplot(PCA_data, aes(PC1, PC2, colour=Stage, shape=Position)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC3", "(",round(PCA_TOT$eig[3,2], 2), "% )"))+
  geom_text_repel(aes(label=substring(colnames(z_dds_map), 12, 15)))+
  theme_bw()
PCA



#Spearman correlation
CorIWIW <- cor(dds_map[,1:46], apply(dds_map[,c(25,26,27)+25], MARGIN = 1, mean), method = c("spearman"))
CorIWAF <- cor(dds_map[,1:46], apply(dds_map[,c(22,23,24)+25], MARGIN = 1, mean), method = c("spearman"))
CorAFAF <- cor(dds_map[,1:46], apply(dds_map[,c(41,42,43)+25], MARGIN = 1, mean), method = c("spearman"))
CorAFIW <- cor(dds_map[,1:46], apply(dds_map[,c(44,45,46)+25], MARGIN = 1, mean), method = c("spearman"))
CorFASI <- cor(dds_map[,1:46], apply(dds_map[,c(31,32,33)+25], MARGIN = 1, mean), method = c("spearman"))
CorF.SI <- cor(dds_map[,1:46], apply(dds_map[,c(34,35,36)+25], MARGIN = 1, mean), method = c("spearman"))
CorDKKn <- cor(dds_map[,1:46], apply(dds_map[,c(28,29,30)+25], MARGIN = 1, mean), method = c("spearman"))

#Cornplot E6.5
y6.5 <- c(1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
x6.5 <- c(0.4, 0.2,0.625, 0.150, 0.650, 0.125, 0.675, 0.1, 0.7 ,0.09, 0.71, 0.08, 0.72, 0.07, 0.73)

E6.5 <- cbind.data.frame(y6.5,x6.5, CorIWIW[,1])
E6.5 <- rbind.data.frame(E6.5, c(2,2,1), c(2,3,0.3))
colnames(E6.5) <- c("Layer", "Region", "Cor")

E6.5_plot <- ggplot(E6.5, aes(x = Region, y= Layer, colour = Cor, width= 25)) + 
  geom_point(size = 12) + 
  xlim(-0.25,1) + 
  ylim(1,8) + 
  scale_colour_gradientn(colours = c("lightgrey", "green", "yellow","red"), 
                         values = c(0,0.1,0.3,1)) +
  labs(title = paste("EpiSC-IwAF" ,"correlation to E6.5"))+
    xlab ("A                  P" )+
  guides(size=F, shape=F)+
  theme(axis.title = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        title = element_text(size=15),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
E6.5_plot


#Cornplot E7.0
y7.0 <- c(1, 1, 2,2,3,rep(c(4,5,6,7,8,9,10,11), each=2))
x7.0 <- c(3.5,1.5, 1,4, 4, rep(c(1,4), 8))

E7.0 <- cbind.data.frame(y7.0,x7.0, CorAFAF[,1])
#E7.0 <- rbind.data.frame(E7.0, c(2,5,1), c(1,5,0.3))
colnames(E7.0) <- c("Layer", "Region", "Cor")


E7.0_plot <- ggplot(E7.0, aes(x = Region, y= Layer, colour = Cor)) + 
  geom_point(size = 12) + 
  xlim(0.5,4.5) + 
  ylim(0,12)+
  scale_colour_gradientn(colours = c("lightgrey", "green", "yellow","red"), 
                         values = c(0,0.1,0.5,1))+
  labs(title = paste("EpiSC-IWIW"))+
  labs(x= "E7.0", ylab="")+
  xlab ("A        P" )+
  guides(size=F, shape=F)+
  theme(axis.title = element_text(size = 20),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        title = element_text(size=15),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
E7.0_plot


#CornPlot E7.5
y7.5 <- as.numeric(substring(colnames(dds_map[1:46]),19,19))
x7.5 <- c(c(2.5,4.5), c(1,6,2.5,4.5), c(rep(c(1,6,2,4,3,5), 6)), c(1,6,2.5,4.5))

E7.5 <- cbind.data.frame(y7.5,x7.5, CorDKKn[,1])
#E7.5 <- rbind.data.frame(E7.5, c(2,5,1), c(1,5,0.3))
colnames(E7.5) <- c("Layer", "Region", "Cor")


E7.5_plot <- ggplot(E7.5, aes(x = Region, y= Layer, colour = Cor)) + 
  geom_point(size = 12, stroke=1, coulour="black", stroke=1) + 
  xlim(0.5,6.5) + 
  ylim(0,10)+
  scale_colour_gradientn(colours = c("lightgrey", "green", "yellow","red"), 
                         values = c(0,0.1,0.5,1))+
  labs(title = paste("EpiSC-DKKn"))+
  labs(x= "E7.5", ylab="")+
  scale_y_continuous("", breaks = y7.5)+
  xlab ("A   L1   R1   L2   R2   P" )+
  guides(size=F, shape=F)+
  theme(axis.title = element_text(size = 18),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        title = element_text(size=15),
        legend.title=element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
E7.5_plot



####################################################################
##DEGs for different EpiSCs 
Lab_Cond <- substring(TABLE_name, 8, 15)
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
dds2 <- dds2[,-c(293,294)]


resIwIwvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_IwIw", "Ost_AFAF"))
resDKKvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_DKKn", "Ost_AFAF"))
resFSIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_F-SI", "Ost_AFAF"))
resFASIvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_FASI", "Ost_AFAF"))
resAFIwvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_AFIw", "Ost_AFAF"))
resIwAFvsAFAF <- results(dds2, contrast = c("Lab_Cond", "Ost_IwAF", "Ost_AFAF"))

#extracting DEGs
#IwIw versus AFAF
setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/DEGs vs AFAF")
DEGIwIwvsAFAF <- cbind.data.frame(rownames(resIwIwvsAFAF), resIwIwvsAFAF$log2FoldChange, resIwIwvsAFAF$padj)
colnames(DEGIwIwvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGIwIwvsAFAF <- DEGIwIwvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 1 )
nrow(DEGIwIwvsAFAF)

DEGIwIwvsAFAF <- DEGIwIwvsAFAF[order(DEGIwIwvsAFAF[,2]), ]
DEGIwIwvsAFAF <- DEGIwIwvsAFAF[-c(grep("Gm", DEGIwIwvsAFAF[,1]),
                                  grep("-ps", DEGIwIwvsAFAF[,1]),
                                  grep("Olfr", DEGIwIwvsAFAF[,1])), ]

nrow(DEGIwIwvsAFAF)

IwIwgene <- DEGIwIwvsAFAF %>% filter( log2FoldChange > 0)
write.csv(IwIwgene, "IwIw_up.csv")
AFAFgene <- DEGIwIwvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "IwIw_down.csv")


#DKKnull versus AFAF
DEGDKKvsAFAF <- cbind.data.frame(rownames(resDKKvsAFAF), resDKKvsAFAF$log2FoldChange, resDKKvsAFAF$padj)
colnames(DEGDKKvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGDKKvsAFAF <- DEGDKKvsAFAF %>% filter ( padj < 0.001 , abs(log2FoldChange) > 1 )
nrow(DEGDKKvsAFAF)

DEGDKKvsAFAF <- DEGDKKvsAFAF[order(DEGDKKvsAFAF[,2]), ]
DEGDKKvsAFAF <- DEGDKKvsAFAF[-c(grep("Gm", DEGDKKvsAFAF[,1]),
                                grep("-ps", DEGDKKvsAFAF[,1])), ]
nrow(DEGDKKvsAFAF)

DKKgene <- DEGDKKvsAFAF %>% filter( log2FoldChange > 0)
write.csv(DKKgene, "DKK_up_p001.csv")
AFAFgene <- DEGDKKvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "DKK_down_p001.csv")


#FSI versus AFAF
DEGFSIvsAFAF <- cbind.data.frame(rownames(resFSIvsAFAF), resFSIvsAFAF$log2FoldChange, resFSIvsAFAF$padj)
colnames(DEGFSIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFSIvsAFAF <- DEGFSIvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 2 )
nrow(DEGFSIvsAFAF)

DEGFSIvsAFAF <- DEGFSIvsAFAF[order(DEGFSIvsAFAF[,2]), ]
DEGFSIvsAFAF <- DEGFSIvsAFAF[-c(grep("Gm", DEGFSIvsAFAF[,1]),
                                grep("-ps", DEGFSIvsAFAF[,1])), ]
nrow(DEGFSIvsAFAF)

FSIgene <- DEGFSIvsAFAF %>% filter( log2FoldChange > 0)
write.csv(FSIgene, "FSI_up.csv")
AFAFgene <- DEGFSIvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "FSI_down.csv")


#FASI versus AFAF
DEGFASIvsAFAF <- cbind.data.frame(rownames(resFASIvsAFAF), resFASIvsAFAF$log2FoldChange, resFASIvsAFAF$padj)
colnames(DEGFASIvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGFASIvsAFAF <- DEGFASIvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 2 )
nrow(DEGFASIvsAFAF)

DEGFASIvsAFAF <- DEGFASIvsAFAF[order(DEGFASIvsAFAF[,2]), ]
DEGFASIvsAFAF <- DEGFASIvsAFAF[-c(grep("Gm", DEGFASIvsAFAF[,1]),
                                grep("-ps", DEGFASIvsAFAF[,1])), ]
nrow(DEGFASIvsAFAF)

FASIgene <- DEGFASIvsAFAF %>% filter( log2FoldChange > 0)
write.csv(FASIgene, "FASI_up.csv")
AFAFgene <- DEGFASIvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "FASI_down.csv")


#IwAF versus AFAF
DEGIwAFvsAFAF <- cbind.data.frame(rownames(resIwAFvsAFAF), resIwAFvsAFAF$log2FoldChange, resIwAFvsAFAF$padj)
colnames(DEGIwAFvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGIwAFvsAFAF <- DEGIwAFvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 2 )
nrow(DEGIwAFvsAFAF)

DEGIwAFvsAFAF <- DEGIwAFvsAFAF[order(DEGIwAFvsAFAF[,2]), ]
DEGIwAFvsAFAF <- DEGIwAFvsAFAF[-c(grep("Gm", DEGIwAFvsAFAF[,1]),
                                  grep("-ps", DEGIwAFvsAFAF[,1]))
                                   ]
nrow(DEGIwAFvsAFAF)

IwAFgene <- DEGIwAFvsAFAF %>% filter( log2FoldChange > 0)
write.csv(IwAFgene, "IwAF_up.csv")
AFAFgene <- DEGIwAFvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "IwAF_down.csv")


#AFIw versus AFAF
DEGAFIwvsAFAF <- cbind.data.frame(rownames(resAFIwvsAFAF), resAFIwvsAFAF$log2FoldChange, resAFIwvsAFAF$padj)
colnames(DEGAFIwvsAFAF) <- c("Gene", "log2FoldChange", "padj")
DEGAFIwvsAFAF <- DEGAFIwvsAFAF %>% filter ( padj < 0.05 , abs(log2FoldChange) > 2 )
nrow(DEGAFIwvsAFAF)

DEGAFIwvsAFAF <- DEGAFIwvsAFAF[order(DEGAFIwvsAFAF[,2]), ]
DEGAFIwvsAFAF <- DEGAFIwvsAFAF[-c(grep("Gm", DEGAFIwvsAFAF[,1]),
                                  grep("-ps", DEGAFIwvsAFAF[,1])), ]
nrow(DEGAFIwvsAFAF)

AFIwgene <- DEGAFIwvsAFAF %>% filter( log2FoldChange > 0)
write.csv(AFIwgene, "AFIw_up.csv")
AFAFgene <- DEGAFIwvsAFAF %>% filter( log2FoldChange < 0)
write.csv(AFAFgene, "AFIw_down.csv")

