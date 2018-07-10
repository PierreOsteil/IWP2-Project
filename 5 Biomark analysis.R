# Analysis of 5 Biomark run for IWP2 project
library(tidyverse)
library(ggplot2)
library(FactoMineR)
library(rafalib)

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/Biomark")



##################################
#function to clean Biomark data
sanitize_biomark <- function (biomark_csv_file) {
  path_to_qPCR_data <- biomark_csv_file
  
  biomark_raw <- read_csv(file=path_to_qPCR_data, col_names = TRUE, skip = 11) #read_csv is the tidyverse function
  biomark_named <- biomark_raw %>% rename( Gene = Name_1 ) %>% dplyr::select(ID, Name, Gene, Value, Call, Threshold)
  
  
  # Filtered out Failed calls and high Ct values
  biomark_filtered <- biomark_named %>% filter( Call != "Fail" & Value < 24 & !grepl("H2O",Name))
  
  
  # Filter out genes and samples with less than 30 good Ct values
  
  biomark_filtered <- biomark_filtered %>% filter(!(Gene %in% (biomark_filtered %>% count(Gene) %>% filter (n < 30) %>% .$Gene)))
  
  biomark_filtered <- biomark_filtered %>% filter(!(Name %in% (biomark_filtered %>% count(Name) %>% filter (n < 30) %>% .$Name)))
  
  # Calculate the housekeeping gene control using the mean of Gadph, Tbp, bActin, B2m and Gapdh
  
  biomark_hk <- biomark_filtered %>% filter( Gene %in% c("Gapdh", "Hsp90ab1", "Actb","B2m", "Gusb")) %>% group_by(Name) %>% summarise( hk_mean = (mean(Value)))
  biomark_dCt <- inner_join(biomark_filtered, biomark_hk, by = c("Name" = "Name")) %>% mutate(dCt = (-1*(Value-hk_mean)))
  biomark_dCt <- biomark_dCt %>% filter(!(Gene %in% c("PPC", "RTC","MGDC" )))
  
  # Spread the data into a matrix
  biomark_spread <- biomark_dCt %>% dplyr::select(Name, Gene, dCt) %>% spread(key = Gene, value= dCt) %>% replace(.,is.na(.), -15) 
  
  # Remove the three housekeeping genes
  biomark_spread_nohk <- biomark_spread %>% dplyr::select(-one_of(c("Gapdh", "Hsp90ab1", "Actb","B2m", "Gusb")))
  
}


###################################
#Open dataset
CL1 <- sanitize_biomark("1362020254 Cell Lineage 1.csv")
CL1 <- CL1[-c(grep("D7", CL1$Name)),]
CL2 <- sanitize_biomark("1362171189 Cell Lineage 2.csv")
CL3 <- sanitize_biomark("1362362584_DKKnull&SB  Cell Lineage 3.csv")
Irene <- CL3[53:80,]
CL3 <- CL3[-(53:80),]
TG1 <- sanitize_biomark("1361935250 TGF 1.csv")
TG1 <- TG1[-c(grep("D7", TG1$Name)),]
WN1 <- sanitize_biomark("1362020249 WNT 1.csv")
WN1 <- WN1[-c(grep("D7", WN1$Name)),]



##################################
#Merge plate
# Normalisation between the two plates

# first Match the names
# CL2 common samples with CL1
CL2$Name[(which(CL2$Name == "HN_DKK-_D0A"))] <- "IwF_D0A"
CL2$Name[(which(CL2$Name == "HN_DKK-_D0B"))] <- "IwF_D0B" #common to three CL plates
CL2$Name[(which(CL2$Name == "HN_DKK-_D0C"))] <- "IwF_D0C" #common to three CL plates

CL2$Name[(which(CL2$Name == "HN_DKK+_D0C"))] <- "IwIw_D0C" #common to three CL plates
CL2$Name[(which(CL2$Name == "HN_DKK+_D4A"))] <- "IwIw_D4A"
CL2$Name[(which(CL2$Name == "HN_DKK+_D4B"))] <- "IwIw_D4B"

CL2$Name[(which(CL2$Name == "HN_FIw_D0#2"))] <- "FIw_D0B"
CL2$Name[(which(CL2$Name == "HN_FIw_D0#3"))] <- "FIw_D0C"


# CL3 common samples with CL1
CL3$Name[(which(CL3$Name == "DKK-_D0_B"))] <- "IwF_D0B" #common to three CL plates
CL3$Name[(which(CL3$Name == "DKK-_D0_C"))] <- "IwF_D0C" #common to three CL plates
CL3$Name[(which(CL3$Name == "DKK-_D1_C"))] <- "IwF_D1C"

CL3$Name[(which(CL3$Name == "DKK+_D0_A"))] <- "IwIw_D0A" 
CL3$Name[(which(CL3$Name == "DKK+_D0_C"))] <- "IwIw_D0C" #common to three CL plates
CL3$Name[(which(CL3$Name == "DKK+_D1_A"))] <- "IwIw_D1A"
CL3$Name[(which(CL3$Name == "DKK+_D1_C"))] <- "IwIw_D1C"

CL3$Name[(which(CL3$Name == "YC8-_D0_B"))] <- "FF_D0B" 
CL3$Name[(which(CL3$Name == "YC8-_D0_C"))] <- "FF_D0C"
CL3$Name[(which(CL3$Name == "YC8-_D1_B"))] <- "FF_D1B"
CL3$Name[(which(CL3$Name == "YC8-_D1_C"))] <- "FF_D1C"


#CL1 and CL3 first
#normalise these using common samples 
normalisation_samples <- intersect(CL1[["Name"]], CL3[["Name"]])
normalisation_samples

normalisation_genes <- intersect(colnames(CL1), colnames(CL3))
normalisation_genes


CL1_norm_gene <- CL1[,c(which(colnames(CL1) %in% normalisation_genes))]
#CL1_norm_gene <- CL1_norm_gene %>% filter(!(Name %in% normalisation_samples))
CL3_norm_gene <- CL3[,c(which(colnames(CL3) %in% normalisation_genes))]
#CL3_norm_gene <- CL3_norm_gene %>% filter(!(Name %in% normalisation_samples))

#Normalisation by gene averages
CL3_norm <- CL3_norm_gene %>% filter(Name %in% normalisation_samples) %>% arrange(Name)
CL1_norm <- CL1_norm_gene %>% filter(Name %in% normalisation_samples) %>% arrange(Name)

delta_normalisation_samples <- map2_df(CL3_norm[2:73], CL1_norm[2:73], `-`)
delta_normalisation_samples[is.na(delta_normalisation_samples)] <- 0
normalisation_factors <- delta_normalisation_samples %>% summarise_all(funs(mean))
CL3_norm_gene[,2:73] <- sweep(CL3_norm_gene[,2:73], 2,unlist(normalisation_factors), FUN="-")
CL3_norm_gene <- CL3_norm_gene %>% filter(!(Name %in% normalisation_samples))


# attach the two data-sets
CL1_CL3_clean <- as.data.frame(bind_rows(CL1_norm_gene, CL3_norm_gene))

#save table
write.csv(CL1_CL3_clean , "CL1_CL3_clean .csv")


#PCAtable
Line1 <- "FF"
Line2 <- "IwF_"
Line3 <- "IwIw_"
Line4 <- "FASI"
Line5 <- "_FSI"
Line6 <- "DKKn"

TABLE <- rbind.data.frame(CL1_CL3_clean[c(grep(Line1, CL1_CL3_clean$Name)),],
                          #CL1_CL3_clean[c(grep(Line2, CL1_CL3_clean$Name)),],
                          #CL1_CL3_clean[c(grep(Line3, CL1_CL3_clean$Name)),]
                          #,CL1_CL3_clean[c(grep(Line4, CL1_CL3_clean$Name)),]
                          #,CL1_CL3_clean[c(grep(Line5, CL1_CL3_clean$Name, ignore.case = T)),],
                          CL1_CL3_clean[c(grep(Line6, CL1_CL3_clean$Name)),]
                          )

CellLine <- c(rep(Line1, 15),
             #rep(Line2, 15),
             #rep(Line3, 15)
             #,rep(Line4, 15)
             #,rep(Line5, 15),
             rep(Line6, 15)
             )
Day <- c(rep(c(rep(c("D0", "D1", "D2", "D3", "D4" ), each =3)), nrow(TABLE)/15))

TABLE <- cbind.data.frame(TABLE, CellLine, Day)

# PCAplot
mypar(1,2)
PCA_TOT=PCA(TABLE[,c(2:(length(colnames(TABLE))-2))] , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)
PCAcoord12 <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])

PCA_data <- cbind.data.frame(PCAcoord12, TABLE$CellLine, TABLE$Day)
colnames(PCA_data) <- c("PC1", "PC2", "CellLine", "Day")

centroids <- aggregate(cbind(PCAcoord[,1], PCAcoord[,2])~Day+CellLine, PCA_data,mean)
colnames(centroids) <- c("Day", "CellLine", "PC1", "PC2")

PCA <- ggplot(PCA_data, aes(PC1,PC2, colour=CellLine, shape = Day)) +
  geom_point(size=4) + 
  geom_path(data=centroids, aes(PC1, PC2, group = CellLine), size = 1, 
            arrow= arrow(angle = 20, type = "closed"), linetype = "dashed")+
  #scale_shape_manual(values=c())+
  #scale_colour_manual(values=c( "#C4D600" ,"#151F6D"))+
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  theme_bw()

PCA

######################################################
## Fold Chang + ANOVA + Tukey HSD on Day 4 Samples
ANOVA_run <- tibble(Gene = character(),Days = integer(), pval = double(), log2FoldChange = double())
sink('ANOVA+Tukey_Day_4.txt')
list_of_genes <- tibble(Gene=integer())

DOI <- "D0"

  for (i in 2:(length(colnames(TABLE))-2)){
    Genes = colnames(TABLE)[i]
    temp_aov <- TABLE %>% filter(Day == DOI) %>% dplyr::select (GOI = Genes, CellLine)
    try(aov_model <- aov(GOI ~ CellLine, data=temp_aov), silent = T)
    pvalue = summary(aov_model)[[1]][["Pr(>F)"]][1]
    if (!is.null(pvalue)){
      if (pvalue < 0.05) {
        print(Genes)
        print(Day)
        list_of_genes <- add_row(list_of_genes, Gene=i)
        print(TukeyHSD(aov_model))
      }}
    Mean_Dat <- aggregate(temp_aov[,1]~CellLine, temp_aov, mean)
    FC <- (Mean_Dat[1,2] - Mean_Dat[2,2]) #!!!!!!WARNING the table is in alphebetical order
    ANOVA_run <- add_row(ANOVA_run,Gene=Genes, pval = pvalue, Days = DOI, log2FoldChange = FC)
  }
sink()

ANOVA_run$Sig <- c(ANOVA_run$pval<0.05)
(ANOVA_run$Sig)


##################################
#VOlcano plot
###Volcano plot
library(tidyverse)
library(ggrepel)
#need a result table from DESeq2 
Layer <- read.csv("Gene Annotation Biomark.txt", h=T, sep = "\t")

Comp <- ANOVA_run

LayerPlot <- Layer$Layer[Layer$Gene %in% Comp$Gene]

gene_list <- cbind.data.frame(Comp$Gene, Comp$log2FoldChange, Comp$pval, LayerPlot)

colnames(gene_list) <- c("Gene","logFC", "adj.P.Val", "Layer")
rownames(gene_list) <- rownames(Comp)

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$adj.P.Val < 0.05)
gene_list[is.na(gene_list)] <- 0

gene_list$threshold= as.factor(gene_list$threshold)
gene_list <- gene_list %>% filter(gene_list$threshold == TRUE)


##Construct the plot object
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=Layer)) +
  geom_point(alpha=0.4, size=1.75) +
  scale_colour_manual(values= c("blue", "green", "red", "purple"))+
  #xlim(c(-3, 3)) + 
  ylim(c(0, NA)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  geom_text_repel(aes(label = gene_list$Gene, size = 3))+
  geom_abline(slope=0,intercept=0.5)+
  geom_hline(yintercept= 0.5, linetype= "dashed", colour = "darkgrey")+
  geom_vline(xintercept = -1, linetype= "dashed", colour = "darkgrey")+
  geom_vline(xintercept = 1, linetype= "dashed", colour = "darkgrey")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background = element_blank())

g

#pie(c(14,5,2,1), labels = c("14","5","2","1"), col = c("blue", "red", "green","purple"))
#pie(c(4,10,8,5), labels = c("4","10","8","5"), col = c("blue", "red", "green","purple"))


#######################################################################
#Cl1, TG1 and WN1 for heatmap.3
#common samples only
WN1_clean <- WN1 %>% filter(WN1$Name %in% CL1$Name)
CL1_WN1_clean <- as.data.frame(bind_cols(CL1, WN1_clean[,-1]))

CL1_WN1_clean <- CL1_WN1_clean %>% filter(CL1_WN1_clean$Name %in% TG1$Name)
CL1_WN1_clean$Name == TG1_clean$Name

CL1_WN1_TG1_clean <- as.data.frame(bind_cols(CL1_WN1_clean, TG1_clean[,-1]))
