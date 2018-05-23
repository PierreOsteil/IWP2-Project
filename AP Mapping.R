# This code will try to obtain the signature genes for each region
# the idea would be to take the Anterior versus POsterior for each day of development
# First I build a table of the mean level of expression for each layer
library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(corrplot)

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/Salmon")

Mapping_dat <- read.csv("TABLEALL.csv", h=T, row.names=1, sep = ",")
#Mapping_dat <- read.csv("TABLEALL_G400.csv", h=T, row.names=1, sep = ",")
Mapping_dat_Peng <- Mapping_dat[,c(grep("Pen", colnames(Mapping_dat)))]

#sample annotation
APaxis <- substring(colnames(Mapping_dat_Peng), 
                    nchar(colnames(Mapping_dat_Peng))-1, nchar(colnames(Mapping_dat_Peng)))


#collect data
Ant_dat <- Mapping_dat_Peng[,c(grep("A", APaxis))]

Ant_Stage <- substring(colnames(Ant_dat), 12, 15)


Ant5.5_dat_med <- apply(Ant_dat[,grep("E5.5", colnames(Ant_dat))], 1, FUN = mean )
Ant6.0_dat_med <- apply(Ant_dat[,grep("E6.0", colnames(Ant_dat))], 1, FUN = mean )
Ant6.5_dat_med <- apply(Ant_dat[,grep("E6.5", colnames(Ant_dat))], 1, FUN = mean )
Ant7.0_dat_med <- apply(Ant_dat[,grep("E7.0", colnames(Ant_dat))], 1, FUN = mean )
Ant7.5_dat_med <- apply(Ant_dat[,grep("E7.5", colnames(Ant_dat))], 1, FUN = mean )

Pos_dat <- cbind.data.frame(Mapping_dat_Peng[,c(grep("P", APaxis))] ,
                            Mapping_dat_Peng[,c(grep("B", APaxis))])

Pos5.5_dat_med <- apply(Pos_dat[,grep("E5.5", colnames(Pos_dat))], 1, FUN = mean )
Pos6.0_dat_med <- apply(Pos_dat[,grep("E6.0", colnames(Pos_dat))], 1, FUN = mean )
Pos6.5_dat_med <- apply(Pos_dat[,grep("E6.5", colnames(Pos_dat))], 1, FUN = mean )
Pos7.0_dat_med <- apply(Pos_dat[,grep("E7.0", colnames(Pos_dat))], 1, FUN = mean )
Pos7.5_dat_med <- apply(Pos_dat[,grep("E7.5", colnames(Pos_dat))], 1, FUN = mean )


AP_dat <- cbind.data.frame(Ant5.5_dat_med, Ant6.0_dat_med, Ant6.5_dat_med, Ant7.0_dat_med, Ant7.5_dat_med,
                           Pos5.5_dat_med, Pos6.0_dat_med, Pos6.5_dat_med, Pos7.0_dat_med, Pos7.5_dat_med)

#check the data on a PCA
mypar(1,2)

PCA_TOT=PCA(t(AP_dat) , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)

PCA_data <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])

colnames(PCA_data) <- c("PC1", "PC2")

PCA <- ggplot(PCA_data, aes(PC1, PC2)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  geom_text_repel(aes(label=substring(colnames(AP_dat),1,6)))+
  theme_bw()
PCA


##correlation with rest of the data

Cor_to_plot <- Mapping_dat[, c(grep("Ost_IwIw", colnames(Mapping_dat)))]
Cor_to_plot2 <- Mapping_dat[, c(grep("Ost_AFAF", colnames(Mapping_dat)))]

Cor_dat <- cbind(AP_dat, Cor_to_plot, Cor_to_plot2)

SpearCor <- cor(Cor_dat, method = "spearman")

x11()
corrplot(SpearCor, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = SpearCor, col = col, symm = TRUE)


##Guilt by association method
#E5.5

AP_E5.5 <- cbind.data.frame(Ant5.5_dat_med, Pos5.5_dat_med)

AP_dat <- AP_E5.5

PCA_TOT=PCA(t(AP_dat) , scale.unit=T,ncp=5, axes = c(1,2))
PCAcoord <- as.data.frame(PCA_TOT$ind)
PCA_data <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])
colnames(PCA_data) <- c("PC1", "PC2")

PCA <- ggplot(PCA_data, aes(PC1, PC2)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  geom_text_repel(aes(label=substring(colnames(AP_dat),1,6)))+
  theme_bw()
PCA





