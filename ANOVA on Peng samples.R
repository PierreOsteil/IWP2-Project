# This code will perform a ANOVA between Peng 5 stages to obtain the signature genes for each region
# the idea would be to take the Anterior versus POsterior for each day of development
# First I build a table of the median level of expression for each layer


Mapping_dat <- read.csv("TABLEALL_znorm.csv", h=T, row.names=1, sep = ",")
Mapping_dat <- Mapping_dat[,c(grep("Pen", colnames(Mapping_dat)))]

#sample annotation
APaxis <- substring(colnames(Mapping_dat), nchar(colnames(Mapping_dat))-1, nchar(colnames(Mapping_dat)))


#collect data
Ant_dat <- Mapping_dat[,c(grep("A", APaxis))]

Ant_Stage <- substring(colnames(Ant_dat), 12, 15)


Ant5.5_dat_med <- apply(Ant_dat[,grep("E5.5", colnames(Ant_dat))], 1, FUN = median )
Ant6.0_dat_med <- apply(Ant_dat[,grep("E6.0", colnames(Ant_dat))], 1, FUN = median )
Ant6.5_dat_med <- apply(Ant_dat[,grep("E6.5", colnames(Ant_dat))], 1, FUN = median )
Ant7.0_dat_med <- apply(Ant_dat[,grep("E7.0", colnames(Ant_dat))], 1, FUN = median )
Ant7.5_dat_med <- apply(Ant_dat[,grep("E7.5", colnames(Ant_dat))], 1, FUN = median )



Pos_dat <- cbind.data.frame(Mapping_dat[,c(grep("P", APaxis))] ,Mapping_dat[,c(grep("B", APaxis))])

Pos5.5_dat_med <- apply(Pos_dat[,grep("E5.5", colnames(Pos_dat))], 1, FUN = median )
Pos6.0_dat_med <- apply(Pos_dat[,grep("E6.0", colnames(Pos_dat))], 1, FUN = median )
Pos6.5_dat_med <- apply(Pos_dat[,grep("E6.5", colnames(Pos_dat))], 1, FUN = median )
Pos7.0_dat_med <- apply(Pos_dat[,grep("E7.0", colnames(Pos_dat))], 1, FUN = median )
Pos7.5_dat_med <- apply(Pos_dat[,grep("E7.5", colnames(Pos_dat))], 1, FUN = median )


AP_dat <- cbind.data.frame(Ant5.5_dat_med, Ant6.0_dat_med, Ant6.5_dat_med, Ant7.0_dat_med, Ant7.5_dat_med,
                           Pos5.5_dat_med, Pos6.0_dat_med, Pos6.5_dat_med, Pos7.0_dat_med, Pos7.5_dat_med)

mypar(1,2)

PCA_TOT=PCA(t(AP_dat) , scale.unit=T,ncp=5, axes = c(1,2))

PCAcoord <- as.data.frame(PCA_TOT$ind)

PCAcoord12 <- cbind.data.frame(PCAcoord[,2], PCAcoord[,2])

colnames(PCA_data) <- c("PC1", "PC2")

PCA <- ggplot(PCA_data, aes(PC1, PC2)) +
  geom_point(size=6) +
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  geom_text_repel(aes(label=colnames(AP_dat)))+
  theme_bw()
PCA
##### REVISION


ANOVA_dat <- as.matrix(ANOVA_dat)

ANOVA_run <- tibble(Gene = character(),Day = integer(), pval = double())

list_of_genes <- tibble(Gene=integer())

for (Days in 7:7){
  for (i in 1:89){
    Genes = colnames(biomark_clean)[i]
    temp_aov <- biomark_clean%>% filter(Day == Days) %>% select (GOI = Genes, Cell_line)
    try(aov_model <- aov(GOI ~ Cell_line, data=temp_aov), silent = T)
    pvalue = summary(aov_model)[[1]][["Pr(>F)"]][1]
    if (!is.null(pvalue)){
      if (pvalue < 0.05) {
        print(Genes)
        print(Days)
        list_of_genes <- add_row(list_of_genes, Gene=i)
        print(TukeyHSD(aov_model))
      }}
    # ANOVA_run <- add_row(ANOVA_run,Gene=Genes, pval = pvalue, Day = Days)
  }
}
sink()
