#################################

E7.0_zipcode <- read.csv2("E7.0.txt", h=F, sep=",")

###Data for Guangdun Mapping
dds_map <- as.data.frame(assay(dds))

dds_map <- dds_map[,c(grep("Ost", colnames(dds_map)))]

dds_map2 <- dds_map %>% rownames_to_column("Gene") %>% filter (rownames(dds_map) %in% E7.0_zipcode$V1)%>% column_to_rownames("Gene")
write.csv(dds_map2, "Mapping7.0.csv")



#zipcode mapping
#Cornplot E7.0
library(ggplot2)

setwd("C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset")
Cordata <- read.table("Correlation Score Zipcode E7.0 from website_AVERAGE2.txt", h=T, row.names=1)


y7.0 <- c(1, 1, rep(c(2,3,4,5,6,7,8,9,10,11), each=4))
x7.0 <- c(1.5,3.5, rep(c(1,4,2,3), 10))

Cor <- Cordata[4,]

E7.0 <- cbind.data.frame(y7.0,x7.0,t(Cor))
#E7.0 <- rbind.data.frame(E7.0, c(2,5,1), c(1,5,0))
colnames(E7.0) <- c("Layer", "Region", "Cor")

E7.0_plot <- ggplot(E7.0, aes(x = Region, y= Layer, colour = Cor)) + 
  geom_point(size = 12) + 
  xlim(0.5,4.5) + 
  ylim(0,12)+
  scale_colour_gradientn(colours = c("lightgrey", "green", "yellow","red"), 
                         values = c(0,0.1,0.5,1))+
  labs(title = paste(rownames(Cor)), subtitle= paste("           Max Cor = ", max(E7.0$Cor)))+
  labs(x= "E7.0", ylab="")+
  xlab ("A     L     R     P" )+
  scale_y_continuous("", breaks = y7.0)+
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
