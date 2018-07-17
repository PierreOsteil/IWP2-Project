###Volcano plot
library(tidyverse)
#need a result table from DESeq2 
Comp <- resCha_FSIvsAFAF

gene_list <- cbind.data.frame(Comp$log2FoldChange, Comp$pvalue)

colnames(gene_list) <- c("logFC", "adj.P.Val")
rownames(gene_list) <- rownames(Comp)

##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$adj.P.Val < 0.05)
gene_list[is.na(gene_list)] <- 0

gene_list$threshold= as.factor(gene_list$threshold)
gene_list <- gene_list %>% filter(gene_list$threshold == TRUE)


##Construct the plot object
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  #xlim(c(-3, 3)) + 
  #ylim(c(0, 25)) +
  xlab("log2 fold change") + ylab("-log10 p-value")+
  #geom_text(aes(label = rownames(gene_list), size = 3))+
  theme_bw()

g
