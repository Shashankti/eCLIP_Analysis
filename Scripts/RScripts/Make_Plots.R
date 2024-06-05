#MMake plots 
library(EnhancedVolcano)
library(ggplot2)
library(viridis)
library(cowplot)


#
source('https://raw.githubusercontent.com/EMBL-Hentze-group/DEWSeq_analysis_helpers/master/Volcano_plot/volcanoplot.R')

volcanoplot(resultWindows,padj_col = 'padj')

E1 <- EnhancedVolcano(as.data.frame(resultRegions),
                      lab = resultRegions$gene_name,
                      x='log2FoldChange_mean',
                      y='padj_mean',
                      title = 'Norm vs Lo',
                      pCutoff = 0.1,
                      FCcutoff = 0.5)
E1
EnhancedVolcano(resultWindows,
                lab=resultWindows$gene_name,
                x='log2FoldChange',
                y='padj',
                title = "Norm vs Input",
                pCutoff = 0.01,
                FCcutoff = 0.5)
EnhancedVolcano(as.data.frame(resultRegions),
                lab = resultRegions$gene_name,
                x='log2FoldChange_mean',
                y='padj_mean',
                pCutoff = 0.01,
                FCcutoff = 0.5)

E2 <- EnhancedVolcano(as.data.frame(resultRegions),
                      lab = resultRegions$gene_name,
                      x='log2FoldChange_mean',
                      y='padj_mean',
                      title = 'Norm vs Gal')
E2

E3 <- EnhancedVolcano(as.data.frame(resultRegions),
                      lab = resultRegions$gene_name,
                      x='log2FoldChange_mean',
                      y='padj_mean',
                      title = 'Gal vs Lo')

  coord_cartesian(xlim = c(-24,17),ylim = c(0,12.5))
  
  
ggplot(resultWindows[resultWindows$significant==TRUE & resultWindows$log2FoldChange>0.5,],
       aes(x=gene_type,fill=gene_type))+
  geom_bar()+
  theme_cowplot()+
  scale_fill_viridis(discrete = T)+
  theme(axis.text.x = element_text(angle=45,hjust=1),
        legend.position = "none")+
  xlab("") 
#make bed files for lfc>20


#percentage of significant intercations
for (i in length(unique(resultWindows$gene_type))){
  tot_counts[i] = dim(resultWindows[resultWindows$gene_type==unique(resultWindows$gene_type)[i],])[1]
}

for ( i in unique(resultWindows$gene_type)){
  tot_counts[i] = dim(resultWindows[resultWindows$gene_type==i,])[1]
  sig_counts[i] = dim(resultWindows[resultWindows$gene_type==i & resultWindows$significant==TRUE,])[1]
}
prop_counts <- sig_counts/tot_counts
prop_counts <- as.data.frame(prop_counts) %>% rownames_to_column("type")
colnames(prop_counts)[2] <- "prop"

ggplot(prop_counts,aes(x=type,fill=type,y=prop))+
  geom_bar(stat = "Identity")+
  theme_cowplot()+
  scale_fill_viridis(discrete = T)+
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        legend.position = "none")+
  xlab("")+
  ylab("Proportion of significant")
#------------------------------------------------------------------------------------------------#

# 



