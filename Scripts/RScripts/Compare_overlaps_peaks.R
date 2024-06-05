#comparing gene names for annotated peaks
library(biomaRt)
library(data.table)
library(dplyr)
library(purrr)

#read in the annotated peaks files
pir_peaks <- fread("Data/Peaks/Piranha/Bed/Lo_R1_annotated.bed")
pir_peaks <- pir_peaks[,c(1,2,3,5)]
colnames(pir_peaks) <- c("chr","start","end","ensembl_id")
pure_peaks <- fread("Data/Peaks/PureClip/Merged_Lo_R1_annotated.bed")
pure_peaks <- pure_peaks[,c(1,2,3,5)]
colnames(pure_peaks) <- c("chr","start","end","ensembl_id")

pir_peaks2 <- fread("Data/Peaks/Piranha/Bed/Lo_R2_annotated.bed")
pir_peaks2 <- pir_peaks2[,c(1,2,3,5)]
colnames(pir_peaks2) <- c("chr","start","end","ensembl_id")
pure_peaks2 <- fread("Data/Peaks/PureClip/Merged_Lo_R2_annotated.bed")
pure_peaks2 <- pure_peaks2[,c(1,2,3,5)]
colnames(pure_peaks2) <- c("chr","start","end","ensembl_id")


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = pir_peaks$ensembl_id, mart= mart)
pir_peaks$gene_name <- gene_IDs[match(pir_peaks$ensembl_id,gene_IDs[,1]),2]

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = pure_peaks$ensembl_id, mart= mart)
pure_peaks$gene_name <- gene_IDs[match(pure_peaks$ensembl_id,gene_IDs[,1]),2]

#list of genes present in both replicates
pir_gene_list <- pir_peaks[which(pir_peaks$ensembl_id %in% pir_peaks2$ensembl_id),]
pure_gene_list <- pure_peaks[which(pure_peaks$ensembl_id %in% pure_peaks2$ensembl_id),]

length(which(unique(pure_gene_list$ensembl_id) %in% unique(pir_gene_list$ensembl_id)))
length(which(unique(gene_list$gene_id) %in% unique(pir_gene_list$ensembl_id)))
length(which(unique(gene_list$gene_id) %in% unique(pure_gene_list$ensembl_id)))

length(which(unique(pure_gene_list$ensembl_id) %in% unique(gene_list$gene_id)))

#write output
overlap1 <- na.omit(pure_gene_list[match(pir_gene_list$ensembl_id,pure_gene_list$ensembl_id),])
overlap1 <- as.data.frame(unique(overlap1$ensembl_id))
colnames(overlap1) <- c("ensembl_id")
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = overlap1$ensembl_id, mart= mart)
overlap1$gene_name <- gene_IDs[match(overlap1$ensembl_id,gene_IDs[,1]),2]
write.csv(overlap1,file = "Data/Peaks/PureClip_vs_Piranha_Gal.csv",row.names = F,quote = F)

mart <- useEnsembl("ensembl","hsapiens_gene_ensembl",mirror = "useast")


overlap2 <- na.omit(pir_gene_list[match(gene_list$gene_id,pir_gene_list$ensembl_id),])
overlap2 <- as.data.frame(unique(overlap2$ensembl_id))
colnames(overlap2) <- c("ensembl_id")
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = overlap2$ensembl_id, mart= mart)
overlap2$gene_name <- gene_IDs[match(overlap2$ensembl_id,gene_IDs[,1]),2]
write.csv(overlap2,file = "Data/Peaks/Piranha_vs_DEWSeq_Lo.csv",row.names = F,quote = F)

overlap3 <- na.omit(pure_gene_list[match(gene_list$gene_id,pure_gene_list$ensembl_id),])
overlap3 <- as.data.frame(unique(overlap3$ensembl_id))
colnames(overlap3) <- c("ensembl_id")
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = overlap3$ensembl_id, mart= mart)
overlap3$gene_name <- gene_IDs[match(overlap3$ensembl_id,gene_IDs[,1]),2]
write.csv(overlap3,file = "Data/Peaks/PureClip_vs_DEWSeq_Lo.csv",row.names = F,quote = F)

merged <- list(gene_list,pure_gene_list$ensembl_id,pir_gene_list$ensembl_id)

common_genes <- as.data.frame(Reduce(intersect,list(gene_list$gene_id,pir_gene_list$ensembl_id,pure_gene_list$ensembl_id)))
colnames(common_genes) <- "ensembl_id"
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = common_genes$ensembl_id, mart= mart)
common_genes$gene_name <- gene_IDs[match(common_genes$ensembl_id,gene_IDs[,1]),2]
write.csv(common_genes,file = "Data/Peaks/All_common_Norm.csv",quote = F,row.names = F)

library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
