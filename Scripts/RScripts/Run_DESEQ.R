#runnning DESeq2

library(DESeq2)
library(tidyverse)
count_matrix$unique_id <- NULL
head(count_matrix)

dds_2f <- DESeqDataSetFromMatrix (countData = round(count_matrix),
                                  colData = col_data,
                                  design = ~ type)


#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)



#releveal condition
dds_2f$condition <- relevel(dds_2f$condition, ref = "Norm")

dds_2f <- DESeq(dds_2f)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f 

res_norm_vs_gal <- results(dds_2f,contrast = c("type","Norm","Gal"))
res_norm_vs_lo <- results(dds_2f,contrast = c("type","Norm","Lo"))
  

#Info regarding metadata
metadata(res_norm_vs_gal)$filterThreshold

#Plot for independent filtering
plot(metadata(res_norm_vs_gal)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_norm_vs_gal)$lo.fit, col="red")
abline(v=metadata(res_norm_vs_gal)$filterTheta)

# Check the filtering criteria

plot(res_norm_vs_gal$baseMean+1, -log10(res_norm_vs_gal$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
use <- res_norm_vs_gal$baseMean > metadata(res_norm_vs_gal)$filterThreshold
h1 <- hist(res_norm_vs_gal$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_norm_vs_gal$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))



##h

hgnc_list <- c("RBM39","ZFP36L2","DGCR8","FXR2","HNRNPK","HNRNPU","AGGF1","AKAP1","AQR","BUD13","CSTF2T","DDX3X","DDX52","DDX55","DDX6","DHX30","DROSHA","EFTUD2","EXOSC5","FAM120A","FASTKD2","FMR1","FTO","FUS","FXR1","GAPDH","GRWD1","GTF2F1","HLTF","HNRNPA1","HNRNPC","HNRNPL","HNRNPM","HNRNPUL1","IGF2BP1","ILF3","KHSRP","LARP4","LARP7","LIN28B","LSM11","MATR3","NCBP2","NOLC1","NSUN2","PCBP1","PRPF8","PTBP1","QKI","RBFOX2","RBM15","RBM22","RPS3","SAFB","SDAD1","SF3B4","SLTM","SMNDC1","SND1","SRSF1","SRSF7","SRSF9SSB","SUPV3L1","TAF15","TARDBP","TBRG4","TIA1","TRA2A","TROVE2","U2AF1","U2AF2","UCHL5","UPF1","UTP18","WDR43","XRCC6","XRN2","YBX3","ZC3H11A","ZNF800","AARS","AATF","ABCF1","ADAT1","AKAP8L","APEX1","APOBEC3C","BCCIP","BCLAF1","CDC40","CPEB4","CPSF6","CSTF2","DDX21","DDX24","DDX42","DDX43","DDX47","DDX51","DDX59","DKC1","EEF2","EIF3D","EIF3G","EIF3H","EIF4E","EIF4G2","ELAC2","ELAVL1","EWSR1","EXOSC10","FKBP4","FUBP3","G3BP1","GARS","GEMIN5","GNL3","GPKOW","GRSF1","IGF2BP2","IGF2BP3","KHDRBS1","MBNL1","METAP2","METTL1","MORC2","MTPAP","NIP7","NIPBL","NKRF","NOL12","NONO","NPM1","PABPC4","PABPN1","PCBP2","PHF6","POLR2G","PPIG","PPIL4","PRPF4","PUM1","PUM2","PUS1","RBM5","RNF187","RPS10","RPS11","RPS6","RYBP","SAFB2","SBDS","SERBP1","SF3A3","SF3B1","SFPQ","SLBP","STAU2","SUB1","SUGP2","TIAL1","UTP3","WDR3","WRN","XPO5","YWHAG","ZC3H8","ZNF622","ZRANB2")
human <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol", values=hgnc_list, mart=human)
gene_coords$size=gene_coords$end_position - gene_coords$start_position
gene_coords