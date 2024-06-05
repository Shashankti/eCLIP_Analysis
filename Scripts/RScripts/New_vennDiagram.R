#Make venndiagrams for eclipseq peaks

library(ChIPpeakAnno)
library(data.table)
library(GenomicRanges)
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

#read in the bed files

Norm_overlap <- fread("Data_Strict_Align/Peaks/PureClip/Norm_R1_fin.bed")
Norm_overlap <- makeGRangesFromDataFrame(Norm_overlap, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Norm_overlap <- Norm_overlap[!duplicated(Norm_overlap)]


Gal_overlap <- fread("Data_Strict_Align/Peaks/PureClip/Norm_R2_fin.bed")
Gal_overlap <- makeGRangesFromDataFrame(Gal_overlap, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Gal_overlap <- Gal_overlap[!duplicated(Gal_overlap)]

Lo_overlap <- fread("Data_Strict_Align/Peaks/PureClip/Merged_Lo_2.bed")
Lo_overlap <- makeGRangesFromDataFrame(Lo_overlap, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Lo_overlap <- Lo_overlap[!duplicated(Lo_overlap)]

#make the overlap object
ol <- findOverlapsOfPeaks(Gal_overlap,Norm_overlap)
averagePeakWidth <- mean(width(unlist(GRangesList(ol$peaklist))))

#read in combined files

Norm_combined <- fread("Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Norm_peaks.bed")
Norm_combined <- makeGRangesFromDataFrame(Norm_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Norm_combined <- Norm_combined[!duplicated(Norm_combined)]


Gal_combined <- fread("Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Gal_peaks.bed")
Gal_combined <- makeGRangesFromDataFrame(Gal_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Gal_combined <- Gal_combined[!duplicated(Gal_combined)]

Lo_combined <- fread("Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Lo_peaks.bed")
Lo_combined <- makeGRangesFromDataFrame(Lo_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Lo_combined <- Lo_combined[!duplicated(Lo_combined)]

ol2 <- findOverlapsOfPeaks(Gal_combined,Lo_combined,Norm_combined)
png("Figures/Venn_diagram_overlap_of_joined_peaks.png",width = 1080,height = 1080,res = 300)
combined_venn <- makeVennDiagram(ol2,
                                 NameOfPeaks = c("Galacatose","Low","Normal"),
                                 # Circles
                                 lwd = 2,
                                 lty = 'blank',
                                 fill = myCol,
                                 
                                 # Numbers
                                 cex = 0.8,
                                 fontfamily = "sans",
                                 
                                 # Set names
                                 cat.cex = 1,
                                 cat.fontface = "bold",
                                 cat.default.pos = "outer",
                                 cat.pos = c(-7, 7, 180),
                                 cat.dist = c(0.075, 0.075, 0.085),
                                 cat.fontfamily = "sans",
                                 rotation = 1,
                                 scaled=FALSE)
dev.off()






#read in intersected files

Norm_combined <- fread("Data_Strict_Align/Peaks/PureClip/Intersection/Merged_Norm_R1_int_R2.bed")
Norm_combined <- makeGRangesFromDataFrame(Norm_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Norm_combined <- Norm_combined[!duplicated(Norm_combined)]


Gal_combined <- fread("Data_Strict_Align/Peaks/PureClip/Intersection/Merged_Gal_R1_int_R2.bed")
Gal_combined <- makeGRangesFromDataFrame(Gal_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Gal_combined <- Gal_combined[!duplicated(Gal_combined)]

Lo_combined <- fread("Data_Strict_Align/Peaks/PureClip/Intersection/Merged_Lo_R1_int_R2.bed")
Lo_combined <- makeGRangesFromDataFrame(Lo_combined, start.field = "V2", end.field = "V3", seqnames.field = "V1")
Lo_combined <- Lo_combined[!duplicated(Lo_combined)]

ol2 <- findOverlapsOfPeaks(Gal_combined,Lo_combined,Norm_combined)
tiff("Figures/Venn_diagram_overlap_of_intersected_peaks.tiff",width = 1080,height = 1080,res = 300)
combined_venn <- makeVennDiagram(ol2,
                                 NameOfPeaks = c("Galacatose","Low","Normal"),
                                 # Circles
                                 lwd = 2,
                                 lty = 'blank',
                                 fill = myCol,
                                 
                                 # Numbers
                                 cex = 0.8,
                                 fontfamily = "sans",
                                 
                                 # Set names
                                 cat.cex = 1,
                                 cat.fontface = "bold",
                                 cat.default.pos = "outer",
                                 cat.pos = c(-7, 7, 180),
                                 cat.dist = c(0.075, 0.075, 0.085),
                                 cat.fontfamily = "sans",
                                 rotation = 1,
                                 scaled=FALSE)
dev.off()








#visualise and compare the gene names


Norm_annotated <- fread("Data_Strict_Align/Peaks/PureClip/Merged_Norm_annotated.bed")
Gal_annotated <- fread("Data_Strict_Align/Peaks/PureClip/Merged_Gal_annotated.bed")
Lo_annotated <- fread("Data_Strict_Align/Peaks/PureClip/Merged_Lo_annotated.bed")

Combined_norm <- fread("Data_Strict_Align/Peaks/PureClip/Union/Norm_joined_annotated.bed")
Combined_Gal <- fread("Data_Strict_Align/Peaks/PureClip/Union/Gal_joined_annotated.bed")
Combined_Lo <- fread("Data_Strict_Align/Peaks/PureClip/Union/Lo_joined_annotated.bed")


ovNorm_genes <- unique(Norm_annotated$V5)
ovGal_genes <- unique(Gal_annotated$V5)
ovLo_genes <- unique(Lo_annotated$V5)

#make venndiagram
venn.diagram(x = list(ol2$peaklist$Gal_combined,ol2$peaklist$Lo_combined,ol2$peaklist$Norm_combined),
             category.names = c("Gal","Lo","Norm"),
             filename = "ol2.png",
             output=T,
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
             rotation = 1)

#combined genes

CBNorm_genes <- unique(Combined_norm$V5)
CBGal_genes <- unique(Combined_Gal$V5)
CBLo_genes <- unique(Combined_Lo$V5)

#make venndiagram
venn.diagram(x = list(CBGal_genes,CBLo_genes,CBNorm_genes),
             category.names = c("Gal","Lo","Norm"),
             filename = "Combined genes_all.png",
             output=T,
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
             rotation = 1)