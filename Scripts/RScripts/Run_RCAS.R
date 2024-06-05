#work with CLIP-Seq peaks

#load the libraries
library(RCAS)
library(data.table)
library(ggplot2)
library(ggtext)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(plyranges)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#load the annotation gff file
gff <- importGtf(filePath="Data/Genome/Homo_sapiens.GRCh38.108.gtf")
#make txdb gene object 
txdbFeatures <- getTxdbFeaturesFromGRanges(gff)
# Bed_file <- fread("Data/Peaks/Piranha/Merged_Gal_R1_dedup.bam.bed_bin100.txt")
# options(scipen = 999)
# write.table(Bed_file,file = "Data/Peaks/Piranha/bed_file.bed",sep = "\t", col.names = F,row.names = F)
#select only the peaks with pvalue non zero pvalues
#read in the bed files
#queryRegions <- importBed("Data/Peaks/Piranha/Bed/Norm_R1_bin100.bed")
queryRegions <- importBed("Data_Strict_Align/Peaks/PureClip/Union/Lo_joined.bed")
#queryRegions <- importBed("Data/Peaks/Piranha/Bed/100_bin/Merged_Norm_R2_peaks.bin100.bed")
#union R1 and R2
queryRegions <- importBed("Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_all_three_peaks.bed")
#intersected r1 and r2
queryRegions <- importBed("Data_Strict_Align/Peaks/PureClip/Intersection/RCAS_any_two_intersection.bed")




#annotation of the input bed with gff gfile
overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))

#plot the distribution of query regions across gene types
biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
df <- overlaps[,length(unique(queryIndex)), by = biotype_col]
colnames(df) <- c("feature", "count")
df$percent <- round(df$count / length(queryRegions) * 100, 1)
df <- df[order(count, decreasing = TRUE)]

plot1 <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_richtext(aes(y = percent + 0.5), label = df$count,angle=0,size=2) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))

tiff("Figures/RCAS_intersection/trans-feat-1_peaks_2.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()
rm(plot1,df)
gc()









summary <- summarizeQueryRegions(queryRegions = queryRegions, 
                                 txdbFeatures = txdbFeatures)

df <- data.frame(summary)
df$percent <- round((df$count / length(queryRegions)), 3) * 100
df$feature <- rownames(df)
plot2 <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 3), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))
tiff("Figures/RCAS_intersection/trans-feat-2_peaks_2.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot2)
dev.off()
rm(plot2,df)
gc()
dt <- getTargetedGenesTable(queryRegions = queryRegions, 
                            txdbFeatures = txdbFeatures)
dt <- dt[order(transcripts, decreasing = TRUE)]

knitr::kable(dt[1:10,])


cvgF <- getFeatureBoundaryCoverage(queryRegions = queryRegions, 
                                   featureCoords = txdbFeatures$transcripts, 
                                   flankSize = 1000, 
                                   boundaryType = 'fiveprime', 
                                   sampleN = 10000)
cvgT <- getFeatureBoundaryCoverage(queryRegions = queryRegions, 
                                   featureCoords = txdbFeatures$transcripts, 
                                   flankSize = 1000, 
                                   boundaryType = 'threeprime', 
                                   sampleN = 10000)

cvgF$boundary <- 'fiveprime'
cvgT$boundary <- 'threeprime'

df <- rbind(cvgF, cvgT)

plot3 <-ggplot2::ggplot(df, aes(x = bases, y = meanCoverage)) + 
  geom_ribbon(fill = 'lightgreen', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + 
  facet_grid( ~ boundary) + theme_bw(base_size = 14) 

tiff("Figures/RCAS_intersection/Two_5prime-3prime_2.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot3)
dev.off()
rm(plot3,df,cvgF,cvgT)
gc()

# coverage profile of query regions for all transcript features

cvgList <- calculateCoverageProfileList(queryRegions = queryRegions, 
                                        targetRegionsList = txdbFeatures,
                                        sampleN=10000)

#custom plot
# p <- plotly::plot_ly(data = cvgList, type = 'scatter', mode = 'lines')
# for (f in unique(cvgList$feature)){
#   data <- cvgList[cvgList$feature == f,]
#   p <- plotly::add_trace(p = p, data = data, x = ~bins, y = ~meanCoverage, 
#                          legendgroup = f, showlegend = FALSE, opacity = 1, color = f)
#   p <- plotly::add_ribbons(p = p, data = data, x = ~bins, 
#                            ymin = data$meanCoverage - data$standardError*1.96,
#                            ymax = data$meanCoverage + data$standardError*1.96, 
#                            legendgroup = f, 
#                            name = f, color = f
#   )
# }
# plotly::layout(p, font = list(size = 14))
# 
# if(params$printProcessedTables == TRUE) {
#   write.table(x = cvgList, file=paste0('Figure',figureCount,'.coverageprofilelist.data.tsv'), quote = FALSE, sep = '\t', row.names = TRUE)
# } 
# figureCount <- figureCount + 1
# 


plot4 <- ggplot2::ggplot(cvgList, aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(fill = 'lightgreen', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + theme_bw(base_size = 14) +
  facet_wrap( ~ feature, ncol = 3) 
tiff("Figures/RCAS_intersection/Two_coverage_2.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot4)
dev.off()
rm(plot4,cvgList)
gc()

#motif analysis
motifResults <- runMotifDiscovery(queryRegions = queryRegions,
                                  resizeN = 15, sampleN = 10000,
                                  genomeVersion = 'hg38', motifWidth = 10,
                                  motifN = 2, nCores = 6)

plot5 <- ggseqlogo::ggseqlogo(motifResults$matches_query)
tiff("Figures/RCAS_Gal_Overlap/motif_Norm_2.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot5)
dev.off()


summary <- getMotifSummaryTable(motifResults)
knitr::kable(summary)
write.table(summary,file = "Results/RCAS_Norm_Overlap/motif_Norm.csv",row.names = T,col.names = T)

targetedGenes <- unique(overlaps$gene_id)

res <- RCAS::findEnrichedFunctions(targetGenes = targetedGenes, species = 'hsapiens')
res <- res[order(res$p_value),]
resGO <- res[grep('GO:BP', res$source),]
knitr::kable(subset(resGO[1:10,], select = c('p_value', 'term_name', 'source')))
fwrite(resGO,file = "Results/RCAS_Norm_Combined/resGo.csv")



#----------------------------------------------------------------------------------------------------#

runReport(queryFilePath = "Data/Peaks/Piranha/Bed/Norm_R1_bin100.bed",
          gffFilePath = "Data/Genome/Homo_sapiens.GRCh38.108.gtf",
          motifAnalysis = T,
          goAnalysis = T,
          printProcessedTables = T,
          genomeVersion = "GRCh38")

#extract the fiveUTR regions
#hits <- findOverlaps(query = queryRegions, subject = txdbFeatures$promoters, ignore.strand=T)
hits <- (subsetByOverlaps(queryRegions, txdbFeatures$fiveUTRs))

genes <- as.data.frame(txdbFeatures$fiveUTRs)
as_granges(genes[,c(1,2,3,4,5,7)])
Overlap_region <- join_overlap_left(hits, as_granges(genes))
#fiveUTRs
regions.df <- data.frame(unique(as.data.frame(Overlap_region)[,8]))
regions.df$type <- rep("fiveUTR", dim(regions.df)[1])
colnames(regions.df)[1] <- "gene_name"
#promoters
regions2.df <- data.frame(unique(as.data.frame(Overlap_region)[,9]))
regions2.df$type <- rep("Promoters", dim(regions2.df)[1]) 
colnames(regions2.df)[1] <- "gene_name"
#3UTRs
regions3.df <- data.frame(unique(as.data.frame(Overlap_region)[,8]))
regions3.df$type <- rep("threeUTR", dim(regions3.df)[1]) 
colnames(regions3.df)[1] <- "gene_name"

regions1.df <- rbind(regions.df,regions2.df,regions3.df)
write.table(regions1.df,file = "Data_Strict_Align/Peaks/Norm_peakRegions.tsv", row.names = F, quote = F, sep = "\t")


#---------------------------------------------------------------------#
#custom runReport function

runReport2 <- function(queryFilePath = 'testdata',
                      gffFilePath = 'testdata',
                      annotationSummary = TRUE,
                      goAnalysis = TRUE,
                      motifAnalysis = TRUE,
                      genomeVersion = 'hg19',
                      outDir = getwd(),
                      printProcessedTables = FALSE,
                      sampleN = 0,
                      quiet = FALSE,
                      selfContained = TRUE) {
  
  db <- checkSeqDb(genomeVersion)
  # get species name 
  # this is needed for gprofiler functional enrichment 
  fields <- unlist(strsplit(db@metadata$organism, ' '))
  species <- tolower(paste0(unlist(strsplit(fields[1], ''))[1], 
                            fields[2]))
  
  if(queryFilePath != 'testdata') {
    queryFilePath <- normalizePath(queryFilePath)
  }
  
  if(gffFilePath != 'testdata') {
    gffFilePath <- normalizePath(gffFilePath)
  }
  
  reportFile <- system.file("reporting_scripts", "report.Rmd", package='RCAS')
  headerFile <- system.file("reporting_scripts", "header.html", package='RCAS')
  footerFile <- system.file("reporting_scripts", "footer.html", package='RCAS')
  
  outFile <- paste0(basename(queryFilePath), '.RCAS.report.html')
  
  rmarkdown::render(
    input = reportFile, 
    output_dir = outDir,
    intermediates_dir = outDir,
    output_file = outFile,
    output_format = rmarkdown::html_document(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'simplex',
      number_sections = TRUE,
      includes = rmarkdown::includes(in_header = headerFile, 
                                     after_body = footerFile), 
      self_contained = selfContained
    ),
    params = list(query = queryFilePath,
                  gff = gffFilePath,
                  annotationSummary = annotationSummary,
                  goAnalysis = goAnalysis,
                  motifAnalysis = motifAnalysis,
                  genomeVersion = genomeVersion,
                  species = species,
                  printProcessedTables = printProcessedTables,
                  sampleN = sampleN,
                  workdir = outDir),
    quiet = quiet
  )
}

