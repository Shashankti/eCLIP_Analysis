library(RCAS)

Gal_R1_path <- "Data/Peaks/Piranha/Bed/Gal_R1_bin.bed"
Gal_R2_path <- "Data/Peaks/Piranha/Bed/Gal_R2_bin.bed"
Norm_R1_path <- "Data/Peaks/Piranha/Bed/Norm_R1_bin.bed"
Norm_R2_path <- "Data/Peaks/Piranha/Bed/Norm_R2_bin.bed"
Lo_R1_path <- "Data/Peaks/Piranha/Bed/Lo_R1_bin.bed"
Lo_R2_path <- "Data/Peaks/Piranha/Bed/Lo_R2_bin.bed"

projData <- data.frame('sampleName' = c ('Gal_R1','Gal_R2','Norm_R1','Norm_R2','Lo_R1','Lo_R2'),
                       'bedFilePath' = c(Gal_R1_path,Gal_R2_path,
                                         Norm_R1_path,Norm_R2_path,
                                         Lo_R1_path,Lo_R2_path),
                       stringsAsFactors = F)
                       
projDataFile <- file.path(getwd(), 'myProjDataFile.tsv') 
write.table(projData, projDataFile, sep = '\t', quote =FALSE, row.names = FALSE) 

gtfFilePath <- "Data/Genome/Homo_sapiens.GRCh38.108.gtf"

databasePath <- file.path(getwd(), 'myProject.sqlite')
invisible(createDB(dbPath = databasePath, projDataFile = projDataFile, 
                   gtfFilePath = gtfFilePath, genomeVersion = 'hg38'))

createDB(dbPath = databasePath, projDataFile = projDataFile, gtfFilePath = gtfFilePath, 
         genomeVersion = 'hg38', update = TRUE,motifAnalysis = F,coverageProfiles = F)


mydb <- RSQLite::dbConnect(RSQLite::SQLite(), dbPath)
projData <- validateProjDataFile(projDataFile, mydb)
gtfData <- insertTableGtfData(conn = mydb, name = "gtfData", 
                              gtfFilePath = gtfFilePath)
message("Parsing transcript features")
txdbFeatures <- getTxdbFeaturesFromGRanges(gtfData)
bedData <- insertTableBedData(conn = mydb, projData = projData)
if (annotationSummary == TRUE) {
  insertTableAnnotationSummaries(conn = mydb, bedData = bedData, 
                                 txdbFeatures = txdbFeatures, nodeN = nodeN)
  geneRanges <- unlist(range(split(gtfData, gtfData$gene_name)))
  insertTableOverlapMatrix(conn = mydb, name = "geneOverlaps", 
                           bedData = bedData, targetRegions = geneRanges, targetRegionNames = unique(names(geneRanges)), 
                           nodeN = 1)