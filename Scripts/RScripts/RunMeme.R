#running meme motif finding via R

library(memes)

options(meme_bin = "~/tools/meme/bin")

runDreme(input = "~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/Norm_joined.fasta",
         dna=TRUE,mink=3,maxk=15,outdir = "~/CLIP_Seq/Data_Strict_Align/Meme/",
         control = "~/CLIP_Seq/Data/Genome/grch38.fa")

## scrapping TCGA data

library(TCGAretriever)
library(TCGA2STAT)


# Obtain a list of cancer studies from cBio
all_studies <- get_cancer_studies()

# Find published TCGA datasets
keep <- grepl("tcga_pub$", all_studies[,1])
tcga_studies <- all_studies[keep, ]

# Show results
head(tcga_studies[, 1:2])

rnaseq.aml <- getTCGA(disease = "LAML",data.type = "RNASeq2",type = "count")

#3copyt function from githib

#annotate results from FIMO
library(dplyr)
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)

found_peaks <- data.table::fread("Data_Strict_Align/fimo_out/fimo.tsv")
stringr::str_split_fixed(found_peaks$sequence_name, ':',2)
found_peaks <- found_peaks[,3]

found_peaks <- found_peaks %>% separate(sequence_name, c('chr','start','end'))
#found_peaks$chr <- as.numeric(found_peaks$chr)
data.table::setorder(found_peaks, chr, start, end)

write.table(found_peaks, file = "Data_Strict_Align/fimo_out/peaks_loc.bed", col.names = F,
            quote = F, row.names = F, sep = "\t")