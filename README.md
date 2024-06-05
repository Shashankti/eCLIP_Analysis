# eCLIP_Analysis

### The repo contains list of all tools which were utilised in the Imig lab for analysis of eCLIP-seq data.

## Preprocessing tools
1. Demultiplexing
2. Trimming and quality checks
3. Alignment
4. Calling Peaks and motifs

## Structure
1. Demux
- Script for demultiplexing files

2. Analysis_Pipeline
- Scripts for generating peaks from fastq file as input
- Requires STAR,samtools,bedtools,pureclip,MEME,cutadapt
- Merge_annotate.sh is used to annotate the peak files and prepare them for RCAS
- Make_bw.sh is used for generating files for visualizing in IGV
- 

3. RScipts
- Scripts used for post processing and generating plots using RCAS
