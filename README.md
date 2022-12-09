# eCLIP_Analysis

### The repo contains list of all tools which were utilised in the Imig lab for analysis of eCLIP-seq data.

## Preprocessing tools
1. Demultiplexing
2. Trimming and quality checks
3. Alignment
4. Peaks

## Structure
1. Demux
- Script for demultiplexing files

2. Analysis_Pipeline
- Scripts for generating peaks from fastq file as input
- Requires STAR,samtools,bedtools,piranha,MEME,cutadapt

3. DEWSeq_Scripts
- Script to prepare input file for DEWSeq
- Run DEWSeq and generate plots
- Requires DEWSeq and htseq-clip

