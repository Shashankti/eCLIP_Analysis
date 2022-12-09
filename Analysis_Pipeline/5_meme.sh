#!/bin/bash

module load samtools 

export MEME=Data/Meme
export MEME1=~/tools/meme/bin
export PEAKS=Data/Peaks/Piranha
export GENOME=Data/Genome/
export SORTED=Data/Processed/Sorted

for i in Data/Processed/Sorted/*.bam
do
SAMPLE=$(echo ${i} | sed "s/_sorted.bam//")
F1=$(basename "${SAMPLE}_sorted.fasta")
F2=$(basename "${SAMPLE}_dedup.bam.bed_bin100.txt")
F3=$(basename "${SAMPLE}_peaks.bin100")
F4=$(basename "${SAMPLE}_meme")

samtools bam2fq ${i} | sed -n '1~4s/^@/>/p;2~4p' > $SORTED/${F1}
 
awk 'NF{NF-=1};1' $PEAK/${F2} | tr ' ' '\t' < - > $PEAK/${F3}_.bed

bedtools getfasta -fi  $GENOME/hg38.fa -bed $PEAK/${F3} -fo $PEAK/${F3}_.fasta

$MEME1/meme $F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_align_long
$MEME1/meme $F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_align_short -minw 15 -maxw 15

$MEME1/meme $F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_peak_long
$MEME1/meme $F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_peak_short -minw 5 -maxw 15
done


