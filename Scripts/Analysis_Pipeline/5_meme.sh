#!/bin/bash

export MEME=Data/Meme
export MEME1=~/tools/meme/bin
export PEAKS=Data/Peaks/Pureclip
export GENOME=Data/Genome/
export SORTED=Data/Processed/Sorted

for i in Data/Processed/Sorted/*.bam
do
SAMPLE=$(echo ${i} | sed "s/_sorted.bam//")
F1=$(basename "${SAMPLE}_sorted.fasta")
F2=$(basename "${SAMPLE}_dedup.bam.bed_bin100.txt")
F3=$(basename "${SAMPLE}_peaks.bin100")
F4=$(basename "${SAMPLE}_meme")

#for the entire bam file 
## can be skipped
samtools bam2fq ${i} | sed -n '1~4s/^@/>/p;2~4p' > $SORTED/${F1}
echo 'samtools done' 
awk 'NF{NF-=1};1' $PEAKS/${F2}| tr ' ' '\t' > $PEAKS/${F3}_.bed
echo 'awk done'
#running meme for the output peaks file
#prepare the input file
awk  '{print "chr"$0}' $PEAKS/${F3}_.bed > $PEAKS/${F3}.bed
bedtools getfasta -fi  $GENOME/hg38.fa -bed $PEAKS/${F3}.bed -fo $PEAKS/${F3}_.fasta
echo 'bedtools fasta done'
#$MEME1/meme $SORTED/$F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_align_long
$MEME1/meme $SORTED/$F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_align_short -minw 5 -maxw 20
echo 'running meme 1 done'
#$MEME1/meme $SORTED/$F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_peak_long
$MEME1/meme $SORTED/$F1 -dna -nmotifs 10 -maxsize 0 -o $MEME/${F4}_peak_short -minw 5 -maxw 20
done


