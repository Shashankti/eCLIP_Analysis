#!/bin/bash

export INPUT=~/CLIP_Seq/Data/Processed
export GENOME=~/CLIP_Seq/Data/Genome
export COUNTS=~/CLIP_Seq/Data/Counts

conda activate htseq-clip

for i in $INPUT/*_dedup.bam
do
SAMPLE1=$(echo ${i} | sed "s/_dedup.bam//")
SAMPLE=$(basename "${SAMPLE1}")
F1=$(basename "${SAMPLE}_counts")
htseq-clip extract -i ${i} -e 2 -s s -g -1 --primary -o $INPUT/${F1}_sites.bed
htseq-clip count -i $INPUT/${F1}_sites.bed -a $GENOME/SLWD_w100s20.txt -o $COUNTS/${F1}.csv
done


