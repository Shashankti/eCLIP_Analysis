#!/bin/bash

for i in ~/CLIP_Seq/Data_Strict_Align/Peaks/Piranha/*dedup_bin.bed
do
SAMPLE=$(echo ${i} | sed "s/_dedup_bin.bed//")
F1=$(basename "${SAMPLE}_annotated.bed")
bedtools intersect -a ${i} -b ~/CLIP_Seq/Data/Genome/gene_loc_file.bed -wa -wb | cut -f1,2,3,10,11 | sed "s/;//g" | tr ' ' '\t'  > ~/CLIP_Seq/Data_Strict_Align/Peaks/Piranha/${F1}
done
