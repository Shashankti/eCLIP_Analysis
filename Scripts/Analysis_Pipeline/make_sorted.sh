#!/bin/bash


for i in /home/shashank.tiwari/CLIP_Seq/Data_wo_rep/Processed/*dedup.bam
do
SAMPLE=$(echo ${i} | sed "s/_dedup.bam//")
F0=$(basename "${SAMPLE}_sorted.bam")
F1=$(basename "${F0}.bw")
echo ${i}
echo ${F0}
samtools sort -@ 8 ${i} > ~/CLIP_Seq/Data_wo_rep/Processsed/${F0}
samtools index ~/CLIP_Seq/Data_wo_rep/Processed/${F0}
bamCoverage -b ~/CLIP_Seq/Data_wo_rep/Processed/${F0} -p 6 -o ~/CLIP_Seq/Data_wo_rep/Processed/${F1}
done

