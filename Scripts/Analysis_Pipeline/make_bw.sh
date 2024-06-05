#!/bin/bash

for i in *.bam
do
SAMPLE=$(echo ${i} | sed "s/.bam//")
F1=$(basename "${SAMPLE}.bw")
bamCoverage -b ${i} -p 6 -o ${F1}
done
