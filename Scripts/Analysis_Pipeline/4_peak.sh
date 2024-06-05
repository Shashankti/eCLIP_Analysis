#!/bin/bash -l 

module load samtools
module load bedtools/2.26.0

export INPUT=/u/samas/GAPDH/Lo-R2/INPUT

echo 'preparing reads for peak calling..'

samtools sort $INPUT/GAPDH_Lo-R2.aligned.merged.sorted.INPUT.bam > $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam && \
samtools index -b $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam.bai && \
bamToBed -i $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam > $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bed && \

echo 'reads are ready for peak calling'

find Sorted/*.bam | parallel 'pureclip -i {} -bai {}.bai -g ../Genome/Homo_sapiens.GRCh38.dna.toplevel.fa -o {}.peaks.bed -or {}_bins.bed -nt 16 -vv'
echo 'peaks are there!'

