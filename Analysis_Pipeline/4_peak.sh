#!/bin/bash -l 

module load samtools
module load bedtools/2.26.0

export INPUT=/u/samas/GAPDH/Lo-R2/INPUT

echo 'preparing reads for peak calling..'

samtools sort $INPUT/GAPDH_Lo-R2.aligned.merged.sorted.INPUT.bam > $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam && \
samtools index -b $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam.bai && \
bamToBed -i $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bam > $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bed && \

echo 'reads are ready for peak calling'

module load gsl/1.16

export CPPFLAGS=-I/mpcdf/soft/SLE_12/packages/haswell/gsl/gcc_9/1.16/include/ 
export LDFLAGS=-L/mpcdf/soft/SLE_12/packages/haswell/gsl/gcc_9/1.16/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mpcdf/soft/SLE_12/packages/haswell/gsl/gcc_9/1.16/lib/
export LD_LIBRARY_PATH

export PIRANHA=/u/samas/piranha/piranha-1.2.1/bin
export OUTPUT_PEAK=/u/samas/GAPDH/Lo-R2/peak/

echo 'starts peak calling with Piranha..'

$PIRANHA/Piranha -s -b 100 $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bed -o $OUTPUT_PEAK/GAPDH_Lo-R2_R2_peaks.bin100.txt && \
$PIRANHA/Piranha -s -b 50 $INPUT/GAPDH_Lo-R2.aligned.merged.sorted2.INPUT.bed -o $OUTPUT_PEAK/GAPDH_Lo-R2_R2_peaks.bin50.txt && \

echo 'peaks are there!'

