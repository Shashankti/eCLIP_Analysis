#!/bin/bash -l 
 
module load gcc

export FASTQ=/scratch1/users/shashank.tiwari/Final_fastq_files
export UMI=/scratch1/users/shashank.tiwari/Final_fastq_files/Umi
export ADAPT1=/scratch1/users/shashank.tiwari/Final_fastq_files/Adapt1
export ADAPT2=/scratch1/users/shashank.tiwari/Final_fastq_files/Adapt2

echo 'adaptor removal from L001, first round' 

for i in $UMI/*_R1.umi.fastq*;
do
SAMPLE=$(echo ${i} | sed "s/_R1.umi.fastq//")
echo ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq
F1=$(basename "${SAMPLE%umi.fastq}_R1.adapt1.fastq") F2=$(basename "${SAMPLE%umi.fastq}_R2.adapt1.fastq")
  cutadapt -f fastq \
  --match-read-wildcards \
  --times 1 \
  -e 0.1 \
  -O 1 \
  --quality-cutoff 6 \
  -m 15 \
  -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -g CTTCCGATCTACAAGTT \
  -g CTTCCGATCTTGGTCCT \
  -A AACTTGTAGATCGGA \
  -A AGGACCAAGATCGGA \
  -A ACTTGTAGATCGGAA \
  -A GGACCAAGATCGGAA \
  -A CTTGTAGATCGGAAG \
  -A GACCAAGATCGGAAG \
  -A TTGTAGATCGGAAGA \
  -A ACCAAGATCGGAAGA \
  -A TGTAGATCGGAAGAG \
  -A CCAAGATCGGAAGAG \
  -A GTAGATCGGAAGAGC \
  -A CAAGATCGGAAGAGC \
  -A TAGATCGGAAGAGCG \
  -A AAGATCGGAAGAGCG \
  -A AGATCGGAAGAGCGT \
  -A GATCGGAAGAGCGTC \
  -A ATCGGAAGAGCGTCG \
  -A TCGGAAGAGCGTCGT \
  -A CGGAAGAGCGTCGTG \
  -A GGAAGAGCGTCGTGT \
  -o $ADAPT1/$F1 \
  -p $ADAPT1/$F2 \
  ${i} \
  ${SAMPLE}_R2.umi.fastq \
> $ADAPT1/MOLM13_Lo_R2_L001_round1.adapt.metrics 
done





echo 'adaptor removal from L001, second round' 


for i in $ADAPT1/*_R1.adapt1.fastq*;
do
SAMPLE=$(echo ${i} | sed "s/_R1.adapt1.fastq//")
#echo ${SAMPLE}_R1.fastq ${SAMPLE}_R2.fastq
F1=$(basename "${SAMPLE%umi.fastq}_R1.adapt2.fastq") F2=$(basename "${SAMPLE%umi.fastq}_R2.adapt2.fastq")

cutadapt -f fastq \
--match-read-wildcards \
--times 1 \
-e 0.1 \
-O 1 \
--quality-cutoff 6 \
-m 15 \
-A AACTTGTAGATCGGA \
-A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA \
-A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG \
-A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA \
-A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG \
-A GTAGATCGGAAGAGC \
-A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG \
-A AGATCGGAAGAGCGT \
-A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG \
-A TCGGAAGAGCGTCGT \
-A CGGAAGAGCGTCGTG \
-A GGAAGAGCGTCGTGT \
-o $ADAPT2/$F1 \
-p $ADAPT2/$F2 \
${i} \
${SAMPLE}_R2.adapt1.fastq \
> $ADAPT2/MOLM13_Lo_R2_L001_round2.adapt.metrics 
done
########################
echo 'pre-alignmnet finished sucessfuly'
