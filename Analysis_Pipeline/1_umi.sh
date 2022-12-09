#!/bin/bash -l
#
#
#
export FASTQ=~/Final_fastq_files
export UMI=/scratch1/users/shashank.tiwari/Final_fastq_files/Umi
#export ADAPT1=~/Final_fastq_files/Lo_R2/Adapt1
#export ADAPT2=~/Final_fastq_files/Lo_R2/Adapt2
#
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L001_R1_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L001_R1_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L001_R2_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L001_R2_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L002_R1_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L002_R1_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L002_R2_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L002_R2_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L003_R1_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L003_R1_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L003_R2_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L003_R2_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L004_R1_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L004_R1_001.fastq && \
gunzip < $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L004_R2_001.fastq.gz > $FASTQ/MOLM13_GAPDH_Lo_R2_S1_L004_R2_001.fastq && \

for f in ~/Final_fastq_files/*.gz; do
  STEM=$(basename "${f}" .gz)
  gunzip -c "${f}" > /scratch1/users/shashank.tiwari/Final_fastq_files/${STEM}
done


echo 'umi-extract from L001 started' 

for i in /scratch1/users/shashank.tiwari/Final_fastq_files/*_R1.fastq
do
SAMPLE1=$(echo ${i} | sed "s/_R1.fastq//")
F1=$(basename "${SAMPLE1%fastq}_R1.umi.fastq") F2=$(basename "${SAMPLE1}_R2.umi.fastq")
echo "Start" ${i}
umi_tools extract \
--random-seed 1 \
--bc-pattern=NNNNNNNNNNNN \
--bc-pattern2=NNNNNNNNNN \
--log $UMI/L001-umi.metrics \
-I ${i} \
-S $UMI/$F1 \
--read2-in=${SAMPLE1}_R2.fastq \
--read2-out=$UMI/$F2
echo "end" $i
done

echo 'umi_extract finished' \
