#!/bin/bash -l 
# star alignment of all samples

export PRE=/scratch1/users/shashank.tiwari/Final_fastq_files/Processed
export GENOME=/scratch1/users/shashank.tiwari/Final_fastq_files/Genome/
export ADAPT2=/scratch1/users/shashank.tiwari/Final_fastq_files/Adapt2
export ALIGN=/scratch1/users/shashank.tiwari/Final_fastq_files/Aligned
export DEDUP=/scratch1/users/shashank.tiwari/Final_fastq_files/Dedup


##################################################3
for i in $ADAPT2/*_R1.adapt2.fastq
do
SAMPLE=$(echo ${i} | sed "s/_R1.adapt2.fastq//")
F1=$(basename "${SAMPLE%adapt2.fastq}")
echo 'L001 alignmnet starts..' && \
STAR --runMode alignReads --runThreadN 16 \
--genomeDir $GENOME \
--readFilesIn ${i} \
${SAMPLE}_R2.adapt2.fastq \
--outFilterMultimapNmax 2 \
--outFilterMatchNmin 0 \
--outFilterMatchNminOverLread 0 \
--outFilterScoreMinOverLread 0 \
--outFilterMismatchNoverLmax 0.05 \
--outSAMtype BAM Unsorted \
--outReadsUnmapped Fastx \
--outFileNamePrefix $ALIGN/$F1
done
echo 'Done alignment'
####################################################
# merge aligned reads

for i in $ALIGN/*_L001Aligned.out.bam
do
SAMPLE1=$(echo ${i} | sed "s/_L001Aligned.out.bam//")
SAMPLE=$(basename "${SAMPLE1}")
F1=$(basename "${SAMPLE1}_merged.aligned.bam")
samtools merge $PRE/${SAMPLE}_merged.aligned.bam \
$ALIGN/${SAMPLE}_L001Aligned.out.bam $ALIGN/${SAMPLE}_L002Aligned.out.bam $ALIGN/${SAMPLE}_L003Aligned.out.bam $ALIGN/${SAMPLE}_L004Aligned.out.bam && \

#pre-dedup: sort and index bam file, make bed file from aligned merged file ###############################

samtools sort $PRE/${SAMPLE}_merged.aligned.bam > $PRE/${SAMPLE}_merged.aligned.sorted.bam && \
samtools index -b $PRE/${SAMPLE}_merged.aligned.sorted.bam $PRE/${SAMPLE}_merged.aligned.sorted.bam.bai && \
bamToBed -i $PRE/${SAMPLE}_merged.aligned.sorted.bam > $PRE/${SAMPLE}_merged.aligned.sorted.bed && \

# deduplication ########################
echo 'deduplication is started..' 
echo $PRE/${SAMPLE}_merged.aligned.sorted.bam
umi_tools dedup \
--random-seed 1 \
--method unique \
-I $PRE/${SAMPLE}_merged.aligned.sorted.bam \
--output-stats=cd deduplicated \
--paired \
-S $DEDUP/${SAMPLE}_merged.dedup.bam
echo 'deduplication finished sucessfully'
done



