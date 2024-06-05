for i in *_bin.bed            
do
F1=$(basename "${i}" | sed "s/_bin.bed//")
#awk 'BEGIN{FS=OFS="\t"} {print $1, ($2-49), ($3+50), $4, $5, $6}' ${i} | tr ' ' '\t' | cut -f1,2,3,4,10,11 | sed "s/;//g" | sort -k1,1V -k2,2g  > ${F1}_fin.bed

done

#bedtools intersect -a - -b ~/CLIP_Seq/Data/Genome/gene_loc_file.bed -wa -wb 
#❯ bedtools intersect -a Norm_R1_fin.bed -b Norm_R2_fin.bed |  bedtools intersect -a - -b ~/CLIP_Seq/Data/Genome/gene_loc_file.bed -wa -wb | cut -f1,2,3,4,8,9 > Merged_Norm_annotated.bed
