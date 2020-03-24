#!/bin/bash
module load bedops/2.4.35
module load crossmap/0.3.3
module load bedtools/2.27.1
for ((i=1;i<=22;i++));
do
wig2bed < chr${i}.wig > chr${i}.bed
# Filter out the fragments with 0 reads
awk '{OFS="\t"}{ if($5>0){print $1, $2, $3, $5}}' chr${i}.bed > chr${i}_filtered.bed
crossmap bed hg18ToHg19.over.chain chr${i}.bed chr${i}_hg19.bed 
## intersect the dyad position with positioning socre
bedtools intersect -a chr${i}_hg19.bed -b mnase_mids_stringency_hg19_sort.bed > chr${i}_hg19_score.bed
done
