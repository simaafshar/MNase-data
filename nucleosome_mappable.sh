#!/bin/bash
module load bedtools/2.29.0
# Convert the mappability file to bed
bigWigToWig wgEncodeCrgMapabilityAlign36mer.bigWig wgEncodeCrgMapabilityAlign36mer.wig
wig2bed --zero-indexed < wgEncodeCrgMapabilityAlign36mer.wig | gzip > wgEncodeCrgMapabilityAlign36mer.bed.gz
# Filter unmappable regions and get the regions with Alignalibility score =1
zcat < wgEncodeCrgMapabilityAlign36mer.bed.gz | awk '$5==1'  > wgEncodeCrgMapabilityAlign36mer_score1.bed


cat mnase_mids_stringency_hg19_sort.bed | sort -k1,1 -k2,2n | \
    awk '{OFS="\t";}{print $1, $2-73, $3+73, $4, $5, $6}'  | \
    intersectBed -a stdin -b wgEncodeCrgMapabilityAlign36mer_score1.bed -wo | \
    awk '{OFS="\t";}{print $1, $2+73, $3-73, $4, $5, $6, $7, $8, $9, $10, $11,$12}' |  \
    sort -u > dyads_mappable.bed

## find the nuclesome with at least 75% of the sequence wityh score =1
python3 Nucleosome_mappable.py
awk 'BENGIN{FS='\t'}{if($5>=0.75) print $0}' dyads_mappability_score1_ratio_all.bed > dyads_mappability_score1_ratio_0.75.bed
