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
