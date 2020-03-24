%%bash

source activate env_nucperiod
scripts=${PWD}/scripts

wig2bed < mnase_mids_combined_147.wig | gzip > mnase_mids_combined_147.bed.gz
# Filter the fragments with 0 reads
zcat < mnase_mids_combined_147.bed.gz | awk '{OFS="\t"}{ if($5>0){print $1, $2, $3, $5}}' | \
    gzip > mnase_mids_filtered.bed.gz
# Generate a bed file for each chromosome
zcat < mnase_mids_filtered.bed.gz | awk '{print $0 >> "mnase_mids_filtered."$1".bed"}'
gzip -f mnase_mids_filtered.*.bed

# Compute the triweight kernel smoothing for each chromosome
# This part can be easily parallelized by launching one process per file
for file in mnase_mids_filtered.*.bed.gz
do
    f="${file/mnase_mids_filtered/mnase_mids_kernel_smoothing}"
    smooth_file="${f/.bed./.tsv.}"
    normalize_file="${smooth_file/mnase_mids_kernel_smoothing/mnase_mids_kernel_smoothing_normalized}"
    python ${scripts}/kernel_smoothing.py ${file} ${smooth_file}
    python ${scripts}/normalize_kernel.py ${smooth_file} ${normalize_file}
done


#. Concatenate all the chromosome split bed files into a single bed file
zcat < mnase_mids_kernel_smoothing_normalized.chr1.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr10.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr11.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr12.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr13.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr14.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr15.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr16.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr17.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr18.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr19.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr2.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr20.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr21.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr22.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr3.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr4.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr5.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr6.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr7.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr8.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chr9.tsv.gz \
    mnase_mids_kernel_smoothing_normalized.chrX.tsv.gz \
    > mnase_mids_kernel_smoothing_normalized.bed

# Get the positions of the local maxima
bedGraphToBigWig mnase_mids_kernel_smoothing_normalized.bed hg18.chrom.sizes \
    mnase_mids_kernel_smoothing_normalized.bw
gzip -f mnase_mids_kernel_smoothing_normalized.bed
bwtool find local-extrema -maxima -min-sep=150 mnase_mids_kernel_smoothing_normalized.bw \
    mnase_mids_local_maxima.bed
gzip -f mnase_mids_local_maxima.bed
# The above command rounds the scores, and we need to get those back
intersectBed -a mnase_mids_local_maxima.bed.gz -b mnase_mids_kernel_smoothing_normalized.bed.gz -wo -sorted | \
    awk '{OFS="\t"}{print $1, $2, $3, $10, $5, $6}' |gzip > mnase_mids_local_maxima.scores.bed.gz

# Extend the local maxima 30 and intersect with mid-fragments
zcat < mnase_mids_local_maxima.scores.bed.gz |
    awk '{OFS="\t";}{print $1, $2-30, $3+30, $1 "_" $2 "_" $3, $4}' | \
    intersectBed -a mnase_mids_filtered.bed.gz -b stdin -wo | \
    gzip > mnase_mids_local_maxima_intersect.bed.gz

# Get the midpoints
python ${scripts}/get_dyads.py mnase_mids_local_maxima_intersect.bed.gz \
    mnase_mids_stringency.bed.gz

# Lift over the coordinates from hg18 to hg19
source deactivate
source activate env_crossmap
CrossMap.py bed hg18ToHg19.over.chain.gz mnase_mids_stringency.bed.gz mnase_mids_stringency_hg19.bed
gzip -f mnase_mids_stringency_hg19.bed
source deactivate
source activate env_nucperiod

# Compute a regions file containing mappable non-genic regions
cat gencode.v19.annotation.gtf |grep protein_coding | grep -w gene | \
    awk '{OFS="\t";} {print $1, $4-501, $5+500}' |gzip > genes_extended.bed.gz
# Convert the mappability file to bed
bigWigToWig wgEncodeCrgMapabilityAlign36mer.bigWig wgEncodeCrgMapabilityAlign36mer.wig
wig2bed --zero-indexed < wgEncodeCrgMapabilityAlign36mer.wig | gzip > wgEncodeCrgMapabilityAlign36mer.bed.gz
# Filter unmappable regions
zcat < wgEncodeCrgMapabilityAlign36mer.bed.gz | awk '$5==1' | gzip > wgEncodeCrgMapabilityAlign36mer_score1.bed.gz
# Generate a mappability file without genic regions
subtractBed -a wgEncodeCrgMapabilityAlign36mer_score1.bed.gz -b genes_extended.bed.gz |\
    gzip > coverage.bed.gz

# Get only nucleosomes in mappable non-genic regions
zcat < mnase_mids_stringency_hg19.bed.gz | sort -k1,1 -k2,2n | \
    awk '{OFS="\t";}{print $1, $2-73, $3+73, $4, $5, $6}'  | \
    intersectBed -a stdin -b coverage.bed.gz -f 1 -sorted | \
    awk '{OFS="\t";}{print $1, $2+73, $3-73, $4, $5, $6}' |  \
    gzip > dyads.bed.gz
    
# Get the coordinates of the exons
grep -w exon gencode.v19.annotation.gtf | grep protein_coding |\
    gzip > gencode.v19.exon_protein_coding.gtf.gz
python ${scripts}/exons.py gencode.v19.exon_protein_coding.gtf.gz \
    cancer_gene_census.csv exons.bed.gz
mergeBed -i exons.bed.gz | gzip > exons.merged.bed.gz
    
# Get only nucleosomes falling in genic regions
zcat <  mnase_mids_stringency_hg19.bed.gz  | \
    awk '{OFS="\t"}{print $1, $2-58, $3+58}' | \
    intersectBed -a stdin -b exons.merged.bed.gz -f 1 | \
    awk '{OFS="\t"}{print $1, $2+58, $3-58}' | \
    sort -k1,1 -k2,2n | gzip > dyads_genic.bed.gz

