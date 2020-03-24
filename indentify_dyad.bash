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
