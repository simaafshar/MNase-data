## software and packages required
python3
crossmap
bedops
bedtools
bigwigtowig
bwtool 

# Idenfity the representative dyad position from MNase-data: indentify_dyad.bash
The Methods and scripts for identifying the representative nucleosme dyads. The scripts are originally taken from "Somatic and Germline Mutation Periodicity Follow the Orientation of the DNA Minor Groove around Nucleosomes,DOI: 10.1016/j.cell.2018.10.004" and you can find the original scripts from: https://nbviewer.jupyter.org/urls/bitbucket.org/bbglab/nucleosome-periodicity/raw/master/nucleosomes/nucleosomes.ipynb


#How to get the files to run the scripts:

The 147 bp length mid-fragments MNAse-seq reads which has been downloaded from http://eqtl.uchicago.edu/nucleosomes/midpoints/mnase_mids_combined_147.wig.gz
Uncompress this file before running the notebook

The file with the hg18 chromosome sizes downloaded from https://github.com/igvteam/igv/blob/master/genomes/sizes/hg18.chrom.sizes

The hg18 to hg19 chain file downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

The CRG mappability track, that can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig

The positions of the genes downloaded from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
Uncompress this file before running the notebook

The Cancer gene census from COSMIC v81, You can find information on how to get this version in https://cancer.sanger.ac.uk/cosmic/help/download#cmd


## Classify the nucleosme by positioning score: position_score.bash
The representative dyad position identified from MNase data:
mnase_mids_stringency_hg19_sort.bed
The Positioning score was taekn from the Prirchard's Plos Genetices paper and you can get it from:
http://eqtl.uchicago.edu/nucleosomes/positioning_scores/mnase_mids_combined_142_to_152_scores/

We intersect the scores with the indentified representative dyad position and classfied as three score ranges:
0 - 0.3: weak position  positioning_score_hg19_0-0.3.bed
0.3 - 0.5: inermediated position positioning_score_hg19_0.3-0.5.bed
0.5 - 1: strong position positioning_score_hg19_0.5.bed

## Filtering nucleosme by mapability score: Nucleosome_mappable.sh and Nucleosome_mappable.py
The mappability score from the CRG36 alignability track in hg19:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign36mer.bigWig
wgEncodeCrgMapabilityAlign36mer_score1.bed.gz:
region in hg19 with mapability score =1, indicating high uniquenss in mapping

