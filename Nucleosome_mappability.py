#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:42:36 2020

@author: pengy10
"""
import pandas as pd 
import numpy as np 
import re
import csv
import statistics 
import seaborn as sns
###select nucleosome with at least 75% sequence have mappability score =1

data = pd.read_csv("dyads_mappable.bed",  sep = "\t", header=None) ## nucleosmal DNA fragments with score1=1
#split by chr
hg_files = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chr10", "chr11", "chr12", "chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", 
             "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

mapability_chr = {}
for chr_name in hg_files:
    mapability_chr[chr_name] = data[data.iloc[:,0] == chr_name]


with open("dyads_mapibility_score1_ratio_all.bed",'w') as file:
    for chr_name in hg_files:
        mapability = mapability_chr[chr_name]
        for i in range(1,len(mapability) -1):
            pos0 = mapability.iloc[i-1][1]
            pos1 = mapability.iloc[i][1]
            pos2 = mapability.iloc[i+1][1]
            sum_len = mapability.iloc[i][11]
            if pos0 != pos1:
                while pos2 == pos1 and i<=len(mapability) -3:
                    i = i +1
                    sum_len = sum_len + mapability.iloc[i][11] ## add up the framents with socre =1 for each nucleosome
                    pos1 = mapability.iloc[i][1]
                    pos2 = mapability.iloc[i +1 ][1]
                score1_ratio = sum_len/147 ## calculate the percentage of sequence with score =1 for each nucleosme
                file.write(str(mapability.iloc[i][0]) +"\t" + str(mapability.iloc[i][1]) + "\t" + str(mapability.iloc[i][2]) + "\t" \
                           + str(mapability.iloc[i][3]) + "\t" + str(score1_ratio) +"\n")
                
 

"""                             
#Find the linker regions:
##Length between 0 to 500
data = pd.read_csv("mnase_mids_stringency_hg19_sort.bed",  sep = "\t", header=None)
#split by chr
hg_files = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chr10", "chr11", "chr12", "chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", 
             "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

dyads_chr = {}
for chr_name in hg_files:
    dyads_chr[chr_name] = data[data.iloc[:,0] == chr_name]
#hg_files = ["chrX"]

with open("dyads_linker.bed",'w') as file:
    for chr_name in hg_files:
        dyads = dyads_chr[chr_name]
        for i in range(0,len(dyads) -1):
 #       for i in range(0,1000):
            pos0 = dyads.iloc[i][1]
            pos1 = dyads.iloc[i+1][1]
            linker_len = pos1-pos0 -147  ## length of linker DNA
            if linker_len >=0 and linker_len <=500:
                file.write(str(dyads.iloc[i][0]) +"\t" + str(dyads.iloc[i][1]+74) + "\t" + str(dyads.iloc[i+1][1]-73) + "\t" \
                           + str(dyads.iloc[i][3]) + "\t" + str(linker_len) +"\n")


###select linker regions with at least 75% sequence have mappability score =1

data = pd.read_csv("dyads_linker_score1.bed",  sep = "\t", header=None)
#split by chr
hg_files = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
             "chr10", "chr11", "chr12", "chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", 
             "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
mapability_chr = {}
for chr_name in hg_files:
    mapability_chr[chr_name] = data[data.iloc[:,0] == chr_name]


with open("dyads_linker_score1_ratio_all.bed",'w') as file:
    for chr_name in hg_files:
        mapability = mapability_chr[chr_name]
        for i in range(1,len(mapability) -1):
 #       for i in range(0,1000):
            pos0 = mapability.iloc[i-1][1]
            pos1 = mapability.iloc[i][1]
            pos2 = mapability.iloc[i+1][1]
            sum_len = mapability.iloc[i][10]
            if pos0 != pos1:
                while pos2 == pos1 and i<=len(mapability) -3:
                    i = i +1
                    sum_len = sum_len + mapability.iloc[i][10]
                    pos1 = mapability.iloc[i][1]
                    pos2 = mapability.iloc[i +1 ][1]
                score1_ratio = sum_len/mapability.iloc[i][4]
                file.write(str(mapability.iloc[i][0]) +"\t" + str(mapability.iloc[i][1]) + "\t" + str(mapability.iloc[i][2]) + "\t" \
                           + str(mapability.iloc[i][3]) + "\t" + str(mapability.iloc[i][4]) + "\t" + str(score1_ratio) +"\n")
"""
"""## Box plot of the linker length
df_linker = pd.read_csv("dyads_linker_score1_ratio_0.75.bed",  sep = "\t", header=None,names=["chr","start","end","reads","linker_length","ratio"])
## box plot to show the distribitions
bplot1=sns.boxplot(y='linker_length', x='chr', 
                 data=df_linker, 
                 width=0.5,
                 palette="colorblind")
 
# add stripplot to boxplot with Seaborn
#bplot1=sns.stripplot(y='linker_length', x='chr', 
#                   data=df_linker, 
#                   jitter=True, 
#                   marker='o', 
#                   alpha=0.5,
#                   color='black')
bplot1.set(xlabel = "chr", ylabel = "linker_length")
bplot1.get_figure().savefig('linker_length.png',  dpi = 300 )


# Density Plot and Histogram of all linker length
bplot2=sns.distplot(df_linker['linker_length'], hist=True, kde=True, 
             bins=int(500/5), color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
bplot2.get_figure().savefig('linker_length_density.png',  dpi = 300 )
"""
                
