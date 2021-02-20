# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: design_plot
#     language: python
#     name: design_plot
# ---

# ### Verify if ATAC-seq peaks from COAD TCGA data are more enriched in the open chromHMM states for each individual PDO
#
# Author: Federica Gervasoni and Lorenzo Drufuca

# #### Import packages

import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as plt
import matplotlib.pyplot as plt
import glob
import scipy as sp
import numpy as np

# +
### To obtain the files used for this analysis you need to intersect dense files produced by chromHMM with TCGA data

##bash

# > zcat $chromHMM_dense | sed -n '1d;p' | awk '{print NR,$0} ' |  awk '$1=("chromState_" $1)' | awk -v OFS='\\t' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' > ${name}_chromHMM_mod.bed 

# > bedtools intersect -a ${name}_chromHMM_mod.bed -b ${params.genome.COAD_ATACseq_TCGA} -wo > ${name}_chromHMM_acces_COAD_TCGA.bed

# > zcat $chromHMM_dense | sed -n '1d;p' | awk '{print NR,$0} ' |  awk '\$1=("chromState_" $1)' | awk -v OFS='\\t' '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | bedtools intersect -a - -b ${COAD_ATACseq_TCGA} -c > ${name}_chromHMM_acces_COAD_TCGA_count.bed
# -

path_chromHMM='/accessibility_TCGA/K_org/'
path_save='/accessibility_TCGA/'

# Prepare the function to import chromHMM overlapped with ATAC-seq data
names=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE' ]
suffix='_chromHMM_acces_COAD_TCGA_count.bed'

# ### Count overlap between chromHMM and ATAC-seq

# +
# Import the chromHMM overlapped with ATAC-seq data + change index and columns names
list_=[]
for name in names:
    df = pd.read_csv(path_chromHMM+'/'+name+suffix,  skiprows=[0], header=None,sep='\t', 
                     names=['chr', 'start', 'end', 'states', 'unkn', 'unkn1', 'unkn2', 'unkn3', 'color', 'numState', 'count'])
    df1=df.groupby('states')['count'].agg('sum')
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)

frame = pd.concat(list_, axis=1)
frame=frame.rename({'7_ActTSS':'1_ActTSS', '8_FlkTSS':'2_FlkTSS', '6_FlkActEnh':'3_FlkActEnh', '5_ActEnh':'4_ActEnh', '4_WkEnh':'5_WkEnh', '3_Elong':'6_Elong','2_Quies':'8_Quies', '1_Repr': '7_Repr' }, axis='index')
frame.columns=names

# +
# Import the chromHMM files and measure the length (bp) of each state for each PDO + change index and columns names
list_=[]
for name in names:
    df = pd.read_csv(path_chromHMM+'/'+name+suffix,  skiprows=[0], header=None,sep='\t', 
                     names=['chr', 'start', 'end', 'states', 'unkn1', 'unkn2', 'unkn3', 'unkn4', 'color', 'numState', 'count'])
    df['length']=abs(df['end']-df['start'])
    df1=df.groupby('states')['length'].agg('sum')
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)

frame_length = pd.concat(list_, axis=1)
frame_length=frame_length.rename({'7_ActTSS':'1_ActTSS', '8_FlkTSS':'2_FlkTSS', '6_FlkActEnh':'3_FlkActEnh', '5_ActEnh':'4_ActEnh', '4_WkEnh':'5_WkEnh', '3_Elong':'6_Elong','2_Quies':'8_Quies', '1_Repr': '7_Repr' }, axis='index')
frame_length.columns=names


# -

# ### Fig 2d - spider plot

# #### p(A|B)
#
# p(A) = ATAC-seq tot (1bp) / len(genome)
#
# p(B) = len(chromHMM state) / len(genome)
#
# p(A intersect B) = ATAC-seq overlap with chromHMM state / len(genome)
#
# p(A|B) = (ATAC-seq overlap with chromHMM state / len(genome)) / (len(chromHMM state) / len(genome))
#
#        = ATAC-seq overlap with chromHMM state / len(chromHMM state)

def my_func(n,x):
    t = (n/x)
    return(-np.log10(t))
vfunc = np.vectorize(my_func)

list_=[]
for c in frame.columns.values:
    out=pd.DataFrame(frame[c]/frame_length[c])
    list_.append(out)
    df2 = pd.concat(list_, axis=1)

# Scaling: sum the probability of each PDO + divide the probability of each state in that PDO for the sum the probability of each PDO
df2 = df2.divide(df2.sum(axis=0), axis=1)
order = ["1_ActTSS", "2_FlkTSS", "3_FlkActEnh", "4_ActEnh", "5_WkEnh", "6_Elong", "7_Repr","8_Quies"]
df2=df2.loc[order]

palette_PDO=['#999999', '#4DAF4A','#E41A1C','#377EB8', '#FF62BC','#984EA3','#FF7F00',  '#FFFF33', '#A65628','#65B7F3']

# ### Web plot

# +
from math import pi
df = df2.T

# number of variable
categories=df.columns.values.tolist()
N = len(categories)

# number of panels (subjects or clusters)
subjects = df.index.values.tolist()
npanels = len(subjects)

fig = plt.figure()
ax = plt.subplot(1,1,1, polar=True)

for p in range(npanels):
    
    # We are going to plot the first line of the data frame.
    # But we need to repeat the first value to close the circular graph:
    values=df.iloc[p,:].values.flatten().tolist()
    values += values[:1]
    values

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Draw ylabels
    ax.set_rlabel_position(0)
 
    ax.set_theta_direction(-1)
    ax.set_theta_offset(pi / 2)

    # Plot data
    ax.plot(angles, values, linewidth=3, linestyle='solid', label=subjects[p])
    
# Draw one axe per variable + add labels labels yet
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_yticks([])
 
plt.legend(loc='best', bbox_to_anchor=(1, 0., 1, 0.5))
plt.suptitle('Yessa')
plt.show()
