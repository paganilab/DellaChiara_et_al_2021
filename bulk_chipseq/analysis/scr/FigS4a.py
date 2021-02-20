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

# ### FigS3 a - pie chart genomic features

# Author: Federica Gervasoni

# See the associated pbs file to extract genomic features from gtf and subset bed file.

import pandas as pd
import numpy as np
import glob 
import os
import seaborn as sns
import matplotlib.pyplot as plt
import math
import scipy.stats as stats

path='/annotation_genomicFeatures/'

df=pd.read_table(path+"K_org_UP_log2FC2_padj0.01_genomicFeatures.bed", sep="\t",  names=['chr', 'start', 'end', 'reg_ID', 'chr_gf', 'start_gf', 'end_gf', 'feature', 'olap'])

df_count=df.groupby('feature').feature.size()

# +
explode = (0.05, 0.05, 0.05, 0.05)
colors = ['#F2CD13','#5DD0FF', '#BFA7F2',  '#8EC63F',]

fig1, ax = plt.subplots()
df_count.plot.pie(subplots=True,figsize=(8, 3), explode=explode, autopct='%1.1f%%', startangle=45, shadow=True,
                            colors=colors, ax=ax)
plt.show()
