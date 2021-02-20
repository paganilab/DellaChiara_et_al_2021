# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: design_plot_new
#     language: python
#     name: design_plot_new
# ---

# ### Study of enhancers shared by different PDOs
#
# Author: Federica Gervsoni

import pandas as pd
import numpy as np
import glob 
import os
import seaborn as sns
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Set3_10
from matplotlib.colors import ListedColormap
import math
import scipy.stats as stats
import scipy
from scipy.stats import mannwhitneyu
import scikits.bootstrap as bootstrap

# ### Recurrence in PDO

path='/recurrence_Korg/'

ref_geneID="NULL"
ref_geneID=pd.read_table("/PDO/PDO_enh_reference.bed", names=['chr', 'start', 'end', 'Ref_peak'])
ref_geneID['diff']=ref_geneID['end']-ref_geneID['start']

# +
allFiles = glob.glob("*_mtx.bed")
frame = pd.DataFrame()
list_ = []
for fileName in allFiles:
     df = pd.read_table(fileName,  names=['chr', 'start', 'end', 'Ref_peak'])
     parts = fileName.split('_chromHMM_enh_mtx.bed')[0]
     ref_geneID[parts] = ref_geneID.Ref_peak.isin(df.Ref_peak).astype(np.int8)

ref_geneID['recurrence'] = ref_geneID.isin({1}).sum(1)
ref_sorted=ref_geneID.sort_values('recurrence', ascending=False, axis=0)
#46508
# -

# ### Import Differential enhancer in PDOs

path_DE='/DEA/DBA_Korg_Ntissue_Input_Jup/'
recur_DE=pd.read_table(path_DE + "recurrence_Korg/all_recurrence_table.txt",  names=['chr_DE', 'start_DE', 'end_DE', 'region_ID', 'chr_rec', 'start_rec', 'end_rec', 'Ref_peak', 'olap'],index_col=False)


# +
def func1(r):
    return 'Yes'

def func2(r):
    return 'No'


# -

ref_sorted_DE=pd.merge(recur_DE, ref_sorted, on='Ref_peak', how='right')
ref_sorted_DE['Ref_peak_group']=ref_sorted_DE.apply(lambda row: func1(row) if row.notnull().all() else func2(row), axis=1)

ref_sorted_DE_rec=ref_sorted_DE.dropna()
ref_sorted_DE_rec=ref_sorted_DE_rec[ref_sorted_DE_rec.recurrence != 1]

# +
groups=["04_KE","08_KE","10_KE","11_KW","13_KE","18_KE","22_KE","24_KE","36_KE","41_KE"]

list_=[]
for group in groups:
    df=ref_sorted_DE_rec.groupby(["recurrence", group]).size().reset_index(name='occ_'+group)
    df1=df.loc[(1,3,5,7,9,11,13,15,16), 'occ_'+group]
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)
    
frame2 = pd.concat(list_, axis=1)
frame2.index=([2,3,4,5,6,7,8,9,10])
# -

# ### Fig3 c - Piechart to plot the recurrent which are DE, dividing them in different class.

freq_df = ref_sorted_DE.groupby(['recurrence', 'Ref_peak_group']).size().unstack()
freq_df=freq_df[['Yes']]
freq_df['rowToSum']=['1To4','1To4','1To4','1To4','5to7','5to7','5to7','8to10','8to10','10to10']
freq_df_pie=freq_df.groupby(['rowToSum'])['Yes'].sum().reset_index()
freq_df_pie.set_index('rowToSum', inplace=True)

# +
explode = (0.2, 0.03, 0.03, 0.03)
colors = ['#590c29','#FF9E11','#FF5733','#C70039', ]

fig1, ax = plt.subplots()
freq_df_pie.plot.pie(subplots=True,figsize=(8, 3), autopct='%1.1f%%', startangle=45, shadow=True, ax=ax,
                    explode=explode, colors=colors)
plt.show()
# -

# ### FigS4 b - Check the distribution among patient of the DE

# +
#cmap
cmap_rec = ListedColormap(Set3_10.mpl_colors)

#plot
fig, ax = plt.subplots()

frame2.T.divide(frame2.T.sum(axis=1), axis=0).plot(kind='bar', stacked=True, ax=ax, cmap=cmap_rec)
ax.set_prop_cycle('color', Set3_10.mpl_colors)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_title('Distribution of the DE conserved regions across patients')
ax.set(xlabel='Patients')
# -

# ### FigS4 c - raw count of recurrent enhancer in 8 to 10 out of 10

recur_DE="NULL"
recur_DE=pd.read_table(path + "486_recurrence.txt", names=['chr_DE', 'start_DE', 'end_DE', 'region_ID', 'chr_rec', 'start_rec', 'end_rec', 'Ref_peak', 'olap'],index_col=None)

df_rec_DE=pd.merge(ref_sorted,recur_DE)

# +
groups=["04_KE","08_KE","10_KE","11_KW","13_KE","18_KE","22_KE","24_KE","36_KE","41_KE"]

list_=[]
for group in groups:
    df=df_rec_DE.groupby(["recurrence", group]).size().reset_index(name='occ_'+group)
    df1=df.loc[(1,3,4), 'occ_'+group]
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)
    
frame = pd.concat(list_, axis=1)
frame.index=([8,9,10])


# +
groups=["04_KE","08_KE","10_KE","11_KW","13_KE","18_KE","22_KE","24_KE","36_KE","41_KE"]

list_=[]
for group in groups:
    df=df_rec_DE.groupby(["recurrence", group]).size().reset_index(name='occ_'+group)
    df1=df.loc[(1,3,4), 'occ_'+group]
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)
    
frame = pd.concat(list_, axis=1)
frame.index=([8,9,10])

# Figure plot
fig, ax = plt.subplots()

frame.T.plot(kind='bar', stacked=True, ax=ax)
ax.set_prop_cycle('color', Set3_10.mpl_colors)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_title('Distribution of the conserved regions across patients')
ax.set(xlabel='Patients')
plt.show()
# -

# ### Fig4 f - WWTR1 olap with chromHMM (each single state)

# +
# Olap of wwtr1 with all the chromHMM states
# bedtools intersect -a chromhmm_dense.bed -b wwtr1_noBL_noUN.narrowPeak > chromhmm_dense_wwtr1.bed
# -

chromhmm_WWTR1=pd.read_table("/chromHMM/chromhmm_dense_wwtr1.bed", sep="\t",  names=['chr', 'start', 'end', 'chromHMM', 'score', 
                                                                                     'position', 'chrEnh', 'startEnh', 'color'])

chromhmm_WWTR1['states']=chromhmm_WWTR1['chromHMM']
chromhmm_WWTR1['states']=chromhmm_WWTR1.states.replace({'8_FlkTSS': 'Promoter', '2_Quies': 'Quiescent', '7_ActTSS': 'Promoter', 
                                                        '4_WkEnh': 'Enhancers', '6_FlkActEnh': 'Enhancers', '5_ActEnh': 'Enhancers', '3_Elong': 'Elongation', '1_Repr': 'Repression' })

chromhmm_WWTR1_df=chromhmm_WWTR1.groupby('states').chromHMM.size()

# +
explode = (0.05, 0.05, 0.05, 0.05, 0.05)
colors = ['#33a02c','#ff7f00', '#e31a1c', '#f6f6f6', '#696969']

fig1, ax = plt.subplots()
chromhmm_WWTR1_df.plot.pie(subplots=True,figsize=(8, 3), explode=explode,  autopct='%1.1f%%', startangle=45, shadow=True,
                            colors=colors, ax=ax)
plt.show()
# -

# ### Fig4 g - Randomic distribution:

# +
path='/DEA/DBA_Korg_Ntissue_Input_Jup/taz_shuffle/'

name="_wwtr1_uniq"

df1=pd.read_table(path+"K_org_UP_log2FC2_padj0.01_randomOlap" + name +".txt",  names=['NumOlap'])
df1["type"]="3_K_org_DE_random"
df1["group"]="K_org_DE"
df1["group_2"]="Tumor"
df1["group_3"]="Random"
df1["NumOlap_n"]=df1.loc[:,"NumOlap"].div(2419, axis=0)*100

df2=pd.read_table(path+"recurrence_10-5_K_org_UP_log2FC2_padj0.01_randomOlap" + name +".txt",  names=['NumOlap'])
df2["type"]="5_rec_10-5_DE_random"
df2["group"]="rec_10-5_DE"
df2["group_2"]="Tumor"
df2["group_3"]="Random"
df2["NumOlap_n"]=df2.loc[:,"NumOlap"].div(1216, axis=0)*100

df3=pd.read_table(path+"recurrence_10-8_K_org_UP_log2FC2_padj0.01_randomOlap"+ name +".txt",  names=['NumOlap'])
df3["type"]="7_rec_10-8_DE_random"
df3["group"]="rec_10-8_DE"
df3["group_2"]="Tumor"
df3["group_3"]="Random"
df3["NumOlap_n"]=df3.loc[:,"NumOlap"].div(486, axis=0)*100

# +
#wwtr1 unqiue counts
d = {'NumOlap': [446,339, 195], 'NumOlap_n': [(446/2419)*100,(339/1216)*100, (195/486)*100], 
     'type': ["4_K_org_DE", "6_rec_10-5_DE", "8_rec_10-8_DE"], 'group':["K_org_DE","rec_10-5_DE", "rec_10-8_DE"] }

df4 = pd.DataFrame(data=d)
df4["group_2"]="Tumor"
df4["group_3"]="Real"
# -

result = df.append([df1,df2, df3,df4])

# +
# https://stackoverflow.com/questions/38892385/scypy-wilcoxon-test-compare-distribution-with-a-single-value
#If the value is within the CIs, then you can't reject that your value is the real mean (or median).
# -

cluster2color={'colorRec': {'3_K_org_DE_random': '#f2964b',
  '5_rec_10-5_DE_random': '#f47d7a',
  '7_rec_10-8_DE_random': '#ce0c4d'}}
cluster2color['colorRec']['3_K_org_DE_random']

# +
types = ['3_K_org_DE_random',
 '5_rec_10-5_DE_random',
 '7_rec_10-8_DE_random']

fig = plt.figure()

# Iterate through the five airlines
for i in types:
    # Subset to the airline
    subset = result[result['type'] == i]

    # Draw the density plot
    ax=sns.distplot(subset['NumOlap'], hist = False, kde = True,
                 kde_kws = {'linewidth': 3},
                 label = i, color=cluster2color['colorRec'][i])
    plt.axvline(result.loc[result['type'] == '4_K_org_DE']['NumOlap'][1], 0, 1, linestyle='--', color='#e26b09')
    plt.axvline(result.loc[result['type'] == '6_rec_10-5_DE']['NumOlap'][2], 0, 1, linestyle='--', color='#EC3A35')
    plt.axvline(result.loc[result['type'] == '8_rec_10-8_DE']['NumOlap'][3], 0, 1, linestyle='--', color='#8e0734')
    
    CIs = bootstrap.ci(subset['NumOlap'], statfunction=np.mean, n_samples=10000)  
    print(i, '-ci: ', CIs)
    print(i, '-mean: ', subset['NumOlap'].mean() )
    print(i, '-max: ', subset['NumOlap'].max() )
    
ax.set_xlim(0,470)
# -

# ### Fig4 g - Barplot with the number of real occurence and the randomic one

result_df=result.groupby('type').mean().round()
result_df=result_df.iloc[[2,3,4,5,6,7]]

neworder=result_df.index.to_list()
newcolors=['#f2964b','#e26b09','#f47d7a','#EC3A35','#ce0c4d','#8e0734']

# +
fig = plt.figure()

ax = sns.barplot(x=result_df.index, y="NumOlap_n", data=result_df, palette =newcolors, order=neworder)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_title('Overlap YAP on randomic and real data normalized on the total number of regions')
ax.set_ylabel('#Overlap Normalized')
ax.set_ylim(0,50)
# -

# ### FigS5 c - Check intensity of YAP/TAZ in different recurrence category

# +
# #%%bash
#bedtools intersect -a $path/recurrence_Korg/all_recurrence_table.bed -b $dataset -wo | cut -f 1,2,3,4,8,9,14,17,18,19 > $path/all_recurrence_table_wwtr1.bed
# -

path='/DEA/DBA_Korg_Ntissue_Input_Jup/'

df=pd.read_csv(path + 'all_recurrence_table_wwtr1.bed', sep='\t', header=None)
df.columns=['chr', 'start', 'end', 'region_ID', 'chromHMM_id','rec', 'taz_peak_id',  'signalValue', 'pValue', 'qValue']
print('df-shape: ', df.shape)

#keep unique values
df=df.sort_values('pValue', ascending=False).drop_duplicates('region_ID').sort_index()

# +
df['g_rec']=df['rec']
df['g_rec']=df['g_rec'].replace({10: '>5', 9:'>5', 8:'>5', 7:'>5', 6:'>5', 5:'<5', 4:'<5', 3:'<5', 2:'<5'})

df['g_rec_8']=df['rec']
df['g_rec_8']=df['g_rec_8'].replace({10: '>8', 9:'>8', 8:'>8', 7:'na', 6:'na', 5:'na', 4: 'na', 3:'na', 2:'na'})

list_val=['signalValue']

# +
fig, axes = plt.subplots(1, 1, figsize=(3, 6))

##Add levels for regions with >5 recurrent    
y_val = 0

x=df.loc[df['g_rec']=='<5']['signalValue']

y_val+=0.5
y=[y_val]*len(x)

plt.boxplot(x, positions=[y_val], showfliers=False, notch=True, widths=0.3) 

##Add levels for regions with >8 recurrent    
y_val = 0.5

x=df.loc[df['g_rec_8']=='>8']['signalValue']

y_val+=0.5
y=[y_val]*len(x)

plt.boxplot(x, positions=[y_val], showfliers=False, notch=True, widths=0.3) 
    
labels=['<5', '>8']

#title
plt.title('TAZ levels across recurrent and DE enhancers')
plt.xlabel('Recurrent DE regions')
plt.ylabel('TAZ levels')

#tick info
plt.xticks(ticks=[0.5, 1], labels=labels)
plt.xlim(0,1.5)
plt.show()
# -

x_less5=df.loc[df['g_rec']=='<5']['signalValue']
x_more8=df.loc[df['g_rec_8']=='>8']['signalValue']

print('comparison less5 and more8: ', (mannwhitneyu(x_less5, x_more8)))
