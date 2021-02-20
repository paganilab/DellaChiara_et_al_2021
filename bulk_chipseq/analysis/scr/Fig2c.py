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

# ### Prouce plots associated to chromHMM analysis
#
# Author: Federica Gervasoni

# #### Import packages

import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as plt
import matplotlib.pyplot as plt
import glob

path='/reordered_states/'

# Create the colormap
#colors = [(0,(1,1,1)), (0.30,(0,0.2,0.8)),(0.4,(0,0.3,0.8)), (0.90,(0,0,0.8)), (1,(0,0,0.8))]  
colors = [(0,(1,1,1)),  (1,(0,0,0.8))]
n_bins = 250  # Discretizes the interpolation into bins
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins) #colormap

# ### Fig2 c - ChromHMM Emission plot

emission=pd.read_table(path+'emissions_8.txt', index_col=0)
emission=emission[['H3K4me3', 'H3K27Ac', 'H3K4me1', 'H3K36me3', 'H3K27me3']]

emission=emission.reindex([7, 8,  6, 5, 4, 3, 1, 2])

emission=emission.rename(index=str, columns={'8':'FlnkTSS', '7':'ActTSS', '6':'FlkActEnh', '5':'ActEnh', '4':'WkEnh', '3':'Elong', '1':'Repress', '2': 'Quiesc'})
emission=emission.rename({'8':'FlnkTSS', '7':'ActTSS', '6':'FlkActEnh', '5':'ActEnh', '4':'WkEnh', '3':'Elong', '1':'Repress', '2': 'Quiesc'}, axis='index')

# +
sns.set_style("white")
fig, ax = plt.subplots() 
ax=sns.heatmap(emission, annot=False, cmap=cm, square=True, cbar=False)

#Add annotation
from matplotlib.patches import Rectangle
rect7 = Rectangle( (-0.6,0.1), 0.4, 0.8, facecolor = '#E80808', clip_on=False)
rect6 = Rectangle( (-0.6,1.1), 0.4, 0.8, facecolor = '#FF3333', clip_on=False)
rect5 = Rectangle( (-0.6,2.1), 0.4, 0.8, facecolor = '#FFA143', clip_on=False)
rect4 = Rectangle( (-0.6,3.1), 0.4, 0.8, facecolor = '#FFD717', clip_on=False)
rect3 = Rectangle( (-0.6,4.1), 0.4, 0.8, facecolor = '#ffff00', clip_on=False)
rect2 = Rectangle( (-0.6,5.1), 0.4, 0.8, facecolor = '#008000', clip_on=False)
rect1 = Rectangle( (-0.6,6.1), 0.4, 0.8, facecolor = '#808080', clip_on=False)
rect0 = Rectangle( (-0.6,7.1), 0.4, 0.8, facecolor = '#E5E5E5', clip_on=False)

ax.add_patch(rect7)
ax.add_patch(rect6)
ax.add_patch(rect5)
ax.add_patch(rect4)
ax.add_patch(rect3)
ax.add_patch(rect2)
ax.add_patch(rect1)
ax.add_patch(rect0)

#Rotate axis
loc, labels = plt.xticks()
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.tick_params(axis='y', pad=17,  labelsize=12) #increase the size of the tick and the distance from the axis

ax=sns.despine();

ax=plt.gcf().subplots_adjust(bottom=0.2)
# -

# ### FigS2 b - Transition plot

transition=pd.read_table(path+'transitions_8.txt', index_col=0)

transition=transition.reindex([7, 8,  6, 5, 4, 3, 1, 2])
transition=transition[['7','8',  '6', '5', '4', '3', '1', '2']]

transition=transition.rename(index=str, columns={'8':'FlnkTSS', '7':'ActTSS', '6':'FlkActEnh', '5':'ActEnh', '4':'WkEnh', '3':'Elong', '1':'Repress', '2': 'Quiesc'})
transition=transition.rename({'8':'FlnkTSS', '7':'ActTSS', '6':'FlkActEnh', '5':'ActEnh', '4':'WkEnh', '3':'Elong', '1':'Repress', '2': 'Quiesc'}, axis='index')

# Create the colormap
#colors_t = [(0,(1,1,1)), (0.50,(0,0.2,0.8)),(0.8,(0,0.3,0.8)), (0.95,(0,0,0.8)), (1,(0,0,0.8))]  
colors_t = [(0,(1,1,1)), (0.02,(1,1,1)),(1,(0,0,0.8))]
n_bins = 250  # Discretizes the interpolation into bins
cmap_name = 'my_list'
cm_t = LinearSegmentedColormap.from_list(cmap_name, colors_t, N=n_bins) #colormap

# +
sns.set_style("white")
ax2 = fig.add_subplot(1,4,2) 
ax2=sns.heatmap(transition, annot=False, cmap=cm_t, square=True, cbar=False)

#Add annotation
from matplotlib.patches import Rectangle
rect7 = Rectangle( (-0.6,0.1), 0.4, 0.8, facecolor = '#E80808', clip_on=False)
rect6 = Rectangle( (-0.6,1.1), 0.4, 0.8, facecolor = '#FF3333', clip_on=False)
rect5 = Rectangle( (-0.6,2.1), 0.4, 0.8, facecolor = '#FFA143', clip_on=False)
rect4 = Rectangle( (-0.6,3.1), 0.4, 0.8, facecolor = '#FFD717', clip_on=False)
rect3 = Rectangle( (-0.6,4.1), 0.4, 0.8, facecolor = '#ffff00', clip_on=False)
rect2 = Rectangle( (-0.6,5.1), 0.4, 0.8, facecolor = '#008000', clip_on=False)
rect1 = Rectangle( (-0.6,6.1), 0.4, 0.8, facecolor = '#808080', clip_on=False)
rect0 = Rectangle( (-0.6,7.1), 0.4, 0.8, facecolor = '#E5E5E5', clip_on=False)

ax2.add_patch(rect7)
ax2.add_patch(rect6)
ax2.add_patch(rect5)
ax2.add_patch(rect4)
ax2.add_patch(rect3)
ax2.add_patch(rect2)
ax2.add_patch(rect1)
ax2.add_patch(rect0)

#Borders
ax2.axhline(y=0, color='k',linewidth=1)
ax2.axhline(y=transition.shape[1], color='k',linewidth=2)
ax2.axvline(x=0, color='k',linewidth=1)
ax2.axvline(x=transition.shape[0], color='k',linewidth=2)

#Tick
loc, labels = plt.xticks()
ax2.set_xticklabels(labels, rotation=45, ha='right')
ax2.tick_params(axis='y', pad=18,  labelsize=12)

fig.savefig(path+'transition_8_ordered.pdf', dpi=200)
# -

# ### FigS2 c - PieChart for each state collapsed

path_chromHMM='/chromHMM_enhancer/percentage_states/K_org/'
path_save='/chromHMM_states/'

names=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE' ]

# +
list_=[]
for name in names:
    df = pd.read_csv(path_chromHMM+name+ '_chromHMM_collapsed.bed', header=None,sep='\t', 
                 names=['chr', 'start', 'end', 'states'])
    df1=df.groupby('states').size()
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)

frame = pd.concat(list_, axis=1)
# -

frame.columns=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE']
frame=frame.rename({'promoter':'1_promoter', 'actEnh': '2_actEnh', 'wkEnh':'3_wkEnh', '3_Elong':'4_Elong','2_Quies':'6_Quiesc', '1_Repr': '5_Repr' }, axis='index')

df_pie=frame.mean(axis=1)

# +
explode = (0.05, 0.05, 0.05, 0.05, 0.05, 0.05)
colors = ['#808080','#E5E5E5', '#008000','#ffa500','#ff0000','#ffff00',   ]

fig1, ax = plt.subplots()
df_pie.plot.pie(subplots=True,figsize=(8, 3), explode=explode,  autopct='%1.1f%%', startangle=180, shadow=True,
                            colors=colors, ax=ax)
plt.show()
# -

# ### FigS2 d - Stacked plot chromHMM, check the ration of each state in each PDO

path_chromHMM='/chromHMM_states/'
path_save='/chromHMM_states/'

names=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE' ]

# +
list_=[]
for name in names:
    df = pd.read_csv(path_chromHMM+name+'/'+name+'_dense.bed.gz', compression='gzip', skiprows=[0], header=None,sep='\t', 
                 names=['chr', 'start', 'end', 'states', 'unkn', 'unkn1', 'start1', 'end1', 'color'])
    df1=df.groupby('states').size()
    list_.append(df1)
    df2 = pd.concat(list_, axis=1)

frame = pd.concat(list_, axis=1)
# -

frame=frame.rename({'7_ActTSS':'1_ActTSS', '8_FlkTSS':'2_FlkTSS', '6_FlkActEnh':'3_FlkActEnh', '5_ActEnh':'4_ActEnh', '4_WkEnh':'5_WkEnh', '3_Elong':'6_Elong','2_Quies':'8_Quies', '1_Repr': '7_Repr' }, axis='index')
frame.columns=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE']

## Calculate percentage
list_=[]
for i in range(10):
    df_per=frame.iloc[:,i].div(frame.iloc[:,i].sum())
    list_.append(df_per)
    df_per2 = pd.concat(list_, axis=1)
df_per2.columns=['04_KE','08_KE', '10_KE','11_KW', '13_KE','18_KE', '22_KE','24_KE','36_KE','41_KE']

df_per2.sort_index(ascending=False).T.plot(kind='bar', stacked=True, color=[ '#E5E5E5','#808080','#008000', '#ffff00', '#FFD717','#FFA143','#FF3333','#E80808'])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
