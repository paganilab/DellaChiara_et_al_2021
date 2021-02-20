# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jupytext//py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Scanpy 1.4.2
#     language: python
#     name: scanpy142
# ---

# ### Single cell RNA-seq of whole CRC tumor tissue from 23 patients - only Epithelial cells - Lee, H. O. et al. 
#
# Author: Federica Gervasoni

# +
### IMPORT LIBRARIES ###

# +
# Key packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import matplotlib as mpl
import scanpy as sc
import seaborn as sns
# import scvelo as scv # for Velocyto
import hashlib
import sqlalchemy as db
from scipy import stats

# Other packages
import csv as cs
from matplotlib.colors import LinearSegmentedColormap # For color palette
import os
from functools import reduce
from scipy.sparse import csr_matrix
import scutils as scu
import re
import itertools
from urllib.parse import quote_plus as urlquote
import networkx as nx
#import matplotlib.colors as clrs

# For jupyter chunk size
from IPython.core.display import display, HTML

# +
### SETTINGS ###
# -

# Scanpy settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi = 200)  # low dpi (dots per inch) yields small inline figures
# sc.logging.print_versions()

# ### Palette

# Colormap definition called cm_2 (gray background with red dots for Umap)
colors = [(0,('#d6d6d6')), (0.2,('#F1F1F1')), (0.7,('#991199')), (1, ('#6e086e'))]  
n_bins = 250  # Discretizes the interpolation into bins
cmap_name = 'new_list'
cm_2 = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bins)

# +
# Set min-center-max colors
colors = [(0,('#243302')), (0.1,('#3E5902')), (0.4,('#838C3B')), (0.5,('#ffffff')), (0.70,('#AB66D9')), (0.8, ('#9337de')), (0.9,('#6F1BB3')),(1,('#4e0e82'))]  
# Discretizes the interpolation into bins
n_bins = 250
# Set Colormap Name
cmap_name = 'diff_list'

# Build colormap
olive2Mel  = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
# -
blue2red=clr.LinearSegmentedColormap.from_list('custom blue', ['#302e82','#110AEA','white',"red","red", "#b20000"], N=256)

# Set chunk size of jupyter notebook
display(HTML("<style>.container { width:90% !important; }</style>"))

# +
### Output

# +
samplesName = ['GSE132465']
SAMPLE= 'GSE132465'

tissue_type='Epithelila_Cells'

raw_mtx='/public_availabe_scRNAseq/crc/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt'

outdir = os.path.join('/public_availabe_scRNAseq/crc/', SAMPLE, tissue_type, hashlib.sha256(
   ','.join(samplesName).encode()).hexdigest())
if not os.path.exists(outdir):
   os.makedirs(outdir)

raw_data_h5ad = os.path.join('/public_availabe_scRNAseq/crc/', SAMPLE, hashlib.sha256(
   ','.join(samplesName).encode()).hexdigest(), 'GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.h5ad') 

results = os.path.join('/public_availabe_scRNAseq/crc/', SAMPLE, tissue_type, hashlib.sha256(
   ','.join(samplesName).encode()).hexdigest(), 'scanpy_analysis_final.h5ad') 

#Set wd
os.chdir(outdir)

print("The result of this analysis are in {}".format(outdir+'/figures'))
print("The h5ad file is located at {}".format(results))
print("The wd is located ad {}".format(outdir))

sc.settings.figdir = outdir + '/' + tissue_type +'/figures'
if not os.path.exists(outdir + '/figures'):
   os.makedirs(outdir + '/figures')
# -

# #### Import data from mtx and annotation (GEO repo) and save as h5ad

# +
# ##read dense matrix, trasform in a sparse matrix and save it
# #before the mtx loading remove the first row and column: 
# #zcat GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz | tail -n +2 | cut -f 2- > GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt

# _X= np.loadtxt(raw_mtx, delimiter='\t')
# X = scipy.sparse.csr_matrix(_X)
# tX=scipy.sparse.csr_matrix.transpose(X)
# del(_X)
# # scipy.io.mmwrite(outdir+'/GSE132465_raw_count.mtx', tX)

# +
# ### import obs and vars
# df_obs=pd.read_csv('/public_availabe_scRNAseq/crc/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz',
#                   sep='\t', index_col=0)
# df_obs['batch']=SAMPLE
# df_vars=pd.read_csv('/public_availabe_scRNAseq/crc/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_count_var.txt',
#                   sep='\t', index_col=0)

# +
# ### create AnnData and save h5ad
# data=anndata.AnnData(X=tX, obs=df_obs, var=df_vars)
# data.write(raw_data_h5ad)

# +
### LOAD H5AD POST ###
# -

# put load == True to load POST analyisis h5ad after running all the chunks above of this
# If load is True you must skip "LOAD AND AGGREGATE SAMPLES" and "PRE ANALYSIS SETTING" below and 
# pay attention to run chunk carefully because data could be overwritten as for example images.
load = False
if load == True:
    # Check results path
    print(results)
    # Reload data
    data = sc.read_h5ad(filename=f"{results}")
    # Restore VISION
    vision = data.uns['VISION']
    # Create shortcuts for access VISION data
    goi = vision['genes_of_interests'] # store different set of genes
    gseas = vision['gseas'] # store Gsea genesets and metadata
    cluster2color = vision['cluster2color'] # store various color combination for clusters
    res_coi = vision['res_clusters_of_interests'] # store data about resolution clusters of interest
    annotations = vision['annotations'] # store various annotations
    ann_coi = vision['ann_clusters_of_interests'] # store data about annotation clusters of interest
    savedImg = vision['saved_images'] # store paths for all saved images
    params = vision['parameters'] # store parameters for all the pipeline
else:
    print('Not load POST data')

# +
### PRE ANALYSIS SETTINGS ###
# -

# Create unstructured data for VISION containing various information about the analysis
vision = {'name': SAMPLE,
         'sha': hashlib.sha256(','.join(SAMPLE).encode()).hexdigest(),
         'user_path': outdir,
         'storage_path': f"/{hashlib.sha256(','.join(SAMPLE).encode()).hexdigest()}",
         'ID': list(SAMPLE),
         'mean_counts_per_cell': '',
         'mean_genes_per_cell': '',
         'HVG': '',
         'genes_of_interests': {},
         'resolution': '',
         'res_differential_genes': '',
         'res_unique_markers': '',
         'gseas': {},
         'res_clusters_of_interests': {},
         'cluster2color': {},
         'annotations': {},
         'annotation': '',
         'ann_differential_genes': '',
         'ann_unique_markers': '',
         'ann_clusters_of_interests': {},
         'saved_images': {},
         'parameters': {}
}

# Create shortcuts for access VISION data
goi = vision['genes_of_interests'] # store different set of genes
gseas = vision['gseas'] # store Gsea genesets and metadata
cluster2color = vision['cluster2color'] # store various color combination for clusters
res_coi = vision['res_clusters_of_interests'] # store data about resolution clusters of interest
annotations = vision['annotations'] # store various annotations
ann_coi = vision['ann_clusters_of_interests'] # store data about annotation clusters of interest
savedImg = vision['saved_images'] # store paths for all saved images
params = vision['parameters'] # store parameters for all the pipeline

# ### Upload data

data=sc.read(raw_data_h5ad)

data = data[data.obs['Cell_type'] == 'Epithelial cells',:]

### Visualize the cells were already filtered
data

SamplesInfo=data.obs
SamplesInfo['id']='GSE132465'
SamplesInfo.head()

# +
### QC ###

# +
import math

min_genes = 200
mc = math.ceil(data.shape[0]*0.001)
print("Minimum genes: {}".format(min_genes))
print("Minimum cells: {}".format(mc))
# Filtering Genes and Barcodes
print("Before \n{}".format(data))
sc.pp.filter_cells(data, min_genes=min_genes)
sc.pp.filter_genes(data, min_cells=mc)
# -

# Plot data distribution for each sample
scu.plot_distributions(data)

# Plot again violin plots for mito-ribo percents and number of genes-counts. You expect your data to be already filtered from the PRE analysis. 
# In other words the dots must be cut at the thresholds level set in PRE analysis. For this reason there is no need to filter again.
thresholds = [5,2,2,3]
scu.plot_metrics_mad(data, 
                     thresholds = thresholds, 
                     details=True, 
                     samples=None)

# Store in Vision and visualize the mean number of counts per cell
vision.update({'mean_counts_per_cell': data.obs['n_counts'].sum()/len(data.obs)})
print(f"The mean number of counts per cell is {vision['mean_counts_per_cell']}")

# Store in Vision and visualize the mean number of genes per cell
vision.update({'mean_genes_per_cell': data.obs['n_genes'].sum()/len(data.obs)})
print(f"The mean number of genes per cell is {vision['mean_genes_per_cell']}")

# +
# Normalize and log data on 10000 counts per cell and save raw data (10000 should change basing on mean number of counts per cell).

# Store parameter for VISION
params.update({'sc.pp.normalize_per_cell': {'counts_per_cell_after': 1e4}})

sc.pp.normalize_per_cell(data, 
                         counts_per_cell_after = params['sc.pp.normalize_per_cell']['counts_per_cell_after'])
sc.pp.log1p(data)
data.raw = data


# +
### ANALYZE DATASET ###

# +
# Define and plot HVG, run mutliple times to select the right number of HVG

# Store parameters for VISION
params.update({'sc.pp.highly_variable_genes': {'min_mean': 0.05,
                                               'max_mean': 4,
                                               'min_disp': 0.7}})

if len(set(SamplesInfo.id)) == 1:
    print('No batch analysis')
    sc.pp.highly_variable_genes(data, 
                                min_mean = params['sc.pp.highly_variable_genes']['min_mean'], 
                                max_mean = params['sc.pp.highly_variable_genes']['max_mean'],
                                min_disp = params['sc.pp.highly_variable_genes']['min_disp'])
else:
    print('Batch analysis')
    sc.pp.highly_variable_genes(data, 
                                batch_key = 'batch', 
                                min_mean = params['sc.pp.highly_variable_genes']['min_mean'],
                                max_mean = params['sc.pp.highly_variable_genes']['max_mean'], 
                                min_disp = params['sc.pp.highly_variable_genes']['min_disp'])
    data.var['highly_variable'] = data.var['highly_variable_intersection']

scu.remove_mito_ribo_from_hvg(data)
sc.pl.highly_variable_genes(data)
print("Selected {} HVG".format(np.sum(data.var['highly_variable'])))

# Add HVG number to VISION
vision.update({'HVG': np.sum(data.var['highly_variable'])})
# -

# Compute cell cycle data
data = scu.define_cell_cycle(data)

# +
# Scale data, compute PCA and visualize it for mito-ribo percents, number of genes-counts and cell cycle

# Store parameters for VISION
params.update({'sc.pp.scale': {'max_value': 10}})
params.update({'sc.tl.pca': {'svd_solver': 'arpack'}})
params.update({'sc.pl.pca_scatter': {'color': ['percent_mito','percent_ribo','n_genes','n_counts','cc_differences', 'S_score', 'G2M_score']}})

sc.pp.scale(data, 
            max_value = params['sc.pp.scale']['max_value'])
sc.tl.pca(data, 
          svd_solver = params['sc.tl.pca']['svd_solver'])
sc.pl.pca_scatter(data, 
                  color = params['sc.pl.pca_scatter']['color'],
                  color_map = 'Purples', 
                  ncols = 4)

# +
# Regressing out for the needed components

# Store parameters for VISION
params.update({'sc.pp.regress_out': {'keys': ['percent_ribo', 'percent_mito','n_counts','n_genes','cc_differences', 'S_score', 'G2M_score']}})

sc.pp.regress_out(data,
                  keys = params['sc.pp.regress_out']['keys'],
                  n_jobs = 10)
# -

# Scale data and perform again PCA to check for changes after regress-out.
sc.pp.scale(data, 
            max_value = params['sc.pp.scale']['max_value'])
sc.tl.pca(data, 
          svd_solver = params['sc.tl.pca']['svd_solver'])
sc.pl.pca_scatter(data, 
                  color = params['sc.pl.pca_scatter']['color'], 
                  color_map = 'Purples',
                  ncols = 4)

# Perform batch effect correction if a batch is present and compute again scaling and PCA to check for changes.
if len(set(SamplesInfo.id)) == 1:
    print("Batch correction NOT performed")
else:
    print("Batch correction")
    sc.pl.pca_scatter(data, 
                      color = ['Sample'], 
                      color_map = 'Purples', 
                      ncols = 1)
    sc.pp.combat(data, 
                 key = 'Patient')
    sc.pp.scale(data, 
                max_value = params['sc.pp.scale']['max_value'])
    sc.tl.pca(data, 
              svd_solver = params['sc.tl.pca']['svd_solver'])
    sc.pl.pca_scatter(data, 
                      color = ['Sample'],
                      color_map = 'Purples', 
                      ncols = 1)

# Plot variance ratio to decide how many components to use for neighbors graph computation
sc.pl.pca_variance_ratio(data, 
                         log = True, 
                         show = False)
sc.pl.pca_variance_ratio(data, 
                         show = False)

# +
# Compute the nearest neighbours network

# Store parameters for VISION
params.update({'sc.pp.neighbors': {'n_neighbors': 15,
                                   'n_pcs': 20}})
sc.pp.neighbors(data, 
                n_neighbors = params['sc.pp.neighbors']['n_neighbors'], 
                n_pcs = params['sc.pp.neighbors']['n_pcs'])
# -

# Compute tsne based on the previous network
sc.tl.tsne(data)

# +
categs=['Cell_type', 'Cell_subtype','Class', 'Sample', 'Patient']

for categ in categs:
    sc.pl.tsne(data, color=categ, frameon=False)
# -

# Plot UMAP
if len(set(SamplesInfo.id)) == 1:
    sc.pl.tsne(data, 
               color = ['n_genes'])
else:
    sc.pl.tsne(data, 
               color = ['n_genes'])
    sc.pl.tsne(data, 
               color=['batch'], 
               title="Batch Correction")

# Plot UMAP for each sample (if batch exists)
if len(set(SamplesInfo.id)) == 1:
    print('No batch available')
else:
    print('Batch available')
    ncol = data.obs.batch.cat.categories.values.size
    fig,ax = plt.subplots(1, ncol, sharex = True, sharey = False, figsize = (4*ncol,4))
    c = 0
    for g, group in data.obs.groupby("batch"):
        idx = group.index
        sc.pl.umap(data[idx,], show=False, ax=ax[c], size=10)
        ax[c].set_title(g)
        c+=1
    plt.tight_layout()
    plt.show()

# +
# Plot useful markers on tsne

# Store parameters for VISION
# It is possible to create any amount of custom genes lists to plot on UMAP.
# It is possible to change the 'list_X' name withh your own custom name.
goi.update({'technical': ['percent_ribo', 'percent_mito','n_counts','n_genes', 'phase'],
            'epithelial_cells': ['KRT18', 'KRT19', 'EPCAM', 'CDH1']})

for i in goi.keys(): 
    print(i)
    sc.pl.tsne(data,
               color = goi[i],
               ncols = 4, 
               color_map = cm_2)
# -

# Compute and plot clusters-batch composition and cluster cells composition
if len(set(SamplesInfo.id)) == 1:
    scu.cluster_origin(data, 'phase', 
                       index = 'Sample', 
                       type = 'percentage', save= './figures/patients_phase.pdf')
else:
    scu.cluster_origin(data, 'batch', 
                       index = vision['resolution'] )
    scu.cluster_origin(data, 'batch', 
                       index = vision['resolution'], 
                       type = 'raw')

# ### Metagene

# +
# Geneset to obtain stemness and crb score
gset_path='./'
gsets=scu.read_gmt(gset_path+'dellaChiara_etAl.gmt')

#add lists in tagged_geneset 
for glist_name in gsets.keys():
    scu.tag_geneset(data, glist_name, gsets[glist_name])
# -

# Stem cells genes
meta_df=scu.metagene_score(data, scu.get_geneset(data.raw, 'stemness_score'), use_raw=True)
data.obs['stemness_score_scaled']=meta_df['combined_score']
data.obs['stemness_score_scaled']=data.obs['stemness_score_scaled']/data.obs['stemness_score_scaled'].max()

## add combined score to obs
meta_df=scu.metagene_score(data, scu.get_geneset(data.raw, 'CRB_score'), use_raw=True)
data.obs['CRB_score_scaled']=meta_df['combined_score']
data.obs['CRB_score_scaled']=data.obs['CRB_score_scaled']/data.obs['CRB_score_scaled'].max()

sc.pl.tsne(data, color='Class', frameon=True)
sc.pl.tsne(data, color='Patient', frameon=True)
sc.pl.tsne(data, color='CRB_score_scaled', cmap=blue2red, frameon=True)
sc.pl.tsne(data, color='stemness_score_scaled', cmap=blue2red, frameon=True)

# +
treshold=0.3
tmp =  (np.where(data.obs['CRB_score_scaled'] > treshold, data.obs['CRB_score_scaled'], 0)) - np.where(data.obs['stemness_score_scaled'] > treshold, data.obs['stemness_score_scaled'], 0)

data.obs['diff_crb_stemness_score']=tmp/tmp.max()

fig, ax = plt.subplots(1,1, frameon=False)

sc=ax.scatter(data.obsm['X_tsne'][:,0], data.obsm['X_tsne'][:,1], 
            c=tmp/tmp.max(), cmap=olive2Mel, 
            edgecolors='black',linewidth=0.1,
            s=40*data.obs['CRB_score_scaled'], 
            alpha=0.7)

ax.grid(b=False)
ax.axis('off')

#Add colorbar
plt.colorbar(sc, ax=ax)

# Add Arrows
ax.annotate("", xy=(0, 0.2), xytext=(0, 0), xycoords='axes fraction', arrowprops=dict(arrowstyle="->"))
ax.annotate("", xy=(0.2, 0), xytext=(0, 0), xycoords="axes fraction", arrowprops=dict(arrowstyle="->"))

plt.show()
# -

ann_coi.keys()

# Save VISION into data object inside unstructured data
data.uns['VISION'] = vision

# visualize VISION into data (unstructured)
data.uns['VISION']

# +
### SAVE H5AD POST ### 
# -

# Check results path
print(results)

# Save the POST result h5ad file.
data.write_h5ad(results)
print(f"Results saved in: {results}")
