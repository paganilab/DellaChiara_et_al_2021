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
#     display_name: gseapy_0.9.9
#     language: python
#     name: gseapy_0.9.9
# ---

# ### GSEA of genes expressed in primary tumors but not in PDOs (Fig. 1e, Venn diagram, n = 3,412) against Isella et al. stromal signature.
#
# Author: Federica Gervasoni

# ### Library

# %matplotlib inline
# %config InlineBackend.figure_format='retina' # mac
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# #### GSEA on the genes different in Ktissue vs PDO: ordered according to -log10(pavlue) x sign FC - Isella et al signatures

path='/DESeq/All.matched.Ksamples.ChIP.Dgroup_1organoid/clean_KcrKorg_genes/'

rnk = pd.read_table(path+'All.matched.Ksamples.ChIP.Dgroup_1organoid.res.all.pval.sign.FC.rnk', header=0)

## Figure 1E
# Run prerank GSEA
# enrichr libraries are supported by prerank module. Just provide the name
pre_res = gp.prerank(rnk=rnk, gene_sets=('stromal_signatures_isella.gmt'),
                     processes=4,
                     permutation_num=1000, # reduce number to speed up test
                     outdir=(path+"GSEApy_3412.diffGenes_Kcr_immune_related_1organoid_log2FC_isella_ifom"), format='pdf', graph_num=20,
                     weighted_score_type=1, figsize=[3,3])

pre_res.res2d
