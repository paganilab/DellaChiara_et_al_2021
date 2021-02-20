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
#     display_name: gseapy_0.9.9
#     language: python
#     name: gseapy_0.9.9
# ---

# ### GSEA of DE genes in PDO vs Ntissue against GRADE genesets
#
# Author: Federica Gervasoni

# ### Library to import

# %matplotlib inline
# %config InlineBackend.figure_format='retina' # mac
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

path_gset=''

path='/DESeq/DGE/KorgvsNcr/KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient/'

rnk = pd.read_table(path+'KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient.res.all.pval.sign.FC.rnk', header=None)

# Run prerank GSEA
# enrichr libraries are supported by prerank module. Just provide the name
pre_res = gp.prerank(rnk=rnk, gene_sets=(path_gset+'GRADE_genesets.gmt'),
                     processes=4,
                     permutation_num=1000, # reduce number to speed up test
                     outdir=(path+"GSEApy_KorgvsNcr"), format='pdf', graph_num=20,
                     weighted_score_type=1, figsize=[3,3])

# ### Fig1g - GSEA 

# Gene signature reported in the paper:
#
# GRADE_COLON_AND_RECTAL_CANCER_UP and GRADE_COLON_CANCER_DN_GRADE_COLON_AND_RECTAL_CANCER_DN

from gseapy.plot import gseaplot
gseaplot(rank_metric=pre_res.ranking, term=terms[4], **pre_res.results[terms[4]], pheno_pos= " 'na_pos' ",pheno_neg=" 'na_neg' ")
gseaplot(rank_metric=pre_res.ranking, term=terms[6], **pre_res.results[terms[6]], pheno_pos= " 'na_pos' ",pheno_neg=" 'na_neg' ")
