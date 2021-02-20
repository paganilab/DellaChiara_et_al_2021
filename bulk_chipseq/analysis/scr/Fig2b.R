# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: R34_ChIPseq_v1
#     language: R
#     name: r34_chipseq_v1
# ---

# ### Correlation plot of all the samples and histone modifications (H3K4me3, H3K27ac, H3K4me1, H3K36me3, H3K27me3)
#
# Author: Federica Gervasoni

source("https://bioconductor.org/biocLite.R")
library(DiffBind)
library(gplots)
library(GenomicFeatures)
library(limma)
library(knitr)
library(data.table)
library(BiocParallel)
library(tidyr)
library(stringr)
library(gdata)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(genefilter)
library(magrittr)

# ### All ChIP-seq samples

samples_table <- read.csv(file = "Diffbind_sample_table_peaks.csv", header = T)

dba <- dba(sampleSheet=samples_table,
         minOverlap = 2, 
         config = data.frame(RunParallel = TRUE,
                             reportInit = "DBA",
                             DataType = DBA_DATA_FRAME))
dba

dba=dba.count(dba, score=DBA_SCORE_TMM_MINUS_FULL,  bParallel=T)

corr=dba.plotHeatmap(dba, colScheme = "YlGnBu", attributes = c("Factor"), ColAttributes=("ID"))

# ### Fig2 b - Correlation plot

matrix_corr = as.matrix(corr)
breaks = seq(0,max(matrix_corr),length.out=1000)
gradient1 = colorpanel( sum( breaks[-1]<=0.1 ), "white", "#FFD1F0" )
gradient2 = colorpanel( sum( breaks[-1]>0.1 ), "#FFD1F0", "#990066" )
hm.colors = c(gradient1,gradient2)

corr=dba.plotHeatmap(dba, colScheme = hm.colors, attributes = c("Factor"), ColAttributes=c("Factor","Condition"), colSideCols=list(c("#999999","#4daf4a","#e41a1c","#377eb8","#FF62BC","#984ea3","#ff7f00","#ffff33","#a65628", "#65B7F3"), c('#E60000', '#F31558', '#F3B619', '#009900', '#666666')))


# ### FigS2 a - PCA plot

# pdf(paste0(path, "PCA_plot_allHMs.pdf"), height=5, width=5)
dba.plotPCA(dba, attributes=DBA_CONDITION, vColors=c("#E60000", "#F31558", "#F3B619", "#009900", "#666666"))
# dev.off()

# ### Save and Load data

# +
## All HMs
name='correlation_plot_peaks'

dba.save(dba, dir=path, file=name, ext='RData', pre='dba_')
# -

dba=dba.load(file=name, path, pre='dba_', ext='RData')
