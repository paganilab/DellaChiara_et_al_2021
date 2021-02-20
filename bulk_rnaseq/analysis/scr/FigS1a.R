# -*- coding: utf-8 -*-
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
#     display_name: R35_RNAseq
#     language: R
#     name: r35_rnaseq
# ---

# ### Analysis of PDOs samples: DE between middle and late stages - matched patients with ChIP-seq
#
# Author: Federica Gervasoni
#
# design = ~ patient + passage

# ### Library

library(knitr)
library(reshape)
library(data.table)
library(DESeq2)
library(BiocParallel)
library(tidyr)
library(stringr)
library(gdata)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(genefilter)
library(magrittr)
library(dplyr)
library(stats)
library(ggpubr)

## set universal plot size:
options(repr.plot.width=6, repr.plot.height=6)

# Rds with all the clinical information about the patients (please see Supplementary Table 1)
sampleinfo=readRDS("sampleinfo.rds")

# ### Read annotation

# +
# Read annotation (excludes PAR_Y genes)
annotation=readRDS("annotation_file_from_gtf.rds")
length(unique(annotation$gene_id))

# Modify the gene_id
annotation$gene_id=gsub("[.].*", "", annotation$gene_id)
## Ensure that the number of unique gene_ids remains the same
length(unique(annotation$gene_id))

# Extract MT genes
MT.genes=annotation[annotation$chromosome=="chrM", "gene_id"]
# -

# ### Read raw counts and remove MT genes

# +
# Read raw counts
counts=readRDS("org_fcounts_all.rds")
rownames(counts)=gsub("[.].*", "", rownames(counts))

# Remove MT genes (PAR_Y are already removed)
    #counts.noMT=counts[grep("PAR_Y",rownames(counts), invert=TRUE),]
    counts.noMT=counts[!rownames(counts)%in%MT.genes,]
    
length(unique(rownames(counts)))
length(unique(rownames(counts.noMT)))

identical(rownames(sampleinfo), colnames(counts.noMT))
# -

# ## Analysis
#
# ### Set comparison

path="/DESeq/DGE/Korg/Korg.MiddlevsLate/"

# +
name="Korg.MiddlevsLate.Dpassage"

cond=c("K_late","K_middle")

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

samples=droplevels(sampleinfo[sampleinfo$condition2%in%cond & sampleinfo$patient%in%matched,])
samples$group=factor(samples$passage, levels=c("middle", "late"))
samples$condition2=factor(samples$condition2, levels=c("K_middle", "K_late"))

file.dds = paste0(path, "org_dds.", name, ".noMT.rds")
file.rld = paste0(path, "org_rld.", name, ".noMT.rds")

contrast=c("passage", "middle", "late")
# -

# ### Create dds and rld objects

# +
if (!file.exists(file.dds)){
    
    dds <- DESeqDataSetFromMatrix(countData = counts.noMT[,colnames(counts.noMT)%in%samples$sample],
                              colData = samples,
                              design = ~ patient + passage)

    # Removing non-expressed genes
    dds <- dds[rowSums(counts(dds)) >= 10,]
    
    #Create dds and rld
    dds <- DESeq(dds)
    rld <- rlog(dds, blind = FALSE)

    saveRDS(object = dds, file = file.dds)
    saveRDS(object = rld, file = file.rld)
        
} else {

    dds=readRDS(file = file.dds)
    rld=readRDS(file = file.rld)
    
}


dds
dds@design
colData(dds)
table(dds$group)
table(dds$condition2)
table(dds$tissue)
table(dds$condition2, dds$patient)
table(dds$passage, dds$patient)
# -

# ### DE

# +
table(colData(dds)[,"condition2"])
contrast
res <- results(dds, contrast=contrast, alpha=0.05)

# Summary
res[order(res$padj),]
summary(res)
sum(res$padj <= 0.05, na.rm=TRUE)
# -

# ### Metadata

# +
metadata(res)$alpha

# Percentage of genes set to NA
# BaseMean threshold of genes set to NA 
metadata(res)$filterThreshold
# -

# ### FigS1. a - MAplot  

# +
# Fig supp 1B
res2=as.data.frame(res)

# Default plot
ggmaplot(res2, 
   fdr = 0.01, fc = 1, size = 2, alpha=0.6,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   legend = "top", top = 20,
   font.label = c("bold", 11),
   font.legend = "bold",
   font.main = "bold",
   ylim=c(-5,5),
   ggtheme = ggplot2::theme_classic()) + 
   geom_hline(yintercept=1, linetype="dashed", size = 2, color = "dodgerblue") +
   geom_hline(yintercept=-1, linetype="dashed", size = 2, color = "dodgerblue") +
   geom_point(alpha=0.6, size=2) +
   geom_hline(yintercept=0,  color = "red", size = 2, alpha=0.6) 
