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

# ### DE analysis between Tumor and Normal tissue - matched patients with ChIP-seq
#
# Author: Federica Gervasoni
#
# design = ~ group + patient

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
library(ggrepel)
library(VennDiagram)
library(genefilter)

# ### Read sampleinfo

## set universal plot size:
options(repr.plot.width=6, repr.plot.height=6)

#Rds with all the clinical information about the patients (please see Supplementary Table 1)
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

counts.noMT=counts[!rownames(counts)%in%MT.genes,]
rownames(counts.noMT)=gsub("[.].*","",rownames(counts.noMT))

length(unique(rownames(counts)))
length(unique(rownames(counts.noMT))) 
identical(rownames(sampleinfo), colnames(counts.noMT))
# -

# ## Analysis
#
# ### Set comparison

path="/DESeq/DGE/KvsNcr/KvsNcr.matched.ChIP.Dgroup.patient/"

# +
name="KvsNcr.matched.ChIP.Dgroup.patient"

cond=c("K_tissue","N_tissue")

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

samples=droplevels(sampleinfo[sampleinfo$condition2%in%cond & sampleinfo$patient%in%matched,])
samples$group=factor(samples$group, levels=c("CRC_N", "CRC_K"))
samples$condition2=factor(samples$condition2, levels=c("N_tissue", "K_tissue"))

file.dds = paste0(path, "org_dds.", name, ".noMT.rds")
file.rld = paste0(path, "org_rld.", name, ".noMT.rds")

contrast=c("group", "CRC_K", "CRC_N")

table(samples$group)
table(samples$condition2)
table(samples$tissue)
table(samples$condition2, samples$patient)
# -

# ### Create dds and rld objects
#
# With moderation of log2 fold changes. betaPrior=TRUE by default in DESeq2 versions < 1.6

# +
if (!file.exists(file.dds)){
    
    dds <- DESeqDataSetFromMatrix(countData = counts.noMT[,colnames(counts.noMT)%in%samples$sample],
                              colData = samples,
                              design = ~ patient +group)

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
sum(res$padj < 0.05, na.rm=TRUE)
# -

# ### Metadata

# +
metadata(res)$alpha

# Percentage of genes set to NA
# BaseMean threshold of genes set to NA 
metadata(res)$filterThreshold
# -

# ### Extract results

# +
# Extract, annotate and save results for all genes 
res2=as.data.frame(res)
res2$gene_id=rownames(res2)
res2$padjxFC=round(-log10(res2$padj)*sign(res2$log2FoldChange), digit=3)
res2=merge(res2, annotation, by="gene_id")
write.table(res2, paste0(path, name, ".res.all.txt"), row.names = F, col.names = T, quote=F, sep="\t")

## Remove all records with padj=NA and with missing or duplicated gene_names
### This files are for GSEA
omit.ids=res2[duplicated(res2$gene_name) | res2$gene_name=="" | !complete.cases(res2[,"padj"]), ]
omit.ids_list=omit.ids$gene_id

write.table(res2[!res2$gene_id%in%omit.ids_list, c("gene_name","log2FoldChange")], paste0(path, name, ".res.all.log2FC.rnk"), row.names = F, col.names = T, quote=F, sep="\t")

write.table(res2[!res2$gene_id%in%omit.ids_list,c("gene_name","padjxFC")], paste0(path, name, ".res.all.pval.sign.FC.rnk"), row.names = F, col.names = T, quote=F, sep="\t")
# -

# ### DGEgenes logFC > 1 & p-value < 0.05

# +
DGEgenes <- rownames(subset(res[order(res$padj),], padj < 0.05  & abs(log2FoldChange)>1))
table(DGEgenes%in%resSig$gene_id)

# Upregulated
DGEgenes.up <- rownames(subset(res[order(res$padj),], padj < 0.05  & log2FoldChange>1))
table(DGEgenes.up%in%DGEgenes)

# Downregulated
DGEgenes.down <- rownames(subset(res[order(res$padj),], padj < 0.05  & log2FoldChange<(-1)))
table(DGEgenes.down%in%DGEgenes)

# +
dds.all=readRDS(paste0(path, "/All.matched.ChIP.Dgroup_1organoid/org_dds.All.matched.ChIP.Dgroup_1organoid.noMT.rds"))

# Modify the gene_id
rownames(dds.all)=gsub("[.].*", "", rownames(dds.all))

norm.all=counts(dds.all, normalized=TRUE)

norm.all.log2=log2(norm.all + 1)

# +
padj_filter=0.01
log2FC=0

DGEgenes <- rownames(subset(res[order(res$padj),], padj <= padj_filter  & abs(log2FoldChange)>=log2FC))
table(DGEgenes%in%resSig$gene_id)

# Upregulated
DGEgenes.up <- rownames(subset(res[order(res$padj),], padj <= padj_filter  & log2FoldChange>=log2FC))
table(DGEgenes.up%in%DGEgenes)

# Downregulated
DGEgenes.down <- rownames(subset(res[order(res$padj),], padj <= padj_filter  & log2FoldChange<=(-log2FC)))
table(DGEgenes.down%in%DGEgenes)
# -

# ### Fig1 f - Heatmap

# +
# Extract common UP and common Down reg genes
DGEgenes_list=c(as.character(DGEgenes.up),as.character(DGEgenes.down))

## Subset the log2 normalised dds.all dataset for the common genes
norm.all.log2.con  <- norm.all.log2[match(DGEgenes_list, row.names(norm.all.log2)), ]

table(row.names(norm.all.log2.con)%in%DGEgenes_list)
table(colnames(norm.all.log2.con)==colData(dds.all)[,"sample"])

### Scale manually
norm.all.log2.con.scale.man=apply(norm.all.log2.con,2, function(x){(x-rowMeans(norm.all.log2.con))/rowSds(norm.all.log2.con)})

### Create the annotation table
anno <- as.data.frame(colData(dds.all)[, c("patient", "condition")])

### reorder columns according to groups
anno=anno[order(anno$condition, decreasing = TRUE),]
anno_columns=row.names(anno)

### Assign specific colors
ann_colors = list(condition = c(K_organoids="#9b1477", K_tissue="#590242", N_tissue="#6A6B06"), 
                  patient=c(CRC10="#e41a1c", CRC11='#377eb8', CRC13='#FF62BC', CRC18='#984ea3', CRC22='#ff7f00', 
                            CRC24='#ffff33', CRC36='#a65628',CRC4='#999999', CRC41="#65B7F3", CRC8='#4daf4a' ))

norm.all.log2.con.scale.man=norm.all.log2.con.scale.man[,anno_columns]

pdf=pheatmap(norm.all.log2.con.scale.man, scale="none",
         cluster_cols =TRUE,
         cluster_rows=F,         
         annotation_col = anno,
         annotation_colors = ann_colors,
         show_rownames=FALSE,
         color = colorRampPalette(rev(c("#8A0807", "#D51214", "#FFFEFD", "#3E44CC", "#0A2A87")))(1250), 
         breaks = c(seq(-5,-2.60,length=250), seq(-1,0.5,length=250),seq(0.6,0.8,length=250),seq(0.9,1.7,length=250), seq(3.1,5,length=250)),
         fontsize=8)

grid::grid.newpage()
grid::grid.draw(pdf)
