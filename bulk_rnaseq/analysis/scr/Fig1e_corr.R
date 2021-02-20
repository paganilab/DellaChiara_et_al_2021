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

# ### Analysis of all tumor RNA-seq samples - matched patients with ChIP-seq
#
# Author: Federica Gervasoni
#
# design = ~ group

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
library(ggcorrplot)

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

#Extract Y genes
Y.genes=annotation[annotation$chromosome=="chrY", "gene_id"]
# -

# ### Read raw counts and remove MT and Y genes

# +
# Read raw counts
counts=readRDS("org_fcounts_all.rds")
rownames(counts)=gsub("[.].*", "", rownames(counts))

counts.noMT=counts[!rownames(counts)%in%MT.genes & !rownames(counts)%in%Y.genes,]
    
length(unique(rownames(counts)))
length(unique(rownames(counts.noMT)))

identical(rownames(sampleinfo), colnames(counts.noMT))
# -

# ## Analysis 
#
# ### Set dds and rld objects

# +
name="All.matched.Ksamples.ChIP.Dgroup_1organoid"

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")
samples$group=factor(samples$tissue, levels=c("tissue", "organoids"))

file.dds = paste0(path, "org_dds.", name, ".noMT.rds")
file.rld = paste0(path, "org_rld.", name, ".noMT.rds")
# -

# ### Remove genes present only in K tissue or K organoids

path_list='/All.matched.Ksamples.ChIP.Dgroup_1organoid/'
intersect_Kcr_anno=read.table(paste0(path_list, "intersect_Kcr_anno.txt"), header=TRUE)
intersect_Korg_anno=read.table(paste0(path_list, "intersect_Korg_anno.txt"), header=TRUE)
intersect_Korg_Kcr_anno=c(as.character(intersect_Kcr_anno$intersect_Kcr), as.character(intersect_Korg_anno$intersect_Korg))

# ### Create dds and rld objects

# +
if (!file.exists(file.dds)){
    
    dds <- DESeqDataSetFromMatrix(countData = counts.noMT[,colnames(counts.noMT)%in%samples$sample],
                              colData = samples,
                              design = ~ group)

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
table(dds$condition)
table(dds$tissue)
table(dds$condition, dds$patient)
# -

## log2 normalized count per correlation matrix.
log2_norm_matx_K=log2(counts(dds,normalized=TRUE)+1)
write.table(log2_norm_matx_K, paste0(path, name, "log2_ncount", ".csv"), sep=',', quote=F, col.names = T, row.names=T)

# #### Correlation

# baseMean according each group
baseMeanPerLvl <- as.data.frame(sapply( levels(dds$group), function(lvl) rowMeans( log2(counts(dds,normalized=TRUE)[,dds$group == lvl] +1) ) ))

# ### Fig1 e - Correlation plot

p=ggplot(data = baseMeanPerLvl, aes(x = tissue, y = organoids)) + 
  geom_point(data=baseMeanPerLvl[as.vector(intersect_Kcr_anno$intersect_Kcr),], aes(x=tissue, y=organoids), colour="#05C7F2", size=0.5) +
  geom_point(data=baseMeanPerLvl[as.vector(intersect_Korg_anno$intersect_Korg),], aes(x=tissue, y=organoids), colour="#8EFFAA", size=0.5) +
  geom_point(alpha=0.1, size=2, colour = "#110AEA") +
  xlab("log2 mean normalized tumor tissue") + ylab("log2 mean normalized PDO") +
  theme_classic() +
  theme(axis.title = element_text(size=10, family="Helvetica"),
        axis.text = element_text(size=10, family="Helvetica"),
        plot.title = element_text(size=10, family="Helvetica"),
        legend.position = "none") +
  geom_smooth(method = "lm", se = TRUE) 
p

cor(baseMeanPerLvl$tissue, baseMeanPerLvl$organoids)

cov(baseMeanPerLvl$tissue, baseMeanPerLvl$organoids)

# summary(lm( baseMeanPerLvl$organoids~baseMeanPerLvl$tissue))
rsq <- function(x, y) summary(lm( baseMeanPerLvl$organoids~baseMeanPerLvl$tissue))$r.squared
rsq(obs, mod) # Rsquared

# ### Sigular correlation plot

pairs=as.list(c('SQ_1957,SQ_1958', 'SQ_1951,SQ_1955', 'SQ_1961,SQ_1964', 'SQ_1968,SQ_2018', 'SQ_1972,SQ_1975', 'SQ_1985,SQ_2016', 'SQ_2020,SQ_2021', 'SQ_2011,SQ_2015', 'SQ_2251,SQ_2253', 'SQ_2263,SQ_2264'))

# baseMean according each group
for(i in 1:length(pairs)) {
    
    tmp <- as.data.frame(sapply( strsplit(pairs[[i]][1], ",", fixed = TRUE)[[1]], function(lvl) log2(counts(dds,normalized=TRUE)[,dds$sample == lvl] +1)))
    write.table(baseMeanPerLvl, paste0(path,"Single_cor_plot_", colnames(tmp)[1],"_", colnames(tmp)[2], '_raw_data',".csv"), sep=',', quote=F, col.names = T, row.names=T)                            
 
    file <- paste0(path,"Single_cor_plot_", colnames(tmp)[1],"_", colnames(tmp)[2],".tiff")
                                
      p=ggplot(data = tmp, aes(x = tmp[,1], y = tmp[,2])) + 
      geom_point(alpha=0.1, size=2, colour = "#110AEA") +
      xlab(paste0("log2 normalized count tumor tissue: ", colnames(tmp)[1])) + ylab(paste0("log2 normalized count PDO: ", colnames(tmp)[2])) +
      theme_classic() +
      theme(axis.title = element_text(size=10, family="Helvetica"),
            axis.text = element_text(size=10, family="Helvetica"),
            plot.title = element_text(size=10, family="Helvetica"),
            legend.position = "none") +
      geom_smooth(method = "lm", se = TRUE) 
      print(p)

      print(paste0('Samples used: ', pairs[[i]][1]))
      
      rsq <- function(x, y) summary(lm( tmp[,2]~tmp[,1]))$r.squared
      print(rsq(obs, mod))                             
}

# ### Extract the 3412 genes and rank them

path='/All.matched.Ksamples.ChIP.Dgroup_1organoid/clean_KcrKorg_genes/'

# +
if (!file.exists(file.dds)){
    
    dds <- DESeqDataSetFromMatrix(countData = counts.noMT[,colnames(counts.noMT)%in%samples$sample],
                              colData = samples,
                              design = ~ group)

    # Removing non-expressed genes
    dds <- dds[rowSums(counts(dds)) >= 10,]
    
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
table(dds$condition)
table(dds$tissue)
table(dds$condition, dds$patient)
# -

contrast=c("group", "tissue", "organoids")

# +
table(colData(dds)[,"condition"])
contrast
res <- results(dds, contrast=contrast, alpha=0.05)

# Summary
res[order(res$padj),]
summary(res)
sum(res$padj <= 0.05, na.rm=TRUE)

# +
# Extract, annotate and save results for all genes 
res2=as.data.frame(results(dds, contrast=contrast))
res2$gene_id=rownames(res2)
res2$padjxFC=round(-log10(res2$padj)*sign(res2$log2FoldChange), digit=3)
res2=merge(res2, annotation, by="gene_id")
rownames(res2)=res2$gene_id

gene_excluded=res2[as.character(intersect_Kcr_anno$intersect_Kcr),]

## Remove all records with padj=NA and with missing or duplicated gene_names
### This files are for GSEA
omit.ids=gene_excluded[duplicated(gene_excluded$gene_name) | gene_excluded$gene_name=="" | !complete.cases(gene_excluded[,"padj"]), ]
omit.ids_list=omit.ids$gene_id

# removing genes with padj= NA
write.table(gene_excluded[!gene_excluded$gene_id%in%omit.ids_list,c("gene_name","padjxFC")], paste0(path, name, ".res.all.pval.sign.FC.rnk"), row.names = F, col.names = T, quote=F, sep="\t")

# considering all of them ranked for log2FC
write.table(gene_excluded[, c("gene_name","log2FoldChange")], paste0(path, name, ".res.all.log2FC.rnk"), row.names = F, col.names = T, quote=F, sep="\t")
