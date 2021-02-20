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

# ### Analysis of all RNA-seq samples - matched patients with ChIP-seq
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

path='/analysis/DESeq/All.matched.ChIP.Dgroup/'

# +
name="All.matched.ChIP.Dgroup"

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

samples=droplevels(sampleinfo)
samples$group=factor(samples$group, levels=c("CRC_N", "CRC_K"))

file.dds = paste0(path, "org_dds.", name, ".noMT.rds")
file.rld = paste0(path,"org_rld.", name, ".noMT.rds")
# -

# ### Create dds and rld objects
#
# With moderation of log2 fold changes. betaPrior=TRUE by default in DESeq2 versions < 1.6

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

    #Save dds and rld
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

# ## PCA 
#
# ### rld

# +
get_prcomp=function(rld, PCA, PCB, color, shape=NULL) {
  
  ntop=500
  rv <- matrixStats::rowVars(assay(rld))
  
  A=as.numeric(gsub("PC", "", PCA))
  B=as.numeric(gsub("PC", "", PCB))

  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(rld)[select,]), center=TRUE, scale=TRUE)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) 
  
  pcaData = pca$x
  pcaData = merge(pcaData, colData(rld), by="row.names")
 
 if(is.null(shape) & color=="condition") {
  color=pcaData[,color]
  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], fill=color)) + 
        geom_point(size=4, pch=21) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
        theme_classic() +
        theme(axis.title=element_text(size=14), axis.text = element_text(size=12), legend.title = element_blank(), legend.key.size = unit(0.6,"cm"), legend.text = element_text(size=12)) +
               scale_fill_manual(values=c("#99258A","#DA1781","#26DB32","green4")) +
     guides(fill = guide_legend(override.aes=list(shape=21)))
  
 } else if(is.null(shape)) {
  color=pcaData[,color]

  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], fill=color)) + 
        geom_point(size=4, pch=21) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
        theme_classic() +
        theme(axis.title=element_text(size=14), axis.text = element_text(size=12), legend.title = element_blank(), legend.key.size = unit(0.6,"cm"), legend.text = element_text(size=12)) +
             geom_point(size=4, colour="black",pch=21)
  } else {
  color=pcaData[,color]
  shape=pcaData[,shape]

  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], fill=color, shape=shape)) + 
              geom_point(size=4) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
        theme_classic() +
        theme(axis.title=element_text(size=14), axis.text = element_text(size=12), legend.title = element_blank(), legend.key.size = unit(0.6,"cm"), legend.text = element_text(size=12)) +          
    scale_shape_manual(values=c(21,24,23,22)) +
    scale_fill_manual(values=c("#e41a1c", "#377eb8", "#FF62BC", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#999999", "#65B7F3", "#4daf4a"))+
     guides(fill = guide_legend(override.aes=list(shape=21)))

  }

  
}

plotPCA(rld, intgroup=c("group"))
# -

# ## Loadings characterization 

# +
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
  
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld)[select,]))
#get loadings from pca
loadings=as.data.frame(pca$rotation)

#HVG
HVG=row.names(loadings)
# -

### Save 500 HVG
HVG_pc1=merge(as.data.frame(HVG), annotation[,c("gene_id", "gene_name" )], by.x="HVG", by.y="gene_id")
head(HVG_pc1)
write.table(HVG_pc1[,2], paste0(path, "HVG_pc1.txt"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# +
### Import genes filtered for Immune related genes (See Methods and Supplementary Table 2)
HVG_filtered=read.table(paste0(path, "immune_related_genes_uniq_manuallyCurated.txt"))
HVG_filtered_eng=merge(HVG_filtered, annotation[,c("gene_id", "gene_name" )], by.x="V1", by.y="gene_name")
HVG_filtered_eng_list=HVG_filtered_eng$gene_id

#Select HVG to keep
HVG_filtered_ens=(setdiff(HVG,HVG_filtered_eng_list))

#Filter the rld
rld_filtered=(assay(rld)[row.names(assay(rld))%in%HVG_filtered_ens,])
# -

# ### Fig1 d

# +
get_prcomp=function(rld_filtered, rld, PCA, PCB, color, shape=NULL) {
    
  pca <- prcomp(t(rld_filtered), center=TRUE, scale=TRUE)
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) 

  A=as.numeric(gsub("PC", "", PCA))
  B=as.numeric(gsub("PC", "", PCB))
    
  pcaData = pca$x
  pcaData = merge(pcaData, colData(rld), by="row.names")
  
  if(is.null(shape) & color=="condition") {
  color=pcaData[,color]

  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], fill=gsub("_", " ", color))) + 
       geom_point(size = 8, pch=21) + theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24), legend.key.size = unit(0.6,"cm")) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
        scale_fill_manual(values=c("#99258A","#DA1781","#26DB32","green4")) +
        guides(fill = guide_legend(override.aes=list(shape=21)))

  } else if(is.null(shape)) {
  color=pcaData[,color]

  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], color=gsub("_", " ", color))) + 
geom_point(size = 8, pch=21) + theme_classic() +
 theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24), legend.key.size = unit(0.6,"cm")) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
        guides(fill = guide_legend(override.aes=list(shape=21)))


  } else {
  color=pcaData[,color]
  shape=pcaData[,shape]

  ggplot(data=pcaData, aes(x=pcaData[,PCA], y=pcaData[,PCB], fill=color, shape=gsub("_", " ", shape))) + 
  geom_point(size = 10, aes(colour=color), stroke = 1) + 
  theme_classic() +
  theme(axis.title = element_text(size=28),
        axis.text = element_text(size=24),
        legend.title = element_blank(),
        legend.text = element_text(size=24), legend.key.size = unit(0.6,"cm")) +
        xlab(paste0(PCA, ": ", round(percentVar[A] * 100),"% variance")) +
        ylab(paste0(PCB, ": ", round(percentVar[B] * 100),"% variance")) +
    scale_shape_manual(values=c(21,24,23,22)) +
    scale_fill_manual(values=c("#e41a1c", "#377eb8", "#FF62BC", "#984ea3", "#ff7f00", "#ffff33", "#a65628","#999999", "#65B7F3", "#4daf4a"))+
    scale_colour_manual(values=c("#930708", "#073c68", "#891055", "#491051", "#b7610c", "#7c7c07", "#662d0c","#635f5f", "#1b5077", "#126610"))+
    guides(fill = guide_legend(override.aes=list(shape=21)))

  }

}

get_prcomp(rld_filtered,rld,"PC1","PC2","patient", "condition")
# -

# ##### 500 HVG
