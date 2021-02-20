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

# ### Analysis of all tumor tissue RNA-seq samples - matched patients with ChIP-seq
#
# Author: Federica Gervasoni

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
library(VennDiagram)

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

path='/All.matched.Ksamples.ChIP.Dgroup_1organoid/'

# +
matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

# K crypts
samples_Kcr=droplevels(sampleinfo[sampleinfo$condition != "K_organoids" & sampleinfo$group !="CRC_N" & sampleinfo$patient%in%matched,])

# K organoids
excl=c("SQ_1959", "SQ_1953", "SQ_1954","SQ_1956", "SQ_1963", "SQ_1970", "SQ_1973", "SQ_1974", "SQ_1976", "SQ_1987", "SQ_2013", "SQ_2014", "SQ_2017", "SQ_2254", "SQ_2265")

samples_Korg=droplevels(sampleinfo[sampleinfo$condition != "K_tissue" & sampleinfo$group !="CRC_N" & sampleinfo$patient%in%matched & !sampleinfo$sample%in%excl,])

# +
nsamples=10
countsMin=5

samples_Kcr_mtx = counts.noMT[,colnames(counts.noMT)%in%samples_Kcr$sample]
samples_Kcr_mtx=samples_Kcr_mtx[rowSums(samples_Kcr_mtx) >= 10,] # The total rowmean must be higher than 10

samples_Korg_mtx = counts.noMT[,colnames(counts.noMT)%in%samples_Korg$sample]
samples_Korg_mtx=samples_Korg_mtx[rowSums(samples_Korg_mtx) >= 10,]
# -

# ### Fig1e - Venny plot

grid.newpage()
draw.pairwise.venn(length(row.names(samples_Kcr_mtx)), length(row.names(samples_Korg_mtx)), length(intersect(row.names(samples_Kcr_mtx),row.names(samples_Korg_mtx))), category = c("Kcr", "Korg"), fill = c("#62D254", "#62D254"), alpha = rep(0.4, 2), cex=3, col="black", lwd=1, cat.dist=0.1)

# ### Check the numbers of overlapping genes specific for each patient

pairTochose=c("CRC4", "SQ_1958")
pairTochose=c("CRC8", "SQ_1955")
pairTochose=c("CRC10", "SQ_1964")
pairTochose=c("CRC11", "SQ_2018")
pairTochose=c("CRC13", "SQ_1975")
pairTochose=c("CRC18", "SQ_2016")
pairTochose=c("CRC22", "SQ_2021")
pairTochose=c("CRC24", "SQ_2015")
pairTochose=c("CRC36", "SQ_2253")
pairTochose=c("CRC41", "SQ_2264")

df_pairs=data.frame("patient"= c("CRC4","CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41"))
df_pairs["sample"]=c("SQ_1958","SQ_1955", "SQ_1964", "SQ_2018", "SQ_1975", "SQ_2016", "SQ_2021", "SQ_2015", "SQ_2253", "SQ_2264")

pairs=as.list(c('CRC4-SQ_1958', 'CRC8-SQ_1955', 'CRC10-SQ_1964', 'CRC11-SQ_2018', 'CRC13-SQ_1975', 'CRC18-SQ_2016', 'CRC22-SQ_2021', 'CRC24-SQ_2015', 'CRC36-SQ_2253', 'CRC41-SQ_2264'))

# +
save <- function(pairs){ 

    Num_count=5

    # K crypts
    samples_Kcr=droplevels(sampleinfo[sampleinfo$condition != "K_organoids" & sampleinfo$tissue !="biopsy" & sampleinfo$group !="CRC_N" & sampleinfo$passage !="early" &  sampleinfo$patient%in%strsplit(pairs[[1]][1], "-", fixed = TRUE)[[1]][1],])
    # K organoids
    samples_Korg=droplevels(sampleinfo[sampleinfo$sample%in%strsplit(pairs[[1]][1], "-", fixed = TRUE)[[1]][2],])

    samples_Kcr_mtx=as.data.frame(rownames(counts.noMT))
    colnames(samples_Kcr_mtx)='sample'
    rownames(samples_Kcr_mtx)=samples_Kcr_mtx$sample
    samples_Kcr_mtx['sample'] = counts.noMT[,colnames(counts.noMT)%in%samples_Kcr$sample]
    samples_Kcr_mtx=samples_Kcr_mtx[rowSums(samples_Kcr_mtx) >= Num_count, , drop=FALSE]
    dim(samples_Kcr_mtx)

    samples_Korg_mtx=as.data.frame(rownames(counts.noMT))
    colnames(samples_Korg_mtx)='sample'
    rownames(samples_Korg_mtx)=samples_Korg_mtx$sample
    samples_Korg_mtx['sample'] = counts.noMT[,colnames(counts.noMT)%in%samples_Korg$sample]
    samples_Korg_mtx=samples_Korg_mtx[rowSums(samples_Korg_mtx) >= Num_count, , drop=FALSE]
    dim(samples_Korg_mtx)

    table_patient_int=data.frame("tot" = c(length(intersect(row.names(samples_Kcr_mtx),row.names(samples_Korg_mtx))) + length(setdiff(row.names(samples_Kcr_mtx),row.names(samples_Korg_mtx))) + length(setdiff(row.names(samples_Korg_mtx),row.names(samples_Kcr_mtx)))))
    table_patient_int['int']=length(intersect(row.names(samples_Kcr_mtx),row.names(samples_Korg_mtx)))  
    table_patient_int['Kcr']=length(setdiff(row.names(samples_Kcr_mtx),row.names(samples_Korg_mtx)))
    table_patient_int['Korg']=length(setdiff(row.names(samples_Korg_mtx),row.names(samples_Kcr_mtx)))
    rownames(table_patient_int)=as.character(strsplit(pairs[[1]][1], "-", fixed = TRUE)[[1]][1])
    return(table_patient_int)
    }

#Execute the function above
res = sapply(pairs, save) 
# -

df=as.data.frame(t(res))
rownames(df)=c("CRC4","CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

colnames(df)=c('tot','3_int','2_Kcr','1_Korg')

df = df  %>% 
    tibble::rownames_to_column('Sample')  %>% 
    gather(key='Param', value = 'Value', - 'Sample')

# ### FigS1 b - Stacked parplot

ggplot(df %>% filter(Param != 'tot'))+
    geom_bar(aes(x=Sample, y=Value, fill=Param), stat = 'identity', position='fill') +
    ylab("Percentage shared genes") +
    xlab("Patients") +
    scale_fill_manual(values=c("#8EFFAA","#05C7F2","#110AEA")) +
    theme_classic()
