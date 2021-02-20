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

# ### DE analysis between PDO and Normal tissue - matched patients with ChIP-seq
#
# Author: Federica Gervasoni
#
# design = ~ Patient class

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
library(ggpubr)
library(VennDiagram)

## set universal plot size:
options(repr.plot.width=5, repr.plot.height=5)

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
rownames(counts)=gsub("[.].*", "", rownames(counts))

counts.noMT=counts[!rownames(counts)%in%MT.genes,]
    
length(unique(rownames(counts)))
length(unique(rownames(counts.noMT)))

identical(rownames(sampleinfo), colnames(counts.noMT))
# -

# ## Analysis 
#
# ### Set comparison

path='/DGE/KorgvsNcr/KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient/'

# +
name="KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient"

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")

samples=droplevels(sampleinfo[sampleinfo$condition !="K_tissue" & sampleinfo$patient%in%matched, ])
samples$group=factor(samples$group, levels=c("CRC_N", "CRC_K"))
samples$condition2=factor(samples$condition2, levels=c("N_tissue", "K_organoids"))

file.dds = paste0(path, "org_dds.", name, ".noMT.rds")
file.rld = paste0(path, "org_rld.", name, ".noMT.rds")

contrast=c("group", "CRC_K", "CRC_N")

table(samples$group)
table(samples$condition)
table(samples$tissue)
table(samples$condition, samples$patient)
# -

# ### Create dds and rld objects
#
# With moderation of log2 fold changes. betaPrior=TRUE by default in DESeq2 versions < 1.6

# +
if (!file.exists(file.dds)){
    
    dds <- DESeqDataSetFromMatrix(countData = counts.noMT[,colnames(counts.noMT)%in%samples$sample],
                              colData = samples,
                              design = ~ patient + group)

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
table(dds$passage, dds$patient)
# -

# ### DGE

# +
table(colData(dds)[,"condition"])
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

# Lowest logFC
resSig <- subset(res2, res2$padj <= 0.05)
min(resSig$log2FoldChange[resSig$log2FoldChange > 0])
max(resSig$log2FoldChange[resSig$log2FoldChange < 0])
# -

# ### DGEgenes logFC >= 1 & p-value <= 0.05

# +
DGEgenes <- rownames(subset(res[order(res$padj),], padj <= 0.05  & abs(log2FoldChange)>=0))
table(DGEgenes%in%resSig$gene_id)

# Upregulated
DGEgenes.up <- rownames(subset(res[order(res$padj),], padj <= 0.05  & log2FoldChange>=0))
table(DGEgenes.up%in%DGEgenes)

# Downregulated
DGEgenes.down <- rownames(subset(res[order(res$padj),], padj <= 0.05  & log2FoldChange<=(0)))
table(DGEgenes.down%in%DGEgenes)
# -
# ### Check if the genes associated to Differentially active regions (H3K27Ac- log2FC >2 and padj<0.01) in chip-seq Korg vs N tissue are also DE

# +
### Check if the genes annotated to DA in Korg and Ntissue are DE in K tissue.
path_save='/DEA/DBA_Korg_Ntissue_Input_Jup/DESeq2/annotation/'

genes_list=read.table("/DEA/DBA_Korg_Ntissue_Input_Jup/DESeq2/annotation/K_org_UP_log2FC2_padj0.01_annoHiC_actProm_sym.txt")
genes=merge(genes_list, annotation[,c("gene_id", "gene_name")], by.x="V1", by.y="gene_name") 
intersect_KorgUP_genes=as.data.frame(intersect(DGEgenes.up,genes$gene_id))
colnames(intersect_KorgUP_genes)="gene_id"

#add gene symbols and save
intersect_KorgUP_genes_df=merge(intersect_KorgUP_genes, annotation[,c("gene_id", "gene_name")], by="gene_id")
write.table(intersect_KorgUP_genes_df, paste0(path_save, "K_org_UP_log2FC2_padj0.01_annoHiC_actProm_sym_diffGenes_Korg.txt"), row.names = F, col.names = T, quote=F, sep="\t")

grid.newpage()
draw.pairwise.venn(length(DGEgenes.up), length(genes$gene_id), length(intersect(DGEgenes.up,genes$gene_id)), category = c("KorgvsNcr", "enh_DEKorg"), fill = c("skyblue", "pink1"), cex=1, col="black", lwd=1, cat.dist=0.1)
# -


# ### FigS4 - Boxplot plot

# +
### Heatmap on the specific genes 

dds.all=readRDS("/All.matched.ChIP.Dgroup_1organoid/org_dds.All.matched.ChIP.Dgroup_1organoid.noMT.rds")

# Modify the gene_id
rownames(dds.all)=gsub("[.].*", "", rownames(dds.all))

norm.all=counts(dds.all, normalized=TRUE)

norm.all.log2=log2(norm.all + 1)

# +
gene_list=as.character(intersect(DGEgenes.up,genes$gene_id))

crc.expr.normalized_genes.check  <- norm.all.log2[row.names(norm.all.log2)%in%gene_list,]

table(row.names(crc.expr.normalized_genes.check)%in%gene_list)

crc.expr.normalized_genes.check.melt=as.data.frame(crc.expr.normalized_genes.check)
crc.expr.normalized_genes.check.melt$gene_id=row.names(crc.expr.normalized_genes.check.melt)
crc.expr.normalized_genes.check.melt=setDT(crc.expr.normalized_genes.check.melt)
crc.expr.normalized_genes.check.melt=melt(crc.expr.normalized_genes.check.melt, id.vars="gene_id", variable.name="sample")
crc.expr.normalized_genes.check.melt=merge(crc.expr.normalized_genes.check.melt, sampleinfo, by="sample")

# pdf(paste0(path, "boxplot_expression495genes_Funcional_NKtis_PDO", ".pdf"))
p <- ggplot(crc.expr.normalized_genes.check.melt, aes(x=condition, y=value, color=condition, fill=condition)) + 
  geom_boxplot(alpha=.5, notch=TRUE) +
  scale_fill_manual(values=c("#C76E24", "#AF0416", "#033E8C")) +
  scale_color_manual(values=c("#C76E24", "#AF0416","#033E8C")) +
  labs(x="", y="Normalised counts log2 transformed") +
  theme_classic() +
  theme(axis.title=element_text(size=14), axis.text = element_text(size=12), axis.text.x = element_text(angle = 60, vjust=0.8, hjust = 1), legend.title = element_blank()) +
  scale_x_discrete(limits=c("N_tissue", "K_tissue","K_organoids"))

# # # # Statistics
my_comparisons <- list( c("K_organoids", "N_tissue"), c("K_tissue", "N_tissue"))
p + stat_compare_means(comparisons = my_comparisons, label.y = c(17, 16), method = "wilcox.test", aes(label = ..p.signif..),
                      paired = FALSE)

compare_means(value ~ condition,  data = crc.expr.normalized_genes.check.melt, ref.group = "N_tissue",
              method = "wilcox.test")

# dev.off()
# -

# ### Save the lists of 495 DE genes in PDOvsNtissue associated with enhancer DEA in PDOvsNtissue
#
# This list was used to run the functional enrichment analysis - Fig3f and Fig3g

# +
#Save 495 genes 
name='K_org_UP_log2FC2_padj0.01_annoHiC_actProm_sym_diffGenes_Korg'

intersect_KorgUP_genes_df_density=merge(res2[,c('chromosome','start_position','end_position', 'gene_id')], intersect_KorgUP_genes, by="gene_id")
intersect_KorgUP_genes_df_density=intersect_KorgUP_genes_df_density[,c(2,3,4,1)]
write.table(intersect_KorgUP_genes_df_density, paste0(path_save, name,"_density.bed"), row.names = F, col.names = F, quote=F, sep="\t")

paste0("Biotype of 495 genes")
table(intersect_KorgUP_genes_df_density$biotype)

# +
### Check up-regulated genes in ChIP-seq (DE and recurrent) and RNA-seq - List used in FigS6 a
chip=read.table("/DEA/DBA_Korg_Ntissue_Input_Jup/recurrence_Korg/annotation/recurrence_10-8_K_org_UP_log2FC2_padj0.01_SQ_2385_toAnno_annoHiC_actProm_sym.txt", header=F, sep="\t", quote=NULL)


KorgNcr=read.table("/DESeq/DGE/KorgvsNcr/KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient/KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient.res.all.txt", header=T, sep="\t", quote=NULL)
KorgNcr=subset(KorgNcr, padj <= 0.05 & log2FoldChange > 0 )

intersect(KorgNcr$gene_id,chip$V1)

grid.newpage()
draw.pairwise.venn(length(KorgNcr$gene_name), length(chip$V1), length(intersect(KorgNcr$gene_name,chip$V1)), category = c("KorgvsNcr", "ChIP"), fill = c("skyblue", "pink1"), cex=3, col="black", lwd=1, cat.dist=0.1)

# +
intersection_ChIPseq_PDO_Ntis=as.character(intersect(KorgNcr$gene_name,chip$V1))

write.table(intersection_ChIPseq_PDO_Ntis, paste0(path,"intersection_ChIPseq_UP_genes_POD_Ntissue.txt"), row.names = F, col.names = T, quote=F, sep="\t")
