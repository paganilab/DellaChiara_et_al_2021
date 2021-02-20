# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: R35_RNAseq
#     language: R
#     name: r35_rnaseq
# ---

# ### Gene expression level of YAP and TAZ
#
# Author: Federica Gervasoni

# ### Library

library(DESeq2)
library(reshape)
library(data.table)
library(ggpubr)

# Read sampleifo
sampleinfo=readRDS("sampleinfo.rds")

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

# +
dds.all=readRDS("/org_dds.All.matched.ChIP.Dgroup_1organoid.noMT.rds")

# Modify the gene_id
rownames(dds.all)=gsub("[.].*", "", rownames(dds.all))

norm.all=counts(dds.all, normalized=TRUE)

norm.all.log2=log2(norm.all + 1)

# +
YAP_TAZ_gene=c( "WWTR1")

norm.all.log2_sym=merge(norm.all.log2, annotation[,c("gene_id", "gene_name")], by.x=0, by.y="gene_id")
unique=unique(norm.all.log2_sym$gene_name)
norm.all.log2_sym=norm.all.log2_sym[norm.all.log2_sym$gene_name==unique, ]
rownames(norm.all.log2_sym)=norm.all.log2_sym$gene_name
norm.all.log2_sym=norm.all.log2_sym[, -c(1, 31)]

norm.all.log2.YT  <- norm.all.log2_sym[row.names(norm.all.log2_sym)%in%YAP_TAZ_gene,]

table(row.names(norm.all.log2.YT)%in%YAP_TAZ_gene)
table(colnames(norm.all.log2.YT)==colData(dds.all)[,"sample"])

norm.all.log2.YT.melt=as.data.frame(norm.all.log2.YT)
norm.all.log2.YT.melt$gene_id=row.names(norm.all.log2.YT.melt)
norm.all.log2.YT.melt=setDT(norm.all.log2.YT.melt)
norm.all.log2.YT.melt=melt(norm.all.log2.YT.melt, id.vars="gene_id", variable.name="sample")
norm.all.log2.YT.melt=merge(norm.all.log2.YT.melt, sampleinfo, by="sample")

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")
norm.all.log2.YT.melt=norm.all.log2.YT.melt[norm.all.log2.YT.melt$patient%in%matched & norm.all.log2.YT.melt$condition!="K_biopsy" & norm.all.log2.YT.melt$condition!="N_biopsy" & norm.all.log2.YT.melt$condition!="N_organoids",]

p <- ggplot(norm.all.log2.YT.melt, aes(x=condition, y=value, color=condition, fill=condition)) + 
  geom_violin(alpha=.2) +
  scale_fill_manual(values=c("#DE4A40", "#D90416","#033E8C")) +
  scale_color_manual(values=c("#DE4A40", "#D90416","#033E8C")) +
  geom_violin(trim=TRUE) +
  labs(x="", y="Normalised counts log2 transformed") +
  theme_classic() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="black", fill="black") + 
  theme(axis.title=element_text(size=14), axis.text = element_text(size=12), axis.text.x = element_text(angle = 60, vjust=0.8, hjust = 1), legend.title = element_blank()) +
  scale_x_discrete(limits=c("N_tissue", "K_tissue", "K_organoids"))

# Statistics
my_comparisons <- list( c("K_organoids", "N_tissue"), c("K_tissue", "N_tissue"))
p + stat_compare_means(comparisons = my_comparisons, label.y = c(10.3, 8.8), method = "t.test")  

# +
print("WWTR1")

print("wilcoxon")
compare_means(value ~ condition,  data = norm.all.log2.YT.melt,ref.group='N_tissue', method = "wilcox.test", p.adjust.method = "fdr")

# +
YAP_TAZ_gene=c( "YAP1")

norm.all.log2_sym=merge(norm.all.log2, annotation[,c("gene_id", "gene_name")], by.x=0, by.y="gene_id")
unique=unique(norm.all.log2_sym$gene_name)
norm.all.log2_sym=norm.all.log2_sym[norm.all.log2_sym$gene_name==unique, ]
rownames(norm.all.log2_sym)=norm.all.log2_sym$gene_name
norm.all.log2_sym=norm.all.log2_sym[, -c(1, 31)]

norm.all.log2.YT  <- norm.all.log2_sym[row.names(norm.all.log2_sym)%in%YAP_TAZ_gene,]

table(row.names(norm.all.log2.YT)%in%YAP_TAZ_gene)
table(colnames(norm.all.log2.YT)==colData(dds.all)[,"sample"])

norm.all.log2.YT.melt=as.data.frame(norm.all.log2.YT)
norm.all.log2.YT.melt$gene_id=row.names(norm.all.log2.YT.melt)
norm.all.log2.YT.melt=melt(norm.all.log2.YT.melt, id.vars="gene_id", variable.name="sample")
norm.all.log2.YT.melt=merge(norm.all.log2.YT.melt, sampleinfo, by="sample")

matched=c("CRC4", "CRC8", "CRC10", "CRC11", "CRC13", "CRC18", "CRC22", "CRC24", "CRC36", "CRC41")
norm.all.log2.YT.melt=norm.all.log2.YT.melt[norm.all.log2.YT.melt$patient%in%matched & norm.all.log2.YT.melt$condition!="K_biopsy" & norm.all.log2.YT.melt$condition!="N_biopsy" & norm.all.log2.YT.melt$condition!="N_organoids",]

p <- ggplot(norm.all.log2.YT.melt, aes(x=condition, y=value, color=condition, fill=condition)) + 
  geom_violin(alpha=.2) +
  scale_fill_manual(values=c("#DE4A40", "#D90416","#033E8C")) +
  scale_color_manual(values=c("#DE4A40", "#D90416","#033E8C")) +
  geom_violin(trim=TRUE) +
  labs(x="", y="Normalised counts log2 transformed") +
  theme_classic() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, color="black", fill="black") + 
  theme(axis.title=element_text(size=14), axis.text = element_text(size=12), axis.text.x = element_text(angle = 60, vjust=0.8, hjust = 1), legend.title = element_blank())  +
  scale_x_discrete(limits=c("N_tissue", "K_tissue", "K_organoids"))

# Statistics: show p value
my_comparisons <- list( c("K_organoids", "N_tissue"), c("K_tissue", "N_tissue"))
p + stat_compare_means(comparisons = my_comparisons, label.y = c(12.3, 12), method = "t.test") 

# +
print("YAP1")

print("wilcoxon")
compare_means(value ~ condition,  data = norm.all.log2.YT.melt,ref.group='N_tissue', method = "wilcox.test", p.adjust.method = "fdr")
