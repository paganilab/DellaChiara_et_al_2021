# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# ### Analysis of H3K27Ac between PDOs and normal tissue - considering the input
#
# Author: Federica Gervasoni

# #### Library

source("https://bioconductor.org/biocLite.R")
library(DiffBind)
library(gplots)
library(data.table)
library(DESeq2)
library(pheatmap)
library(hyperSpec)
library(forcats)

## set universal plot size:
options(repr.plot.width=6, repr.plot.height=6)

# ### Read sampleinfo and extract the information

# +
samples_table <- read.csv(file = "/DEA/DBA_Korg_Ntissue_Input_Jup/DiffBind_H3K27Ac_table.csv", header = T)

sampleSheet=samples_table[samples_table$Tissue=="4KE" | samples_table$Tissue=="8KE" | samples_table$Tissue=="10KE" | samples_table$Tissue=="11KW" | samples_table$Tissue=="13KE" | samples_table$Tissue=="18KE" | samples_table$Tissue=="22KE" | samples_table$Tissue=="24KE" | samples_table$Tissue=="36KE" | samples_table$Tissue=="41KE" | samples_table$Tissue=="CR_28" | samples_table$Tissue=="CR_29" | samples_table$Tissue=="CR_37" | samples_table$Tissue=="CR_41_mp"| samples_table$Tissue=="CR_28_mp",]

K_org=samples_table[samples_table$Tissue=="4KE" | samples_table$Tissue=="8KE" | samples_table$Tissue=="10KE" | samples_table$Tissue=="11KW" | samples_table$Tissue=="13KE" | samples_table$Tissue=="18KE" | samples_table$Tissue=="22KE" | samples_table$Tissue=="24KE" | samples_table$Tissue=="36KE" | samples_table$Tissue=="41KE",]

N_tissue=samples_table[samples_table$Tissue=="CR_28" | samples_table$Tissue=="CR_29" | samples_table$Tissue=="CR_37" | samples_table$Tissue=="CR_41_mp"| samples_table$Tissue=="CR_28_mp",]
# -

# ### Create dba

dba <- dba(sampleSheet=sampleSheet,
         minOverlap = 2, 
         config = data.frame(RunParallel = TRUE,
                             reportInit = "DBA",
                             DataType = DBA_DATA_FRAME))
dba
#33131 sites in matrix (61994 total) 

# ### Count reads
# ##### DBA_SCORE_TMM_MINUS_FULL
# TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size

dba=dba.count(prova, score=DBA_SCORE_TMM_MINUS_FULL,  bParallel=T)

# ### DEA

# +
# Add K_org vs N_crypts contrast
dba = dba.contrast(dba, group1 =dba$masks$K_org, group2 = dba$masks$N_crypts, name1 = "K_org", name2 = "N_crypts" , minMembers=2)

dba = dba.analyze(dba, bParallel=T)

# Report
db_K_org_vs_N_crypts = dba.report(dba, DataType = DBA_DATA_GRANGES, contrast=1)
# -

# ### Save dba and regions object

# +
path="/DEA/DBA_Korg_Ntissue_Input_Jup/"
name= "Korg_Ntissue"

# Save the Peakset
dba_peakSet <- dba.peakset(dba, bRetrieve=TRUE, DataType = DBA_DATA_GRANGES)
dba_rg <- data.frame(seqnames=seqnames(dba_peakSet),
starts=start(dba_peakSet)-1,
ends=end(dba_peakSet))
dba_rg$region_ID <- paste0("reg_", (1:nrow(dba_rg)))

write.table(dba_rg, file=paste0(path, "dba_", name, ".bed"), sep='\t', quote=F, col.names = F, row.names = F)

#Save dba object
dba.save(dba, dir=path, file=name, ext='RData', pre='dba_')


# +
path="/DEA/DBA_Korg_Ntissue_Input_Jup/"
name= "Korg_Ntissue"

dba=dba.load(file=name, path, pre='dba_', ext='RData')
# -

# ### Calculate the perchentage of the genome coverade by the enhancerome (33K total peaks considered)

dba_rg$length=abs(dba_rg$starts-dba_rg$ends)
paste0('Percentage of the genome covered by the enhancerome considere for the analysis: ', round((sum(dba_rg$length)/3049315783)*100), '%')

# ### DESeq2 analysis
# #### Obtain the dba set for DESEQ2 analysis

dba = dba.analyze(dba, bParallel=T, bReduceObjects=F)

# ### K_organoids vs N crypts

# +
dds= dba$contrasts[[1]]$DESeq2$DEdata

colnames=c("SQ_2157", "SQ_1990", "SQ_2010", "SQ_2163", "SQ_2204", "SQ_2212", "SQ_2216", "SQ_2222", "SQ_2288", "SQ_2303", "SQ_2298", "SQ_2145", "GSM2058021", "GSM2058022", "GSM2058023")
colnames(dds)=colnames

sampleinfo=samples_table
samples=droplevels(sampleinfo[sampleinfo$SampleID%in%colnames,])
samples$Treatment=factor(samples$Treatment, levels=c("K_org", "N_crypts"))
samples$Tissue=factor(samples$Tissue, levels=c("10KE", "11KW", "13KE" ,"18KE", "22KE" ,"24KE", "36KE" ,"41KE", "4KE" ,"8KE","CR_41_mp", "CR_28_mp", "CR_28", "CR_29", "CR_37"))
samples$SampleID=factor(samples$SampleID, levels=c("SQ_2157", "SQ_1990", "SQ_2010", "SQ_2163", "SQ_2204", "SQ_2212", "SQ_2216", "SQ_2222", "SQ_2288", "SQ_2303", "SQ_2298", "SQ_2145", "GSM2058021", "GSM2058022", "GSM2058023"))

dds$Tissue=samples$Tissue
dds$SampleID=samples$SampleID
dds$Treatment=samples$Treatment
# -

# ### Add info to dds

# +
rld <- rlog(dds, blind = FALSE)

saveRDS(object = dds, file=paste0(path, "DESeq2/", "dds_", name, ".rds")) 
saveRDS(object = rld, file=paste0(path, "DESeq2/", "rld_", name, ".rds"))
# -

dds=readRDS(file=paste0(path, "DESeq2/", "dds_", name, ".rds"))
rld=readRDS(file=paste0(path, "DESeq2/", "rld_", name, ".rds"))

# ### Fig3 a - Correlation matrix of H3K27Ac for all the available samples

# +
# Select the 500 HVG from dds
ntop=500
rv <- matrixStats::rowVars(assay(dds))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# +
#matrix of sample-to-sample distances
sampleDist <- pearson.dist(t(assay(dds)[select,]))
sampleDistMat <- as.matrix(sampleDist)

matrix_corr = as.matrix(sampleDistMat)
breaks = seq(0,max(sampleDistMat),length.out=1000)
gradient1 = colorpanel( sum( breaks[-1]<=0.27 ), "#666A01" , "#dddebd" )
gradient2 = colorpanel( sum( breaks[-1]>0.27 ), "#dddebd", "white" )
hm.colors = c(gradient1,gradient2)

# #annotation
anno <- as.data.frame(colData(dds)[, c("groups", "Tissue")])

print('Pearson distance')

#heatmap across all groups
pdf=pheatmap(sampleDistMat,
         clustering_distance_rows=sampleDist,
         clustering_distance_cols=sampleDist,
         col=hm.colors,
         annotation_col = anno,
         annotation_names_row=F,
         border_color = NA)
# -

# ### Results Differential binding analysis
# lfcThreshold = 0

# +
table(colData(dds)[,"groups"])
contrast=c("groups", "K_org", "N_crypts")
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
# Create an annotation to get the GRANGE of the peaks and a unique identifier
annotation_GRANGE=dba_rg
annotation_GRANGE$region_ID <- paste0("reg_", (1:nrow(annotation_GRANGE)))
rownames(annotation_GRANGE)=annotation_GRANGE$region_ID

# Extract, annotate and save results for all genes 
res2=as.data.frame(res)
res2$region_ID =paste0("reg_", (1:nrow(annotation_GRANGE))) #add id to each regions
rownames(res2)=res2$region_ID

res2$padjxFC=round(-log10(res2$padj)*sign(res2$log2FoldChange), digit=3)
res2=merge(res2, annotation_GRANGE, by="region_ID", sort = F)
# -

# ### DGEgenes logFC >= 2 & p-value <= 0.01

# +
DGEgenes <- rownames(subset(res[order(res$padj),], padj < 0.01  & abs(log2FoldChange)>2))

# Upregulated
DGEgenes.up <- rownames(subset(res[order(res$padj),], padj <= 0.01  & log2FoldChange>=2))
table(DGEgenes.up%in%DGEgenes)

# Downregulated
DGEgenes.down <- rownames(subset(res[order(res$padj),], padj <= 0.01  & log2FoldChange<=(-2)))
table(DGEgenes.down%in%DGEgenes)
# -

# ### Fig3 b - Volcano plot

volcano= as.data.frame(res2)
volcano$threshold=ifelse(volcano$log2FoldChange >= 2 & volcano$padj < 0.01,"A", ifelse(volcano$log2FoldChange <= -2 & volcano$padj < 0.01, "B", "C"))


ggplot(volcano, aes(x=log2FoldChange, y =-log10(padj), color=threshold)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c( "A"="#D90416","B"="#033E8C", "C"="grey")) +
  xlab("log2 fold change") + ylab("-log10 adj p-value") +
  theme_classic() +
  theme(legend.position="none") +
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed") +
    scale_y_continuous(trans = "log1p")

# ### Pathway analysis


gprof=read.table("KEGG_Korg_log2FC2_padj0.01_DERNAseq.txt", sep="\t", header=T)
gprof_df=gprof[0:10, ]

# ### Fig3 f - Bubble plot 

p=ggplot(gprof_df, aes(x = -log10(FDR), y = fct_reorder(Description, -log10(FDR)))) + 
               geom_point(aes(size = GeneRatio, color = -log10(FDR))) +
               scale_size(range = c(3, 9)) +
               theme_bw(base_size = 14) +
        scale_colour_gradient2(limits=c(-4, 4), low="#033E8C", high="#D90416", midpoint = -1) +
        geom_line() +
        geom_vline(xintercept = -log10(0.05), linetype="dotted", 
                color = "grey", size=1) +
        ylab(NULL) +
        ggtitle("GO pathway enrichment PDO") +
        theme(legend.position="none") +
        theme_classic() +
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank())
p
