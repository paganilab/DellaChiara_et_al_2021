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

# ### Analysis of H3K27Ac from Public available data coming from differe types of cancers - considering the input
#
# Author: Federica Gervasoni

# ##### Library

# +
source("https://bioconductor.org/biocLite.R")
library(DiffBind)
library(GenomicFeatures)

##RNAseq
library(reshape2)

##graph
library(pheatmap)
library(ggplot2)

##format
library(tidyr)
library(dplyr)

##stat
library(ggpubr)

##batch
library("sva")
# -

## set universal plot size:
options(repr.plot.width=6, repr.plot.height=6)

# ### Read sampleinfo and extract the information

# +
###import table and select the samples
samples_table <- read.csv(file = "/diffbind/DiffBind_H3K27Ac_table_tumor-normal_primary_tissue_PDOs.csv", header = T)

target=c('tumor_breast', 'tumor_endometrium')
samplesID=c('SQ_2298', 'SQ_2145', 'GSM2058021', 'GSM2058022', 'GSM2058023', #our normal colon
            'GSM1252278', 'GSM1252286', 'GSM1252294', 'GSM1252312', 'GSM1252320', #gastric
                            'GSM1127173', 'GSM2534707', #liver tissue roadmap
                            'GSM896163', 'GSM1013126', #adrenal gland tissue roadmap
                            'GSM2698612', 'GSM2527562', #adrenal gland encode        
                            'GSM1013129', 'GSM906397' ,  #pancreas tissue roadmap
                            'GSM2700597', 'GSM2700502', 'GSM4247380', 'GSM2534105', #body of pancreas encode 
                            'GSM4247568', 'GSM2700584', 'GSM2527594', 'GSM2534198', #thiroid encode
                            'GSM2698835', 'GSM4247275', #uterus encode
                            #tumor
                            'GSM1252274', 'GSM1252282', 'GSM1252290', 'GSM1252316', #gastric
                            'GSM3149068', 'GSM3149076', 'GSM3149084', 'GSM3149092', #breast
                            'GSM2870620', 'GSM2870626', 'GSM2870638', 'GSM2870644', #osteosarcoma
                            'GSM3149096', 'GSM3149100', 'GSM3149104') #endometrial

samples_table=samples_table[samples_table$Treatment!="GSE75898"  ,]
sampleSheet=filter(samples_table, Condition %in% target | SampleID %in% samplesID) 
# -

# ### Prepare Samples info file

###prepare a table with all the info
samples=sampleSheet[,c('SampleID','Tissue','Factor','Replicate','Condition','Treatment','ControlID')]
samples=droplevels(samples)
samples$group=factor(samples$Treatment, levels=c('Encode', 'Roadmap', 'GSE114737', 'GSE51776', 'GSE74230', 'pagani'))

# ### Create dba

###create a dba object with the 33131 regions (epigenome of PDOs and normal colon tissues)
dba <- dba(sampleSheet=sampleSheet,
         minOverlap = 1, 
         config = data.frame(RunParallel = TRUE,
                             reportInit = "DBA",
                             DataType = DBA_DATA_FRAME))
dba
#33131 sites in matrix 

# ### Count reads
#
# ##### DBA_SCORE_TMM_MINUS_FULL
# TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size

# +
dba=dba.count(dba, score=DBA_SCORE_TMM_MINUS_FULL,  bParallel=T)
dba_READS=dba.count(dba, score=DBA_SCORE_RPKM_FOLD, peak=NULL, bParallel=T)

dba=dba_READS
# -

# # Save

# +
path="/public_primary_tissue/diffbind/"
name= "H3K27Ac_primary_tissue_Input"

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
path="/public_primary_tissue/diffbind"
name= "H3K27Ac_primary_tissue_Input"

dba=dba.load(file=name, path, pre='dba_', ext='RData')
# -

# ### Extract rownames of the 33K regions

reg_names=read.table("/DEA/DBA_Korg_Ntissue_Input_Jup/dba_Korg_Ntissue_primary_tissues.bed", header=T)

# ### Create a dds from dba

sset <- dba(dba, bSummarizedExperiment=TRUE)
rownames(sset)=reg_names$reg_id
se=sset

# ### Remove batch effect

# #### Combat

library(stringr)
samples_d=as.data.frame(samples)
samples_d=cbind(as.data.frame(str_split_fixed(samples_d$Condition, "_", 2)), samples_d)
colnames(samples_d)[1]="type"
colnames(samples_d)[2]="tissue_orig"

mod = model.matrix(~as.factor(type), data=samples_d)

scoresCombat_Ref_pagani=ComBat(dat=log2(assays(se)$scores+1), batch=se$Treatment, BPPARAM=SerialParam(), mod=mod, par.prior=TRUE,  ref.batch="pagani") #
assays(se)$scoresCombat_Ref_pagani <- scoresCombat_Ref_pagani

# ### Heatmap

annot195reg=read.table('/TCGA_ATACseq/TCGA_PanCancer_peakSet/1accesibleRegionxEnhancer_TCGAseg_SQ_2358_w2000.bed', head=T)

#Import the order of the heatmap in Fig5 (see  bulk_chipseq/scripts/analyses/TCGA_ATACseq/ATAC-seq.Rmd
df_ord=read.table('/recurrence_10-8_K_org_UP_log2FC2_padj0.01_ATAC_195reg_ref.bed', header=F)
head(df_ord)

# +
#Prepare annotations
samples_anno=as.data.frame(samples)
rownames(samples_anno)=samples_anno$SampleID
samples_anno=samples_anno[,c('Condition','Treatment')]

df_ord_anno=as.data.frame(df_ord)
rownames(df_ord_anno)=df_ord_anno$V4
df_ord_anno=as.data.frame(df_ord_anno[,c('V4','V6')])
colnames(df_ord_anno)=c('enh_reg_id', 'cluster_all')
# -

colNames=c('SQ_2298', 'SQ_2145', 'GSM2058021', 'GSM2058022', 'GSM2058023', #our normal colon
                            'GSM1252278', 'GSM1252286', 'GSM1252294', 'GSM1252312', 'GSM1252320', #gastric
                            'GSM1127173', 'GSM2534707', #liver tissue roadmap
                            'GSM896163', 'GSM1013126', #adrenal gland tissue roadmap
                            'GSM2698612', 'GSM2527562', #adrenal gland encode        
                            'GSM1013129', 'GSM906397' ,  #pancreas tissue roadmap
                            'GSM2700597', 'GSM2700502', 'GSM4247380', 'GSM2534105', #body of pancreas encode 
                            'GSM4247568', 'GSM2700584', 'GSM2527594', 'GSM2534198', #thiroid encode
                            'GSM2698835', 'GSM4247275', #uterus encode
                            #tumor
                            'GSM1252274', 'GSM1252282', 'GSM1252290', 'GSM1252316', #gastric
                            'GSM3149068', 'GSM3149076', 'GSM3149084', 'GSM3149092', #breast
                            'GSM2870620', 'GSM2870626', 'GSM2870638', 'GSM2870644', #osteosarcoma
                            'GSM3149096', 'GSM3149100', 'GSM3149104') #endometirum

# +
#Extract the 195 regions from the total matrix
mtx_2use=assays(se)$scoresCombat_Ref_pagani

counts_n_reg195=mtx_2use[rownames(mtx_2use) %in% annot195reg$enh_reg_id,]
# -

### order the table with the column order we decided above and the rows according to the TCGA-ATAC-seq clustering in fig5
counts_n_reg195_ord=counts_n_reg195[order(match(rownames(counts_n_reg195),df_ord$V4)),]
counts_n_reg195_ord=counts_n_reg195_ord[,order(match(colnames(counts_n_reg195_ord),colNames))]

# ### Aggregate information

##aggregate clusters - Put together cluster 5 and 3 (46 regions) "conserved" and all the others "not conserved"
df_ord_anno$cluster_all_aggr=paste0('clust_', df_ord_anno$cluster_all)
df_ord_anno$cluster_all_aggr=gsub('clust_5', 'clust_cons', df_ord_anno$cluster_all_aggr)
df_ord_anno$cluster_all_aggr=gsub('clust_3', 'clust_cons', df_ord_anno$cluster_all_aggr)
df_ord_anno$cluster_all_aggr=gsub('clust_1', 'clust_notcons', df_ord_anno$cluster_all_aggr)
df_ord_anno$cluster_all_aggr=gsub('clust_2', 'clust_notcons', df_ord_anno$cluster_all_aggr)
df_ord_anno$cluster_all_aggr=gsub('clust_4', 'clust_notcons', df_ord_anno$cluster_all_aggr)
head(df_ord_anno)

###create a tidy object melting the matrix with the samples information
counts_n_melt=cbind(gene_id=rownames(counts_n_reg195_ord), counts_n_reg195_ord)
counts_n_melt=melt(as.data.frame(counts_n_melt), id.vars="gene_id", variable.name="sample")
counts_n_melt=merge(counts_n_melt, samples_d, by.x="sample", by.y="SampleID")
counts_n_melt=merge(counts_n_melt, df_ord_anno, by.x="gene_id", by.y="enh_reg_id")
counts_n_melt$cluster_all=paste0('clust_', counts_n_melt$cluster_all)
head(counts_n_melt)

# ### Boxplot with N and K divided by tissue

# +
## 46 regions
counts_n_melt_clust=counts_n_melt[counts_n_melt$cluster_all=="clust_5" | counts_n_melt$cluster_all=="clust_3" ,]

dfnorm_nk=counts_n_melt_clust %>%
group_by( type, Condition, sample) %>% 
summarise(type_orig_mean=mean(as.numeric(value))) 

dfnorm_nk$Condition<-factor(dfnorm_nk$Condition, levels=c('tumor_endometrium', 'tumor_osteosarcoma', 'tumor_breast','tumor_gastric', 'normal_colon', 'normal_gastric', 'normal_adrenal_gland', 'normal_thyroid_gland', 'normal_pancreas',  'normal_liver',  
                 'normal_uterus'))

p1=ggplot(dfnorm_nk, aes(x=Condition, y=type_orig_mean, fill=type)) + 
         geom_boxplot(alpha=0.6,outlier.size = 0.7) +
         ggtitle('TAZ-bound enhancers') +
         ylim(0.95, 1.45) +
         scale_fill_manual(values=c('#475D25','#9F68A9')) +
         theme_classic() +
         theme(axis.title=element_text(size=14), axis.text = element_text(size=12), 
               axis.text.x = element_text(angle = 90, vjust=0.8, hjust = 1), 
               legend.title = element_blank()) 
p1

##Compute two-samples Wilcoxon test 
wilcox.test(type_orig_mean ~ type, data = dfnorm_nk,
                   exact = FALSE)

print("wilcoxon")
compare_means(type_orig_mean ~ type, data = dfnorm_nk,
                   exact = FALSE)
# -


# <!-- # K organoids -->
# <!-- ```{r} -->
# <!-- dds=readRDS(file = "/mnt/rnd/projects/organoids/scratch/bulk_chipseq/analyses/diffbind/chromHMM_enhancer_org_pA/DEA/DBA_Korg_Ntissue_Input/DESeq2/dds_Korg_Ncrypts.rds") -->
#
# <!-- dds_Korg_norm=as.data.frame(counts(dds, normalized=TRUE)) -->
# <!-- dds_Korg_norm$meanKorg=rowMeans(dds_Korg_norm[,1:10]) -->
# <!-- dds_Korg_norm$meanNtissue=rowMeans(dds_Korg_norm[,11:13]) -->
#
# <!-- dds_Korg_norm$FC=dds_Korg_norm$meanKorg/dds_Korg_norm$meanNtissue -->
#
# <!-- # Add identifier -->
# <!-- # Create an annotation to get the GRANGE of the peaks and a unique identifier -->
# <!-- annotation_GRANGE=dba_rg -->
# <!-- annotation_GRANGE$region_ID <- paste0("reg_", (1:nrow(annotation_GRANGE))) -->
# <!-- rownames(annotation_GRANGE)=annotation_GRANGE$region_ID -->
#
# <!-- # Extract, annotate and save results for all genes  -->
# <!-- dds_Korg_norm$region_ID =paste0("reg_", (1:nrow(dds_Korg_norm))) #add id to each regions -->
# <!-- rownames(dds_Korg_norm)=dds_Korg_norm$region_ID -->
#
# <!-- dds_Korg_norm=merge(dds_Korg_norm, annotation_GRANGE, by="region_ID", sort = F) -->
#
# <!-- # Divide dataset in decile -->
# <!-- dds_Korg_norm= dds_Korg_norm %>% -->
# <!--      mutate(quantile = ntile(meanKorg, 10)) -->
#
# <!-- # Extract only the 10 -->
# <!-- top_basemean=as.data.frame(dds_Korg_norm[dds_Korg_norm$quantile==10, ]) -->
#
# <!-- # Extract the regions with a FOLD CHANGE between 1 and 2 -->
# <!-- top_basemean_select=top_basemean[which(top_basemean$FC>1 & top_basemean$FC<2),] -->

# <!-- top_basemean_select_org.cr= merge(top_basemean_select, org.cr, by="region_ID", sort = F) -->
# <!-- ``` -->
