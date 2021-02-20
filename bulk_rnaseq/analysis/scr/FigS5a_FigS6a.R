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

# ### Analysis using COAD-TCGA data
#
# Author: Federica Gervasoni and Raoul Bonnal

# +
library("CMSclassifier")
library("biomaRt")
library("pheatmap")
library("genefilter")
library("TCGAbiolinks")
library("dplyr")
library("DESeq2")
library("data.table")
library("ggplot2")
library("ggpubr")

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
library(magrittr)
library(dplyr)
library(VennDiagram)
library(dplyr)
# -

name="coad_tcga"

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

# ### Get read counts

# Import Raoul Bonnal script with formatting functions
source("/datasets/TCGA_data/TCGABiolinks_fixes.R")

wd <- "/TCGA_data"
datadir <- paste(wd, "Data", sep='/')

# +
coad.gbm <- GDCquery(project =  "TCGA-COAD",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "HTSeq - Counts")

coad.clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
coad.biospecimen <- GDCquery_clinic(project = "TCGA-COAD", type = "biospecimen")
# -

# INCASE OF DOWNLOAD ROW DATA: START
GDCdownload(coad.gbm, method = "api", directory = paste0(datadir, 'TCGA', sep='/'), files.per.chunk = 50)
# INCASE OF DOWNLOAD ROW DATA: END

coad.data <- GDCprepare(coad.gbm, directory = paste0(datadir, 'TCGA', sep='/'))
# dim: 6 521

# Data Cleanup and split N vs T
coad.data<-coad.data[, !coad.data$is_ffpe]
# dim: 6 508

# +
file.coad = paste0(datadir, 'TCGA/', "TCGA-COAD/", "coad_data", ".rds")

if (!file.exists(file.coad)){

    saveRDS(object = coad.data, file =file.coad)
    
} else {

    coad.data=readRDS(file = file.coad)
}

# +
# Split NT vs TP data and Remove duplicated
## NB: Using Raoul's function, NOT a function from TCGAbiolinks package
coad.data.tp=TCGAanalyze_RemoveReplicateSamples(coad.data, "TP")

coad.data.nt=TCGAanalyze_RemoveReplicateSamples(coad.data, "NT")
# -

counts = cbind(assay(coad.data.tp),assay(coad.data.nt))

# ### Sample annotation

# +
# Prepare N and K annotation
samplesNT <- TCGAquery_SampleTypes(colnames(counts), typesample = c("NT")) # Returning a vector of length(0)
samplesTP <- TCGAquery_SampleTypes(colnames(counts), typesample = c("TP"))

samples=as.data.frame(c(samplesTP, samplesNT))
samples$type=c(rep("TP", length(samplesTP)), rep("NT", length(samplesNT)))
colnames(samples)=c("sample", "type")
rownames(samples)=samples$sample
samples$group=samples$type

samples$group=factor(samples$group, levels=c("NT", "TP"))

file.dds = paste0(datadir, 'TCGA/', "TCGA-COAD/", "dds.", name, ".noMT.rds")
file.rld = paste0(datadir, 'TCGA/', "TCGA-COAD/", "rld.", name, ".noMT.rds")
# -

# ### Read raw counts and remove MT and Y genes

# +
# Read raw counts
rownames(counts)=gsub("[.].*", "", rownames(counts))

counts.noMT=counts[!rownames(counts)%in%MT.genes & !rownames(counts)%in%Y.genes,]

identical(rownames(samples), colnames(counts.noMT))
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
    
    dds <- DESeq(dds)

    saveRDS(object = dds, file = file.dds)
    
} else {

    dds=readRDS(file = file.dds)
    
}

dds
dds@design
colData(dds)
table(dds$group)

# +
colData(dds)$patient=c(coad.data.tp$patient, coad.data.nt$patient)
colData(dds)$ajcc_pathologic_stage=c(coad.data.tp$ajcc_pathologic_stage, coad.data.nt$ajcc_pathologic_stage)
colData(dds)$tumor_stage=c(coad.data.tp$tumor_stage, coad.data.nt$tumor_stage)

#tumor stage
colData(dds)$simplify_tumor_stage=c(SimplifyTumorStage(coad.data.tp$tumor_stage), SimplifyTumorStage(coad.data.nt$tumor_stage))
colData(dds)$simplify_tumor_stage[is.na(colData(dds)$simplify_tumor_stage)] <- "not reported"

# +
##Count mtx
norm.all=counts(dds, normalized=TRUE)

norm.all.log2=log2(norm.all + 1)
# -

# ### Deconvolution of CRC TCGA data using scRNA-seq from 23 CRC samples

## BisqueRNA package https://cran.r-project.org/web/packages/BisqueRNA/vignettes/bisque.html
library(BisqueRNA)

# ### Bulk RNA-seq

##extract count mtx for TCGA data
mtx_raw=counts(dds, normalized=TRUE)
head(mtx_raw)

phenoData <- new("AnnotatedDataFrame", data=as.data.frame(colData(dds)))

## create eset for bulk RNA-seq
bulk.eset <- ExpressionSet(assayData=mtx_raw, phenoData=phenoData)

# ### scRNA-seq data from the 23 patients (Lee et al.)

path_dec='/sc_RNAseq_data/'

# +
features <- data.table::fread(paste0(path_dec, 'features.tsv'), header=F)

barcodes <- data.table::fread(paste0(path_dec, 'barcodes_obs.tsv'), header=T)
# -

matrix_file = paste0(path_dec, 'raw_count_mtx.mtx')

# Load dgTmatrix
matrix <- Matrix::readMM(gzfile(matrix_file))

# +
# add dimnames #
colnames(matrix) <- barcodes$index
rownames(matrix) <- features$V1

matrix <- as(matrix, 'matrix')
print(dim(matrix))
idx = rowSums(matrix)>0
matrix <- matrix[idx,]
print(dim(matrix))
# -

rownames(barcodes)=barcodes$index
colnames(barcodes)=c('index', 'cellType', 'Sample', 'SubjectName')
barcodes=barcodes[,c('index', 'cellType', 'SubjectName')]

# +
sample.ids <- colnames(matrix)
# individual.ids and cell.types should be in the same order as in sample.ids
sc.pheno <- data.frame(check.names=F, check.rows=F,
                       stringsAsFactors=F,
                       row.names=sample.ids,
                       SubjectName=barcodes$SubjectName,
                       cellType=barcodes$cellType)

sc.meta <- data.frame(labelDescription=c("SubjectName",
                                         "cellType"),
                      row.names=c("SubjectName",
                                  "cellType"))

sc.pdata <- new("AnnotatedDataFrame",
                data=sc.pheno,
                varMetadata=sc.meta)
# -

## create eset for scRNA-seq with Samples and Celltype information
sc.eset <- Biobase::ExpressionSet(assayData=matrix,
                                  phenoData=sc.pdata)

###  Run the deconvolution (marker=NULL so it can keep all the genes and use.overlap=FALSE because we don't have matched samples)
res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)

### Calculate a coefficient to get the results for epithelial cells: Fraction of epithelial/(Fraction of epithelial + Fraction of stromal)
dec_melt=as.data.frame(t(as.data.frame(res$bulk.props[c('Epithelial cells'),])))
rownames(dec_melt)=c('epith')
dec_melt$gene_id=row.names(dec_melt)
dec_melt=setDT(dec_melt)
dec_melt=melt(dec_melt, id.vars="gene_id", variable.name="sample")
dec_melt=as.data.frame(merge(dec_melt, as.data.frame(colData(dds)), by="sample"))
rownames(dec_melt)=dec_melt$sample
head(dec_melt)

# ### Fig S5a - YAP/TAZ levels in COAD

# +
gene_list=c( 'YAP1','WWTR1')

for (i in gene_list){
    print(i)

    ## Calculate mean value of gene level in all the dataset 
    gene_list=annotation[annotation$gene_name==i , "gene_id"]

    norm.all.log2_merged=as.data.frame(t(norm.all.log2[match(gene_list, row.names(norm.all.log2)), ]))
    
    ##Level of gene multiplied for the coefficient of epithelial cells
    norm.all.log2_merged_norm=as.data.table(plyr::aaply(norm.all.log2_merged, 1, "*", dec_melt$value))
    
    rownames(norm.all.log2_merged_norm)=i

    norm.all.log2_melt=as.data.frame(norm.all.log2_merged_norm)

    norm.all.log2_melt$gene_id=i
    norm.all.log2_melt=setDT(norm.all.log2_melt)
    norm.all.log2_melt=melt(norm.all.log2_melt, id.vars="gene_id", variable.name="sample")
    norm.all.log2_melt=merge(norm.all.log2_melt, as.data.frame(colData(dds)), by="sample")
    head(norm.all.log2_melt)

    p <- ggplot(norm.all.log2_melt, aes(x=type, y=as.numeric(value), color=type, fill=type)) + 
      labs(x="", y="log2normCounts adj freq Epith") +
      geom_violin(aes(alpha=0.5), trim=FALSE, scale='width')+
      ggtitle(i)+
      scale_fill_manual(values=c("#033E8C", "#AF0416")) +
      scale_color_manual(values=c("#033E8C", "#AF0416")) +
      ggtitle(i)+
      theme_classic() +
      ylim(1,5) +
      theme(axis.title=element_text(size=14), axis.text = element_text(size=12), axis.text.x = element_text(angle = 60, vjust=0.8, hjust = 1), legend.title = element_blank()) 

    #Statistics
    my_comparisons <- list( c("TP", "NT"))
    print(p + stat_compare_means(comparisons = my_comparisons, label.y = c(4.7,4.7), method = "wilcox.test", label = "p.signif", paired = FALSE,
                                method.args = list(alternative = "greater")) )

    #stat
    print("wilcoxon")
    print(compare_means(value ~ type,  data = norm.all.log2_melt, method = "wilcox.test", p.adjust.method = "fdr", method.args = list(alternative = "greater")))
}
# -

# ### FigS6a - Genes associated to DE enhancer + UP regulated in PDO vs Ntissue in COAD dataset

genes_de=read.table('/DESeq/DGE/KorgvsNcr/KorgvsNcr.matched.ChIP.Dgroup_1organoid.patient/intersection_ChIPseq_UP_genes_POD_Ntissue.txt',
                    sep='\t', header=T)
genes_de=merge(genes_de, annotation[,c("gene_name", "gene_id")], by.x="x", by.y="gene_name")

# +
gene_list=as.character(genes_de$gene_id)

#Adjust normalized count for deconvolution res
norm.all.log2_dec=norm.all.log2

##Normalize the whole matrix
norm.all.log2_dec=as.data.frame(plyr::aaply(norm.all.log2_dec, 1, "*", dec_melt$value))  
rownames(norm.all.log2_dec)=rownames(norm.all.log2)

## Subset the log2 normalised dds.all dataset for the common genes
norm.all.log2.con=norm.all.log2_dec[match(gene_list, row.names(norm.all.log2_dec)), ]
norm.all.log2.con=norm.all.log2.con[apply(norm.all.log2.con[,-1], 1, function(x) !all(x==0)),]

table(row.names(norm.all.log2.con)%in%gene_list)

#Annotable
anno <- as.data.frame(colData(dds)[, c("group","simplify_tumor_stage")])
anno$stage_group=anno$simplify_tumor_stage
   
anno$stage_group=ifelse(anno$group=='NT', "NT", anno$stage_group)                                         
                                          
# ### reorder columns according to groups
anno=anno[order(anno$stage_group, decreasing = TRUE),]
anno_columns=row.names(anno)
                                          
# ### Assign specific colors
ann_colors = list(condition = c(TP="#9b1477", NT="#6A6B06"),
                  simplify_tumor_stage=c(i="#fdff00", ii="#ffdb00", iii="#ff5a00", iv="#ff1a00", "not reported"="grey"),
                  stage_group=c(i="#fdff00", ii="#ffdb00", iii="#ff5a00", iv="#ff1a00", "not reported"="grey", "NT"="white")
                 )

norm.all.log2.con=norm.all.log2.con[,anno_columns]
                                          
matrix_rownames=annotation[annotation$gene_id %in% gene_list, c("gene_name", "gene_id")]

rownames(matrix_rownames)=matrix_rownames$gene_id
                                         
norm.all.log2.con=as.matrix(norm.all.log2.con[match(rownames(matrix_rownames), rownames(norm.all.log2.con)),])                                                                                                                 

pdf=pheatmap(norm.all.log2.con, scale="row",
         cluster_cols=F,
         cluster_rows=F,         
         annotation_col = anno[,c("group","stage_group")],
         annotation_colors = ann_colors,
         show_rownames=T,
         show_colnames=FALSE,
         labels_row = matrix_rownames$gene_name,
         color = colorRampPalette(rev(c("#8A0807", "#D51214", "#FFFEFD", "#3E44CC", "#0A2A87")))(1250), 
         breaks = c(seq(-5,-2.60,length=250), seq(-1,-0.1,length=250),seq(0,0.1,length=250),seq(0.3,1.2,length=250), seq(2.9,4,length=250)),
         fontsize=8)                                       
