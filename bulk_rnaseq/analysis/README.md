## Instruction to install the conda environment used to run bulk RNA-seq downstream analysis.

**Install R 3.5**

```
# Activate the environment
conda activate R35_RNAseq

# Install the basic R packages
conda install -y -c r r-base==3.5.1 r-irkernel==0.8.12
conda install -y -c r r-devtools==1.13.6

R --vanilla <<code
chooseCRANmirror(graphics=FALSE, ind=88) #Italy GARR
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install()

# Install this environment as kernel in Jupyter
IRkernel::installspec(name = 'R35_RNAseq', displayname = 'R35_RNAseq')
code
```

**Install main packages for the analysis**

```
R --vanilla <<code
#BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#requirement for Deseq2
install.packages("RColorBrewer")

packageurl <- "http://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

#DESeq2
BiocManager::install("DESeq2", version = "3.8")

#dplyr
install.packages("dplyr")

#Reshape
install.packages("reshape")

#Pheatmap
install.packages("pheatmap")

#gdata
install.packages("gdata")

#tidyr
install.packages("tidyr")

#VennDiagram
install.packages("VennDiagram")

#ggrepel
install.packages("ggrepel")

#ggpubr
packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-7.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

#ggpubr
install.packages("ggpubr")

#ggcorrplot
install.packages("ggcorrplot")

#TCGAbiolink
BiocManager::install("TCGAbiolinks")

q()
```

`sessionInfo()` after `library(DESeq2)`
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Ubuntu 18.04.1 LTS

Matrix products: default
BLAS/LAPACK: ~/miniconda3univa/envs/R35_RNAseq/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] DESeq2_1.22.2               SummarizedExperiment_1.12.0
 [3] DelayedArray_0.8.0          BiocParallel_1.16.5
 [5] matrixStats_0.54.0          Biobase_2.42.0
 [7] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1
 [9] IRanges_2.16.0              S4Vectors_0.20.1
[11] BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] bit64_0.9-7            splines_3.5.1          Formula_1.2-3
 [4] assertthat_0.2.0       latticeExtra_0.6-28    blob_1.1.1
 [7] GenomeInfoDbData_1.2.0 pillar_1.3.1           RSQLite_2.1.1
[10] backports_1.1.3        lattice_0.20-38        glue_1.3.0
[13] digest_0.6.18          RColorBrewer_1.1-2     XVector_0.22.0
[16] checkmate_1.8.5        colorspace_1.3-2       htmltools_0.3.6
[19] Matrix_1.2-15          plyr_1.8.4             XML_3.98-1.16
[22] pkgconfig_2.0.2        genefilter_1.64.0      zlibbioc_1.28.0
[25] purrr_0.2.5            xtable_1.8-3           scales_1.0.0
[28] htmlTable_1.13         tibble_2.0.0           annotate_1.60.0
[31] ggplot2_3.1.0          nnet_7.3-12            lazyeval_0.2.1
[34] survival_2.43-3        magrittr_1.5           crayon_1.3.4
[37] memoise_1.1.0          foreign_0.8-71         tools_3.5.1
[40] data.table_1.11.8      stringr_1.3.1          locfit_1.5-9.1
[43] munsell_0.5.0          cluster_2.0.7-1        AnnotationDbi_1.44.0
[46] bindrcpp_0.2.2         compiler_3.5.1         rlang_0.3.0.1
[49] grid_3.5.1             RCurl_1.95-4.11        rstudioapi_0.8
[52] htmlwidgets_1.3        bitops_1.0-6           base64enc_0.1-3
[55] gtable_0.2.0           DBI_1.0.0              R6_2.3.0
[58] gridExtra_2.3          knitr_1.21             dplyr_0.7.8
[61] bit_1.1-14             bindr_0.1.1            Hmisc_4.1-1
[64] stringi_1.2.4          Rcpp_1.0.0             geneplotter_1.60.0
[67] rpart_4.1-13           acepack_1.4.1          tidyselect_0.2.5
[70] xfun_0.4
```

