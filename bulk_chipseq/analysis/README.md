### Instruction to recreate the conda environment used to run bulk ChIP-seq downstream analysis.

#### ChIP-seq downstream analysis (R)

**Install R**

```
# Create a new conda environment
conda create -n R34_ChIPseq_v1
# Activate the environment
conda activate R34_ChIPseq_v1

conda install -y -c r r-base==3.4.3 r-irkernel==0.8.11
conda install -c conda-forge r-devtools==1.13.4
conda install -c conda-forge udunits2==2.2.27.6
conda install -c conda-forge zlib==1.2.11

git clone https://github.com/ropensci/git2r.git

R CMD INSTALL git2r --configure-vars='LIBS=-L~/miniconda3univa/envs/R34_ChIPseq_v1/lib CPPFLAGS=-I~/miniconda3univa/envs/R34_ChIPseq_v1/lib'
conda install -c bioconda mysqlclient 
```

**Install main packages for the analysis**

```
R --vanilla <<code
chooseCRANmirror(graphics=FALSE, ind=41) #Italy GARR
# Install Bioconductor http://bioconductor.org/install/
source("https://bioconductor.org/biocLite.R")
# Install this environment as kernel in Jupyter
IRkernel::installspec(name = 'R34_ChIPseq_v1', displayname = 'R34_ChIPseq_v1')

#Install packages

##The following are requirements for DiffBind (You might not need this)
install.packages("lattice")
install.packages("RColorBrewer")
install.packages("bitops")
packageurl <- "http://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-71.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "http://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "http://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "http://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "http://cran.r-project.org/src/contrib/Archive/amap/amap_0.8-16.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
biocLite("ShortRead", configure.vars=paste0("LDFLAGS=-L", Sys.getenv("CONDA_PREFIX"),"/lib"))

biocLite("DiffBind")

##Annotation
biocLite("ChIPpeakAnno")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.7-4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.7-4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
biocLite("ChIPseeker")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
install.packages("tidyr")
biocLite("regioneR")

##Correlation matrix
biocLite("bioDist")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/hyperSpec/hyperSpec_0.99-20180627.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

#xlsx
install.packages("xlsx")

##stats
install.packages("minqa")
install.packages("nloptr")
install.packages("RcppEigen")
install.packages("quantreg")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/lme4/lme4_1.1-14.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.4-7.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <- "https://cran.r-project.org/src/contrib/Archive/car/car_2.1-3.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

packageurl <-"https://cran.r-project.org/src/contrib/Archive/cowplot/cowplot_0.9.4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

install.packages("ggpubr")

code
```

`sessionInfo()` after `library(DiffBind)`

```
R version 3.4.3 (2017-11-30)
Platform: x86_64-conda_cos6-linux-gnu (64-bit)
Running under: Ubuntu 18.04.1 LTS

Matrix products: default
BLAS/LAPACK: ~/miniconda3univa/envs/R34_ChIPseq_v1/lib/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] magrittr_1.5               genefilter_1.60.0         
 [3] pheatmap_1.0.12            RColorBrewer_1.1-2        
 [5] ggplot2_3.1.0              gdata_2.18.0              
 [7] stringr_1.4.0              tidyr_0.8.3               
 [9] BiocParallel_1.12.0        DESeq2_1.18.1             
[11] data.table_1.12.0          reshape2_1.4.3            
[13] knitr_1.21                 limma_3.34.9              
[15] GenomicFeatures_1.30.3     AnnotationDbi_1.40.0      
[17] ChIPpeakAnno_3.12.7        VennDiagram_1.6.20        
[19] futile.logger_1.4.3        Biostrings_2.46.0         
[21] XVector_0.18.0             gplots_3.0.1.1            
[23] ChIPseeker_1.14.2          DiffBind_2.6.6            
[25] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
[27] matrixStats_0.54.0         Biobase_2.38.0            
[29] GenomicRanges_1.30.3       GenomeInfoDb_1.14.0       
[31] IRanges_2.12.0             S4Vectors_0.16.0          
[33] BiocGenerics_0.24.0        BiocInstaller_1.28.0  
```

#### ChIP-seq downstream analysis (python)

```
# Create a new conda environment
conda create -n design_plot
# Activate the environment
conda activate design_plot

# Install python and the basic packages
conda install -y ipython ipykernel
conda install -y seaborn=0.10.0 scikit-learn=0.22.1 numba=0.48.0 numpy=1.18.1 scipy=1.4.1 matplotlib=3.1.3 pandas=1.0.2 cython=0.29.15

# Create kernell link
ipython kernel install --user --name design_plot --display-name "design_plot"
```
