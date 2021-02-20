# DellaChiara_et_al_2021

In this repository you can find the code to reproduce the work and figures of **“Epigenetic landscape of human colorectal cancer unveils an aberrant core of pan-cancer enhancers orchestrated by YAP/TAZ”**

-----

The code is divided based on technology used.

To see the code for:

- [Bulk RNA-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/bulk_rnaseq/)

- [Bulk ChIP-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/bulk_chipseq/)

- [Single cell RNA-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/sc_rnaseq/)

Inside each folder you can find the *code* and the instructions to reproduce the *environments* ([Docker](https://www.docker.com/) or [Conda](https://docs.conda.io/en/latest/)) used to run the analyses.

The main folder is subdivided in:
- **nf-pipeline**: contains the custom pipelines of RNA-seq and ChIP-seq managed by [Nextflow](https://www.nextflow.io/).
- **analysis**: contains the scripts of the downstream analysis and the code to produce the figures.

To run the downstreamn analysis we used [JupyterLab](https://jupyter.org/). If you want to convert the scripts and run them as notebook using Jupyter, please see [Jupytext](https://github.com/mwouts/jupytext)

-----

Below a more detailed link to the script for each figure present in the paper (figure legend descriptions from the paper are reported for a better understanding).

Main figures are reported first, then the supplementary figures.

**Figure 1**

[Fig. 1d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1d.R): PCA on normalised gene counts from RNA-seq data distinguished normal colon mucosa, primary tumor and PDOs.

[Fig. 1e - Venny](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_venny.R): Venn diagram showing the number of concordant expressed genes between tumors and PDOs. 

[Fig 1e - Correlation](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_corr.R): Mean log2 normalized gene counts between primary tumors and PDOs were well correlated. 

[Fig. 1f](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1f.R): Hierarchical clustering analysis using differentially expressed genes (DEG, adjusted P-value ≤ 0.01) between tumor and normal colon tissues clustered PDOs together with parental tumors. Tissue populations and patients are represented by color-coded bars above the heatmap.

[Fig. 1g](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1g_GSEA.py): PDOs are enriched in gene signatures of CRC clinical specimen. GSEA on the ranked list of  genes from the [comparison between PDOs and normal colon tissue](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1g_DE.R) with the normalized enrichment score (NES) and P-value reported.

**Figure 2**
Fig2. a:

Fig2. b:

Fig2. c:

Fig2. d:

Fig2. e-f: See Hepic website

-----

**Supplementary Figure 1**

[Supplementary Fig1. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1a.R): MA plot of log2 mean gene expression over log2 fold-change showing the lack of differentially expressed genes between early and late passages of organoids.

[Supplementary Fig1. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_venny.R): Concordance of expressed genes detected in PDOs and corresponding tumors. Bar graph represents the proportion of expressed genes (gene count > 5) that is shared between each PDO and its corresponding tumor, and those detected only in the PDO or parental tumor.

[Supplementary Fig1. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1c.py): Genes expressed in primary tumors but not in PDOs (Fig. 1e, Venn diagram, n = 3,412), are enriched for gene signatures of stromal cells.

[Supplementary Fig1. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_corr.R): Correlation of gene expression between matching pairs of tumor tissue and derived PDO.

**Supplementary Figure 2**

FigS2. a:

FigS2. b:

FigS2. c:

FigS2. d:
