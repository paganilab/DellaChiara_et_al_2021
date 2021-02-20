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

[Fig2. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2a.pbs): Histone modifications localization in relation to the gene body as well as regions surrounding +/- 3 Kb of the transcription start (TSS) and end (TES) sites. Representative density plots of average intensity (top) and corresponding heatmaps (bottom) display the relative distribution of H3K4me3 (red), H3K27ac (pink), H3K4me1 (yellow), H3K36me3 (green), and H3K27me3 (grey) signals for all the genes present in the GENCODEv25 annotation. 

[Fig2. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2b.R): Pearson correlation heatmap of ChIP- seq data for the complete set of five histone modifications across all PDOs

[Fig2. c]:Combinatorial pattern of histone marks in an 8-state model using ChromHMM. The heatmap (Emission plot) displays the frequency of the histone modifications 10 found in each state. 

Fig2. d: The probability of each ChromHMM- defined chromatin state overlapping ATAC-seq regions for TCGA colon adenocarcinoma samples is shown across PDOs using a spider plot. 

Fig2. e-f see Hepic website: Representative tracks of ChromHMM states for the FABP1 and LAMA5 genomic loci in all PDOs. The expanded regions show H3K4me3, H3K27ac, H3K4me1, H3K36me3 and H3K27me3 profiles, along with RNA-seq signal and ChromHMM states for PDOs of different molecular subtypes as indicated.

**Figura 4**

[Fig4. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2a.pbs): Signal density plot (top) and corresponding heatmaps (bottom) displaying the relative distribution of TAZ peaks around ChromHMM-defined active enhancers (n = 33,131) and promoters. Data for enhancers and promoters are depicted in yellow and red, respectively.

-----

**Supplementary Figure 1**

[Supplementary Fig1. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1a.R): MA plot of log2 mean gene expression over log2 fold-change showing the lack of differentially expressed genes between early and late passages of organoids.

[Supplementary Fig1. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_venny.R): Concordance of expressed genes detected in PDOs and corresponding tumors. Bar graph represents the proportion of expressed genes (gene count > 5) that is shared between each PDO and its corresponding tumor, and those detected only in the PDO or parental tumor.

[Supplementary Fig1. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1c.py): Genes expressed in primary tumors but not in PDOs (Fig. 1e, Venn diagram, n = 3,412), are enriched for gene signatures of stromal cells.

[Supplementary Fig1. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_corr.R): Correlation of gene expression between matching pairs of tumor tissue and derived PDO.

**Supplementary Figure 2**

[Supplementary Fig2. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2b.R): PCA of input normalized ChIP-seq signals for the five histone modifications used to build the ChromHMM 8- state model.

Supplementary Fig2. b: Heatmaps showing the annotation of the ChromHMM 8-states with known genomic features (Overlap) and the probability that a state is found in the proximity of another state (Transition).

Supplementary Fig2. c: Average proportion of each chromatin state over all PDOs. The chromatin segments for active/flanking TSS and active/flanking enhancer states are merged into the promoter and enhancer functional elements, respectively.

Supplementary Fig2. d:  Distribution of the eight ChromHMM 4 states for each PDO.

**Supplementary Figure 5**

[Supplementary Fig5. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS5a_FigS6a.R): YAP and TAZ are transcriptionally upregulated in primary tumors compared to normal colon tissues in the TCGA colon adenocarcinoma dataset. Violin plots show the distribution of RNA-seq log2 normalized gene counts adjusted for epithelial cell frequency (see Methods for details). P < 0.0001, Wilcoxon rank sum test. 


**Supplementary Figure 6**

[Supplementary Fig6. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS5a_FigS6a.R): Heatmap of RNA-seq log2 normalized counts for target genes of YAP/TAZ-controlled enhancers in the TCGA COAD dataset. Tissue populations and clinical stages are represented by color-coded bars above the heatmap. Expression values of TCGA bulk RNA-seq data were adjusted for epithelial cell frequency (see Methods for details).


