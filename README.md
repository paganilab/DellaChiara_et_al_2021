# DellaChiara_et_al_2021

In this repository you can find the code to reproduce the work and figures of **“Epigenetic landscape of human colorectal cancer unveils an aberrant core of pan-cancer enhancers orchestrated by YAP/TAZ”**

-----

The code is divided based on technology used.

To see the code for:

- [Bulk RNA-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/bulk_rnaseq/)

- [Bulk ChIP-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/bulk_chipseq/)

- [Single cell RNA-seq data](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/sc_RNAseq/)

Inside each folder you can find the *code* and the instructions to reproduce the *environments* ([Docker](https://www.docker.com/) or [Conda](https://docs.conda.io/en/latest/)) used to run the analyses.

The main folder is subdivided in:
- **nf-pipeline**: contains the custom pipelines of RNA-seq and ChIP-seq managed by [Nextflow](https://www.nextflow.io/).
- **analysis**: contains the scripts of the downstream analysis and the code to produce the figures.

To run the downstreamn analysis we used [JupyterLab](https://jupyter.org/). If you want to convert the scripts and run them as notebook using Jupyter, please see [Jupytext](https://github.com/mwouts/jupytext)

The **RNA-seq and ChIP-seq data** generated during this study are **available** at the European Nucleotide Archive under these links:

- [RNA-seq data](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8448/)

- [ChIP-seq data](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8416/)

All the data used to produced the genomic overview present in the paper are available at our genomic web browser:

> https://ifom.eu/bioinformatics/hepic

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

[Fig2. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2c.py):Combinatorial pattern of histone marks in an 8-state model using ChromHMM (see [nf-chromhmm pipeline](https://github.com/paganilab/DellaChiara_et_al_2021/tree/main/bulk_chipseq/nf-chromhmm)). The heatmap (Emission plot) displays the frequency of the histone modifications 10 found in each state. 

[Fig2. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2d.py): The probability of each ChromHMM- defined chromatin state overlapping ATAC-seq regions for TCGA colon adenocarcinoma samples is shown across PDOs using a spider plot. 

**Figure 3**

[Fig3. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3a_b_f.R): Unsupervised clustering analysis and pearson correlation heatmap of H3K27ac ChIP-seq data for the 33,131 ChromHMM-defined active enhancers clearly distinguish PDOs from normal colon tissues.

[Fig3. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3a_b_f.R): Volcano plot of differentially enriched enhancer regions between PDOs and normal colon mucosa. Dotted lines indicate thresholds for adjusted P-value < 0.01 and |log2 fold-change| > 2.

[Fig3. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): Percentage of gained enhancers shared by different PDOs.

[Fig3. f](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3a_b_f.R): The Hippo signaling pathway is the most significantly enriched pathway related to the gained enhancer-associated genes that are upregulated in PDOs compared to normal tissues. The size of the circles corresponds to the number of gained-enhancer associated genes present in the geneset of a particular pathway (Gene Ratio). The dotted line indicates the threshold for significantly enriched pathways (false discovery rate < 0.05).

**Figura 4**

[Fig4. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig4a.R): YAP and TAZ are transcriptionally upregulated in primary tumors and PDOs compared to normal colon tissues. Violin plots show the distribution of normalized gene counts for samples in the three groups. ** P < 0.01, *** P < 0.001, Wilcoxon rank sum test. 

[Fig4. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2a.pbs): Signal density plot (top) and corresponding heatmaps (bottom) displaying the relative distribution of TAZ peaks around ChromHMM-defined active enhancers (n = 33,131) and promoters. Data for enhancers and promoters are depicted in yellow and red, respectively.

[Fig4. f](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): Distribution of TAZ peaks across ChromHMM-defined functional elements for active and inactive genomic regions.

[Fig4. g](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): TAZ enrichment in human CRC gained enhancers (G.E.). The percentage of TAZ-bound enhancers increases with the level of conservation across PDOs. The barplots show the percentage of enhancers in each G.E. subset that overlap a TAZ peak or the percentage (mean ± s.d.) of TAZ-bound regions in 1000 random sets generated for each of the G.E. subsets: i) all G.E., ii) G.E. conserved in at least 5 patients, and iii) G.E. conserved in at least 8 patients (see Methods for details). P < 0.001, empirical p-value. 

**Figure 5**

[Fig5. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig5a.Rmd): Chromatin accessibility profiles of the 195 conserved gained enhancers in 23 diverse primary human cancer types reveals a signature of 46 pan-CRC enhancers with active chromatin profiles across cancer types (cluster in purple). The heatmap shows hierarchical clustering of log2 normalized insertion counts of ATAC-seq data derived from TCGA (Corces et al., 2018). Colon adenocarcinoma (COAD) samples are the first cancer type reported on the left of the heatmap.

[Fig5. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig5b.R): Differences in the H3K27ac intensities of the pan-cancer enhancers in primary tumors relative to normal tissues from public available ChIP-seq data.

**Figure 6**

[Fig6. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6.py): t-distributed stochastic neighbor embedding (t-SNE) visualization depicting the major cell types identified in scRNA-seq data of primary tumor and normal tissues from 23 CRC patients.

[Fig6. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6.py): t-SNE visualization of the CRB score across all cell populations. 

[Fig6. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6.py): Malignant cells display significantly higher CRB score. Distribution of the CRB score in the malignant (grey) and non-malignant epithelial clusters (green), and in all other major cell populations. Boxplots describe the median and interquartile range with whiskers denoting the 1.5 x interquartile range.  P < 0.0001, Wilcoxon rank sum test

[Fig6. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6-epith.py): t-SNE representation of the CRB scores across 18,539 epithelial cells. Contour lines denote normal epithelial cells. 

[Fig6. e](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6-epith.py): t-SNE representation of the stemness scores across 18,539 epithelial cells. Contour lines denote normal epithelial cells. 

[Fig6. f](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6-epith.py): The cancer regulatory blueprint does not relate to stemness. t-SNE representation of the difference between the CRB and stemness score across all epithelial cells. Cells depicted in grey display similar levels of the two scores, whereas cells in purple and green display mutually exclusive high levels of CRB or stemness, respectively. Contour lines denote normal epithelial cells.

-----

**Supplementary Figure 1**

[Supplementary Fig1. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1a.R): MA plot of log2 mean gene expression over log2 fold-change showing the lack of differentially expressed genes between early and late passages of organoids.

[Supplementary Fig1. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_venny.R): Concordance of expressed genes detected in PDOs and corresponding tumors. Bar graph represents the proportion of expressed genes (gene count > 5) that is shared between each PDO and its corresponding tumor, and those detected only in the PDO or parental tumor.

[Supplementary Fig1. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS1c.py): Genes expressed in primary tumors but not in PDOs (Fig. 1e, Venn diagram, n = 3,412), are enriched for gene signatures of stromal cells.

[Supplementary Fig1. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1e_corr.R): Correlation of gene expression between matching pairs of tumor tissue and derived PDO.

**Supplementary Figure 2**

[Supplementary Fig2. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2b.R): PCA of input normalized ChIP-seq signals for the five histone modifications used to build the ChromHMM 8- state model.

[Supplementary Fig2. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2c.py): Heatmaps showing the annotation of the ChromHMM 8-states with known genomic features (Overlap) and the probability that a state is found in the proximity of another state (Transition).

[Supplementary Fig2. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2c.py): Average proportion of each chromatin state over all PDOs. The chromatin segments for active/flanking TSS and active/flanking enhancer states are merged into the promoter and enhancer functional elements, respectively.

[Supplementary Fig2. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig2c.py):  Distribution of the eight ChromHMM 4 states for each PDO.

**Supplementary Figure 4**

[Supplementary Fig4. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/FigS4a.py): Pie chart showing the localization of the 2,419 ChromHMM-defined gained enhancers within functional features of the genome. UTR, untranslated region. see [this](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/FigS4a.pbs) for the gtf pre-processing.

[Supplementary Fig4. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): Percentage of concordant gained enhancers across PDOs.

[Supplementary Fig4. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): Number of conserved gained enhancers across PDOs. d, Boxplots of RNA-seq log2 normalized counts showing the expression distribution of genes that are annotated to gained active enhancers and upregulated in PDOs across normal colon tissues, primary tumors, and PDOs. P < 0.0001, Wilcoxon rank sum test.

[Supplementary Fig4. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/Fig1g_DE.R): Boxplots of RNA-seq log2 normalized counts showing the expression  distribution of genes that are annotated to gained active enhancers and upregulated in PDOs across normal colon tissues, primary tumors, and PDOs. P < 0.0001, Wilcoxon rank sum test.

**Supplementary Figure 5**

[Supplementary Fig5. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS5a_FigS6a.R): YAP and TAZ are transcriptionally upregulated in primary tumors compared to normal colon tissues in the TCGA colon adenocarcinoma dataset. Violin plots show the distribution of RNA-seq log2 normalized gene counts adjusted for epithelial cell frequency (see Methods for details). P < 0.0001, Wilcoxon rank sum test. 

[Supplementary Fig5. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_chipseq/analysis/scr/Fig3c_4b_4c_4f_4g_s5c.py): TAZ signal intensity is higher in conserved (shared by 8-10 PDOs) compared to non-conserved gained enhancers (shared by 1-5 PDOs). * P = 0.028, Mann Whitney U test exact P-value.  

**Supplementary Figure 6**

[Supplementary Fig6. a](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/bulk_rnaseq/analysis/scr/FigS5a_FigS6a.R): Heatmap of RNA-seq log2 normalized counts for target genes of YAP/TAZ-controlled enhancers in the TCGA COAD dataset. Tissue populations and clinical stages are represented by color-coded bars above the heatmap. Expression values of TCGA bulk RNA-seq data were adjusted for epithelial cell frequency (see Methods for details).

[Supplementary Fig6. b](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6.py): Detailed t-SNE visualization of all cell populations identified in the dataset. Cell type designations are based on the original clustering.

[Supplementary Fig6. c](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6-epith.py):  Sub-clustering analysis of 18,539 tumor and normal epithelial cells showing patient-specific clustering of most tumor cells. 

[Supplementary Fig6. d](https://github.com/paganilab/DellaChiara_et_al_2021/blob/main/sc_RNAseq/scr/Fig6.py): t-SNE plot depicting the distribution of the stemness score across all cell populations.

