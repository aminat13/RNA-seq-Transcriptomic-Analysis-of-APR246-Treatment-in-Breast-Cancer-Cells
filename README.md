# RNA-seq Transcriptomic Analysis of APR246 Treatment in Breast Cancer Cells

This repository contains an RNA-seq transcriptomic analysis investigating gene expression changes in BT459 breast cancer cells treated with APR246 compared with DMSO controls.

The project was completed as part of the Computational Biology module in the MSc Technologies and Analytics in Precision Medicine.

The analysis combines RNA-seq preprocessing on an HPC cluster with downstream differential expression and pathway enrichment analysis in R.

---

## Project Overview

RNA sequencing (RNA-seq) enables genome-wide quantification of transcript abundance and allows the identification of genes that are differentially expressed between biological conditions.

In this project, RNA-seq data from BT459 breast cancer cells were analysed to evaluate the transcriptional effects of APR246 treatment relative to DMSO controls.

The report describes the following workflow:

- FastQC assessment of raw reads
- adapter trimming with BBMap
- repeat FastQC on trimmed reads
- index generation
- read alignment
- gene counting and annotation
- differential expression analysis
- pathway enrichment analysis

The aim of the study was to identify genes and biological pathways affected by APR246 treatment.

---

### Research Question

Which genes and biological pathways are differentially expressed in BT459 breast cancer cells following APR246 treatment compared with DMSO controls?

---

## Dataset

The dataset consisted of paired-end RNA-seq data generated from:

- BT459 breast cancer cells
- APR246-treated samples
- DMSO control samples

Each condition contained three biological replicates, giving a total of six samples.

---

## RNA-seq Preprocessing

According to the report, the raw sequencing data were processed using several Linux-based workflows on an HPC cluster.

The preprocessing steps included:

1. FastQC quality assessment of raw reads  
2. Adapter trimming using BBMap  
3. FastQC on trimmed reads  
4. Index generation
5. Read alignment 
6. Gene count generation
7. Gene annotation

The report also notes that commands such as `sacct -j <jobID>` were used to monitor job status and `more slurm-<jobID>.out` was used to inspect batch job output.

---

## Exploratory Quality Assessment

The report includes interpretation of the following RNA-seq QC outputs:

- Per Base Sequence Quality: All sequences were described as showing good per-base sequence quality, suggesting low sequencing error rates.

- Per Sequence GC Content: Most sequences showed moderate GC-content deviation, interpreted as mild but not problematic technical deviation.

- Sequence Duplication Levels: All sequences showed poor duplication reports, suggesting low library complexity likely due to PCR amplification or highly abundant transcripts.

A summary table of FastQC results for raw and trimmed sequences is included in the report.

---

## Summary Plots

The report presents three main summary plots.

- Mapped vs Unmapped Reads: All six samples showed approximately **23–30 million reads** with a very small unmapped fraction, indicating high mapping efficiency.
-  Library Size per Sample: Library sizes ranged from approximately **17–23 million counts**, with no clear outliers across APR246 and DMSO replicates.

- Multidimensional Scaling (MDS) Plot: APR246 and DMSO samples formed two distinct clusters along principal component 1, indicating that treatment was a major source of variation.

---

## Differential Expression Analysis

Differential expression analysis was performed in R using the limma-voom workflow.

The report shows code for:

- constructing a design matrix
- defining the contrast APR246 - DMSO
- normalising expression data using voom
- fitting a linear model
- extracting differentially expressed genes
- subsetting the top 5 upregulated and top 5 downregulated genes
- generating a volcano plot

Positive log fold-change values represent genes upregulated in APR246, while negative values indicate relative downregulation in APR246 or higher expression in DMSO.

---

## Differential Expression Results

The report highlights the following top upregulated genes in APR246-treated cells:

- WISP3
- TAS2R5
- LOC44887
- FBXO39
- TCEB3CL2

The top downregulated genes reported were:

- MYCL
- PCDH18
- HIST1H1A
- TET1
- SNORD73A

A volcano plot was used to visualise the differentially expressed genes between APR246 and DMSO.

---

## Pathway Enrichment Analysis

The report includes KEGG pathway enrichment analysis of the differentially expressed genes.

The most enriched pathways reported were:

- Cell cycle
- MicroRNAs in cancer
- FoxO signalling pathway

These pathways were discussed in the context of breast cancer biology and the proposed anti-tumour activity of APR246.

## Running the Analysis

The statistical analysis can be reproduced using the R script located in:

scripts/transcriptomics_analysis.R

The script performs:

- exploratory transcriptomic analysis
- differential expression testing
- visualisation of gene expression patterns
- pathway enrichment analysis

Note that the original raw sequencing files are not included in this repository.

## Skills Demonstrated

This project demonstrates practical skills in:

- RNA-seq data analysis
- high-throughput sequencing quality control
- gene expression normalisation
- differential expression analysis using limma-voom
- pathway enrichment analysis
- biological interpretation of transcriptomic data
- data visualisation in R

## Reproducibility Note

This repository represents an archival portfolio version of the project.
The original raw sequencing files and intermediate alignment files are not included because access to the source data is no longer available.

As a result:

- figures included in this repository were extracted from the final report
- the R script documents the downstream statistical analysis performed during the project
