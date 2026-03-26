# RNA-seq-Transcriptomic-Analysis-of-APR246-Treatment-in-Breast-Cancer-Cells

This repository contains an RNA-seq transcriptomic analysis investigating gene expression changes in BT459 breast cancer cells treated with APR246 compared with DMSO controls.
The project was completed as part of the Computational Biology module in the MSc Technologies and Analytics in Precision Medicine.
The analysis integrates Linux-based preprocessing pipelines with statistical analysis in R to identify differentially expressed genes and biological pathways associated with APR246 treatment.

## Project Overview

RNA sequencing (RNA-seq) allows genome-wide measurement of transcript abundance and enables the identification of genes that are differentially expressed between experimental conditions.
In this project, RNA-seq data from BT459 breast cancer cells were analysed to evaluate the transcriptional effects of APR246 treatment, a compound known to restore mutant p53 function and induce tumour cell apoptosis.

The analysis involved:

- preprocessing and quality control of raw sequencing reads
- read alignment and gene quantification
- statistical analysis of differential gene expression
- pathway enrichment analysis to interpret biological effects

The goal of the study was to identify genes and pathways affected by APR246 treatment.

### Research Question

Which genes and biological pathways are differentially expressed in BT459 breast cancer cells following APR246 treatment compared with DMSO controls?

## Dataset

The dataset consisted of RNA-seq reads generated from:

- BT459 breast cancer cells
- APR246-treated samples
- DMSO control samples

Each experimental condition contained three biological replicates, producing a total of six RNA-seq samples.

## RNA-seq Preprocessing Workflow

Raw sequencing data were processed on a high-performance computing (HPC) cluster using a standard RNA-seq preprocessing workflow provided during the course.

The preprocessing workflow included the following steps:

### 1. Read Quality Control

Initial sequencing quality was assessed using FastQC, which evaluates metrics such as:

- per-base sequence quality
- GC content distribution
- sequence duplication levels
- adapter contamination

These metrics help identify potential issues with sequencing data before further processing.

### 2. Adapter Trimming and Filtering

Low-quality bases and adapter sequences were removed using BBMap trimming tools.
Trimming improves alignment performance by removing sequencing artifacts and low-confidence bases.
Quality control was repeated after trimming to confirm improvements in read quality.

### 3. Genome Alignment

Trimmed reads were aligned to the reference genome using the Rsubread alignment algorithm.
This step generates BAM alignment files containing the genomic locations of sequencing reads.

### 4. Gene Quantification

Gene-level read counts were generated using featureCounts, which assigns aligned reads to annotated genomic features.
The resulting count matrix contains the number of reads mapping to each gene for each sample and forms the basis of downstream statistical analysis.

# Differential Expression Analysis

Downstream transcriptomic analysis was performed in R using the limma-voom framework.

Key steps included:

- Creating a DGEList object from the gene count matrix
- Filtering low-expression genes
- Normalising library sizes
- Transforming count data using voom
- Fitting a linear model to estimate gene expression differences between conditions
- Applying treat() to test for biologically meaningful fold-changes

Genes were considered significantly differentially expressed based on adjusted p-values and fold-change thresholds.

## Exploratory Data Analysis

Several plots were used to assess sample quality and relationships between samples.

- Library Size Distribution: Library size plots were used to examine differences in sequencing depth across samples.

- Multidimensional Scaling (MDS): MDS plots visualised similarity between samples based on gene expression profiles.

Samples clustered according to treatment condition, indicating a measurable transcriptional response to APR246 treatment.

## Differential Expression Results

Differential expression analysis identified genes that were significantly up-regulated or down-regulated following APR246 treatment.

Visualisations generated during the analysis include:

- Volcano plots highlighting significantly regulated genes
- Heatmaps of differentially expressed genes
- Sample expression distribution plots

These visualisations help identify transcriptional signatures associated with treatment.

## Pathway Enrichment Analysis

To interpret the biological significance of the differentially expressed genes, KEGG pathway enrichment analysis was performed.

Pathway analysis identified enrichment of biological pathways associated with:

- apoptosis and programmed cell death
- cellular stress response pathways
- cancer signalling pathways

These findings are consistent with the known mechanism of action of APR246, which restores p53 function and promotes tumour cell apoptosis.

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
