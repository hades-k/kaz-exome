# Integrative Genomic Study of Kazakh Exomes in the Context of Global Populations

This repository contains the complete workflow, scripts, and results for a comprehensive population genetics study of Kazakh exomes. The analysis integrates quality control, variant calling, and comparative analysis with data from the 1000 Genomes Project to investigate the genetic structure and ancestry of the Kazakh population.

## Project Overview

This research aims to characterize the genetic landscape of the Kazakh population by:

- Processing and performing rigorous quality control on raw exome sequencing data.
- Calling and annotating genetic variants to identify high-quality SNPs and indels.
- Merging the Kazakh exome data with a global reference panel from the 1000 Genomes Project.
- Conducting detailed population structure analyses, including Principal Component Analysis (PCA), clustering, and admixture analysis, to understand the genetic relationships between the Kazakh population and other global populations.
- Performing in-depth variant-level statistics and assessing data quality through missingness and heterozygosity analysis.

The project workflow is implemented through a series of shell commands, Python scripts, and Jupyter notebooks, with detailed protocols available in markdown files.

## Repository Structure

- **`1000_gen_v2.ipynb`**: The main Jupyter notebook for population genetics analysis, including data merging, PCA, clustering, and admixture analysis.
- **`100_samp_pca.ipynb`**: A supplementary notebook focused on PCA for a subset of 100 samples, providing a more detailed look at population substructure.
- **`variant_calling.md`**: A step-by-step protocol for raw data quality control, read alignment (BWA), variant calling (GATK), and annotation.
- **`1000_gen_v2.md`**: A detailed guide for processing, filtering, and merging genotype data using PLINK, including steps for preparing data for admixture and other population-level analyses.
- **`plink_start.md`**: A foundational document detailing the initial steps of creating and quality-controlling PLINK datasets from VCF files, including missingness, sex checks, and MAF filtering.
- **`shared_vars_and_missingness.ipynb`**: A notebook for analyzing shared variants and missingness statistics across different populations.
- **`vcf_stats.py`**: A Python script for computing summary statistics from VCF files, such as variant types, Ti/Tv ratio, and quality metrics.
- **`count_rs_variants.py`**: A Python script to count and compare rsID variants between different VCF files.
- **`fastqc/`**: A directory containing FastQC reports and quality plots for both raw and trimmed sequencing reads.
- **`pca_plots/`**, **`clustering_plots/`**, **`admixture_plots/`**, **`fst_plots/`**: Directories containing visualization of the analysis results.

## Analysis Workflow

The analysis is divided into several key stages:

### 1. Raw Data Quality Control and Preprocessing
- **Quality Assessment**: Raw sequencing reads are assessed for quality using FastQC.
- **Trimming**: Trimmomatic is used to remove adapters and low-quality reads.
- **Alignment**: Reads are aligned to the hg38 reference genome using BWA.
- **Preprocessing**: SAM/BAM files are processed, sorted, and indexed with SAMtools. GATK is used to add read groups, mark duplicates, and perform Base Quality Score Recalibration (BQSR).

### 2. Variant Calling and Annotation
- **Variant Discovery**: GATK HaplotypeCaller is used for initial variant discovery.
- **Joint Genotyping**: Samples are merged into a cohort VCF through joint genotyping.
- **Filtering and Annotation**: Hard filters are applied to SNPs and indels, and variants are annotated with dbSNP rsIDs using Annovar.

### 3. PLINK Data Preparation and Quality Control
- **VCF to PLINK**: The VCF file is converted to PLINK format.
- **Missingness Filtering**: SNPs and individuals with high rates of missing genotypes are removed.
- **Sex Discrepancy Check**: Individuals with mismatched sex information are identified and removed.
- **MAF and HWE Filtering**: SNPs with low Minor Allele Frequency (MAF) and those deviating from Hardy-Weinberg Equilibrium (HWE) are filtered out.
- **Heterozygosity and Relatedness Filtering**: Individuals with outlier heterozygosity rates and related individuals are removed to ensure a dataset of independent samples.

### 4. Population Genetics and Comparative Analysis
- **Data Merging**: The filtered Kazakh dataset is merged with data from the 1000 Genomes Project using PLINK.
- **Principal Component Analysis (PCA)**: PCA is performed to investigate population structure and identify genetic outliers.
- **Clustering and Admixture**: Clustering algorithms and ADMIXTURE are used to infer ancestry components and genetic affinities between populations.
- **Variant Statistics**: Custom Python scripts are used to compute and compare variant statistics between the Kazakh and 1000 Genomes datasets.

This structured workflow ensures high-quality data and robust downstream analyses, providing a solid foundation for investigating the genetic history of the Kazakh population.