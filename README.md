# Integrative Genomic Study of Kazakh Exomes in the Context of Global Populations

This repository contains the complete workflow, scripts, and results for a population genetics study of Kazakh exomes. The analysis integrates quality control, variant calling, and comparative analysis with data from the 1000 Genomes Project to investigate the genetic structure and ancestry of the Kazakh population.

## Project Overview

This study aims to characterize the genetic landscape of the Kazakh population by:

- Processing and performing quality control on raw exome sequencing data.
- Calling and annotating genetic variants to identify high-quality SNPs and indels.
- Merging the Kazakh exome data with a global reference panel from the 1000 Genomes Project.
- Conducting population structure analyses, including Principal Component Analysis (PCA), ADMIXTURE, Fixation Index (Fst), and Locus-Specific Branch Length (LSBL) to understand the genetic relationships and selective pressures.
- Performing in-depth variant-level statistics and assessing data quality through missingness and heterozygosity analysis.

The project workflow is implemented through a series of shell commands, Python scripts, and Jupyter notebooks, with detailed protocols available in markdown files.

## Repository Structure

- **`1000_gen_v2.ipynb`**: The main Jupyter notebook for population genetics analysis, including data merging, PCA, clustering, and admixture analysis.
- **`100_samp_pca.ipynb`**: A supplementary notebook focused on PCA for a larger subset of ~100 samples per population.
- **`lsbl_calc.ipynb`**: A notebook for calculating and plotting Locus-Specific Branch Length (LSBL) statistics.
- **`variant_calling.md`**: A step-by-step protocol for raw data quality control, read alignment (BWA), variant calling (GATK), and annotation.
- **`1000_gen_v2.md`**: A detailed guide for processing, filtering, and merging genotype data using PLINK, including steps for preparing data for admixture and other population-level analyses.
- **`plink_start.md`**: A foundational document detailing the initial steps of creating and quality-controlling PLINK datasets from VCF files.
- **`shared_vars_and_missingness.ipynb`**: A notebook for analyzing shared variants and missingness statistics.
- **`vcf_stats.py`**: A Python script for computing summary statistics from VCF files.
- **`fastqc/`**: Directory containing FastQC reports for raw and trimmed sequencing reads.
- **`pca_plots/`**, **`clustering_plots/`**, **`admixture_plots/`**, **`fst_plots/`**, **`lsbl_plots/`**: Directories containing visualizations of the analysis results.

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
- **Filtering and Annotation**: Hard filters are applied to SNPs and indels, and variants are annotated with dbSNP rsIDs.

### 3. Data Harmonization with 1000 Genomes
- **QC and Filtering**: Both Kazakh and 1000 Genomes datasets are filtered for missingness (`--mind`, `--geno`) and minor allele frequency (`--maf`).
- **SNP Matching**: Variants are matched between the two datasets to create a common set of SNPs.
- **Harmonization**: The datasets are harmonized by updating chromosome coordinates, setting a uniform reference allele, flipping strands to resolve mismatches, and removing ambiguous or problematic SNPs.
- **Outlier Removal**: PCA is used to identify and remove genetic outliers from the Kazakh cohort.
- **LD Pruning**: Linkage disequilibrium (LD) pruning is performed to ensure SNPs are independent for population structure analyses.

### 4. Population Genetics and Comparative Analysis
- **Dataset Subsetting**: A balanced dataset of ~100 individuals per population is created for robust analysis.
- **Principal Component Analysis (PCA)**: PCA is performed to investigate population structure.
- **ADMIXTURE Analysis**: ADMIXTURE is used to infer ancestry components, with cross-validation to determine the optimal number of ancestral populations (K).
- **Fixation Index (Fst) Calculation**: Pairwise Fst values are calculated using Arlequin and VCFtools to quantify genetic differentiation between populations.
- **Locus-Specific Branch Length (LSBL)**: LSBL is calculated for population triads (e.g., Kazakh, British, Han Chinese) to identify genomic regions under positive selection in the Kazakh lineage.

## Key Findings

- **Genetic Structure**: PCA places the Kazakh population on a clear East-West Eurasian axis, intermediate between European and East Asian clusters. PC2 further separates South Asian populations. ![pca](/pca_plots/pca_100_samp/pca_supergroup_100.png)
- **Admixture Analysis**: The most informative ADMIXTURE model (K=5) reveals that the Kazakh genome is a composite of multiple ancestral components, including East Asian, Southeast Asian, European, and South Asian lineages. Additionally, share a North Eurpean component with Finnish. ![admixture](/admixture_plots/5_ancestry_bs_kaz.png)
- **Population Differentiation (Fst)**: Fst values indicate that the Kazakh population is genetically closest to South Asian groups (e.g., Bengali) and most distant from European and East Asian populations, confirming its admixed nature. All pairwise Fst values were statistically significant. ![fst](/fst_plots/fst_plot_bar.png)
- **Positive Selection (LSBL)**: (ongoing)

