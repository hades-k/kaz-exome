# Integrative Genomic Study of Kazakh Exomes in the Context of Global Populations

This repository contains scripts, notebooks, and documentation for the analysis of Kazakh exome sequencing data, including quality control, variant calling, population genetics, and comparative analysis with 1000 Genomes Project data.

## Project Overview

The project aims to:
- Process and quality-check raw exome sequencing data
- Call and annotate genetic variants
- Compare Kazakh exome data to global populations (1000 Genomes)
- Perform population structure analysis (PCA, clustering, admixture)
- Assess variant statistics and missingness

The workflow is implemented using a combination of shell commands, Python scripts, and Jupyter notebooks, with detailed step-by-step documentation in markdown files.

## Repository Structure

- `1000_gen_v2.ipynb` — Main Jupyter notebook with Python code for population genetics analysis and visualization
- `1000_gen_v2.md` — Step-by-step markdown protocol for processing, filtering, and merging genotype data (PLINK, ADMIXTURE, etc.)
- `variant_calling.md` — Detailed protocol for raw data QC, alignment, variant calling (GATK), and annotation
- `shared_vars_and_missingness.ipynb` — Notebook for analyzing shared variants and missingness statistics
- `vcf_stats.py` — Python script to compute summary statistics from a VCF file (variant types, Ti/Tv, quality, etc.)
- `count_rs_variants.py` — Python script to count and compare rsID variants between VCFs
- `fastqc/` — Directory containing FastQC HTML reports and quality plots for raw and trimmed reads

## Main Analysis Steps

### 1. Quality Control & Trimming
- FastQC is used to assess raw read quality (see `fastqc/` for reports)
- Trimmomatic is used for adapter and quality trimming
- Post-trimming quality is re-assessed with FastQC

### 2. Read Alignment & Preprocessing
- Reads are aligned to the hg38 reference genome using BWA
- SAM/BAM processing, sorting, and indexing with SAMtools
- Read groups are added and duplicates marked (GATK)
- Base quality recalibration (BQSR) is performed

### 3. Variant Calling & Filtering
- GATK HaplotypeCaller is used for variant discovery
- Joint genotyping merges samples into a cohort VCF
- Hard filtering is applied to SNPs and indels
- Variants are annotated with Annovar (dbSNP rsIDs)

### 4. Comparative and Population Analysis
- PLINK is used to filter, match, and merge Kazakh and 1000 Genomes data
- PCA and clustering (scikit-learn) reveal population structure
- ADMIXTURE is used for ancestry estimation
- Outlier and missingness analysis is performed

### 5. Variant Statistics & Overlap
- `vcf_stats.py` computes variant type, quality, and Ti/Tv statistics
- `count_rs_variants.py` compares rsID overlap between VCFs
- `shared_vars_and_missingness.ipynb` analyzes shared variants and missingness per variant

