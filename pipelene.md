#### Initial quality check

Steps were repeated for 2 samples WE001 and WE002 forward and reverse strands using command line FastQC

```bash
fastqc -o fastqc_results WE_run_first_fastq_gz/WE001.first.R1.fastq.gz
```
| -     | Per base quality WE001 | Per base quality WE002 |
| ----- | ------                 | ------                 |
| Forward | ![per_base_quality](https://github.com/user-attachments/assets/43c9b67d-1e22-4a80-9d88-2850e1b14134) | ![per_base_quality](https://github.com/user-attachments/assets/cdd54e31-3bfa-418c-8984-3b262e7c576c) |
| Reverse | ![per_base_quality](https://github.com/user-attachments/assets/e5cc4ae6-554d-4083-bec8-73fa64322d3b) | ![per_base_quality](https://github.com/user-attachments/assets/8ae41847-60d6-4b3b-8e5e-61aa7561f79b) |

#### Adapter and Quality Trimming

After the quality check, the reads were trimmed with Trimmomatic, the tool available to me on the server. 

```bash
trimmomatic PE -threads 16 \
  -trimlog trimming_res/001_log.txt \
  WE_run_first_fastq_gz/WE001.first.R1.fastq.gz WE_run_first_fastq_gz/WE001.first.R2.fastq.gz \
  trimming_res/WE001_R1_paired.fastq.gz trimming_res/WE001_R1_unpaired.fastq.gz \
  trimming_res/WE001_R2_paired.fastq.gz trimming_res/WE001_R2_unpaired.fastq.gz \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30 \
  -summary trimming_res/001_summary.txt
```

Using this command, bases with quality lower than 20 were trimmed from both ends of a read. Additionally, trimmomatic would scan the read using a 4bp sliding window, and would trim bases if the average quality would be < 20. After completion, if the read length would be under 30bp, the read would be discarded. The summary of the process is available below: 

| WE001 | WE002 |
| ---- | ----|
| Input Read Pairs: 63429143 | Input Read Pairs: 48713866 |
| Both Surviving Reads: 61400551 | Both Surviving Reads: 44617486 |
| Surviving Read Percent: 96.80 | Both Surviving Read Percent: 91.59 |
| Forward Only Surviving Reads: 1328997 | Forward Only Surviving Reads: 2639723 |
| Forward Only Surviving Read Percent: 2.10 | Forward Only Surviving Read Percent: 5.42 |
| Reverse Only Surviving Reads: 476903 | Reverse Only Surviving Reads: 878593 |
| Reverse Only Surviving Read Percent: 0.75 | Reverse Only Surviving Read Percent: 1.80 |
| Dropped Reads: 222692 | Dropped Reads: 578064 |
| Dropped Read Percent: 0.35 | Dropped Read Percent: 1.19 |

After trimming, a second quality check was performed to assess the results. 

| -  | Per base quality WE001 | Per base quality WE002 |
| ----- | ------| ------|
| Forward | ![per_base_quality](https://github.com/user-attachments/assets/4893d329-2e2d-4b60-b8c7-869631520ba9) | ![per_base_quality](https://github.com/user-attachments/assets/fea045cb-8f75-45a9-9a8d-5a65b709df09) |
| Reverse | ![per_base_quality](https://github.com/user-attachments/assets/df41bacd-00b4-4e18-8f09-8d3944feb1b5) | ![per_base_quality](https://github.com/user-attachments/assets/a90785b0-76dd-4ea2-8ac1-e8695574a2b2) | 

#### Reference preparation 

Hg38 assembly was retrieved from the UCSC database. In reality, I have made a mistake in choosing the reference fasta file, as I dowloaded the primary assembly, that contains alternative chromosomes and random sequences, and not the analysis set, containing only the canonical chromosomes. More on this later. 

Prior to alignment, the reference genome file was indexed using BWA. This will enable efficient and accurate alignment of reads by reducing the time and memor required to search for matching sequences across the genome. 

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bwa index hg38.fa
```

#### Read alignment

BWA was chosen for read alignment due to its high accuracy, efficiency with short Illumina reads, and support for gapped alignment, and compatibility with downstream variant calling tools. (Another hiccup: did not add read groups in this step. Fixed later)

```bash
bwa mem -t 22 hg38.fa \
  trimming_res/WE001_R1_paired.fastq.gz trimming_res/WE001_R2_paired.fastq.gz \
  > WE001_aligned.sam
```

#### SAM → BAM, Sorting, and Stats

Sam file was converted into bam:

```bash 
samtools view -@ 22 -bo bwa_res/WE001_aligned.bam WE001_aligned.sam
```

To get alignment statistics the following command was run: 

```bash
samtools stats bwa_res/WE001_aligned.samtools statsbam > bwa_res/stats.txt
```
Some of the key stats are: 

| Metric                      | WE001                | WE002               | 
| -------                     | -----                | -----               |
| Total reads (raw sequences) | 122,801,102          | 89,234,972          |
| Mapped reads                | 122,794,992          | 89,225,028          |         
| Properly paired reads       | 122,158,754 (≈99.5%) | 88,513,144 (≈99.2%) | 
| Unmapped reads              | 6,110                | 9,944               |
| Reads MQ=0                  | 4,225,839            | 3,414,963           |
| Error rate                  | 0.158%               | 0.193%              | 
| Average read length         | 98 bp                | 95 bp               |
| Average quality             | 37.4                 | 36.5                |
| Insert size (mean ± SD)     | 212.0 ± 68.1 bp      | 190.5 ± 46.4 bp     |

The sam files were sorted and indexed
```bash
samtools sort -@ 22 -o samtools_res/WE001_sorted.bam bwa_res/WE001_aligned.bam
samtools index samtools_res/WE001_sorted.bam
```

#### Add Read Groups (fixing previous omission)

As said prior, the read groups were not added during alignment step. Therefore, gatk AddOrReplaceReadGroups function was used to add a 'generic' read group.

```bash
gatk AddOrReplaceReadGroups \
  -I samtools_res/WE001_sorted.bam \
  -O samtools_res/WE001_rg.bam \
  -RGID WE001 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM WE001
```


