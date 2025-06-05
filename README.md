#### Initial quality check

Steps were repeated for 2 samples WE001 and WE002 forward and reverse strands using command line FastQC

```bash
fastqc -o fastqc_results WE_run_first_fastq_gz/WE001.first.R1.fastq.gz
```
| -  | Per base quality WE001 | Per base quality WE002 |
| ----- | ------| ------|
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








