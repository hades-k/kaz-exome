### Quality check and trimming

#### Initial quality check

Steps were repeated for 2 samples WE001 and WE002 forward and reverse strands using command line FastQC

```bash
fastqc -o fastqc_results WE_run_first_fastq_gz/WE001.first.R1.fastq.gz
```
| -     | Per base quality WE001 | Per base quality WE002 |
| ----- | ------                 | ------                 |
| Forward | <img src="https://github.com/hades-k/thesis/blob/main/fastqc/WE001.R1.png">| <img src = "https://github.com/hades-k/thesis/blob/main/fastqc/WE002.R1.png"> |
| Reverse | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE001.R2.png'> | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE002.R2.png'> |

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

| WE001                                     | WE002                                     |
| ----                                      | ----                                      |
| Input Read Pairs: 63429143                | Input Read Pairs: 48713866                |
| Both Surviving Reads: 61400551            | Both Surviving Reads: 44617486            |
| Surviving Read Percent: 96.80             | Both Surviving Read Percent: 91.59        |
| Forward Only Surviving Reads: 1328997     | Forward Only Surviving Reads: 2639723     |
| Forward Only Surviving Read Percent: 2.10 | Forward Only Surviving Read Percent: 5.42 |
| Reverse Only Surviving Reads: 476903      | Reverse Only Surviving Reads: 878593      |
| Reverse Only Surviving Read Percent: 0.75 | Reverse Only Surviving Read Percent: 1.80 |
| Dropped Reads: 222692                     | Dropped Reads: 578064                     |
| Dropped Read Percent: 0.35                | Dropped Read Percent: 1.19                |

After trimming, a second quality check was performed to assess the results. 

| -  | Per base quality WE001 | Per base quality WE002 |
| ----- | ------| ------|
| Forward | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE001.R1.paired.png'> | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE002.R1.paired.png'> |
| Reverse | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE001.R2.paired.png'> | <img src = 'https://github.com/hades-k/thesis/blob/main/fastqc/WE002.R2.paired.png'> | 

______

### Read alignment and pre-processing for variant discovery

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

#### Mark Duplicates 

After alignment, I used GATK’s MarkDuplicates to identify and flag PCR or optical duplicate reads—i.e., reads that originate from the same DNA fragment but appear multiple times due to library amplification artifacts. These duplicates do not represent independent observations and can bias downstream analyses, particularly variant allele frequency estimation and variant calling confidence.

```bash
gatk MarkDuplicates \
  -I samtools_res/WE001_rg.bam \
  -O gatk_res/WE001_marked.bam \
  -M gatk_res/WE001_marked.metrics.txt

samtools index gatk_res/WE001_marked.bam
```
Some of the key stats form marked metrics files are: 

| Stat                           | WE001      | WE002      |
| ---                            | ---        | ---        |
| UNPAIRED_READS_EXAMINED        | 3350       | 5748       |
| READ_PAIRS_EXAMINED            | 61395821   | 44609640   |
| SECONDARY_OR_SUPPLEMENTARY_RDS | 26704      | 52170      |
| UNMAPPED_READS                 | 6110       | 9944       |
| UNPAIRED_READ_DUPLICATES       | 982        | 1494       |
| READ_PAIR_DUPLICATES           | 1865459    | 647259     |
| READ_PAIR_OPTICAL_DUPLICATES   | 22709      | 4005       |
| PERCENT_DUPLICATION            | 0.030391   | 0.014525   |
| ESTIMATED_LIBRARY_SIZE         | 1001459378 | 1531656110 |


After this, the marked bam file was indexed again. 

```bash
samtools index gatk_res/WE001_marked.bam
```
#### Base quality recalibration

The next step is Base Quality Recalibration (BQSR). This applies machine learning to correct patterns of systematic errors in the base quality scores. Base quality scores play an important role in weighing the evidence for or against possible variant alleles during the variant discovery process, so it's important to correct any systematic bias observed in the data. To distinguish between true sequencing errors and real variants, it needs a list of "known" trusted variants.

```bash
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz.tbi
```

The chromosome naming is inconsistent beweem NCBI and UCSC, namely '1, 2, 3, ...' and 'chr1, chr2, chr3, ...' To fix this bcftools annontate was used

```bash
bcftools query -f '%CHROM\n' 00-common_all.vcf.gz | uniq > chrom.txt
# Edit chrom.txt manually or script the rename map
bcftools annotate --rename-chrs chrom.txt 00-common_all.vcf.gz -Oz -o renamed_common.vcf.gz
```

After all this GATK BaseRecalibrator can be used to create a table, that is later used by ApplyBQSR

```bash
gatk BaseRecalibrator \
  -I gatk_res/WE001_marked.bam \
  -R hg38.fa \
  --known-sites renamed_common.vcf.gz \
  -O gatk_res/WE001_recal_data.table

gatk ApplyBQSR \
  -R hg38.fa \
  -I gatk_res/WE001_marked.bam \
  --bqsr-recal-file gatk_res/WE001_recal_data.table \
  -O gatk_res/WE001_recalibrated.bam
```

_______

### Variant discovery 

Following GATK best prractices:

HaplotypeCaller wau used for this on the recallibrated set. This tool supports multi-thread on the native pair hidden markov models stage only. 

```bash
gatk HaplotypeCaller \
  -R hg38.fa \
  -I gatk_res/WE001_recalibrated.bam \
  -O gatk_res/WE001.g.vcf.gz \
  --native-pair-hmm-threads 25 \
  -ERC GVCF
```
#### Joint genotyping 

After both samples have been processed, we can join them for further work 

```bash 
echo -e "WE001\t/media/dell_816/aidana/KAZ_WE/gatk_res/WE001.g.vcf.gz\nWE002\t/media/dell_816/aidana/KAZ_WE/gatk_res/WE002.g.vcf.gz" > sample_map.txt
```
``` bash 
nohup gatk GenomicsDBImport \
  --genomicsdb-workspace-path gatk_res/cohort \
  --sample-name-map sample_map.txt \
  -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
  -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16 -L chr17 -L chr18 -L chr19 \
  -L chr20 -L chr21 -L chr22 -L chrX -L chrY \
  --tmp-dir /media/dell_815/aidana/KAZ_WE/gatk_res/tmp \
  > GenomicsDBImport.log 2>&1 &

```

```bash
gatk GenotypeGVCFs \
  -R hg38.fa \
  -V gendb:///media/dell_815/aidana/KAZ_WE/gatk_res/cohort \
  -O /media/dell_815/aidana/KAZ_WE/gatk_res/cohort/joint_genotyped.vcf.gz \
  --tmp-dir /media/dell_815/aidana/KAZ_WE/gatk_res/tmp
```

This way we have joined the two vcf files into one for a population analysis

#### Variant filtration 

We do not have enough samples to perform VQSR, so just hard-filtering for now. We need to split the vcf into SNPs and Indels. 

```bash
gatk SelectVariants \
  -V gatk_res/cohort/joint_genotyped.vcf.gz \
  --select-type-to-include INDEL \
  -O gatk_res/cohort/joint_indels.vcf.gz

gatk SelectVariants \
  -V gatk_res/cohort/joint_genotyped.vcf.gzz \
  --select-type-to-include SNP \
  -O gatk_res/cohort/joint_snps.vcf.gz
```

To visualize the variant quality distribution: 

```bash
gatk VariantsToTable \
  -V joint_snps.vcf.gz \
  -F QD -F FS -F MQ -F MQRankSum -F ReadPosRankSum \
  -O joint_snps_annotations.table
```

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("WE001_snps_annotations.table", sep="\t")
annotations = ['QD', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum']

plt.figure(figsize=(15,10))
for i, col in enumerate(annotations, 1):
    plt.subplot(3, 2, i)
    df[col].hist(bins=50)
    plt.title(col)
    plt.xlabel(col)
    plt.ylabel("Count")

plt.tight_layout()
plt.show()

```
Using the recommended tresholds for hard fitering: 

```bash
gatk VariantFiltration   -V joint_snps.vcf.gz   -O joint_filtered_snps.vcf.gz \
  --filter-name "QD_lt_2" --filter-expression "QD < 2.0" \
  --filter-name "FS_gt_60" --filter-expression "FS > 60.0"  \
  --filter-name "MQ_lt_40" --filter-expression "MQ < 40.0"   \
  --filter-name "MQRankSum_lt_-12.5" --filter-expression "MQRankSum < -12.5"   \
  --filter-name "ReadPosRankSum_lt_-8" --filter-expression "ReadPosRankSum < -8.0"

gatk VariantFiltration   -V joint_indels.vcf.gz  -O joint_filtered_indels.vcf.gz   \
  --filter-name "QD_lt_2"   --filter-expression "QD < 2.0"   \
  --filter-name "FS_gt_200" --filter-expression "FS > 200.0"   \
  --filter-name "ReadPosRankSum_lt_-20" --filter-expression "ReadPosRankSum < -20.0"
```

And Merging the snp vcf and indel vcf back together:

```bash
gatk MergeVcfs   -I joint_filtered_snps.vcf.gz   -I joint_filtered_indels.vcf.gz   -O joint_filtered.vcf.gz
```

### Annotation and comparison to the pre-processed files

#### Annotation 

Annovar was used for filter-based annotation. Using avsnp150 (an abbreviated version of dbSNP 150 with left-normalization by ANNOVAR developers), RS IDs were added to the vcf file:

```bash 
table_annovar.pl joint_filtered.vcf /Applications/annovar/humandb   \
  -buildver hg38   -out annovar_res/annotated   -remove   \
  -protocol avsnp150   -operation f   -vcfinput -otherinfo
```

 This file was compared to the pre-processed data. Additionally, sample 1 (WE001) was compared to the pre-processed vcf containing variants for sample one using python:
 
 [stats](https://github.com/hades-k/thesis/blob/main/vcf_stats.py)
 
 [rs id comparison](https://github.com/hades-k/thesis/blob/main/count_rs_variants.py)

| Statistic        | My VCF full | Gold VCF full | My VCF WE001 | Gold VCF WE001 |
| ---              | ---         | ---           | ---          | ---            | 
| Total variants   | 900195      | 500376        | 520962       | 496454         |
| SNPs             | 793691      | 442486        | 451009       | 439326         |
| INDELs           | 106504      | 57890         | 69953        | 57128          |
| Quality min      | 30.0        | 30.0          | 30.01        | 30.0           |
| Quality max      | 66550.13    | 5149018.88    | 41806.06     | 5149018.88     |
| Quality mean     | 466.43      | 28494.77      | 489.81       | 28698.81       |
| Pass filters     | 885256      | 474413        | 510564       | 470905         |
| Transitions      | 531449      | 310132        | 305999       | 308057         |
| Transversions    | 262242      | 132354        | 145010       | 131269         |
| Ti/Tv ratio      | 2.03        | 2.34          | 2.11         | 2.35           |
| RS variant count | 834048      | 369515        | 485420       | 366361         |
| % dbsnp vars in total |  92.65 | 73.85         | 93.18        | 73.79          |
| Overlapping rs   | 129337      |      <-       | 97240        | <-             |
| % overlap        | 15.5        | 35            | 20           | 26.5           |

Possible explanations of differences: 

- Different reference genome builds were used hg38 v hg19 - would explain a small difference, not such a drastic one
- My VCF uses hard filtering, the gold VCF uses VQSR - does not impact total number
- Different pre processing - need to figure this out
- Different reference dataset of known variants for BQSR - mine was most likely smaller, used NCBI common vars
- Gold VCF has a higher Ti/Tv - likely higher specificity and stricter filtering
- My VCF has more variants in dbsnp - different version of dbsnp or inclusion of lower quality variants
- I used joint calling, so variants were added together even if they are present in only one sample, full gold VCF has roughly the same number of variants as one-sample gold. Some kind of filter was applied, keeping high quality variants present in most samples.
- Most likely gold was manually filtered to exclude variants with a certain percentage of missingness - according to bioinfromatician who performed it
- Gold was processed 10 years ago, tools, algorithms, etc change and evolve

[Calculations here](https://github.com/hades-k/thesis/blob/main/shared_vars_and_missingness.ipynb)

Calculating the number of variant present in both of my samples in my full VCF, I got **392741 (~44%)**. Less than half of the variants were present in both samples. In addition, I looked at the gold vcf missingness per variant. The mean value is **5.58%**, and the median is **0.8%**, meaning that afterwards variants with missingness above a certain treshold were removed. Interestingly, this filter was imperfect, as there are very few variants with very high missingness present, having escaped filters, max value is **99.2%**. 








