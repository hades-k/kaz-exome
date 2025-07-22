## PLINK data creation from a VCF file

MS_KAZ_WE_125_kaz_samples.MSHC.gold.vcf is used for the following analyses.

```bash
plink --vcf ../MS_KAZ_WE_125_kaz_samples.MSHC.gold.vcf --make-bed \
  --out PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
```
Output:
500376 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.

##  Identify missingness per individual and per SNP

Excludes SNPs that are missing in a large proportion of the subjects. In this step, SNPs with low genotype calls are removed.

```bash
plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold --missing
```
Output: Total genotyping rate is 0.939852.

Get the histogram in R
```R
indmiss<-read.table(file="plink.imiss", header=TRUE)
snpmiss<-read.table(file="plink.lmiss", header=TRUE)

pdf("histimiss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness")

pdf("histlmiss.pdf") 
hist(snpmiss[,5],main="Histogram SNP missingness")  
dev.off()
```

[histimiss.pdf](https://github.com/user-attachments/files/21369672/histimiss.pdf)

[histlmiss.pdf](https://github.com/user-attachments/files/21369673/histlmiss.pdf)

First, we filter SNPs and individuals based on a relaxed threshold (0.2; >20%), as this will filter out SNPs and individuals with very high levels of missingness. Then a filter with a more stringent threshold can be applied (0.02).

```bash
plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold --geno 0.2 \
  --make-bed --out PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter
```

Output: Total genotyping rate is 0.939852.
50618 variants removed due to missing genotype data (--geno).
449758 variants and 125 people pass filters and QC.

```bash
plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold \
  --geno 0.02 --make-bed --out PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter
```

Output: Total genotyping rate is 0.939852.
167744 variants removed due to missing genotype data (--geno).
332632 variants and 125 people pass filters and QC.

##  Sex discrepancy

```bash
plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter --check-sex
```
Output: --check-sex: 5666 Xchr and 0 Ychr variant(s) scanned, 125 problems detected.

```R
gender <- read.table("plink.sexcheck", header=T,as.is=T)
pdf("Gender_check.pdf")
hist(gender[,6],main="Gender", xlab="F")
dev.off()
```

<img width="2100" height="2100" alt="Gender_check" src="https://github.com/user-attachments/assets/17bf8bee-d9e7-4d44-ab50-7a9422b806f8" />

```bash
plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter \
  --impute-sex 0.3 0.7 --make-bed --out plink_sex_impute_03_07_filtered
```

Now comparing the expected meta data (sex_info_sorted_exp.csv) and the plink-imputed data

```bash
awk -v OFS="\t" '{print $1,$4}' plink_sex_impute_03_07_filtered.sexcheck > sex_info_obs.csv
sort sex_info_obs.csv > sex_info_obs_sorted.csv
awk 'NR==FNR {a[$1]=$2; next} {if ($1 in a && a[$1] != $2) print $1, a[$1], $2}' \
  sex_info_obs_sorted.csv ../../sex_info_sorted_exp.csv 
```

there are two mismatching files:
*  WE100: female in plink, male in metadata
*  WE104 : female in plink, male in metadata

Due to the discrepancy, these samples will be deleted.

```bash
grep -E 'WE100|WE104' plink_sex_impute_03_07_filtered.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt
plink --bfile plink_sex_impute_03_07_filtered --remove sex_discrepancy.txt \
  --make-bed --out gold_sex_imputed_2
```

Output: 332632 variants loaded from .bim file.
125 people (64 males, 61 females) loaded from .fam.
--remove: 123 people remaining.

##  Fill missing SNP IDs with chromosome:position information

```bash
plink --bfile gold_sex_imputed_2 --set-missing-var-ids @:# \
  --make-bed --out gold_filled_miss_after_imp_3
```

Output: 332632 variants loaded from .bim file.
86580 missing IDs set. Total genotyping rate is 0.997221.

##  Autosomal + MAF filtering

```bash
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' gold_filled_miss_after_imp_3.bim \
  > snp_1_22_new.txt
plink --bfile gold_filled_miss_after_imp_3 --extract snp_1_22_new.txt \
  --make-bed --out gold_autosomal_only_4
```

Output: 332632 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
--extract: 326870 variants remaining.

Generate a plot of MAF distribution

```bash
plink --bfile gold_autosomal_only_4 --freq --out MAF_check
```

```R
maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)
pdf("MAF_distribution.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()
```

[MAF_distribution.pdf](https://github.com/user-attachments/files/21370118/MAF_distribution.pdf)

```bash
plink --bfile gold_autosomal_only_4 --maf 0.01 --make-bed \
  --out gold_maf_filtered_5
```
Output: 149179 variants removed due to minor allele threshold(s)
177691 variants and 123 people pass filters and QC.

##  Hardy-Weinberg equilibrium filtering HWE)
SNPs that are not in HWE will be filtered out. Check the distribution of HWE p-values of all SNPs. p < 0.001

```bash
plink --bfile gold_maf_filtered_5 --hardy
awk '{if ($9 <   0.0001) print $0}' plink.hwe > plinkzoomhwe.hwe
```

```R
hwe<-read.table (file="plink.hwe", header=TRUE)
pdf("histhwe.pdf")
hist(hwe[,9],main="Histogram HWE")
dev.off()

hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE)
pdf("histhwe_below_theshold.pdf")
hist(hwe_zoom[,9],main="Histogram HWE: strongly deviating SNPs only")
dev.off()
```

[histhwe.pdf](https://github.com/user-attachments/files/21370176/histhwe.pdf)

[histhwe_below_theshold.pdf](https://github.com/user-attachments/files/21370175/histhwe_below_theshold.pdf)

```bash
plink --bfile gold_maf_filtered_5 --hwe 0.0001 --make-bed --out gold_hwe_filtered_6
```

Output: --hwe: 2580 variants removed due to Hardy-Weinberg exact test.
175111 variants and 123 people pass filters and QC.

##  Heterozygosity filtering

The heterozygosity rate of an individual is the proportion of heterozygous genotypes. High levels of heterozygosity within an individual might indicate low sample quality whereas low levels of heterozygosity may be due to inbreeding.

The general threshold of Het rate is 3 standard deviations or ~99.7% of the total data.
*Heterozygosity Rate = (Number of Heterozygous Genotypes) / (Total Number of Genotyped Loci)*

###  Generate a plot of the distribution of the heterozygosity rate of your subjects and remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.

Checks for heterozygosity are performed on a set of SNPs that are not highly correlated. Therefore, to generate a list of non-(highly)correlated SNPs, we prune the SNPs with high LD using the command --indep-pairwise. The parameters 50 5 0.2 stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for an SNP being regressed on all other SNPs simultaneously.
Pruning:
```
plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --indep-pairwise 50 5 0.2 --out pruned
```
Before: 175111
After: 65518 (NOT APPLIED YET, will be removed during relatedness filtering)

Note, we don't delete the file indepSNP.prune.in, we will use this file in later steps of the tutorial.

```
plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --extract pruned.prune.in --het --out R_check
```
This file contains our pruned dataset.

Plot of the heterozygosity rate distribution

```R
het <- read.table("R_check.het", head=TRUE)
pdf("heterozygosity.pdf")
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()
```

A heterozygosity rate peaking between 0.15 and 0.20
The next code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.

```R
het <- read.table("R_check.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)
```

Output of the command above: fail-het-qc.txt:

```
"FID" "IID" "O.HOM." "E.HOM." "N.NM." "F" "HET_RATE" "HET_DST"
"WE1209921" "WE1209921" 49222 54460 65497 -0.4742 0.248484663419699 4.69825524126072
"WE1216321" "WE1216321" 45885 54460 65504 -0.7768 0.299508426966292 7.60408858816003
```

Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
```
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
```
Filter these individuals out:
```
plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --remove het_fail_ind.txt --make-bed --out gold_hetero_filtered_7
```

Output: 175111 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
--remove: 121 people remaining.


