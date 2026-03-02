## Missingness and Minor Allele Frequency check

```bash
plink --vcf ../MS_KAZ_WE_125_kaz_samples.MSHC.gold.vcf \
  --make-bed --out import/kaz_wes_125_gold

plink --bfile 0_import/kaz_wes_125_gold --missing
plink --bfile 0_import/kaz_wes_125_gold --freq
```
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
Total genotyping rate is 0.939852.

![lmiss](/qc_plots/lmiss.png)
![imiss](/qc_plots/imiss.png)
![maf](/qc_plots/maf.png)

```bash 
plink --bfile 0_import/kaz_wes_125_gold \
  --geno 0.02 --maf 0.01 --mind 0.2 \
  --make-bed --out 1_missing_maf/missingness_and_maf
```
2 people removed due to missing genotype data (--mind).
155373 variants removed due to missing genotype data (--geno).
156122 variants removed due to minor allele threshold(s)
188881 variants and 123 people pass filters and QC.

WE113	and WE135	were removed due to missing genotypes.

## Sex Imputation
```bash
plink --bfile 1_missing_maf/missingness_and_maf \
  --check-sex --out 2_sexcheck/initial_sexcheck

plink --bfile 1_missing_maf/missingness_and_maf \
  --impute-sex 0.2 0.8 --make-bed --out 2_sexcheck/sex_imputed
```
![sexcheck](/qc_plots/sexcheck.png)

And now compare metadata sex with the imputed one

```python
metadata = pd.read_csv('metadata_sex.csv', sep = ';', header = None)
metadata.columns = ['IID', 'SEX_verb']

sex_map = {"male": 1, "female": 2}
metadata["SEX"] = metadata["SEX_verb"].map(sex_map).fillna(0).astype(int)

metadata[['IID', 'SEX']].to_csv('sex_sample_info.csv', index = False)
```

```bash
sort sex_sample_info.csv > 2_sexcheck/sex_info_sorted_exp.csv

awk -v OFS="\t" '{print $2,$4}' 2_sexcheck/sex_imputed.sexcheck \
  > 2_sexcheck/sex_info_obs.csv

sort 2_sexcheck/sex_info_obs.csv > 2_sexcheck/sex_info_obs_sorted.csv
```

```python
exp = pd.read_csv("2_sexcheck/sex_info_sorted_exp.csv")
obs = pd.read_csv("2_sexcheck/sex_info_obs_sorted.csv", sep = '\t')

combined = pd.merge(exp, obs, on="IID", how="outer")
mismatch = combined[combined.SEX != combined.SNPSEX]
```

there are two mismatching samples:

- WE100: female in plink, male in metadata
- WE104: female in plink, male in metadata

In addition, for one sample sex was not imputed: WE047
  
Due to the discrepancy, these samples will be deleted.

```bash
plink --bfile 2_sexcheck/sex_imputed \
  --remove 2_sexcheck/remove_samples.txt \
  --make-bed --out 2_sexcheck/no_sex_discrepancy
```
Total genotyping rate in remaining samples is 0.997417.
188881 variants and 120 people pass filters and QC.

## Keep autosomal SNPs and fill in missing ids

```bash
plink --bfile 2_sexcheck/no_sex_discrepancy /
  --set-missing-var-ids @:# --make-bed /
  --out 3_autosomal_only/filled

awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' \
  3_autosomal_only/filled.bim \
  > 3_autosomal_only/snp_1_22_new.txt

plink --bfile 3_autosomal_only/filled \
  --extract 3_autosomal_only/snp_1_22_new.txt \
  --make-bed --out 3_autosomal_only/autosomal
```
Total genotyping rate is 0.997442.
185685 variants and 120 people pass filters and QC.


## Hardy - Weinberg equillibrium

```bash
plink --bfile 3_autosomal_only/autosomal \
  --hardy --out 4_hwe/hardy

plink --bfile 3_autosomal_only/autosomal \
  --hwe 0.0001 --make-bed --out 4_hwe/hwe_filtered
```
2673 variants removed due to Hardy-Weinberg exact test.
183012 variants and 120 people pass filters and QC.

![hwe](/qc_plots/hwe.png)

## Heterozygocity check

```bash
plink --bfile 4_hwe/hwe_filtered --het --out 5_het/het_stats
```

```python
het = pd.read_csv('5_het/het_stats.het', delim_whitespace=True)

mean_f = het['F'].mean()
sd_f   = het['F'].std()

upper = mean_f + 3*sd_f
lower = mean_f - 3*sd_f

outliers = het[(het['F'] > upper) | (het['F'] < lower)][['FID','IID','F']]
```
Outliers:
|   |       FID |       IID  |     F|
|- |-|-|-|
|87|  WE1209921|  WE1209921| -0.4858|
|93|  WE1216321 | WE1216321 |-0.8092|

These two samples will be removed 

```bash
plink --bfile 4_hwe/hwe_filtered \
  --remove 5_het/remove_samples.txt \
  --make-bed --out 5_het/het_filtered
```
Total genotyping rate in remaining samples is 0.997415.
183012 variants and 118 people pass filters and QC.

## Runs of Homozygocity 

For this, using default parameters from PLINK manual 

```bash
plink --bfile 5_het/het_filtered \
  --homozyg \
  --homozyg-window-snp 50 \
  --homozyg-snp 100 \
  --homozyg-kb 1000 \
  --homozyg-gap 1000 \
  --homozyg-density 50 \
  --homozyg-window-het 1 \
  --homozyg-window-missing 5 \
  --homozyg-window-threshold 0.05 \
  --out 6_roh/roh_stats
```

## Relatedness check

```bash
plink --bfile 5_het/het_filtered --genome \
  --min 0.2 --out 7_ibd/pihat_0.2

awk 'NR>1 { print $2, $4, $10 }' \
    7_ibd/pihat_0.2.genome > 7_ibd/0.2related_pairs.txt
```
We get the following pairs of related samples:

- WE002 WE091 0.5000
- WE002 WE092 0.5061
- WE015 WE057 0.4552
- WE017 WE019 0.5000
- WE018 WE019 0.5000
- WE033 WE109 0.4971
- WE058 WE060 0.5000
- WE060 WE117 0.2364
- WE063 WE070 0.4845
- WE089 WE090 0.4017
- WE091 WE092 0.4934
- WE124 WE129 0.4426
- WE127 WE128 0.4369

Therefore removing the following samples: WE091, WE092,	WE057,	WE019, WE109,	WE060, WE090,	WE129, WE128, WE070

```bash
plink --bfile 5_het/het_filtered \
  --remove 7_ibd/ids_to_remove.txt \
  --make-bed --out 7_ibd/related_removed
```
183012 variants and 108 people pass filters and QC.

## Extract exome only snps

Download reference bef file https://support.illumina.com/downloads/Illumina-dna-prep-exome-20-bed-files.html

```bash
plink --bfile 7_ibd/related_removed \
  --extract range hg19_Twist_Bioscience_for_Illumina_Exome_2_5.bed \
  --allow-extra-chr --make-bed --out 8_exome_only/exomes
```
128373 variants excluded.
54639 variants and 108 people pass filters and QC.



