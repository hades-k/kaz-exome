(base) Mac:plink_copy aidana$ plink --vcf ../MS_KAZ_WE_125_kaz_samples.MSHC.gold.vcf --make-bed --out PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold.log.
Options in effect:
  --make-bed
  --out PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
  --vcf ../MS_KAZ_WE_125_kaz_samples.MSHC.gold.vcf

24576 MB RAM detected; reserving 12288 MB for main workspace.
--vcf: PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold-temporary.bed +
PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold-temporary.bim +
PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold-temporary.fam written.
500376 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.939852.
500376 variants and 125 people pass filters and QC.
Note: No phenotypes present.
--make-bed to PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold.bed +
PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold.bim +
PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold --missing
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink.log.
Options in effect:
  --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
  --missing

24576 MB RAM detected; reserving 12288 MB for main workspace.
500376 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.939852.
--missing: Sample missing data report written to plink.imiss, and variant-based
missing data report written to plink.lmiss.
(base) Mac:plink_copy aidana$ plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold --geno 0.2 --make-bed --out PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter.log.
Options in effect:
  --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
  --geno 0.2
  --make-bed
  --out PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter

24576 MB RAM detected; reserving 12288 MB for main workspace.
500376 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.939852.
50618 variants removed due to missing genotype data (--geno).
449758 variants and 125 people pass filters and QC.
Note: No phenotypes present.
--make-bed to PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter.bed +
PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter.bim +
PLINK_MS_KAZ_WE_125_kaz_samples_0.2_miss_filter.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold --geno 0.02 --make-bed --out PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter.log.
Options in effect:
  --bfile PLINK_MS_KAZ_WE_125_kaz_samples.MSHC.gold
  --geno 0.02
  --make-bed
  --out PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter

24576 MB RAM detected; reserving 12288 MB for main workspace.
500376 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.939852.
167744 variants removed due to missing genotype data (--geno).
332632 variants and 125 people pass filters and QC.
Note: No phenotypes present.
--make-bed to PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter.bed +
PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter.bim +
PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter --check-sex
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink.log.
Options in effect:
  --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter
  --check-sex

24576 MB RAM detected; reserving 12288 MB for main workspace.
332632 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997261.
332632 variants and 125 people pass filters and QC.
Note: No phenotypes present.
--check-sex: 5666 Xchr and 0 Ychr variant(s) scanned, 125 problems detected.
Report written to plink.sexcheck .
(base) Mac:plink_copy aidana$ 
(base) Mac:plink_copy aidana$ plink --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter --impute-sex 0.3 0.7 --make-bed --out plink_sex_impute_03_07_filtered
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink_sex_impute_03_07_filtered.log.
Options in effect:
  --bfile PLINK_MS_KAZ_WE_125_kaz_samples_0.02_miss_filter
  --impute-sex 0.3 0.7
  --make-bed
  --out plink_sex_impute_03_07_filtered

24576 MB RAM detected; reserving 12288 MB for main workspace.
332632 variants loaded from .bim file.
125 people (0 males, 0 females, 125 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink_sex_impute_03_07_filtered.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 125 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997261.
332632 variants and 125 people pass filters and QC.
Note: No phenotypes present.
--impute-sex: 5666 Xchr and 0 Ychr variant(s) scanned, all sexes imputed.
Report written to plink_sex_impute_03_07_filtered.sexcheck .
--make-bed to plink_sex_impute_03_07_filtered.bed +
plink_sex_impute_03_07_filtered.bim + plink_sex_impute_03_07_filtered.fam ...
done.
(base) Mac:plink_copy aidana$ awk -v OFS="\t" '{print $1,$4}' plink_sex_impute_03_07_filtered.sexcheck > sex_info_obs.csv
(base) Mac:plink_copy aidana$ sort sex_info_obs.csv > sex_info_obs_sorted.csv
(base) Mac:plink_copy aidana$ awk 'NR==FNR {a[$1]=$2; next} {if ($1 in a && a[$1] != $2) print $1, a[$1], $2}' sex_info_obs_sorted.csv ../../sex_info_sorted_exp.csv 
awk: can't open file ../../sex_info_sorted_exp.csv
 input record number 126, file ../../sex_info_sorted_exp.csv
 source line number 1
(base) Mac:plink_copy aidana$ grep -E 'WE100|WE104' plink_sex_impute_03_07_filtered.sexcheck | awk '{print$1,$2}' > sex_discrepancy.txt
(base) Mac:plink_copy aidana$ plink --bfile plink_sex_impute_03_07_filtered --remove sex_discrepancy.txt --make-bed --out gold_sex_imputed_2
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_sex_imputed_2.log.
Options in effect:
  --bfile plink_sex_impute_03_07_filtered
  --make-bed
  --out gold_sex_imputed_2
  --remove sex_discrepancy.txt

24576 MB RAM detected; reserving 12288 MB for main workspace.
332632 variants loaded from .bim file.
125 people (64 males, 61 females) loaded from .fam.
--remove: 123 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3419 het. haploid genotypes present (see gold_sex_imputed_2.hh ); many
commands treat these as missing.
Total genotyping rate in remaining samples is 0.997221.
332632 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_sex_imputed_2.bed + gold_sex_imputed_2.bim +
gold_sex_imputed_2.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile gold_sex_imputed_2 --set-missing-var-ids @:# --make-bed --out gold_filled_miss_after_imp_3
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_filled_miss_after_imp_3.log.
Options in effect:
  --bfile gold_sex_imputed_2
  --make-bed
  --out gold_filled_miss_after_imp_3
  --set-missing-var-ids @:#

24576 MB RAM detected; reserving 12288 MB for main workspace.
332632 variants loaded from .bim file.
86580 missing IDs set.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 3419 het. haploid genotypes present (see
gold_filled_miss_after_imp_3.hh ); many commands treat these as missing.
Total genotyping rate is 0.997221.
332632 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_filled_miss_after_imp_3.bed +
gold_filled_miss_after_imp_3.bim + gold_filled_miss_after_imp_3.fam ... done.
(base) Mac:plink_copy aidana$ awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' gold_filled_miss_after_imp_3.bim > snp_1_22_new.txt
(base) Mac:plink_copy aidana$ plink --bfile gold_filled_miss_after_imp_3 --extract snp_1_22_new.txt --make-bed --out gold_autosomal_only_4
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_autosomal_only_4.log.
Options in effect:
  --bfile gold_filled_miss_after_imp_3
  --extract snp_1_22_new.txt
  --make-bed
  --out gold_autosomal_only_4

24576 MB RAM detected; reserving 12288 MB for main workspace.
332632 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
--extract: 326870 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997252.
326870 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_autosomal_only_4.bed + gold_autosomal_only_4.bim +
gold_autosomal_only_4.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile gold_autosomal_only_4 --freq --out MAF_check
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to MAF_check.log.
Options in effect:
  --bfile gold_autosomal_only_4
  --freq
  --out MAF_check

24576 MB RAM detected; reserving 12288 MB for main workspace.
326870 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997252.
--freq: Allele frequencies (founders only) written to MAF_check.frq .
(base) Mac:plink_copy aidana$ plink --bfile gold_autosomal_only_4 --maf 0.01 --make-bed --out gold_maf_filtered_5
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_maf_filtered_5.log.
Options in effect:
  --bfile gold_autosomal_only_4
  --maf 0.01
  --make-bed
  --out gold_maf_filtered_5

24576 MB RAM detected; reserving 12288 MB for main workspace.
326870 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997252.
149179 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
177691 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_maf_filtered_5.bed + gold_maf_filtered_5.bim +
gold_maf_filtered_5.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile gold_maf_filtered_5 --hardy
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to plink.log.
Options in effect:
  --bfile gold_maf_filtered_5
  --hardy

24576 MB RAM detected; reserving 12288 MB for main workspace.
177691 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997031.
--hardy: Writing Hardy-Weinberg report (founders only) to plink.hwe ... done.
(base) Mac:plink_copy aidana$ awk '{if ($9 <   0.0001) print $0}' plink.hwe > plinkzoomhwe.hwe
(base) Mac:plink_copy aidana$ plink --bfile gold_maf_filtered_5 --hwe 0.0001 --make-bed --out gold_hwe_filtered_6
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_hwe_filtered_6.log.
Options in effect:
  --bfile gold_maf_filtered_5
  --hwe 0.0001
  --make-bed
  --out gold_hwe_filtered_6

24576 MB RAM detected; reserving 12288 MB for main workspace.
177691 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.997031.
--hwe: 2580 variants removed due to Hardy-Weinberg exact test.
175111 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_hwe_filtered_6.bed + gold_hwe_filtered_6.bim +
gold_hwe_filtered_6.fam ... done.
(base) Mac:plink_copy aidana$ plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --indep-pairwise 50 5 0.2 --out pruned
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to pruned.log.
Options in effect:
  --bfile ../6_HWE_filter/gold_hwe_filtered_6
  --indep-pairwise 50 5 0.2
  --out pruned

24576 MB RAM detected; reserving 12288 MB for main workspace.
Error: Failed to open ../6_HWE_filter/gold_hwe_filtered_6.bed.
(base) Mac:plink_copy aidana$ plink --bfile gold_hwe_filtered_6 --indep-pairwise 50 5 0.2 --out pruned
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to pruned.log.
Options in effect:
  --bfile gold_hwe_filtered_6
  --indep-pairwise 50 5 0.2
  --out pruned

24576 MB RAM detected; reserving 12288 MB for main workspace.
175111 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.99704.
175111 variants and 123 people pass filters and QC.
Note: No phenotypes present.
Pruned 12051 variants from chromosome 1, leaving 6646.
Pruned 8533 variants from chromosome 2, leaving 4824.
Pruned 6847 variants from chromosome 3, leaving 4001.
Pruned 4841 variants from chromosome 4, leaving 3036.
Pruned 5923 variants from chromosome 5, leaving 3407.
Pruned 5413 variants from chromosome 6, leaving 3331.
Pruned 5908 variants from chromosome 7, leaving 3668.
Pruned 4560 variants from chromosome 8, leaving 2837.
Pruned 4907 variants from chromosome 9, leaving 3059.
Pruned 5009 variants from chromosome 10, leaving 3107.
Pruned 7071 variants from chromosome 11, leaving 3824.
Pruned 4027 variants from chromosome 12, leaving 2806.
Pruned 1841 variants from chromosome 13, leaving 1456.
Pruned 3019 variants from chromosome 14, leaving 2002.
Pruned 3483 variants from chromosome 15, leaving 2071.
Pruned 4631 variants from chromosome 16, leaving 2601.
Pruned 5551 variants from chromosome 17, leaving 3225.
Pruned 1595 variants from chromosome 18, leaving 1254.
Pruned 6199 variants from chromosome 19, leaving 3725.
Pruned 3213 variants from chromosome 20, leaving 2003.
Pruned 1689 variants from chromosome 21, leaving 959.
Pruned 3282 variants from chromosome 22, leaving 1676.
Pruning complete.  109593 of 175111 variants removed.
Marker lists written to pruned.prune.in and pruned.prune.out .
(base) Mac:plink_copy aidana$ Rscript --no-save check_heterozygosity_rate.R
Fatal error: cannot open file 'check_heterozygosity_rate.R': No such file or directory
(base) Mac:plink_copy aidana$ Rscript --no-save check_heterozygosity_rate.R
Error in file(file, "rt") : cannot open the connection
Calls: read.table -> file
In addition: Warning message:
In file(file, "rt") :
  cannot open file 'R_check.het': No such file or directory
Execution halted
(base) Mac:plink_copy aidana$ plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --extract pruned.prune.in --het --out R_check
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to R_check.log.
Options in effect:
  --bfile ../6_HWE_filter/gold_hwe_filtered_6
  --extract pruned.prune.in
  --het
  --out R_check

24576 MB RAM detected; reserving 12288 MB for main workspace.
Error: Failed to open ../6_HWE_filter/gold_hwe_filtered_6.bed.
(base) Mac:plink_copy aidana$ plink --bfile gold_hwe_filtered_6 --extract pruned.prune.in --het --out R_check
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to R_check.log.
Options in effect:
  --bfile gold_hwe_filtered_6
  --extract pruned.prune.in
  --het
  --out R_check

24576 MB RAM detected; reserving 12288 MB for main workspace.
175111 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
--extract: 65518 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 123 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.99661.
65518 variants and 123 people pass filters and QC.
Note: No phenotypes present.
--het: 65518 variants scanned, report written to R_check.het .
(base) Mac:plink_copy aidana$ Rscript --no-save check_heterozygosity_rate.R
null device 
          1 
(base) Mac:plink_copy aidana$ Rscript --no-save heterozygosity_outliers_list.R
(base) Mac:plink_copy aidana$ sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
(base) Mac:plink_copy aidana$ plink --bfile ../6_HWE_filter/gold_hwe_filtered_6 --remove het_fail_ind.txt --make-bed --out gold_hetero_filtered_7
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_hetero_filtered_7.log.
Options in effect:
  --bfile ../6_HWE_filter/gold_hwe_filtered_6
  --make-bed
  --out gold_hetero_filtered_7
  --remove het_fail_ind.txt

24576 MB RAM detected; reserving 12288 MB for main workspace.
Error: Failed to open ../6_HWE_filter/gold_hwe_filtered_6.bed.
(base) Mac:plink_copy aidana$ plink --bfile gold_hwe_filtered_6 --remove het_fail_ind.txt --make-bed --out gold_hetero_filtered_7
PLINK v1.90p 64-bit (6 Sep 2023)               www.cog-genomics.org/plink/1.9/
(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to gold_hetero_filtered_7.log.
Options in effect:
  --bfile gold_hwe_filtered_6
  --make-bed
  --out gold_hetero_filtered_7
  --remove het_fail_ind.txt

24576 MB RAM detected; reserving 12288 MB for main workspace.
175111 variants loaded from .bim file.
123 people (64 males, 59 females) loaded from .fam.
--remove: 121 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 121 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.996996.
175111 variants and 121 people pass filters and QC.
Note: No phenotypes present.
--make-bed to gold_hetero_filtered_7.bed + gold_hetero_filtered_7.bim +
gold_hetero_filtered_7.fam ... done.
