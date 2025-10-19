## Download and prepare the files

```bash
wget "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1" -O all_phase3.pgen.zst
wget 'https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1' -O all_phase3.pvar.zst
wget 'https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1' -O all_phase3.psam

plink2 --zst-decompress all_phase3.pgen.zst all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst all_phase3.pvar

plink2 --pfile all_phase3 --max-alleles 2 --make-pgen --out all_phase3_biall
```
2504 samples (1270 females, 1233 males, 1 ambiguous; 2497 founders) loaded from
all_phase3.psam.
84358431 out of 84805772 variants loaded from all_phase3.pvar.

Extract snps that are present in kazakh dataset: 

```bash
awk '{ print $2 }' 8_exome_only/exomes.bim > exome_snps.txt

plink2 --pfile all_phase3_biall \
  --extract exome_snps.txt \
  --make-bed --out all_phase3_exome
```
49842 variants remaining after main filters.

```bash
plink --bfile 3_matched/all_phase3_exome --snps-only just-acgt \
	--make-bed --out 4_snp/1000_gen_snps
```
49209 variants and 2504 people pass filters and QC.

now doing the same for kazaks

```bash
awk '{print $2}' 4_snp/1000_gen_snps.bim > matched_snps.txt
plink --bfile ../plink_attempt_2.0/8_exome_only/exomes \
	--extract matched_snps.txt --make-bed \
	--out 4_snp/kazakh_matched
```
49209 variants and 108 people pass filters and QC.

## Matching the genome build

```bash 
awk '{print $2, $4}' 4_snp/kazakh_matched.bim > 5_coord_sync/buildmap.txt
plink --bfile 4_snp/1000_gen_snps --update-map 5_coord_sync/buildmap.txt \
	--make-bed --out 5_coord_sync/1000_gen_buildsync
```
49209 values updated.
Warning: Base-pair positions are now unsorted!
49209 variants and 2504 people pass filters and QC.

Sorting the positions:

```bash
plink2 --bfile 5_coord_sync/1000_gen_buildsync --sort-vars --make-pgen \
	--out 5_coord_sync/1000_gen_sorted_pgen
plink2 --pfile 5_coord_sync/1000_gen_sorted_pgen --make-bed \
	--out 5_coord_sync/1000_gen_sorted
```

## Set reference genome

```bash
awk '{print $2, $5}' 5_coord_sync/1000_gen_sorted.bim > ref_alleles.txt
plink --bfile 4_snp/kazakh_matched --reference-allele ref_alleles.txt  \
	--make-bed --out 6_set_ref/kazakh_refset
```
--a1-allele: 49209 assignments attempted, 49178 made.

### Removing impossible A1 alleles

```bash
grep 'Impossible A1 allele' 6_set_ref/kazakh_refset.log | awk '{print $(NF)}' \
	> 6_set_ref/impossible_a1_snps.txt
```
The rs ids have dots in the end, remove them manually 

```bash
plink --bfile 6_set_ref/kazakh_refset --exclude 6_set_ref/impossible_a1_snps.txt \
	--make-bed --out 6_set_ref/kazakh_noimpossible
```
49178 variants and 108 people pass filters and QC.


Exclude them from the 1000 genomes files as well

```bash
awk '{print $2}' 6_set_ref/kazakh_noimpossible.bim > 6_set_ref/snp_list.txt
plink --bfile 5_coord_sync/1000_gen_sorted --extract 6_set_ref/snp_list.txt \
	--make-bed --out 6_set_ref/1000_genomes_noimpossible
```

49178 variants and 2504 people pass filters and QC.

```bash
plink --bfile 6_set_ref/1000_genomes_noimpossible \
  --remove 6_set_ref/1000_genomes_noimpossible.nosex \
  --make-bed --out 6_set_ref/1000_genomes_no_nosex
```

## Strand Flpping

```bash
awk '{print $2, $5, $6}' 6_set_ref/kazakh_noimpossible.bim | sort > 7_flip/kazakh_vars.txt
awk '{print $2, $5, $6}' 6_set_ref/1000_genomes_no_nosex.bim | sort > 7_flip/1000gen_vars.txt
comm -3 7_flip/kazakh_vars.txt 7_flip/1000gen_vars.txt > 7_flip/differences.txt
wc -l 7_flip/differences.txt
```
    4134 7_flip/differences.txt
	
```bash
awk '{print $1}' 7_flip/differences.txt | sort -u > 7_flip/to_flip.txt
plink --bfile 6_set_ref/kazakh_noimpossible --flip 7_flip/to_flip.txt \
	--make-bed --out 7_flip/kazakh_flipped
```
--flip: 2067 SNPs flipped.

Check for differences again. We get     3640 7_flip/differences_2.txt

Removing them: 

```bash
awk '{print $1}' 7_flip/differences_2.txt | sort -u > 8_exclude_diff/final_exclude.txt

plink --bfile 7_flip/kazakh_flipped --exclude 8_exclude_diff/final_exclude.txt \
	--make-bed --out 8_exclude_diff/kazakh_final

plink --bfile 6_set_ref/1000_genomes_no_nosex \
	--exclude 8_exclude_diff/final_exclude.txt --make-bed \
	--out 8_exclude_diff/1000_gen_final
```


___________


```bash
plink --bfile 4_snp/1000_gen_snps --missing --out 5_missinges_freq/missingness
plink --bfile 4_snp/1000_gen_snps --freq --out 5_missinges_freq/freqs
```
