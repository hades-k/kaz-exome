# Continuation of 1000_gen_v2

Creating a subset with 100 samples from each reference population 

```bash
grep 'BEB' filtered_igsr_samples.tsv | shuf -n 100 \
	> 100_samp_subset/Bengali_subset.tsv
```

We do not have enough samples of some of them, however the differences are not so large: 

```
      86 100_samp_subset/Bengali_subset.tsv
      91 100_samp_subset/British_subset.tsv
      93 100_samp_subset/Dai_Chinese_subset.tsv
      99 100_samp_subset/Finnish_subset.tsv
     100 100_samp_subset/Gujarati_subset.tsv
     100 100_samp_subset/Han_Chinese_subset.tsv
     100 100_samp_subset/Iberian_subset.tsv
     100 100_samp_subset/Japanese_subset.tsv
      99 100_samp_subset/Kinh_Vietnamesee_subset.tsv
      96 100_samp_subset/Punjabi_subset.tsv
     100 100_samp_subset/Tamil_subset.tsv
     100 100_samp_subset/Toscani_subset.tsv
```

```bash 
cat 100_samp_subset/*_subset.tsv | awk '{print "0 "$1}' \
  > 100_samp_subset/ids_to_keep.txt

plink --bfile 7_final_exclude/1000_gen_final \
  --keep 100_samp_subset/ids_to_keep.txt \
	--make-bed --out 100_samp_subset/1000gen_subset

plink --bfile 11_remove_outliers/kazakh_no_outliers \
  --bmerge 00_samp_subset/1000gen_subset --make-bed  \
  --out 100_samp_subset/merged_large

plink --bfile 100_samp_subset/merged_large --pca 10 \
  --out 100_samp_subset/pca
```

Admixture 

```bash 
admixture --cv 100_samp_subset/merged_large.bed 2 | tee 100_samp_subset/log2.out
```
