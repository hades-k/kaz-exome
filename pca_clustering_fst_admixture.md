## Subsetting 

```bash
grep 'BEB' filtered_igsr_samples.tsv | shuf -n 100 \
	> subsets/Bengali_subset.tsv
```
Repeat this for alls populations

```bash
cat subsets/*_subset.tsv | awk '{print "0 "$1}' > ids_to_keep.txt

plink --bfile ../1000_gen_qc/8_exclude_diff/1000_gen_final  \
  --keep subsets/ids_to_keep.txt --make-bed \
  --out 1000gen_subset
```
47358 variants and 1164 people pass filters and QC.

## Merge and LD prune

```bash
plink --bfile ../1000_gen_qc/8_exclude_diff/kazakh_final \
	--bmerge 1000gen_subset --make-bed --out merged

plink --bfile merged --indep-pairwise 50 5 0.2 --out pruned

plink --bfile merged --extract pruned.prune.in --make-bed --out merged_pruned
```
19109 of 47358 variants removed.
28249 variants remaining.

## PCA

```bash
plink --bfile merged_pruned --pca 10 --out pca/pca
```

#### Exact population split:

- Bengali: 86 samples
- British: 91 samples
- Dai Chinese: 93 samples
- Finnish: 99 samples
- Gujarati: 100 samples
- Han Chinese: 100 samples
- Iberian: 100 samples
- Japanese: 100 samples
- Kazakh: 108 samples
- Kinh Vietnamese: 99 samples
- Punjabi: 96 samples
- Tamil: 100 samples
- Toscani: 100 samples

Total samples: 1272

![scree](/NEW_ATTEMPT/pca_plots/scree_plot.png)

#### Variance explained by each PC:
- PC1: 56.3%
- PC2: 16.0%
- PC3: 6.5%
- PC4: 3.8%
- PC5: 3.4%
- PC6: 3.0%
- PC7: 2.9%
- PC8: 2.8%
- PC9: 2.7%
- PC10: 2.7%

#### Cumulative variance explained:
- PC1 to PC1: 56.3%
- PC1 to PC2: 72.2%
- PC1 to PC3: 78.7%
- PC1 to PC4: 82.5%
- PC1 to PC5: 85.9%
- PC1 to PC6: 88.9%
- PC1 to PC7: 91.8%
- PC1 to PC8: 94.6%
- PC1 to PC9: 97.3%
- PC1 to PC10: 100.0%

![3x3](/NEW_ATTEMPT/pca_plots/3x3_pca.png)

![supergroup](/NEW_ATTEMPT/pca_plots/pca_supergroup.png)

## K-means clustering

Based on the scree plot, using PC1-PC3, 78.7% cumulative variance explained.

Using elbow method and silhouette score to determine the right number of clusters

![elbow](/NEW_ATTEMPT/clustering_plots/elbow_plot.png)

Scores:
- 3: 0.6389660469195052
- 4: 0.7671864998612629
- 5: 0.8212778836471257
- 6: 0.7684536786444122

K = 5 has the highest silhouette score.

![cluster](/NEW_ATTEMPT/clustering_plots/5_kmeans_cluster.png)
![3d_cluster](/NEW_ATTEMPT/clustering_plots/3d_clustering.png)

Measuring euclidean distances to the cluster containing most Kazakhs:

Kazakh cluster: 3
- Kazakh → Cluster 0: 0.0892
- Kazakh → Cluster 1: 0.0918
- Kazakh → Cluster 2: 0.0618
- Kazakh → Cluster 4: 0.1131

## Fst calculations

Create a vcf from plink files and write a popmap

```bash
plink --bfile merged_pruned --recode vcf-iid --out merged_pruned

awk -F'\t' '{print $1, $5}' *_subset.tsv > id2pop.txt
```

Fixing the spaces (changing 'Han Chinese' to 'Han_Chinese') to avoid counfusion

```python
input_file = 'id2pop.txt'
output_file = 'id2pop_fix.txt'

with open(input_file, 'r', encoding='utf-8') as f:
    lines = f.readlines()

fixed_lines = [
    line.replace(" Chinese", "_Chinese").replace(" Vietnamese", "_Vietnamese")
    for line in lines
]
with open(output_file, 'w', encoding='utf-8') as f:
    f.writelines(fixed_lines)
```

```bash
bcftools query -l merged_pruned.vcf | \
	awk 'NR==FNR{a[$1];next} !($1 in a){print $1 "\tKazakh"}' \
	id2pop_fix.txt - >> id2pop_fix.txt
```

Use PGDSpider to convert vcf into arp

```bash
java -Xmx10g -Xms512m -jar PGDSpider3.jar
```
Using following parameters in Arlequin: 

        Compute pairwise differences
        No. of permutations for significance = 10000
        No. of permutations for Mantel test  = 1000

       Distance matrix:
            Compute distance matrix
            Molecular distance : Pairwise difference
            Gamma a value      = 0

We get the following results: 

|Label  |	Population name|
|-----  |	---------------|
|  1|	Bengali|
|  2|	Gujarati|
|  3|	Punjabi|
|  4|	Tamil|
|  5|	British|
|  6|	Finnish|
 | 7|	Iberian|
|  8|	Toscani|
 | 9|	Dai_Chinese|
  |10|	Han_Chinese|
|  11|	Japanese|
 | 12|	Kinh_Vietnamese|
 | 13|	Kazakh|

#### Population pairwise FSTs
```
Distance method: Pairwise difference
                     1         2         3         4         5         6         7         8         9        10        11        12        13
           1   0.00000
           2   0.00467   0.00000
           3   0.00396   0.00375   0.00000
           4   0.00243   0.00481   0.00411   0.00000
           5   0.04099   0.03602   0.02996   0.04323   0.00000
           6   0.04023   0.03676   0.03127   0.04331   0.00706   0.00000
           7   0.04081   0.03598   0.03003   0.04267   0.00265   0.01083   0.00000
           8   0.03994   0.03450   0.02874   0.04149   0.00420   0.01253   0.00164   0.00000
           9   0.05802   0.07578   0.07293   0.06994   0.10824   0.10105   0.10743   0.10747   0.00000
          10   0.05579   0.07296   0.07000   0.06759   0.10391   0.09627   0.10340   0.10347   0.00951   0.00000
          11   0.05813   0.07474   0.07218   0.06943   0.10645   0.09852   0.10560   0.10598   0.01800   0.00699   0.00000
          12   0.05301   0.07034   0.06751   0.06474   0.10292   0.09590   0.10234   0.10229   0.00192   0.00659   0.01488   0.00000
          13   0.02208   0.02930   0.02530   0.02899   0.03920   0.03523   0.03978   0.03975   0.03600   0.02648   0.02833   0.03142   0.00000
```

#### FST P values

Number of permutations : 10100
```
                       1                 2                 3                 4                 5                 6                 7                 8                 9                10                11                12                13
           1           *
           2   0.00000+-0.0000           *
           3   0.00000+-0.0000   0.00000+-0.0000           *
           4   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
           5   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
           6   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
           7   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
           8   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
           9   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
          10   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
          11   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
          12   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
          13   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000   0.00000+-0.0000           *
```

![fst](/NEW_ATTEMPT/fst_matrix.png)

## ADMIXTURE

```bash
admixture --cv merged_pruned.bed 2 | tee log2.out
```

CV error rates:
| K | Error rate |
| - | ---------- |
| 2 | 0.31486 |
| 3 | 0.30929 |
| 4 | 0.30832 |
| 5 | 0.30805 |
| 6 | 0.30820 |
| 7 | 0.30909 |
| 8 | 0.30902 |
| 9 | 0.30948 |
| 10 | 0.31008 |

The lowest error rate is at K = 5.

![cv](/NEW_ATTEMPT/admixture_plots/cv_er.png)

Creating population map file for from id2pop_fix.txt

```bash
awk '
  NR==FNR { pop[$1] = $2; next }   
  {
    iid = $2
    p = (iid in pop ? pop[iid] : "UNK")
    printf "%s\t%s\n", p, iid   
  }
' subsets/id2pop_fix.txt merged_pruned.fam > admixture/popmap.txt
```

Plotting 
```python
q_file = "merged_pruned.2.Q"
pop_file = "popmap.txt" 

df_q = pd.read_csv(q_file, sep=r"\s+", header=None)              
df_pop = pd.read_csv(pop_file, sep=r"\s+", header=None, names=["Population", "Individual"])
pop_order = ['Punjabi', 'Gujarati', 'Tamil', 'Bengali', 'Kinh_Vietnamese',
               'Dai_Chinese', 'Han_Chinese', 'Japanese', 'Kazakh',
               'Finnish', 'British', 'Iberian', 'Toscani']

df_q["Population"] = pd.Categorical(
    df_pop["Population"],
    categories=pop_order,  # explicit custom order
    ordered=True)

df_sorted = df_q.sort_values(by="Population").reset_index(drop=True)

df_plot = df_sorted.drop(columns=["Population"])               
N, K = df_plot.shape
x_pos = np.arange(N)
colors = plt.cm.tab10.colors[:K]

fig, ax = plt.subplots(figsize=(17,4))
left = np.zeros(N, dtype=float)
for i in range(K):
    vals = df_plot.iloc[:, i].to_numpy(float)
    ax.bar(x_pos, vals, bottom=left, color=colors[i], width=1.0)
    left += vals

pop_labels = df_sorted["Population"].to_numpy()
groups = pd.Series(pop_labels).groupby(pop_labels, sort=False).size()
starts = groups.cumsum().shift(fill_value=0).to_numpy()
ends = groups.cumsum().to_numpy()

for boundary in ends[:-1]:
    ax.axvline(x=boundary - 0.5, color="black", linestyle="--", linewidth=1)

for pop, s, e in zip(groups.index, starts, ends):
    mid = (s + e) / 2.0
    ax.text(mid - 0.5, -0.05, pop, ha="center", transform=ax.get_xaxis_transform())

ax.set_title(f"Admixture Proportions (K={K})")
ax.set_ylabel("Ancestry Proportion")

legend_patches = [mpatches.Patch(color=colors[i], label=f"Ancestry {i+1}") for i in range(K)]
ax.legend(handles=legend_patches, loc="upper right")

ax.set_xticks([])
plt.tight_layout()
plt.show()
```
![k2](/NEW_ATTEMPT/admixture_plots/k2.png)
![k3](/NEW_ATTEMPT/admixture_plots/k3.png)
![k4](/NEW_ATTEMPT/admixture_plots/k4.png)
![k5](/NEW_ATTEMPT/admixture_plots/k5.png)
![k6](/NEW_ATTEMPT/admixture_plots/k6.png)
![k7](/NEW_ATTEMPT/admixture_plots/k7.png)
![k8](/NEW_ATTEMPT/admixture_plots/k8.png)
![k9](/NEW_ATTEMPT/admixture_plots/k9.png)
![k10](/NEW_ATTEMPT/admixture_plots/k10.png)







