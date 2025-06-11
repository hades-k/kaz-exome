import pandas as pd

mine = 'my_vcf/annovar_res/WE001_annotated.hg38_multianno.txt'
df = pd.read_csv(mine, sep='\t')
my_rsids = set()
for rsid in df['avsnp150']:
    if isinstance(rsid, str) and rsid.startswith('rs'):
        my_rsids.add(rsid)
print(f'Count of rs variants in my vcf : {len(my_rsids)}')

gold = 'MS_KAZ_WE_125_kaz_samples.MSHC.gold.WE001.vcf'
gold_rsids = set ()
with open(gold) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) > 2 and fields[2].startswith('rs'):
            gold_rsids.add(fields[2])

print(f'Count of rs variants in VCF gold: {len(gold_rsids)}')

overlap = my_rsids.intersection(gold_rsids)
print (f'Number of overlap variants between vcfs: {len(overlap)}')
