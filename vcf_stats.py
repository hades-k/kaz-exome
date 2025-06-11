import sys
from collections import defaultdict
from typing import Dict, Any

vcf_file = 'my_vcf/joint_filtered.vcf'

def calculate_vcf_stats(vcf_file: str) -> Dict[str, Any]:
    stats = {
        'total_variants': 0,
        'variant_types': defaultdict(int),
        'chromosomes': defaultdict(int),
        'quality_stats': {
            'sum': 0,
            'count': 0,
            'min': float('inf'),
            'max': float('-inf')
        },
        'filter_stats': defaultdict(int),
        'transition_transversion': {'transitions': 0, 'transversions': 0}
    }
    
    transitions = {('A','G'), ('G','A'), ('C','T'), ('T','C')}
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
                
            chrom = fields[0]
            ref = fields[3]
            alt = fields[4].split(',')[0]
            
            # Basic counts
            stats['total_variants'] += 1
            stats['chromosomes'][chrom] += 1
            
            # Variant type
            if len(ref) == len(alt) == 1:  # SNP
                stats['variant_types']['SNP'] += 1
                if (ref.upper(), alt.upper()) in transitions:
                    stats['transition_transversion']['transitions'] += 1
                else:
                    stats['transition_transversion']['transversions'] += 1
            elif len(ref) != len(alt):
                stats['variant_types']['INDEL'] += 1
            else:
                stats['variant_types']['OTHER'] += 1

            try:
                if len(fields) > 5:
                    qual = float(fields[5])
                    if qual != '.':
                        stats['quality_stats']['sum'] += qual
                        stats['quality_stats']['count'] += 1
                        stats['quality_stats']['min'] = min(stats['quality_stats']['min'], qual)
                        stats['quality_stats']['max'] = max(stats['quality_stats']['max'], qual)
            except (ValueError, IndexError):
                pass

            try:
                if len(fields) > 6:
                    filter_field = fields[6]
                    if filter_field == '.' or filter_field == 'PASS':
                        stats['filter_stats']['PASS'] += 1
                    else:
                        for f in filter_field.split(';'):
                            stats['filter_stats'][f] += 1
            except IndexError:
                pass

    if stats['quality_stats']['count'] > 0:
        stats['quality_stats']['average'] = stats['quality_stats']['sum'] / stats['quality_stats']['count']
        if stats['quality_stats']['min'] == float('inf'):
            stats['quality_stats']['min'] = 'N/A'
        if stats['quality_stats']['max'] == float('-inf'):
            stats['quality_stats']['max'] = 'N/A'

    if stats['transition_transversion']['transversions'] > 0:
        stats['transition_transversion']['ratio'] = (
            stats['transition_transversion']['transitions'] /
            stats['transition_transversion']['transversions']
        )
    
    return stats

def print_stats(stats: Dict[str, Any]) -> None:
    print("\n=== VCF Statistics ===")
    print(f"\nTotal variants: {stats['total_variants']}")
    
    print("\nVariant Types:")
    for vtype, count in stats['variant_types'].items():
        print(f"  {vtype}: {count}")
    
    print("\nChromosome Distribution:")
    for chrom, count in sorted(stats['chromosomes'].items()):
        print(f"  {chrom}: {count}")
    
    print("\nQuality Statistics:")
    print(f"  Minimum: {stats['quality_stats'].get('min', 'N/A')}")
    print(f"  Maximum: {stats['quality_stats'].get('max', 'N/A')}")
    print(f"  Average: {stats['quality_stats'].get('average', 'N/A'):.2f}")
    
    print("\nFilter Statistics:")
    for filt, count in stats['filter_stats'].items():
        print(f"  {filt}: {count}")
    
    print("\nTransition/Transversion Statistics:")
    ti_tv = stats['transition_transversion']
    print(f"  Transitions: {ti_tv['transitions']}")
    print(f"  Transversions: {ti_tv['transversions']}")
    print(f"  Ti/Tv ratio: {ti_tv.get('ratio', 'N/A'):.2f}")

stats = calculate_vcf_stats(vcf_file)
print_stats(stats) 