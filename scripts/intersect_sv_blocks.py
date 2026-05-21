#!/usr/bin/env python3
"""intersect_sv_blocks.py
Simple intersection of SV regions with block orientation data.
For demonstration, it joins SV records with block orientation by chromosome
(and adds a placeholder for genes). Output written to
`genomica_comparativa/r570/tables/sv_block_overlap.tsv`.
"""
import os, csv, pandas as pd

def main():
    # Paths
    sv_path = os.path.abspath('genomica_comparativa/r570/tables/sv_regions.bed')
    block_path = os.path.abspath('data/block_orientation.tsv')  # output from block_orientation.py
    out_path = os.path.abspath('genomica_comparativa/r570/tables/sv_block_overlap.tsv')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # Load SVs
    sv_df = pd.read_csv(sv_path, sep='\t', header=0)
    # Load block orientation (contains Chr1, Chr2 fields)
    block_df = pd.read_csv(block_path, sep='\t', header=0)

    # Simplistic join on chromosome matching (Chr1 or Chr2)
    merged = []
    for _, sv in sv_df.iterrows():
        matches = block_df[(block_df['Chr1'] == sv['CHROM']) | (block_df['Chr2'] == sv['CHROM'])]
        for _, blk in matches.iterrows():
            merged.append({
                'Block_ID': blk['Block_ID'],
                'SVTYPE': sv['SVTYPE'],
                'Chrom': sv['CHROM'],
                'SV_Start': sv['START'],
                'SV_End': sv['END'],
                'Orientation': blk['Orientation'],
                'NumGenes': blk['NumGenes'],
                'OrphanFlag': 'NA'  # placeholder – real implementation would check PAV list
            })
    out_df = pd.DataFrame(merged)
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"SV‑block overlap written to {out_path}")

if __name__ == '__main__':
    main()
