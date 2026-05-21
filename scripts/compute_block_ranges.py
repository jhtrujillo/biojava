#!/usr/bin/env python3
"""compute_block_ranges.py
Generate a BED‑like table with the genomic span of every collinear block.
Output (data/block_ranges.tsv):
Block_ID\tChr1\tStart1\tEnd1\tChr2\tStart2\tEnd2\tOrientation
"""
import os, pandas as pd

REPORT = os.path.abspath('genomica_comparativa/r570/reporte_comparativo.tsv')
OUT_TSV = os.path.abspath('data/block_ranges.tsv')
os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
        'Gene2_ID','Chr2','Start2','End2','Strand2']
rep = pd.read_csv(REPORT, sep='\t', header=0, usecols=range(len(cols)), names=cols)

# Helper to decide orientation (majority of strand concordance)
def orientation(group):
    same = (group['Strand1'] == group['Strand2']).sum()
    return 'direct' if same >= len(group)/2 else 'inverted'

records = []
for block_id, sub in rep.groupby('Block_ID'):
    chr1 = sub['Chr1'].mode()[0]
    chr2 = sub['Chr2'].mode()[0]
    start1 = sub['Start1'].min()
    end1   = sub['End1'].max()
    start2 = sub['Start2'].min()
    end2   = sub['End2'].max()
    orient = orientation(sub)
    records.append([block_id, chr1, start1, end1, chr2, start2, end2, orient])

out_df = pd.DataFrame(records,
    columns=['Block_ID','Chr1','Start1','End1','Chr2','Start2','End2','Orientation'])
out_df.to_csv(OUT_TSV, sep='\t', index=False)
print(f'Wrote {len(out_df)} block ranges to {OUT_TSV}')
