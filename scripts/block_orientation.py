#!/usr/bin/env python3
"""block_orientation.py
Compute block orientation (direct vs inverted) from reporte_comparativo.tsv.
Output: data/block_orientation.tsv with columns:
Block_ID\tOrientation\tNumGenes\tChr1\tChr2
Orientation = "direct" if majority of rows have Strand1 == Strand2, else "inverted".
"""
import os, pandas as pd

REPORT = os.path.abspath('genomica_comparativa/r570/reporte_comparativo.tsv')
OUT_DIR = os.path.abspath('data')
os.makedirs(OUT_DIR, exist_ok=True)
OUT_TSV = os.path.join(OUT_DIR, 'block_orientation.tsv')

# Load report (tab‑separated, first line is header)
cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1','Gene2_ID','Chr2','Start2','End2','Strand2']
df = pd.read_csv(REPORT, sep='\t', usecols=range(len(cols)), names=cols, header=0)

# Group by block
summary = []
for block, sub in df.groupby('Block_ID'):
    # majority orientation
    same = (sub['Strand1'] == sub['Strand2']).sum()
    orientation = 'direct' if same >= len(sub)/2 else 'inverted'
    # representative chromosomes (most common)
    chr1 = sub['Chr1'].mode()[0]
    chr2 = sub['Chr2'].mode()[0]
    summary.append([block, orientation, len(sub), chr1, chr2])

out_df = pd.DataFrame(summary, columns=['Block_ID','Orientation','NumGenes','Chr1','Chr2'])
out_df.to_csv(OUT_TSV, sep='\t', index=False)
print(f"Block orientation written to {OUT_TSV}")
