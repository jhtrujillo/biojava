#!/usr/bin/env python3
"""generate_sugar_list.py
Extract all gene IDs that contain sucrose‑related keywords from the two GFFs
and write them to data/sugar_gene_ids.txt (one ID per line)."""
import os, re

GFF1 = os.path.abspath('benchmarks/genomas/1940/CC-01-1940.gff3')
GFF2 = os.path.abspath('benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3')
OUT  = os.path.abspath('data/sugar_gene_ids.txt')
os.makedirs(os.path.dirname(OUT), exist_ok=True)

keywords = ['sucrose','sps','sut','susy','invertase','sweet',
            'sugar transporter','fructosyltransferase','galactosyltransferase']

def extract_ids(gff_path):
    ids = set()
    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ('gene','mRNA'):
                continue
            attrs = parts[8]
            id_match = re.search(r'ID=([^;]+)', attrs)
            if not id_match:
                continue
            gene_id = id_match.group(1).strip()
            note_match = re.search(r'Note=([^;]+)', attrs, flags=re.I)
            desc = note_match.group(1).lower() if note_match else ''
            if any(kw in desc for kw in keywords):
                ids.add(gene_id)
    return ids

sugar_ids = extract_ids(GFF1) | extract_ids(GFF2)
with open(OUT, 'w') as out:
    for gid in sorted(sugar_ids):
        out.write(gid + '\n')
print(f'Found {len(sugar_ids)} sucrose‑related gene IDs → {OUT}')
