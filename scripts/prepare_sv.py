#!/usr/bin/env python3
"""prepare_sv.py
Extract structural variant (SV) records from the gzipped VCF and write a TSV.
Output columns: Chrom, Start (1‑based), End, SVTYPE, Size, Info.
"""
import os
import sys
import gzip
import csv

VCF_PATH = os.path.abspath('benchmarks/vcfs/1940/cc-01-1940_flye_polishing_allhic_220_standarfiltered.vcf.gz')
OUT_DIR = os.path.abspath('data')
OUT_TSV = os.path.join(OUT_DIR, 'sv_calls.tsv')

os.makedirs(OUT_DIR, exist_ok=True)

def parse_info(info_str):
    info = {}
    for entry in info_str.split(';'):
        if '=' in entry:
            k, v = entry.split('=', 1)
            info[k] = v
        else:
            info[entry] = True
    return info

with gzip.open(VCF_PATH, 'rt') as vcf, open(OUT_TSV, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t')
    writer.writerow(['Chrom', 'Start', 'End', 'SVTYPE', 'Size', 'Info'])
    for line in vcf:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue
        chrom = parts[0]
        pos = int(parts[1])  # 1‑based already
        ref = parts[3]
        alt = parts[4]
        info_str = parts[7]
        info = parse_info(info_str)
        svtype = info.get('SVTYPE')
        if not svtype:
            continue
        end = int(info.get('END', pos + len(ref) - 1))
        size = end - pos + 1
        writer.writerow([chrom, pos, end, svtype, size, info_str])
print(f"SV TSV written to {OUT_TSV}")
