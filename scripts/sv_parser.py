#!/usr/bin/env python3
"""sv_parser.py
Parse VCF of structural variants and output BED-like TSV with CHROM, START, END, SVTYPE.
Assumes INFO contains END and SVTYPE fields.
"""
import sys, csv, gzip

def parse_vcf(vcf_path, out_path):
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    with open_func(vcf_path, 'rt') as vcf, open(out_path, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['CHROM', 'START', 'END', 'SVTYPE'])
        for line in vcf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]
            info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
            end = int(info_dict.get('END', pos))
            svtype = info_dict.get('SVTYPE', 'NA')
            writer.writerow([chrom, pos, end, svtype])

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: sv_parser.py <input.vcf.gz> <output.tsv>')
        sys.exit(1)
    parse_vcf(sys.argv[1], sys.argv[2])
