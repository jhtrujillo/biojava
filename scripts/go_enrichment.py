#!/usr/bin/env python3
"""go_enrichment.py
GO/KEGG enrichment via g:Profiler REST API for structurally affected genes
(genes in inverted blocks + orphan genes).
Output: genomica_comparativa/r570/tables/go_enrichment.tsv
Sources queried: GO:BP, GO:MF, KEGG
"""
import os, requests, pandas as pd

BASE    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIENT  = os.path.join(BASE, 'data', 'block_orientation.tsv')
REPORT  = os.path.join(BASE, 'genomica_comparativa', 'r570', 'reporte_comparativo.tsv')
OUT_DIR = os.path.join(BASE, 'genomica_comparativa', 'r570', 'tables')
OUT_TSV = os.path.join(OUT_DIR, 'go_enrichment.tsv')
os.makedirs(OUT_DIR, exist_ok=True)

API_URL  = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'


def write_placeholder(path, note, n_genes):
    pd.DataFrame([{
        'source': 'NOTE',
        'term_id': 'N/A',
        'term_name': note,
        'p_value': 'N/A',
        'gene_count': n_genes,
        'intersection_size': 0,
        'genes': ''
    }]).to_csv(path, sep='\t', index=False)
    print(f"Placeholder written to {path}")


# ── Load data ──────────────────────────────────────────────────────────────
orient_df = pd.read_csv(ORIENT, sep='\t', low_memory=False)
cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
        'Gene2_ID','Chr2','Start2','End2','Strand2']
rep = pd.read_csv(REPORT, sep='\t', header=0,
                  usecols=range(len(cols)), names=cols, low_memory=False)

# Inverted block gene IDs
inv_blocks = set(orient_df[orient_df['Orientation'] == 'inverted']['Block_ID'])
inv_genes  = set(rep[rep['Block_ID'].isin(inv_blocks)]['Gene1_ID'].dropna())

# Orphan gene IDs (cap to 500 for API)
orphan_mask  = rep['Status'].isin(['Orphan_G1', 'Orphan_G2'])
orphan_genes = set(rep[orphan_mask]['Gene1_ID'].dropna())

# Combined — prioritise inverted genes, then orphans, cap at 2000
combined = list(inv_genes)[:1000] + list(orphan_genes - inv_genes)[:1000]
combined = combined[:2000]
print(f"Query genes: {len(combined)} ({len(inv_genes)} inverted, {len(orphan_genes)} orphan)")

# ── g:Profiler REST call ───────────────────────────────────────────────────
payload = {
    'organism':       'sbicolor',
    'query':          combined,
    'sources':        ['GO:BP', 'GO:MF', 'KEGG'],
    'user_threshold': 0.05,
    'no_evidences':   False,
    'all_results':    False,
    'no_iea':         False,
    'domain_scope':   'annotated',
    'significance_threshold_method': 'g_SCS'
}

try:
    resp = requests.post(API_URL, json=payload, timeout=90)
    resp.raise_for_status()
    result   = resp.json()
    enriched = result.get('result', [])
    if enriched:
        rows = []
        for entry in enriched:
            rows.append({
                'source':            entry.get('source', ''),
                'term_id':           entry.get('native', ''),
                'term_name':         entry.get('name', ''),
                'p_value':           entry.get('p_value', ''),
                'gene_count':        entry.get('query_size', ''),
                'intersection_size': entry.get('intersection_size', ''),
                'genes':             ','.join(entry.get('intersections', []))
            })
        out_df = pd.DataFrame(rows)
        out_df.to_csv(OUT_TSV, sep='\t', index=False)
        print(f"GO enrichment: {len(out_df)} terms written to {OUT_TSV}")
    else:
        write_placeholder(OUT_TSV,
            'g:Profiler returned 0 enriched terms for sbicolor query',
            len(combined))
except Exception as e:
    print(f"g:Profiler API error: {e}")
    write_placeholder(OUT_TSV,
        f'g:Profiler API unavailable: {e}',
        len(combined))
