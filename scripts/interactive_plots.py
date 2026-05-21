#!/usr/bin/env python3
"""interactive_plots.py
Generate 4 standalone Plotly HTML files for structural change analysis.

Outputs (genomica_comparativa/r570/plots/):
  1. blocks_orientation.html   – bar + pie charts of block orientation
  2. dotplot_synteny_sv.html   – dot-plot coloured by orientation / status
  3. gene_sv_snp_network.html  – Sankey: chromosomes → block type → sugar/non-sugar
  4. sv_gene_table.html        – filterable DataTable with all gene records
"""
import os, sys
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

BASE      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ORIENT    = os.path.join(BASE, 'data', 'block_orientation.tsv')
RANGES    = os.path.join(BASE, 'data', 'block_ranges.tsv')
REPORT    = os.path.join(BASE, 'genomica_comparativa', 'r570', 'reporte_comparativo.tsv')
SUGAR_IDS = os.path.join(BASE, 'data', 'sugar_gene_ids.txt')
SNP_OVL   = os.path.join(BASE, 'genomica_comparativa', 'r570', 'tables', 'sugar_snp_overlap.tsv')
GO_TSV    = os.path.join(BASE, 'genomica_comparativa', 'r570', 'tables', 'go_enrichment.tsv')
OUT_DIR   = os.path.join(BASE, 'genomica_comparativa', 'r570', 'plots')
os.makedirs(OUT_DIR, exist_ok=True)

DARK_BG  = '#0f172a'
CARD_BG  = '#1e293b'
TEXT     = '#e2e8f0'
MUTED    = '#94a3b8'
GREEN    = '#22c55e'
ORANGE   = '#f97316'
BLUE     = '#3b82f6'
PURPLE   = '#a855f7'
GOLD     = '#f59e0b'
RED      = '#ef4444'

# ── Load shared data ────────────────────────────────────────────────────────
orient_df = pd.read_csv(ORIENT, sep='\t', low_memory=False)
ranges_df = pd.read_csv(RANGES, sep='\t', low_memory=False)

cols = ['Block_ID','Status','Gene1_ID','Chr1','Start1','End1','Strand1',
        'Gene2_ID','Chr2','Start2','End2','Strand2']
rep = pd.read_csv(REPORT, sep='\t', header=0,
                  usecols=range(len(cols)), names=cols, low_memory=False)

with open(SUGAR_IDS) as f:
    sugar_ids = set(l.strip() for l in f if l.strip())

rep['IsSugarGene'] = rep['Gene1_ID'].isin(sugar_ids)

# Merge orientation into report
orient_map = dict(zip(orient_df['Block_ID'], orient_df['Orientation']))
rep['Orientation'] = rep['Block_ID'].map(orient_map).fillna('unknown')

snp_df = pd.read_csv(SNP_OVL, sep='\t') if os.path.exists(SNP_OVL) else pd.DataFrame()

# ─────────────────────────────────────────────────────────────────────────────
# 1. blocks_orientation.html
# ─────────────────────────────────────────────────────────────────────────────
def plot_orientation():
    # Counts per chr-pair
    orient_df['ChrPair'] = orient_df['Chr1'] + ' ↔ ' + orient_df['Chr2']
    pair_counts = (orient_df.groupby(['ChrPair','Orientation'])
                            .size().reset_index(name='Count'))

    top_pairs = (pair_counts.groupby('ChrPair')['Count'].sum()
                             .nlargest(20).index.tolist())
    pc_top = pair_counts[pair_counts['ChrPair'].isin(top_pairs)]

    direct_c  = orient_df[orient_df['Orientation']=='direct'].shape[0]
    inv_c     = orient_df[orient_df['Orientation']=='inverted'].shape[0]

    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{"type": "bar"}, {"type": "pie"}]],
        subplot_titles=('Orientación por Par Cromosómico (Top 20)',
                        'Distribución Global'),
        column_widths=[0.65, 0.35]
    )

    for orient, color in [('direct', GREEN), ('inverted', ORANGE)]:
        sub = pc_top[pc_top['Orientation'] == orient]
        fig.add_trace(go.Bar(
            name=orient.capitalize(),
            x=sub['Count'], y=sub['ChrPair'],
            orientation='h',
            marker_color=color,
            hovertemplate='<b>%{y}</b><br>%{x} bloques<extra></extra>'
        ), row=1, col=1)

    fig.add_trace(go.Pie(
        labels=['Directos', 'Invertidos'],
        values=[direct_c, inv_c],
        marker_colors=[GREEN, ORANGE],
        hole=0.45,
        textinfo='label+percent',
        hovertemplate='<b>%{label}</b><br>%{value} bloques (%{percent})<extra></extra>',
        textfont=dict(size=13, color=TEXT)
    ), row=1, col=2)

    fig.update_layout(
        title=dict(
            text='Orientación de Bloques Sinténicos — CC 1940 vs R570',
            font=dict(size=20, color=TEXT, family='Inter, sans-serif'),
            x=0.5
        ),
        paper_bgcolor=DARK_BG,
        plot_bgcolor=CARD_BG,
        font=dict(color=TEXT, family='Inter, sans-serif'),
        barmode='stack',
        legend=dict(bgcolor='rgba(30,41,59,0.8)', bordercolor=MUTED,
                    borderwidth=1, font=dict(color=TEXT)),
        margin=dict(l=20, r=20, t=100, b=20),
        annotations=[dict(
            text=f'Total: {direct_c+inv_c:,} bloques',
            x=0.82, y=0.5, xref='paper', yref='paper',
            showarrow=False, font=dict(size=15, color=GOLD)
        )]
    )
    fig.update_xaxes(gridcolor='#334155', zerolinecolor='#334155')
    fig.update_yaxes(gridcolor='#334155', tickfont=dict(size=9))

    out = os.path.join(OUT_DIR, 'blocks_orientation.html')
    fig.write_html(out, include_plotlyjs='cdn')
    print(f'✓ {out}')

# ─────────────────────────────────────────────────────────────────────────────
# 2. dotplot_synteny_sv.html
# ─────────────────────────────────────────────────────────────────────────────
def plot_dotplot():
    # Use block_ranges for coordinates; merge orientation
    ranges_df['Orientation'] = ranges_df['Block_ID'].map(orient_map).fillna('unknown')
    ranges_df['ChrPair'] = ranges_df['Chr1'] + ' ↔ ' + ranges_df['Chr2']

    # Top 12 pairs by gene count
    top_pairs = (ranges_df.groupby('ChrPair')['Block_ID'].count()
                          .nlargest(12).index.tolist())
    sub = ranges_df[ranges_df['ChrPair'].isin(top_pairs)].copy()

    color_map  = {'direct': GREEN, 'inverted': ORANGE, 'unknown': MUTED}
    symbol_map = {'direct': 'circle', 'inverted': 'diamond', 'unknown': 'x'}

    fig = go.Figure()

    for orient in ['direct', 'inverted']:
        s = sub[sub['Orientation'] == orient]
        mid1 = (s['Start1'] + s['End1']) / 2
        mid2 = (s['Start2'] + s['End2']) / 2
        size = (s['End1'] - s['Start1']).clip(lower=0)
        size_scaled = (size / size.max() * 18 + 4).fillna(4)

        fig.add_trace(go.Scatter(
            x=mid1, y=mid2,
            mode='markers',
            name=orient.capitalize(),
            marker=dict(
                color=color_map[orient],
                size=size_scaled,
                symbol=symbol_map[orient],
                opacity=0.7,
                line=dict(width=0.5, color='rgba(255,255,255,0.2)')
            ),
            customdata=s[['Block_ID','Chr1','Chr2','Orientation']].values,
            hovertemplate=(
                '<b>Bloque %{customdata[0]}</b><br>'
                'Chr1: %{customdata[1]}  Pos1: %{x:,.0f}<br>'
                'Chr2: %{customdata[2]}  Pos2: %{y:,.0f}<br>'
                'Orientación: <b>%{customdata[3]}</b><extra></extra>'
            )
        ))

    fig.update_layout(
        title=dict(
            text='Dot-Plot de Colinealidad — CC 1940 vs R570<br>'
                 '<sup>Top 12 pares cromosómicos · tamaño ∝ longitud del bloque</sup>',
            font=dict(size=18, color=TEXT, family='Inter, sans-serif'),
            x=0.5
        ),
        xaxis=dict(title='Posición genómica CC 1940 (bp)',
                   gridcolor='#1e293b', color=MUTED, showgrid=True),
        yaxis=dict(title='Posición genómica R570 (bp)',
                   gridcolor='#1e293b', color=MUTED, showgrid=True),
        paper_bgcolor=DARK_BG,
        plot_bgcolor='#0f1f35',
        font=dict(color=TEXT, family='Inter, sans-serif'),
        legend=dict(bgcolor='rgba(30,41,59,0.85)', bordercolor=MUTED,
                    borderwidth=1, font=dict(color=TEXT)),
        margin=dict(l=60, r=20, t=100, b=60)
    )

    out = os.path.join(OUT_DIR, 'dotplot_synteny_sv.html')
    fig.write_html(out, include_plotlyjs='cdn')
    print(f'✓ {out}')

# ─────────────────────────────────────────────────────────────────────────────
# 3. gene_sv_snp_network.html  (Sankey)
# ─────────────────────────────────────────────────────────────────────────────
def plot_sankey():
    # Aggregate flows: Chr1 → Status → Sugar/NonSugar
    # Build flow counts
    flow = (rep.groupby(['Chr1','Status','IsSugarGene'])
               .size().reset_index(name='count'))

    # Limit to top 15 chromosomes by total gene count
    top_chrs = (rep.groupby('Chr1').size().nlargest(15).index.tolist())
    flow = flow[flow['Chr1'].isin(top_chrs)]

    # Node lists
    chrs     = sorted(flow['Chr1'].unique().tolist())
    statuses = ['Syntenic', 'Orphan_G1', 'Orphan_G2', 'Direct', 'Inverted']
    status_set = sorted(flow['Status'].unique().tolist())
    leaves   = ['Sugar Gene', 'Non-Sugar Gene']

    # Map orientation into status for clarity
    all_nodes = chrs + status_set + leaves
    node_idx  = {n: i for i, n in enumerate(all_nodes)}

    node_colors = []
    for n in all_nodes:
        if n in chrs:          node_colors.append(BLUE)
        elif 'Orphan' in n:    node_colors.append(PURPLE)
        elif n == 'Syntenic':  node_colors.append(GREEN)
        elif n == 'Sugar Gene': node_colors.append(GOLD)
        else:                  node_colors.append(MUTED)

    sources, targets, values, colors = [], [], [], []

    for _, row in flow.iterrows():
        chr_node    = row['Chr1']
        status_node = row['Status']
        leaf_node   = 'Sugar Gene' if row['IsSugarGene'] else 'Non-Sugar Gene'
        cnt         = row['count']

        if chr_node not in node_idx or status_node not in node_idx:
            continue

        # Chr → Status
        sources.append(node_idx[chr_node])
        targets.append(node_idx[status_node])
        values.append(cnt)
        colors.append('rgba(59,130,246,0.4)')

        # Status → leaf
        sources.append(node_idx[status_node])
        targets.append(node_idx[leaf_node])
        values.append(cnt)
        link_color = 'rgba(245,158,11,0.5)' if row['IsSugarGene'] else 'rgba(148,163,184,0.2)'
        colors.append(link_color)

    fig = go.Figure(go.Sankey(
        node=dict(
            pad=15, thickness=20,
            line=dict(color='rgba(255,255,255,0.1)', width=0.5),
            label=all_nodes,
            color=node_colors,
            hovertemplate='<b>%{label}</b><br>%{value} genes<extra></extra>'
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=colors,
            hovertemplate='%{value} genes<extra></extra>'
        )
    ))

    fig.update_layout(
        title=dict(
            text='Flujo Genómico: Cromosomas → Tipo de Bloque → Genes de Azúcar',
            font=dict(size=18, color=TEXT, family='Inter, sans-serif'),
            x=0.5
        ),
        paper_bgcolor=DARK_BG,
        font=dict(color=TEXT, family='Inter, sans-serif', size=11),
        margin=dict(l=20, r=20, t=80, b=20)
    )

    out = os.path.join(OUT_DIR, 'gene_sv_snp_network.html')
    fig.write_html(out, include_plotlyjs='cdn')
    print(f'✓ {out}')

# ─────────────────────────────────────────────────────────────────────────────
# 4. sv_gene_table.html  (full interactive DataTable)
# ─────────────────────────────────────────────────────────────────────────────
def plot_table():
    # Build display dataframe
    disp = rep[['Block_ID','Status','Gene1_ID','Chr1','Start1','End1',
                'Gene2_ID','Chr2','Start2','End2','Orientation','IsSugarGene']].copy()
    disp = disp.fillna('—')
    disp['IsSugarGene'] = disp['IsSugarGene'].map({True: '🍬 Sí', False: 'No', '—': '—'})
    disp_sample = disp.head(10000)  # cap for browser performance

    rows_html = ''
    for _, r in disp_sample.iterrows():
        sugar_cls = 'sugar-row' if r['IsSugarGene'] == '🍬 Sí' else ''
        orient_badge = ''
        if r['Orientation'] == 'direct':
            orient_badge = '<span class="badge-direct">Direct</span>'
        elif r['Orientation'] == 'inverted':
            orient_badge = '<span class="badge-inv">Invertido</span>'
        else:
            orient_badge = f'<span class="badge-other">{r["Orientation"]}</span>'
        status_badge = f'<span class="badge-status">{r["Status"]}</span>'
        rows_html += f'''<tr class="{sugar_cls}">
          <td>{r["Block_ID"]}</td><td>{status_badge}</td>
          <td class="gene-id">{r["Gene1_ID"]}</td><td>{r["Chr1"]}</td>
          <td>{int(r["Start1"]):,}</td><td>{int(r["End1"]):,}</td>
          <td class="gene-id">{r["Gene2_ID"]}</td><td>{r["Chr2"]}</td>
          <td>{int(r["Start2"]):,}</td><td>{int(r["End2"]):,}</td>
          <td>{orient_badge}</td><td>{r["IsSugarGene"]}</td>
        </tr>'''

    html = f'''<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8">
  <title>Genes & Cambios Estructurales — CC 1940 vs R570</title>
  <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: 'Inter', sans-serif; background: #0f172a; color: #e2e8f0; }}
    .header {{ background: linear-gradient(135deg, #1e3a8a 0%, #7c3aed 100%);
               padding: 24px 32px; border-bottom: 1px solid #334155; }}
    .header h1 {{ font-size: 22px; font-weight: 800; letter-spacing: -0.02em; }}
    .header p  {{ font-size: 13px; color: rgba(255,255,255,0.7); margin-top: 4px; }}
    .controls {{ display: flex; gap: 12px; padding: 16px 32px; background: #1e293b;
                  border-bottom: 1px solid #334155; flex-wrap: wrap; align-items: center; }}
    .controls label {{ font-size: 11px; font-weight: 700; color: #94a3b8;
                        text-transform: uppercase; letter-spacing: 0.05em; }}
    .controls input, .controls select {{
        padding: 7px 12px; border-radius: 8px; border: 1px solid #334155;
        background: #0f172a; color: #e2e8f0; font-size: 12px; font-family: inherit; outline: none;
    }}
    .controls input:focus, .controls select:focus {{ border-color: #3b82f6; }}
    .count-badge {{ margin-left: auto; background: #1d4ed8; color: #fff;
                    padding: 5px 14px; border-radius: 20px; font-size: 12px; font-weight: 700; }}
    .table-wrap {{ overflow-x: auto; padding: 0 32px 32px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 12px; margin-top: 16px; }}
    thead tr {{ background: #1e293b; }}
    th {{ padding: 10px 14px; text-align: left; font-size: 10px; font-weight: 700;
           color: #94a3b8; text-transform: uppercase; letter-spacing: 0.05em;
           border-bottom: 2px solid #334155; white-space: nowrap; cursor: pointer; user-select: none; }}
    th:hover {{ color: #e2e8f0; }}
    td {{ padding: 8px 14px; border-bottom: 1px solid #1e293b; vertical-align: middle; }}
    tr:hover td {{ background: #1e293b; }}
    .sugar-row td {{ background: rgba(245,158,11,0.08); }}
    .sugar-row:hover td {{ background: rgba(245,158,11,0.15); }}
    .gene-id {{ font-family: monospace; font-size: 11px; color: #93c5fd; }}
    .badge-direct {{ background: #14532d; color: #4ade80; padding: 2px 8px;
                      border-radius: 12px; font-size: 10px; font-weight: 700; }}
    .badge-inv    {{ background: #7c2d12; color: #fb923c; padding: 2px 8px;
                      border-radius: 12px; font-size: 10px; font-weight: 700; }}
    .badge-other  {{ background: #1e293b; color: #94a3b8; padding: 2px 8px;
                      border-radius: 12px; font-size: 10px; font-weight: 700; }}
    .badge-status {{ background: #1e3a8a; color: #93c5fd; padding: 2px 8px;
                      border-radius: 12px; font-size: 10px; font-weight: 600; }}
    #no-results {{ text-align: center; padding: 40px; color: #64748b; font-size: 14px; display: none; }}
  </style>
</head>
<body>
<div class="header">
  <h1>🧬 Genes &amp; Cambios Estructurales — CC 1940 vs R570</h1>
  <p>Tabla interactiva con orientación de bloques, estado sinténico y anotación de genes de azúcar
     (mostrando hasta 10,000 registros)</p>
</div>
<div class="controls">
  <div style="display:flex;gap:6px;align-items:center;">
    <label>🔍 Buscar gen</label>
    <input type="text" id="gene-search" placeholder="Gene ID..." oninput="filterTable()" style="width:200px;">
  </div>
  <div style="display:flex;gap:6px;align-items:center;">
    <label>Estado</label>
    <select id="status-filter" onchange="filterTable()">
      <option value="">Todos</option>
      <option>Syntenic</option><option>Orphan_G1</option><option>Orphan_G2</option>
    </select>
  </div>
  <div style="display:flex;gap:6px;align-items:center;">
    <label>Orientación</label>
    <select id="orient-filter" onchange="filterTable()">
      <option value="">Todas</option>
      <option value="direct">Directa</option>
      <option value="inverted">Invertida</option>
    </select>
  </div>
  <div style="display:flex;gap:6px;align-items:center;">
    <label>Cromosoma</label>
    <select id="chr-filter" onchange="filterTable()">
      <option value="">Todos</option>
      {"".join(f'<option>{c}</option>' for c in sorted(disp_sample["Chr1"].unique()) if c != '—')}
    </select>
  </div>
  <div style="display:flex;gap:6px;align-items:center;">
    <label>Azúcar</label>
    <select id="sugar-filter" onchange="filterTable()">
      <option value="">Todos</option>
      <option value="Sí">Solo genes azúcar</option>
      <option value="No">No azúcar</option>
    </select>
  </div>
  <span class="count-badge" id="count-badge">Cargando...</span>
</div>
<div class="table-wrap">
  <table id="gene-table">
    <thead>
      <tr>
        <th onclick="sortTable(0)">Bloque ↕</th>
        <th onclick="sortTable(1)">Estado ↕</th>
        <th onclick="sortTable(2)">Gen G1 ↕</th>
        <th onclick="sortTable(3)">Chr G1 ↕</th>
        <th onclick="sortTable(4)">Inicio G1 ↕</th>
        <th onclick="sortTable(5)">Fin G1 ↕</th>
        <th onclick="sortTable(6)">Gen G2 ↕</th>
        <th onclick="sortTable(7)">Chr G2 ↕</th>
        <th onclick="sortTable(8)">Inicio G2 ↕</th>
        <th onclick="sortTable(9)">Fin G2 ↕</th>
        <th onclick="sortTable(10)">Orientación ↕</th>
        <th onclick="sortTable(11)">Azúcar ↕</th>
      </tr>
    </thead>
    <tbody id="table-body">{rows_html}</tbody>
  </table>
  <div id="no-results">No se encontraron resultados con los filtros aplicados.</div>
</div>
<script>
  const tbody = document.getElementById('table-body');
  const allRows = Array.from(tbody.querySelectorAll('tr'));
  let sortDir = {{}};

  function filterTable() {{
    const gene    = document.getElementById('gene-search').value.toLowerCase();
    const status  = document.getElementById('status-filter').value.toLowerCase();
    const orient  = document.getElementById('orient-filter').value.toLowerCase();
    const chr     = document.getElementById('chr-filter').value.toLowerCase();
    const sugar   = document.getElementById('sugar-filter').value.toLowerCase();
    let visible = 0;
    allRows.forEach(row => {{
      const cells = row.querySelectorAll('td');
      const txt   = row.textContent.toLowerCase();
      const show  =
        (!gene   || cells[2].textContent.toLowerCase().includes(gene)) &&
        (!status || cells[1].textContent.toLowerCase().includes(status)) &&
        (!orient || cells[10].textContent.toLowerCase().includes(orient)) &&
        (!chr    || cells[3].textContent.toLowerCase() === chr) &&
        (!sugar  || cells[11].textContent.toLowerCase().includes(sugar));
      row.style.display = show ? '' : 'none';
      if (show) visible++;
    }});
    document.getElementById('count-badge').textContent = visible.toLocaleString() + ' registros';
    document.getElementById('no-results').style.display = visible === 0 ? 'block' : 'none';
  }}

  function sortTable(col) {{
    const dir = sortDir[col] = !(sortDir[col]);
    const sorted = allRows.slice().sort((a, b) => {{
      const va = a.querySelectorAll('td')[col]?.textContent.trim() || '';
      const vb = b.querySelectorAll('td')[col]?.textContent.trim() || '';
      const na = parseFloat(va.replace(/,/g,'')), nb = parseFloat(vb.replace(/,/g,''));
      if (!isNaN(na) && !isNaN(nb)) return dir ? na - nb : nb - na;
      return dir ? va.localeCompare(vb) : vb.localeCompare(va);
    }});
    sorted.forEach(r => tbody.appendChild(r));
    filterTable();
  }}

  // Init count
  filterTable();
</script>
</body>
</html>'''

    out = os.path.join(OUT_DIR, 'sv_gene_table.html')
    with open(out, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'✓ {out}')

# ─────────────────────────────────────────────────────────────────────────────
# Run all
# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('Generating interactive plots...')
    plot_orientation()
    plot_dotplot()
    plot_sankey()
    plot_table()
    print('\nAll visualizations written to:', OUT_DIR)
