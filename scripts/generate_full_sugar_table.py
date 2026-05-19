#!/usr/bin/env python3
import os
import csv
import re
import urllib.parse

gff1_path = "benchmarks/genomas/1940/CC-01-1940.gff3"
gff2_path = "benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3"
report_path = "genomica_comparativa/r570/reporte_comparativo.tsv"
kaks_path = "genomica_comparativa/r570/kaks_1940_vs_r570.tsv"
informe_path = "informe/metodos_y_resultados_fase1.md"

def load_gff_notes(gff_path):
    notes = {}
    if not os.path.exists(gff_path):
        return notes
    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type not in ('gene', 'mRNA'):
                continue
            
            attributes = parts[8]
            id_match = re.search(r'ID=([^;]+)', attributes)
            if not id_match:
                continue
            gene_id = id_match.group(1).strip()
            
            note_match = re.search(r'Note=([^;]+)', attributes)
            if note_match:
                note_val = note_match.group(1).strip()
                notes[gene_id] = urllib.parse.unquote(note_val)
    return notes

def norm_g1(x):
    x = x.strip().replace("gene:", "").replace("transcript:", "").replace("mrna:", "")
    x = re.sub(r"CC(\d+)t(\d+)", r"CC\1g\2", x)
    x = x.split(".")[0]
    return x

def norm_g2(x):
    x = x.strip().replace("gene:", "").replace("transcript:", "").replace("mrna:", "")
    match = re.match(r"(SoffiXsponR570\.[^.]+g\d+)", x)
    if match:
        return match.group(1)
    return x.split(".")[0]

print("Cargando datos...")
notes1 = load_gff_notes(gff1_path)
notes2 = load_gff_notes(gff2_path)

kaks_data = {}
with open(kaks_path, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    for row in reader:
        if not row or len(row) < 5:
            continue
        g1_norm = norm_g1(row[0])
        g2_norm = norm_g2(row[1])
        try:
            ka = float(row[2])
            ks = float(row[3])
            ratio = float(row[4])
            kaks_data[(g1_norm, g2_norm)] = (ka, ks, ratio)
        except ValueError:
            continue

keywords = ['sucrose', 'sps', 'sut', 'susy', 'invertase', 'sweet', 'sugar transporter', 'fructosyltransferase', 'galactosyltransferase']
sugar_genes = []

block_total = {}
block_sugar = {}
rows_data = []

with open(report_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        status = row['Status']
        if status != 'Syntenic':
            continue
        g1_id = row['Gene1_ID']
        g2_id = row['Gene2_ID']
        
        desc1 = notes1.get(g1_id, notes1.get(g1_id.split('.')[0], ""))
        desc2 = notes2.get(g2_id, notes2.get(g2_id.split('.')[0], ""))
        
        is_sugar = False
        for kw in keywords:
            if kw in desc1.lower() or kw in desc2.lower() or kw in g1_id.lower() or kw in g2_id.lower():
                is_sugar = True
                break
        
        block_id = row['Block_ID']
        block_total[block_id] = block_total.get(block_id, 0) + 1
        if is_sugar:
            block_sugar[block_id] = block_sugar.get(block_id, 0) + 1
        
        rows_data.append((row, desc1, desc2, is_sugar))

for row, desc1, desc2, is_sugar in rows_data:
    if is_sugar:
        g1_id = row['Gene1_ID']
        g2_id = row['Gene2_ID']
        g1_norm = norm_g1(g1_id)
        g2_norm = norm_g2(g2_id)
        
        ka, ks, ratio = kaks_data.get((g1_norm, g2_norm), (0.0, 0.0, 0.0))
        block_id = row['Block_ID']
        tot = block_total.get(block_id, 0)
        sug = block_sugar.get(block_id, 0)
        
        desc = desc1 or desc2 or "Sugar-related gene"
        desc = desc.replace("\t", " ").replace("|", "\\|")
        
        if ratio > 1.0:
            pressure = "**Selección Positiva**: Divergencia adaptativa acelerada."
        elif ratio < 0.1:
            pressure = "**Selección Purificadora Extrema**: Restricción evolutiva absoluta."
        else:
            pressure = "**Selección Purificadora**: Conservación funcional estándar."
            
        sugar_genes.append({
            'G1': g1_id,
            'Chr1': row['Chr1'],
            'G2': g2_id,
            'Chr2': row['Chr2'],
            'Block': f"**{block_id}** ({tot}/{sug})",
            'Desc': desc,
            'Ka': f"{ka:.4f}",
            'Ks': f"{ks:.4f}",
            'Ratio': f"**{ratio:.4f}**" if ratio > 1.0 or ratio < 0.1 else f"{ratio:.4f}",
            'Pressure': pressure
        })

# Sort: positive selection first, then by G1 chromosome, then by G1 ID
def sort_key(x):
    r = float(x['Ratio'].replace('**', ''))
    # prioritize positive, then purifying extreme, then moderate
    if r > 1.0:
        cat = 0
    elif r < 0.1:
        cat = 2
    else:
        cat = 1
    # extract chr number
    chr_match = re.search(r'\d+', x['Chr1'])
    chr_num = int(chr_match.group()) if chr_match else 99
    return (cat, chr_num, x['G1'])

sugar_genes.sort(key=sort_key)

print(f"Total sugar genes in table: {len(sugar_genes)}")

# Generate Markdown table
md_rows = []
for g in sugar_genes:
    row_str = f"| **`{g['G1']}`** | `{g['Chr1']}` | `{g['G2']}` | `{g['Chr2']}` | {g['Block']} | {g['Desc']} | {g['Ka']} | {g['Ks']} | {g['Ratio']} | {g['Pressure']} |"
    md_rows.append(row_str)

table_header = """| Gen CC-01-1940 | Chr CC (G1) | Gen R570 | Chr R570 (G2) | Bloque (Total/Azúcar) | Anotación Funcional | Ka | Ks | Ka/Ks | Presión Selectiva e Implicación Fisiológica |
| :--- | :---: | :--- | :---: | :---: | :--- | :---: | :---: | :---: | :--- |"""

full_table_md = table_header + "\n" + "\n".join(md_rows)

# Load existing report and replace table
with open(informe_path, 'r', encoding='utf-8') as f:
    content = f.read()

# Find and replace the table
# We want to replace from "| Gen CC-01-1940 |" down to the next "#### Discusión"
pattern = r"\| Gen CC-01-1940 \|.*?\n\n#### Discusión"
new_section = full_table_md + "\n\n#### Discusión"

content, count = re.subn(pattern, new_section, content, flags=re.DOTALL)
print(f"Substitutions made: {count}")

with open(informe_path, 'w', encoding='utf-8') as f:
    f.write(content)
print("Report updated successfully!")
