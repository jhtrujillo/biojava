#!/usr/bin/env python3
import os
import sys
import re

# Default file paths
CDS1_PATH = "benchmarks/genomas/1940/CC-01-1940.cds.fna"
CDS2_PATH = "benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.cds.fna"
GFF1_PATH = "benchmarks/genomas/1940/CC-01-1940.gff3"
GFF2_PATH = "benchmarks/genomas/r570/Saccharum_hybrid_cultivar_R570.gff3"

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def load_annotation(gff_path, gene_id):
    if not os.path.exists(gff_path):
        return "No GFF found"
    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            attributes = parts[8]
            if gene_id in attributes:
                note_match = re.search(r'Note=([^;]+)', attributes)
                if note_match:
                    return note_match.group(1).strip()
    return "Anotación no encontrada"

def find_sequence(fasta_path, gene_id):
    if not os.path.exists(fasta_path):
        return None
    seq_lines = []
    found = False
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                if gene_id in header:
                    found = True
                elif found:
                    break
            else:
                if found:
                    seq_lines.append(line.strip())
    return "".join(seq_lines) if found else None

def needleman_wunsch(seq1, seq2):
    match = 2
    mismatch = -1
    gap = -2
    n, m = len(seq1), len(seq2)
    
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i * gap
    for j in range(m + 1):
        dp[0][j] = j * gap
        
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )
            
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            a1.append(seq1[i-1])
            a2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            a1.append(seq1[i-1])
            a2.append('-')
            i -= 1
        else:
            a2.append(seq2[j-1])
            a1.append('-')
            j -= 1
            
    return "".join(reversed(a1)), "".join(reversed(a2))

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 compare_genes.py <ID_Gene_CC1940> <ID_Gene_R570>")
        print("Ejemplo: python3 compare_genes.py CC01t042770.1 SoffiXsponR570.01Ag380000.1")
        sys.exit(1)
        
    g1_id = sys.argv[1]
    g2_id = sys.argv[2]
    
    print(f"Buscando secuencias CDS para G1 ({g1_id}) y G2 ({g2_id})...")
    s1 = find_sequence(CDS1_PATH, g1_id)
    s2 = find_sequence(CDS2_PATH, g2_id)
    
    if not s1:
        print(f"Error: No se encontró la secuencia para {g1_id} en {CDS1_PATH}")
        sys.exit(1)
    if not s2:
        print(f"Error: No se encontró la secuencia para {g2_id} en {CDS2_PATH}")
        sys.exit(1)
        
    print(f"Secuencias encontradas.")
    print(f"  G1 Length: {len(s1)} bp")
    print(f"  G2 Length: {len(s2)} bp")
    
    note1 = load_annotation(GFF1_PATH, g1_id)
    note2 = load_annotation(GFF2_PATH, g2_id)
    print(f"Anotaciones:")
    print(f"  G1: {note1}")
    print(f"  G2: {note2}")
    
    print("\nAlineando secuencias (Needleman-Wunsch)...")
    al1, al2 = needleman_wunsch(s1, s2)
    print("Alineamiento finalizado.")
    
    mutations = []
    codon_idx = 1
    for idx in range(0, len(al1), 3):
        c1 = al1[idx:idx+3]
        c2 = al2[idx:idx+3]
        if len(c1) < 3 or len(c2) < 3:
            continue
        
        # Check for gap
        if '-' in c1 or '-' in c2:
            mutations.append({
                'codon': codon_idx,
                'type': 'Indel/Gap',
                'c1': c1, 'c2': c2,
                'aa1': '-', 'aa2': '-'
            })
        elif c1 != c2:
            aa1 = GENETIC_CODE.get(c1.upper(), 'X')
            aa2 = GENETIC_CODE.get(c2.upper(), 'X')
            mtype = 'Sinónima' if aa1 == aa2 else 'No Sinónima'
            mutations.append({
                'codon': codon_idx,
                'type': mtype,
                'c1': c1, 'c2': c2,
                'aa1': aa1, 'aa2': aa2
            })
        codon_idx += 1
        
    syn = [m for m in mutations if m['type'] == 'Sinónima']
    nonsyn = [m for m in mutations if m['type'] == 'No Sinónima']
    indels = [m for m in mutations if m['type'] == 'Indel/Gap']
    
    print("\n==================================================")
    print("           RESULTADO DE LA COMPARACIÓN            ")
    print("==================================================")
    print(f"Mutaciones Sinónimas (silenciosas): {len(syn)}")
    print(f"Mutaciones No Sinónimas (cambian AA): {len(nonsyn)}")
    print(f"Indels (Inserciones/Deleciones): {len(indels)}")
    print("==================================================")
    
    if len(nonsyn) > 0:
        print("\nDetalle de mutaciones no sinónimas (cambio de aminoácido):")
        for m in nonsyn:
            print(f"  Codón {m['codon']:03d}: CC-1940 ({m['c1']} -> {m['aa1']}) vs R570 ({m['c2']} -> {m['aa2']})")
    else:
        print("\nNo se encontraron cambios en la secuencia proteica (100% de identidad en aminoácidos).")

if __name__ == '__main__':
    main()
