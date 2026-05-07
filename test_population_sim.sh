#!/bin/bash
[ignoring loop detection]

mkdir -p results/population_test

# Definición de Secuencias Base (Longitud 100)
# Ref: 'A' x 49 + 'C' + 'A' x 50
ref_seq="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
# SNP All (G at 120, T at 150, G at 180)
snp_all="AAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

# Calidad Phred (100 'I's)
qual="IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"

echo "Generando muestra 1 (Homocigoto Referencia)..."
cat <<EOF > results/population_test/Sample1.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
for i in {1..10}; do
  echo -e "read${i}_s1\t0\tchr1\t101\t60\t100M\t*\t0\t0\t${ref_seq}\t${qual}" >> results/population_test/Sample1.sam
done

echo "Generando muestra 2 (Mutaciones Variadas - Heterocigoto)..."
cat <<EOF > results/population_test/Sample2.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
# 5 lecturas de referencia
for i in {1..5}; do
  echo -e "read${i}_ref_s2\t0\tchr1\t101\t60\t100M\t*\t0\t0\t${ref_seq}\t${qual}" >> results/population_test/Sample2.sam
done
# 3 lecturas forward con mutaciones
for i in {6..8}; do
  echo -e "read${i}_alt_fwd_s2\t0\tchr1\t101\t60\t100M\t*\t0\t0\t${snp_all}\t${qual}" >> results/population_test/Sample2.sam
done
# 2 lecturas reverse con mutaciones
for i in {9..10}; do
  echo -e "read${i}_alt_rev_s2\t16\tchr1\t101\t60\t100M\t*\t0\t0\t${snp_all}\t${qual}" >> results/population_test/Sample2.sam
done

echo "Generando muestra 3 (Homocigoto Alternativo)..."
cat <<EOF > results/population_test/Sample3.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
for i in {1..5}; do
  echo -e "read${i}_fwd_s3\t0\tchr1\t101\t60\t100M\t*\t0\t0\t${snp_all}\t${qual}" >> results/population_test/Sample3.sam
done
for i in {6..10}; do
  echo -e "read${i}_rev_s3\t16\tchr1\t101\t60\t100M\t*\t0\t0\t${snp_all}\t${qual}" >> results/population_test/Sample3.sam
done

echo "Ejecutando BioJava Population Variant Caller..."
java -jar target/biojava.jar call-variants \
  -i results/population_test \
  -r results/simulated_ref.fasta \
  -o results/final_population.vcf \
  -p 4 \
  --preset freebayes \
  -t 2

echo -e "\n=== VCF RESULTANTE ==="
cat results/final_population.vcf | grep -v "##"

echo -e "\n\n🚀 TEST DE ESTIMACIÓN AUTOMÁTICA DE PLOIDÍA (--auto-ploidy)..."
java -jar target/biojava.jar call-variants \
  -i results/population_test \
  -r results/simulated_ref.fasta \
  -o results/final_population_autoploidy.vcf \
  --auto-ploidy \
  --preset freebayes \
  -t 2

echo -e "\n\n🚀 TEST DE FILTRADO POBLACIONAL (Call Rate > 80%): Esperamos 0 SNPs..."
java -jar target/biojava.jar call-variants \
  -i results/population_test \
  -r results/simulated_ref.fasta \
  -o results/final_population_filtered.vcf \
  --min-call-rate 0.8 \
  --preset freebayes \
  -t 2

