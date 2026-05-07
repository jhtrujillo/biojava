#!/bin/bash

# Script para simular una población y probar el Variant Caller Conjunto
mkdir -p results/population_test

echo "Generando muestra 1 (Homocigoto Referencia - 0/0/0/0)..."
cat <<EOF > results/population_test/Sample1.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
for i in {1..10}; do
  echo -e "read${i}_s1\t0\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample1.sam
done

echo "Generando muestra 2 (Heterocigoto - 0/0/1/1)..."
cat <<EOF > results/population_test/Sample2.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
for i in {1..5}; do
  # Reads Reference (C)
  echo -e "read${i}_ref_s2\t0\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample2.sam
done
for i in {6..8}; do
  echo -e "read${i}_alt_fwd_s2\t0\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample2.sam
done
for i in {9..10}; do
  echo -e "read${i}_alt_rev_s2\t16\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample2.sam
done

echo "Generando muestra 3 (Homocigoto Alternativo - 1/1/1/1)..."
cat <<EOF > results/population_test/Sample3.sam
@SQ	SN:chr1	LN:280
@SQ	SN:chr2	LN:280
EOF
for i in {1..5}; do
  echo -e "read${i}_fwd_s3\t0\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample3.sam
done
for i in {6..10}; do
  echo -e "read${i}_rev_s3\t16\tchr1\t101\t60\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> results/population_test/Sample3.sam
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
