#!/bin/bash
set -e

echo "================================================="
echo "Simulando Datos de Segregación de Población (VCF)"
echo "================================================="

mkdir -p results
cat <<EOF > results/simulated_mapping_population.vcf
##fileformat=VCFv4.2
##source=BioJavaGeneticSimulator_1.0
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Integer,Description="Dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Ind1	Ind2	Ind3	Ind4	Ind5	Ind6	Ind7	Ind8	Ind9	Ind10
chr1	100	M1_1	A	G	100	PASS	DP=20	GT:DS	0/1:1	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0
chr1	200	M1_2	C	T	100	PASS	DP=20	GT:DS	0/1:1	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0
chr1	300	M1_3	G	A	100	PASS	DP=20	GT:DS	0/1:1	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0
chr2	100	M2_1	T	C	100	PASS	DP=20	GT:DS	0/0:0	1/1:2	0/1:1	0/1:1	1/1:2	0/0:0	1/1:2	0/1:1	0/0:0	1/1:2
chr2	200	M2_2	A	T	100	PASS	DP=20	GT:DS	0/0:0	1/1:2	0/1:1	0/1:1	1/1:2	0/0:0	1/1:2	0/1:1	0/0:0	1/1:2
chr2	300	M2_3	C	G	100	PASS	DP=20	GT:DS	0/0:0	1/1:2	0/1:1	0/1:1	1/1:2	0/0:0	1/1:2	0/1:1	0/0:0	1/1:2
chr3	100	M3_1	G	C	100	PASS	DP=20	GT:DS	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2
chr3	200	M3_2	T	A	100	PASS	DP=20	GT:DS	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2
chr3	300	M3_3	A	G	100	PASS	DP=20	GT:DS	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2	0/0:0	0/1:1	1/1:2
EOF

echo "🎉 Archivo VCF de segregación simulado en: results/simulated_mapping_population.vcf"

echo -e "\nCompilando paquete con Maven..."
mvn clean package -DskipTests

echo -e "\nEjecutando Constructor de Mapas Genéticos (BioJava genetic-map)..."
java -jar target/biojava.jar genetic-map \
  -i results/simulated_mapping_population.vcf \
  -o results/final_genetic_map.map \
  --lod 2.5 \
  --max-r 0.35 \
  --mapping-function kosambi

echo -e "\n=== MAPA GENÉTICO RESULTANTE ==="
cat results/final_genetic_map.map
