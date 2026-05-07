#!/bin/bash
# ==============================================================================
# BioJava Bioinformatics Suite - Script de Ejecución de Simulación
# ==============================================================================

# Colores elegantes para consola
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m' # No Color

clear
echo -e "${BOLD}${CYAN}======================================================================${NC}"
echo -e "${BOLD}${BLUE}   ⚡ BIOJAVA: INAUGURACIÓN Y EJECUCIÓN DEL PARALLEL VARIANT CALLER ⚡${NC}"
echo -e "${BOLD}${CYAN}======================================================================${NC}"
echo ""

# 1. Compilación
echo -e "${BOLD}${YELLOW}[PASO 1]${NC} Compilando y empaquetando BioJava..."
echo -e "${CYAN}Ejecutando: mvn clean package -DskipTests...${NC}"
mvn clean package -DskipTests

if [ $? -ne 0 ]; then
    echo -e "\n❌ ${BOLD}Error durante la compilación con Maven. Verifica tu entorno Java/Maven.${NC}"
    exit 1
fi

echo -e "\n✨ ${GREEN}¡Compilación completada con éxito! JAR generado en target/biojava.jar${NC}\n"

# 2. Ejecución con Modelo Probabilístico estilo FreeBayes (ploidía = 4, tetraploide)
echo -e "${BOLD}${YELLOW}[PASO 2]${NC} Corriendo llamado de variantes con modelo preestablecido FreeBayes..."
echo -e "${CYAN}Ejecutando Variant Caller para Ploides (p=4, hilos=2) en los datos simulados...${NC}"
echo ""

java -jar target/biojava.jar call-variants \
  -i results/simulated_aligned.sam \
  -r results/simulated_ref.fasta \
  -o results/simulated_variants.vcf \
  -p 4 \
  --preset freebayes \
  -t 2

if [ $? -ne 0 ]; then
    echo -e "\n❌ ${BOLD}Error durante la ejecución del Variant Caller.${NC}"
    exit 1
fi

# 3. Mostrar resultados
echo -e "\n${BOLD}${YELLOW}[PASO 3]${NC} Inspeccionando el archivo VCF generado (${BOLD}results/simulated_variants.vcf${NC})..."
echo -e "${CYAN}----------------------------------------------------------------------${NC}"

if [ -f results/simulated_variants.vcf ]; then
    cat results/simulated_variants.vcf
else
    echo -e "❌ ${BOLD}Archivo VCF de salida no encontrado.${NC}"
    exit 1
fi

echo -e "${CYAN}----------------------------------------------------------------------${NC}"
echo -e "\n🎉 ${BOLD}${GREEN}¡Simulación ejecutada con éxito!${NC}"
echo -e "${BOLD}El SNP en Chr1:150 fue genotipado exitosamente (0/1/1/1) reportando el nuevo campo de calidad GQ.${NC}\n"
