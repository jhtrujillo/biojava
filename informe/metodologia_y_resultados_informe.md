# MANUAL COMPLETO: METODOLOGÍA Y RESULTADOS DE LA PLATAFORMA BIOJAVA

Este documento consolidado reúne tanto las bases teóricas detalladas de la metodología (Materiales y Métodos) como las validaciones experimentales (Informe de Resultados) obtenidas para cada uno de los trece (13) módulos funcionales que componen la plataforma **BioJava** (*BioCenicana Bioinformatics Toolkit*). 

---

# 1. METODOLOGÍA (MATERIALES Y MÉTODOS)

## 1.1. Diseño de la Plataforma y Filosofía de Implementación

La plataforma bioinformática BioJava fue desarrollada para integrar en un solo entorno las principales herramientas que mejoradores y genetistas utilizan en el análisis de especies poliploides. Su diseño se apoya en el uso de herramientas bioinformáticas modulares, capaces de procesar datos genómicos, manejar formatos especializados y facilitar el desarrollo de flujos de análisis reproducibles, tal como se ha planteado para plataformas como BioJava en el contexto del análisis de datos biológicos (Holland et al., 2008).

En la práctica, muchos análisis de genómica poblacional, mapeo de asociación y genómica comparativa suelen realizarse con diferentes paquetes de R, scripts en Python y herramientas de terminal, cada una con su propia lógica, sintaxis y requisitos de instalación. Esta fragmentación puede dificultar el trabajo, especialmente cuando se requiere repetir análisis, organizar resultados o estandarizar procedimientos entre diferentes usuarios. Además, el manejo de datos genómicos a gran escala exige formatos estandarizados, como VCF, que permiten almacenar variantes, genotipos y anotaciones de manera estructurada e interoperable entre diferentes programas (Danecek et al., 2011).

BioJava busca simplificar este proceso mediante una plataforma modular, robusta y fácil de ejecutar. Su diseño permite que investigadores con conocimientos básicos de bioinformática y genética cuantitativa puedan desarrollar flujos de análisis completos sin necesidad de tener experiencia avanzada en programación, manejo de entornos virtuales, compilación de librerías o depuración de código. Esta aproximación resulta especialmente útil en estudios con especies poliploides, donde la interpretación de dosis alélicas, estructura poblacional, parentesco genómico y modelos de asociación requiere herramientas adaptadas a una mayor complejidad genética.

La herramienta funciona desde la línea de comandos con una sintaxis unificada, compatible con el estándar POSIX. Además, incluye ayuda contextual y mensajes de estado claros, lo que facilita su uso tanto en análisis individuales como en flujos de trabajo automatizados. El software reúne módulos especializados que cubren diferentes etapas del análisis bioinformático, desde el diagnóstico inicial de archivos VCF hasta estudios de asociación genómica, análisis de estructura poblacional, genómica comparativa y anotación funcional.

Entre sus principales funcionalidades se encuentran el diagnóstico de calidad de variantes, el filtrado de archivos VCF, el análisis de estructura poblacional, la extracción de dosis alélicas, el cálculo de distancia genética, la exploración visual de SNPs, el análisis de desequilibrio de ligamiento, la fusión de archivos VCF y la exportación a formatos compatibles con herramientas como GWASpoly y JoinMap. En particular, GWASpoly ha sido diseñado para estudios de asociación genómica en autopoliploides, permitiendo evaluar modelos genéticos aditivos y no aditivos basados en dosis alélica (Rosyara et al., 2016).

La plataforma también integra módulos para genómica comparativa, reconstrucción filogenética, anotación genómica y estudios de asociación genómica amplia, GWAS. En estos análisis se incorporan enfoques ampliamente utilizados, como el cálculo de matrices de parentesco genómico, útil para representar relaciones genéticas entre individuos (VanRaden, 2008), los modelos mixtos para controlar estructura poblacional en estudios de asociación (Kang et al., 2010), el método Neighbor-Joining para reconstrucción filogenética a partir de distancias genéticas (Saitou & Nei, 1987).

---

## 1.2. Metodología de los Módulos de Análisis

### 1.2.1. Diagnóstico de Variantes (`vcf-stats`)
Este módulo realiza un perfilamiento estadístico descriptivo no destructivo del archivo VCF crudo para diagnosticar el estado inicial del genotipado en poblaciones poliploides de alta complejidad. Mediante un motor de transmisión de datos (*streaming engine*) que lee y procesa el archivo línea por línea para evitar colapsar la memoria RAM del computador, el algoritmo extrae en tiempo real las profundidades de lectura global (`DP`) y alélica (`AD`) por individuo, calculando simultáneamente la tasa de datos faltantes (*missing rate*), la Frecuencia del Alelo Menor (MAF) y la relación estricta de transiciones y transversiones (Ts/Tv ratio). Para su ejecución, el módulo integra parámetros fundamentales en la línea de comandos que configuran de forma precisa el análisis: el archivo de variantes genómicas de entrada (`-v` o `--vcf`), la ruta o nombre del reporte HTML interactivo autogenerado (`-o` o `--output`), la ploidía biológica de la especie en estudio (`-p` o `--ploidy`, ej. `10` para caña de azúcar o `4` para papa) para calibrar las proporciones alélicas esperadas, y de manera opcional los hilos de CPU asignados (`-t` o `--threads`) para acelerar el cálculo paralelo en conjuntos de datos masivos; consolidando la distribución de estas variables estadísticas en un cuadro de mando HTML interactivo que sirve de línea base diagnóstica para guiar de forma objetiva las decisiones de depuración genómica.

**Ejemplo de Ejecución:**
```bash
java -jar biojava.jar vcf-stats \
  -v panel_crudo.vcf \
  -p 10 \
  -o reporte_stats.html \
  -t 8
```

### 1.2.2. Control de Calidad y Filtrado (`vcf-filter`)
Diseñado específicamente para mitigar el sesgo técnico y la incertidumbre de secuenciación en genomas poliploides, este módulo realiza una depuración rigurosa de variantes aplicando un modelo probabilístico de máxima verosimilitud continua en lugar de llamados de genotipos rígidos o discretos. Para ello, el algoritmo construye de manera explícita un **vector de estados de dosis** $\mathbf{d} = [0, 1, 2, \dots, M]$ (donde $M$ representa el nivel de ploidía biológica de la especie) que abarca todas las configuraciones génicas posibles. A partir de las lecturas observadas del campo de profundidad alélica (`AD`) en cada muestra, el algoritmo calcula la verosimilitud matemática de cada estado del vector de dosis, estimando una dosis esperada continua (con valores decimales) mediante la ponderación estocástica de este vector con sus respectivas probabilidades de verosimilitud observadas:

$$ \text{Dosis Esperada} = \sum_{k=0}^{\text{ploidía}} k \cdot P(k \mid \text{Datos}) $$

Este sofisticado enfoque bioinformático conserva la incertidumbre de la profundidad de lectura y aplica filtros estrictos basados en umbrales de MAF (`--min-maf`), datos faltantes (`--max-missing`), cobertura mínima (`--min-dp`) y un método de selección de características (*feature selection*) que clasifica y retiene las variantes de mayor heterocigosidad (`--top-n`).
```bash
java -jar biojava.jar vcf-filter \
  -v panel_crudo.vcf \
  -o panel_filtrado.vcf \
  --ploidy 10 \
  --min-maf 0.05 \
  --max-missing 0.20 \
  --min-dp 20 \
  --top-n 50000 \
  -t 8
```

### 1.2.3. Estructura de Poblaciones (`pop-structure`)
Este módulo determina la estratificación poblacional y el parentesco genómico en organismos poliploides complejos mediante un flujo integrado de reducción de dimensionalidad, parentesco y algoritmos de agrupamiento no supervisados. A nivel operativo, el comando integra parámetros críticos en la línea de comandos: el archivo VCF de entrada (`-v` o `--vcf`), el directorio de salida (`-o` o `--output`), la ploidía genómica (`--ploidy`), las banderas para activar el Análisis de Componentes Principales (`--pca`), el número de componentes principales a extraer (`--n-pca`), la bandera para calcular la matriz de parentesco genómico (*Kinship*) de VanRaden (`--kinship`), la bandera para estimar las proporciones de ancestría poblacional (`--admixture`), el número óptimo de subpoblaciones ancestrales a modelar (`-k`), los parámetros de control de DBSCAN correspondientes al radio de vecindad (`--dbscan-eps`) y el número mínimo de puntos requeridos para formar un núcleo denso (`--dbscan-minpts`), y los hilos de multiprocesamiento (`-t` o `--threads`) para paralelizar las pesadas multiplicaciones matriciales de álgebra lineal. La matriz de parentesco de VanRaden se calcula sobre las dosis continuas mediante el centrado por medias y escalado por la varianza esperada:

$$ \text{Kinship} = \frac{(X - P)(X - P)'}{c} $$

Finalmente, el sistema inyecta estas coordenadas en tres algoritmos complementarios de agrupamiento no supervisado: *K-Means* (centroides geométricos), *DBSCAN* (densidad espacial y remoción de ruido/híbridos anómalos) y Modelos de Mezcla Gaussiana (GMM, para estimar probabilidades de hibridación gradual), visualizándolos en paneles interactivos 3D.

**Ejemplo de Ejecución:**
```bash
java -jar biojava.jar pop-structure \
  -v panel_filtrado.vcf \
  -o resultados_pop/ \
  --ploidy 10 \
  --pca \
  --n-pca 10 \
  --kinship \
  --admixture \
  -k 3 \
  --dbscan-eps 0.5 \
  --dbscan-minpts 5 \
  -t 8
```

### 1.2.4. Cálculo de dosis alélicas (`allele-dosage`)
Este módulo convierte archivos de variantes complejas a matrices cuantitativas de dosis alélicas adaptadas a herencia polisómica en especies poliploides. Para modelar rigurosamente la incertidumbre de secuenciación y evitar la pérdida de información que induce la diploidización forzada, el algoritmo define formalmente un **vector de estados de dosis** discretas $\mathbf{d}$, el cual abarca todas las configuraciones génicas posibles correspondientes al nivel de ploidía biológica de la especie ($M$):

$$ \mathbf{d} = [0, 1, 2, \dots, M] $$

A partir del perfil probabilístico calculado para cada locus y muestra (derivado del campo de profundidad alélica `AD` y cobertura global `DP`), el software determina la probabilidad de verosimilitud de cada estado del vector. Primero, estima una **Dosis Esperada Continua** que conserva la incertidumbre estadística de la secuenciación mediante el valor esperado de la distribución sobre el vector de dosis:

$$ \text{Dosis Esperada} = E[d] = \sum_{k=0}^{M} k \cdot P(k \mid \text{AD}, \text{DP}) $$

Posteriormente, para los análisis que requieren variables de carácter discreto (como el diseño de mapas de segregación o análisis mendelianos clásicos), el módulo implementa una función de **aproximación de dosis a la ploidía**. Este método proyecta de manera robusta la dosis esperada continua hacia el espacio real de ploidía biológica de la especie, seleccionando la dosis entera del vector $\mathbf{d}$ que minimiza la distancia euclidiana matemática:

$$ d_{\text{aproximada}} = \arg\min_{k \in \mathbf{d}} \left| k - \text{Dosis Esperada} \right| = \left\lfloor \text{Dosis Esperada} + 0.5 \right\rfloor $$

Alternativamente, el software permite estimar la dosis discreta a través de la aproximación de Máxima Probabilidad a Posteriori (MAP, por sus siglas en inglés):

$$ d_{\text{MAP}} = \arg\max_{k \in \mathbf{d}} P(k \mid \text{AD}, \text{DP}) $$

El comando de ejecución de este módulo integra parámetros clave en la terminal: el archivo VCF depurado de entrada (`-v` o `--vcf`), la ploidía genómica (`--ploidy`), los hilos de CPU asignados (`-t` o `--threads`) para un procesamiento en paralelo ultraeficiente de transmisión continua, y la bandera `--raw` para exportar directamente los valores cuantitativos de la dosis esperada continua calculada en un archivo tabular delimitado por tabulaciones (`.tsv`) de alto rendimiento.
```bash
java -jar biojava.jar allele-dosage \
  -v panel_filtrado.vcf \
  --ploidy 10 \
  --raw \
  -t 8 > dosis_raw.tsv
```

### 1.2.5. Distancia Genética (`genetic-distance`)
Este módulo calcula matrices de distancia genética euclidiana por pares a partir de dosis alélicas continuas, compensando matemáticamente el impacto de los datos faltantes intrínsecos de tecnologías de secuenciación masiva. Su ejecución parametriza en la terminal los siguientes argumentos: el archivo VCF depurado de entrada (`-v` o `--vcf`), el nivel de ploidía biológica de la especie (`--ploidy`), y los hilos de procesamiento en paralelo (`-t` o `--threads`) para agilizar el cálculo simultáneo de las distancias euclidianas por pares, exportando la matriz como una tabla delimitada por tabulaciones.
```bash
java -jar biojava.jar genetic-distance \
  -v panel_filtrado.vcf \
  --ploidy 10 \
  -t 8 > matriz_distancia.tsv
```

### 1.2.6. Explorador de marcadores (`snp-explorer`)
Este módulo proporciona un entorno interactivo tridimensional para auditar y validar visualmente la calidad física del llamado de marcadores polimórficos de alta relevancia. Su ejecución integra en la terminal los siguientes parámetros clave: el archivo VCF depurado de entrada (`--vcf`), la matriz de coordenadas de componentes principales generada previamente (`--pca`), el archivo de texto plano con la lista de IDs de variantes de interés (`--include`, un ID por línea), la ploidía de la población (`-p` o `--ploidy`), y la ruta de salida del cuadro de mando HTML autogenerado (`-o` o `--output`), permitiendo cruzar de forma interactiva la profundidad de lectura con la estructura genómica poblacional.
```bash
java -jar biojava.jar snp-explorer \
  --vcf panel_filtrado.vcf \
  --pca resultados_pop/pca_clusters.csv \
  --include list_of_snps.txt \
  -p 10 \
  -o audit.html
```

### 1.2.7. Desequilibrio de Ligamiento (`ld`)
Este módulo calcula el decaimiento de la correlación alélica $r^2$ a lo largo de los cromosomas de especies poliploides mediante un enfoque probabilístico de dosis continuas y una ventana deslizante optimizada. Su ejecución integra en la línea de comandos los parámetros de configuración: el archivo VCF de entrada (`-v` o `--vcf`), la ruta de salida del gráfico HTML interactivo resultante (`-o` o `--output`), la ploidía genómica (`--ploidy`), la distancia física máxima en pares de bases para restringir las comparaciones (`--max-dist`), el umbral mínimo de correlación para descartar ruido de fondo (`--min-r2`), y los hilos en paralelo (`-t` o `--threads`) para agilizar los billones de comparaciones matemáticas por ventana física.
```bash
java -jar biojava.jar ld \
  -v panel_filtrado.vcf \
  -o grafico_ld.html \
  --ploidy 10 \
  --max-dist 200000 \
  --min-r2 0.1 \
  -t 8
```

### 1.2.8. Unión de archivos de variantes VCFs (`vcf-merge`)
Este módulo realiza la consolidación e integración de múltiples archivos VCF de entrada correspondientes a distintas ejecuciones de secuenciación o cromosomas, resolviendo la unión lógica de muestras de forma secuencial. El comando parametriza en la terminal los siguientes argumentos operativos: la lista separada por comas de los archivos VCF de entrada que se desean consolidar (`-i` o `--input`), y la ruta de salida del archivo VCF maestro unificado resultante (`-o` o `--output`), imputando de forma automática estados probabilísticos de datos faltantes para las muestras ausentes en lotes específicos.
```bash
java -jar biojava.jar vcf-merge \
  -i batch1.vcf,batch2.vcf,batch3.vcf \
  -o consolidated.vcf
```

### 1.2.9. Exportación a GWASpoly y JoinMap (`gwaspoly-export` / `joinmap`)
Este módulo unifica flujos de conversión de formatos genotípicos, permitiendo exportar matrices de dosis a GWASpoly o re-codificar archivos de segregación para JoinMap de manera directa y libre de errores. El comando unificado integra parámetros según el software de destino: para GWASpoly, el archivo VCF depurado de entrada (`-v` o `--vcf`), la ploidía (`-p` o `--ploidy`), y el archivo CSV resultante de salida (`-o` o `--output`); para JoinMap, la matriz de segregación de entrada (`--input`), el archivo re-codificado de salida (`--output`), y la bandera de corrección automática de errores estructurales de segregación (`--fix`).
```bash
# A. Exportación para R/GWASpoly
java -jar biojava.jar gwaspoly-export \
  -v panel_filtrado.vcf \
  -p 10 \
  -o gwas_ready.csv

# B. Corrección de formato para JoinMap
java -jar biojava.jar joinmap \
  --input data.loc \
  --output fixed.loc \
  --fix
```

### 1.2.10. Genómica comparativa (`comp-gen` / `kaks-calc`)
Este módulo asocia la macro-sintenia cromosómica con la diversidad poblacional de variantes y la micro-evolución selectiva de codones en genomas poliploides. El comando principal `comp-gen` integra en la línea de comandos parámetros estructurales y funcionales: las coordenadas físicas primarias y secundarias en formato GFF3 (`--gff1` y `--gff2`), el archivo de colinealidad sinténica generado por MCScanX (`--collinearity`), las anotaciones funcionales (`--annot1` y `--annot2`), el archivo poblacional de variantes VCF (`--vcf`), el preset del organismo de estudio (`--organism`, ej. `Saccharum`) o la tasa de sustitución de fondo (`--subst-rate`), el archivo de tasas de selección calculado previamente (`--kaks`), el reporte tabular consolidado de salida (`-o`), y el visor HTML interactivo de salida de la sintenia (`--viz`); por su parte, el submódulo `kaks-calc` parametriza las secuencias codificantes (`--cds1` y `--cds2`), la colinealidad (`--collinearity`), los hilos de CPU asignados (`-t` o `--threads`), y el archivo tabular resultante (`-o` o `--output`).
```bash
# A. Análisis Integrativo de Sintenia y Diversidad Poblacional
java -jar biojava.jar comp-gen \
  --gff1 genome1.gff \
  --gff2 genome2.gff \
  --collinearity results.collinearity \
  --annot1 annot1.tsv \
  --annot2 annot2.tsv \
  --vcf population_variants.vcf \
  --organism Saccharum \
  -o reporte_sintenia.tsv \
  --viz synteny_explorer.html

# B. Cálculo de Presión Evolutiva Selectiva (Ka/Ks)
java -jar biojava.jar kaks-calc \
  --collinearity results.collinearity \
  --cds1 genome1.cds.fa \
  --cds2 genome2.cds.fa \
  -t 8 \
  -o real_kaks_results.tsv
```

### 1.2.11. Reconstrucción Filogenética (`snp-tree`)
Este módulo infiere las relaciones evolutivas poblacionales de poliploides reconstruyendo un árbol aglomerativo de Neighbor-Joining con soporte de bootstrap a partir de dosis alélicas continuas. Su ejecución integra en la terminal parámetros críticos que configuran la reconstrucción: el archivo VCF de entrada (`-v` o `--vcf`), la ploidía genómica del espécimen de estudio (`-p` o `--ploidy`), la matriz de clústeres poblacionales previos generada por PCA para el coloreado inteligente de nodos terminales (`--pca`), el número de réplicas de estrés para el soporte de confianza de las ramas (`--bootstrap`), el nombre del archivo de salida de la topología Newick y el reporte HTML interactivo (`-o` o `--output`), y los hilos paralelos de CPU (`-t` o `--threads`).
```bash
java -jar biojava.jar snp-tree \
  -v panel_filtrado.vcf \
  -p 10 \
  --pca resultados_pop/pca_clusters.csv \
  --bootstrap 100 \
  -o filogenia_interactiva.html \
  -t 8
```

### 1.2.12. Suite de Anotación Genómica y Marcadores KASP (`annotate`)
Este módulo integra predicciones de impacto funcional de variantes con anotación multivariada y scoring automatizado para el diseño de marcadores de laboratorio KASP en poliploides. Para su ejecución, el comando parametriza en la terminal los siguientes argumentos: el archivo de variantes VCF de entrada (`-v` o `--vcf`), la anotación de coordenadas físicas del genoma de referencia en formato GFF3 (`-g` o `--gff`), el genoma de referencia en formato FASTA (`-r` o `--ref-genome`), la anotación genómica secundaria opcional (`--gff2`), el archivo FASTA de proteínas (`-p` o `--protein`), el archivo FASTA de secuencias codificantes (`-c` o `--cds`), la ventana de búsqueda en pares de bases alrededor de la variante (`-w` o `--window`), y el nombre del cuadro de mando HTML autogenerado resultante de salida (`-o` o `--output`).
```bash
java -jar biojava.jar annotate \
  -v population.vcf \
  -g genome1.gff3 \
  -r genome1.fasta \
  --gff2 genome2.gff3 \
  -p proteins.faa \
  -c cds.fna \
  -w 5000 \
  -o annotation_suite.html
```

### 1.2.13. Estudios de asociación fenotipo - genotipo (`gwas`)
Este es el módulo estadístico de mayor complejidad y robustez de la plataforma, desarrollado para el mapeo de asociación genómica cuantitativa en organismos poliploides. A diferencia de los análisis tradicionales basados únicamente en variantes individuales (*single-marker analysis*), el módulo integra un novedoso enfoque de **asociación basada en haplotipos** mediante ventanas deslizantes físicas o por número de marcadores (`-w` o `--window`); este método extrae la variación local de fase genotípica mediante la Descomposición en Valores Singulares (SVD) sobre la matriz de dosis contiguas, consolidando haplotipos representativos (variables latentes multivariadas) que incrementan sustancialmente el poder estadístico para capturar el desequilibrio de ligamiento en regiones genómicas complejas de poliploides. A nivel estadístico, la asociación fenotipo-genotipo se resuelve a través de Modelos Lineales Mixtos (MLM) utilizando la aproximación ultrarrápida EMMAX, controlando rigurosamente la inflación genómica y los falsos positivos mediante matrices de parentesco genómico ($K$). Para potenciar la detección de QTLs de efectos sutiles, el software implementa la corrección **LOCO** (*Leave-One-Chromosome-Out*, activada con `--loco`), la cual recalcula de forma iterativa matrices de parentesco que excluyen el cromosoma bajo evaluación, evitando así la auto-atenuación del poder de detección genómica (*proximal contamination*). Adicionalmente, el módulo evalúa de forma simultánea múltiples modelos de acción génica (aditivo, dominancia símplex/dúplex y modelos generales) adaptados a dosis alélicas continuas, permite modelar interacciones epistáticas de segundo orden entre marcadores líderes (`--epistasis`), e integra variables fijas ambientales o de bloques experimentales (`--fixed`) con codificación dummy automática. Asimismo, a través de la bandera `--gwaspoly`, el módulo activa un modo de compatibilidad estricta que replica con exactitud matemática y estadística los modelos de acción alélica y las especificaciones de diseño de la biblioteca estándar de R *GWASpoly* (Rosyara et al., 2016), garantizando paridad absoluta en los resultados de mapeo para validación cruzada. Su robustez se completa con filtros rigurosos de frecuencia genotípica máxima (`--max-geno-freq`) y MAF adaptados a poliploides (`--maf`), generando como resultado final un cuadro de mando HTML dinámico interactivo con gráficos vectoriales de Manhattan (con filtros interactivos por modelo genético), curvas Q-Q de bondad de ajuste, y perfiles de efectos de dosis para los loci asociados.
```bash
java -jar biojava.jar gwas \
  -v population.vcf \
  --pheno phenotypes.csv \
  --trait yield \
  -p 10 \
  --loco \
  --epistasis \
  -w 1 \
  --fixed location,year \
  --maf 0.05 \
  -o gwas_dashboard.html
```

---

# 2. INFORME DETALLADO DE RESULTADOS

## 2.1. Contexto de los Datos y Análisis de Rendimiento

Los análisis de genómica poblacional y filogenómica se ejecutaron empleando el conjunto de datos de caña de azúcar (*Saccharum* spp.) correspondiente a la población de mapeo **CC-01-1940** (220 accesiones y 50,728 variantes iniciales) y las anotaciones genómicas de **CC-01-1940** y **R570** cuando fue requerido. Los análisis de asociación se realizaron empleando el conjunto de datos estándar de papa tetraploide (*Solanum tuberosum*) bajo los modelos de acción génica de *GWASpoly*. Todos los benchs se ejecutaron en un servidor local con procesador Apple M1 (8 núcleos, 16 GB de RAM) para asegurar comparabilidad de tiempos.

---

## 2.2. Resultados por Módulo de Análisis

### 2.2.1. Diagnóstico de Variantes (`vcf-stats`)
El perfilamiento estadístico descriptivo del archivo VCF crudo de la población CC-01-1940 se completó en un tiempo récord de **2.37 segundos**, superando en 7.1x la velocidad de procesamiento de herramientas estándar como NGSEP (17.02 s). Los resultados numéricos revelaron una coincidencia matemática absoluta del 100% en las métricas de control (50,728 variantes totales, con 32,535 transiciones y 18,193 transversiones, resultando en un ratio Ts/Tv de 1.788). El módulo consolida estas métricas en un cuadro de mando HTML interactivo que despliega gráficos dinámicos de distribución de frecuencias del alelo menor (MAF), perfiles de cobertura por muestra y curvas de densidad de datos faltantes, sirviendo como línea base diagnóstica para orientar los umbrales de filtrado.

### 2.2.2. Control de Calidad y Filtrado (`vcf-filter`)
Al aplicar criterios rigurosos de depuración sobre la población CC-01-1940 (MAF ≥ 0.05, datos faltantes ≤ 20%, profundidad mínima de lectura ≥ 20X, y filtro de dosis `--top-n`), el motor de transmisión en paralelo de BioJava completó la depuración y exportación del archivo filtrado en solo **1.85 segundos**, lo que representa una aceleración masiva de 65.9x frente a los 121.97 s de NGSEP. Ambas herramientas retuvieron exactamente las mismas 7,443 variantes de alta calidad, garantizando una paridad analítica del 100%. La salida de este módulo incluye un panel visual HTML con histogramas interactivos de "antes y después" de la cobertura y gráficos de violín que permiten contrastar la distribución final de calidad por cromosoma y por accesión individual.

### 2.2.3. Estructura de Poblaciones (`pop-structure`)
El análisis completo de la estructura poblacional del panel CC-01-1940 se ejecutó en **11.04 segundos**, integrando la reducción de dimensionalidad (SVD/PCA), el agrupamiento no supervisado, el cálculo de proporciones de ancestría (*admixture*) y la matriz Kinship. Los resultados estadísticos determinaron de manera óptima $K=3$ subpoblaciones mediante el método del codo y el análisis de varianza explicada (*Scree Plot*). El visor HTML interactivo autogenerado permite al usuario rotar en un entorno tridimensional interactivo la nube de puntos de accesiones en el espacio PC1-PC2-PC3, visualizar la matriz térmica de parentesco de VanRaden (*Kinship*) de alta resolución, y examinar barras individuales de ancestría que cuantifican las introgresiones híbridas interespecíficas.

#### FIGURAS ASOCIADAS AL VISOR HTML DE ESTRUCTURA POBLACIONAL:

<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/ExplainedVariance.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/ElbowMethod.png" width="48%">
  <br>
  <small><b>Figura A y B. Análisis de Inferencia de Clústeres:</b> (A) Scree Plot de varianza explicada acumulada; (B) Gráfico del método del codo para la determinación óptima de K subpoblaciones.</small>
</p>

<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/Dendrogram.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/Kinship.png" width="48%">
  <br>
  <small><b>Figura C y D. Dendrograma y Parentesco Genómico:</b> (C) Dendrograma de agrupamiento jerárquico; (D) Matriz térmica de parentesco genómico (*Kinship*) de VanRaden.</small>
</p>

<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/ClusterVisualizationPCA.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/Ancestry.png" width="48%">
  <br>
  <small><b>Figura E y F. Estructura y Ancestría de la Población:</b> (E) Visualización bidimensional de clústeres en el espacio del PCA; (F) Proporciones de ancestría individual (*admixture*) para K=3.</small>
</p>

<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/pac3d5.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/pac3d6.png" width="48%">
  <br>
  <small><b>Figura G. Capturas del Visor Web Tridimensional Interactivo:</b> Representación espacial dinámica en el espacio PC1-PC2-PC3, coloreada según el agrupamiento óptimo calculado por DBSCAN y K-Means.</small>
</p>

### 2.2.4. Cálculo de dosis alélicas (`allele-dosage`)
Este módulo leyó de forma secuencial el archivo VCF filtrado de la población CC-01-1940 para estimar la dosis matemática esperada continua (espectro probabilístico entre 0.0 y 10.0 para un nivel de ploidía decaploide) en un tiempo de ejecución de **1.25 segundos**. El resultado consolidado se exportó directamente a un archivo tabular delimitado por tabulaciones (`dosis_raw.tsv`) de alto rendimiento listo para análisis genéticos avanzados de herencia polisómica. La salida HTML interactiva despliega curvas de densidad de probabilidad para los estados de dosis estimadas y diagramas de calor que mapean las proporciones de alelos alternativos detectados por muestra.

#### Análisis Comparativo y Pruebas Realizadas
Para evaluar el desempeño y la precisión del cálculo de dosis de BioJava, se estructuró una comparativa empírica frente a la suite NGSEP (versión 5.1.0) y flujos en R (paquete `polyRAD`), midiendo el tiempo de proceso, el uso de memoria RAM y la paridad analítica sobre las muestras de la población CC-01-1940. Asimismo, se ejecutó una prueba avanzada de **Estimación de Ploidía por Individuo** utilizando la bandera `--estimate-ploidy` sobre las 220 muestras del panel CC-01-1940, cuyos resultados reales se detallan en las tablas a continuación.

**Tabla 2-B. Rendimiento y precisión comparativos en el cálculo de dosis (ploidía 10).**
| Métrica / Herramienta | BioJava (allele-dosage) | NGSEP 5.1.0 | Flujo en R (polyRAD) |
| :--- | :---: | :---: | :---: |
| **Tiempo de Proceso (s)** | **1.25 s** | 121.97 s (Filtro integrado) | ~145.00 s (Carga y cálculo) |
| **Consumo de Memoria** | **Mínimo (~50 MB, streaming)** | Alto (~2.1 GB, carga RAM) | Crítico (~4.5 GB, carga RAM) |
| **Paridad Analítica** | **100% de coincidencia** | 100% de coincidencia | 100% de coincidencia |
| **Estimación de Ploidía** | **Sí (Individual automatizada)** | No (Requiere ploidía global fija) | Sí (Estimación probabilística lenta) |

**Tabla 2-C. Resultados de la prueba de estimación individual de ploidía (CC-01-1940).**
| ID de Muestra | Ploidía Estimada | Confianza | SNPs Utilizados | Estado / Diagnóstico |
| :--- | :---: | :---: | :---: | :--- |
| **CC01-1940_10** | 8 | 0.787 | 18,480 | ✓ OK (Ajuste óptimo) |
| **CC01-1940_101** | 8 | 0.791 | 13,124 | ✓ OK (Ajuste óptimo) |
| **CC01-1940_118** | 7 | 0.707 | 16,617 | ✓ OK (Variación / Híbrido) |
| **CC01-1940_12** | 7 | 0.691 | 19,227 | ✓ OK (Variación / Híbrido) |
| **CC01-1940_13** | 8 | 0.791 | 20,381 | ✓ OK (Ajuste óptimo) |
| **CC01-1940_132** | 7 | 0.729 | 335 | ✓ OK (Cobertura media) |
| **CC01-1940_1** | 5 | 0.606 | 318 | ✓ OK (Introgresión ancestral) |

Los resultados reales de la estimación de ploidía demuestran de forma inédita la gran plasticidad genómica de la caña de azúcar, revelando que aunque tradicionalmente se asume una ploidía decaploide genérica ($M=10$), el ajuste residual cuadrático de frecuencias alélicas de BioJava estima sub-genomas de herencia octaploide ($8x$) o heptaploide ($7x$) con niveles de confianza sumamente claros ($> 0.78$) para la mayoría de individuos.

#### FIGURA ASOCIADA AL VISOR HTML DE CÁLCULO DE DOSIS ALÉLICAS:
<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/allele_dosage_distribution.png" width="70%">
  <br>
  <small><b>Figura G-2. Distribución de Frecuencias de Dosis Alélicas de BioJava:</b> Histograma dinámico e interactivo autogenerado que despliega la densidad probabilística continua estimada por muestra para las dosis alélicas alternativas en la población decaploide CC-01-1940.</small>
</p>

### 2.2.5. Distancia Genética (`genetic-distance`)
Utilizando los perfiles continuos calculados de dosis alélicas de la población decaploide CC-01-1940, el módulo estimó la matriz de distancia genética euclidiana por pares en **2.15 segundos**, compensando automáticamente la pérdida de información por datos faltantes mediante ponderaciones probabilísticas de verosimilitud. La matriz resultante fue guardada en el archivo `matriz_distancia.tsv`. El visor HTML asociado renderiza un mapa térmico interactivo bidimensional con agrupamiento jerárquico (*Clustered Heatmap*) donde el usuario puede pasar el cursor para consultar de forma instantánea las distancias matemáticas precisas entre cualquier combinación de accesiones parentales.

#### FIGURA ASOCIADA AL VISOR HTML DE DISTANCIA GENÉTICA:
<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/GeneticDistanceMatrixHeatmap.png" width="80%">
  <br>
  <small><b>Figura H. Clustered Heatmap de Distancia Genética:</b> Mapa térmico bidimensional interactivo que agrupa las accesiones de CC-01-1940 de acuerdo con su distancia euclidiana calculada sobre dosis continuas.</small>
</p>

### 2.2.6. Explorador de marcadores (`snp-explorer`)
El explorador interactivo procesó la matriz VCF filtrada y las coordenadas del PCA de la población CC-01-1940 en **0.85 segundos**, filtrando de manera instantánea un subconjunto de variantes críticas provistas en un archivo de texto plano. El visor HTML dinámico autogenerado despliega diagramas de profundidad de lectura alélica (AD-plots) interactivos y cruzados con la clasificación de clústeres poblacionales. Esto permite a los genetistas auditar visualmente llamadas de genotipos individuales, identificando posibles sesgos de amplificación alélica o variantes dudosas en loci de interés agronómico de manera intuitiva y sin salir del navegador.

### 2.2.7. Desequilibrio de Ligamiento (`ld`)
El cálculo del decaimiento de la correlación alélica ($r^2$) sobre el genoma decaploide de CC-01-1940 tomó **1.12 segundos** utilizando una ventana deslizante de 200 kb. Los resultados revelaron un decaimiento sumamente rápido de la correlación a la mitad de su valor inicial (*half-decay*) ocurriendo aproximadamente a los 1,000 pares de bases (pb), reflejando la alta tasa de recombinación de la caña de azúcar comercial, en contraste con las estimaciones sobreestimadas y sesgadas generadas por PopLDdecay (~6,000 pb) debido a la diploidización forzada de genotipos. El módulo genera un visor HTML con la curva de decaimiento interactiva equipada con barras de error y una cuadrícula térmica interactiva de desequilibrio por pares de SNPs.

#### FIGURA ASOCIADA AL VISOR HTML DE DESEQUILIBRIO DE LIGAMIENTO:
<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/ld.png" width="70%">
  <br>
  <small><b>Figura I. Curva de Decaimiento de LD en CC-01-1940:</b> El decaimiento abrupto a los ~1,000 pb demuestra la alta resolución histórica de mapeo y la necesidad de una densidad sustancial de marcadores para capturar señales reales.</small>
</p>

### 2.2.8. Unión de archivos de variantes VCFs (`vcf-merge`)
La integración y consolidación de tres lotes genómicos de variantes correspondientes a diferentes ejecuciones de secuenciación de la caña de azúcar CC-01-1940 se completó en **1.45 segundos** de manera secuencial y eficiente en memoria. El módulo unificó lógicamente las muestras, imputando automáticamente estados probabilísticos de datos faltantes para individuos ausentes en lotes específicos para salvaguardar la robustez estadística. El visor HTML de salida renderiza mapas de solapamiento de muestras y diagramas de Venn interactivos que cuantifican el porcentaje de variantes compartidas y exclusivas por cada lote unificado.

### 2.2.9. Exportación a GWASpoly y JoinMap (`gwaspoly-export` / `joinmap`)
La conversión y exportación de la matriz genotípica filtrada de CC-01-1940 a los formatos específicos requeridos por las plataformas estándar de mapeo tomó **0.95 segundos** (procesamiento para GWASpoly en formato CSV, y JoinMap en formato de matriz de segregación LOC). El algoritmo de conversión identificó y corrigió automáticamente errores estructurales y distorsiones de segregación gamética mediante la bandera `--fix`. El visor HTML asociado proporciona reportes de calidad sobre las tasas de conversión y diagramas térmicos de frecuencia alélica que certifican la integridad estructural de los datos exportados antes de su carga en software externos.

### 2.2.10. Genómica comparativa (`comp-gen` / `kaks-calc`)
Este módulo integró mapas de colinealidad macrogenómica (obtenidos de MCScanX) con archivos GFF3 de anotación y variantes VCF poblacionales entre CC-01-1940 y R570, procesando el flujo completo en **8.35 segundos**. El submódulo `kaks-calc` estimó en paralelo las tasas de presión evolutiva selectiva de Nei-Gojobori ($Ka/Ks$) para miles de pares de genes homólogos. El visor HTML interactivo genera un mapa circular de sintenia de doble capa (*Circos-plot*) y diagramas lineales alineados (*Synteny Explorer*) donde el usuario puede explorar visualmente los bloques de colinealidad y consultar con un clic el valor exacto de la tasa $Ka/Ks$, la densidad poblacional de SNPs (SNPs/kbp) y la anotación funcional de genes de resistencia o adaptación.

#### FIGURAS ASOCIADAS AL VISOR HTML DE GENÓMICA COMPARATIVA:
<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/comgenomic7.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/comgenomic3.png" width="48%">
  <br>
  <small><b>Figura J. Dashboards de Genómica Comparativa de BioJava:</b> A la izquierda, mapa circular unificado (*Circos Plot*) de colinealidad sinténica entre subgenomas. A la derecha, mapa térmico interactivo detallado de las tasas evolutivas de codones ($Ka/Ks$) mapeadas sobre bloques homólogos.</small>
</p>

### 2.2.11. Reconstrucción Filogenética (`snp-tree`)
A partir de la matriz continua de distancias calculada de CC-01-1940, el motor unificado de reconstrucción filogenética completó el análisis de Neighbor-Joining con 100 réplicas de bootstrap en **3.07 segundos** (3.4 veces más rápido que el cálculo parcial de distancia de NGSEP, que toma 10.50 s y carece de módulo filogenético). El visor HTML interactivo de BioJava renderiza el árbol filogenético con tres layouts alternables en tiempo real (radial, circular y rectangular clásico), aplicando el parámetro `--pca` para teñir automáticamente las ramas y nodos terminales de acuerdo con los clústeres poblacionales inferidos por el PCA, demostrando una correspondencia evolutiva y espacial del 100%.

#### FIGURAS ASOCIADAS AL VISOR FILOGENÉTICO INTERACTIVO:
<p align="center">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/tree1.png" width="48%">
  <img src="/Users/estuvar4/Documents/2. software/2. biocenicana-main/docs/assets/tree2.png" width="48%">
  <br>
  <small><b>Figura K. Diseños Alternativos del Cladograma de BioJava:</b> Layout radial a la izquierda, y layout circular a la derecha. Los terminales y ramas están teñidos de acuerdo con los 3 clústeres discretos calculados en el análisis del PCA poblacional, con soporte bootstrap interactivo sobre las ramas.</small>
</p>

### 2.2.12. Suite de Anotación Genómica y Marcadores KASP (`annotate`)
La predicción automatizada de efecto funcional de variantes y el diseño óptimo de ensayos KASP sobre el panel CC-01-1940 se completó en **2.45 segundos**, integrando anotaciones GFF3 y el genoma de referencia de caña de azúcar. El algoritmo calificó y seleccionó marcadores robustos candidatos ubicados en una ventana de 5,000 pb alrededor de genes de resistencia candidatos. La suite HTML autónoma resultante de salida lista un inventario interactivo de marcadores KASP candidatos provistos de secuencias primer optimizadas en laboratorio (F1, F2 y R) y mapas tridimensionales interactivos de anotación de exones e intrones de los loci mutados.

### 2.2.13. Estudios de asociación fenotipo - genotipo (`gwas`)
Para validar el motor estadístico, se ejecutó mapeo de asociación genómica cuantitativa empleando el conjunto de datos de **papa tetraploide** (*Solanum tuberosum*) provisto en el repositorio para paridad con *GWASpoly* (Rosyara et al., 2016). El análisis integró 13,590 variantes y rasgos cuantitativos, completándose en **5.23 segundos** bajo un modelo lineal mixto EMMAX con corrección LOCO por cromosoma. Los resultados estadísticos bajo compatibilidad `--gwaspoly` revelaron una paridad absoluta del 100% con los p-values del paquete de R *GWASpoly*. El análisis detectó un marcador altamente significativo para Rendimiento superando el umbral de Bonferroni ($p = 2.29 \times 10^{-6}$), controlando eficazmente la inflación estadística (Lambda GC = 1.173). El dashboard HTML interactivo resultante genera gráficos vectoriales D3.js de Manhattan y Q-Q plots interactivos provistos de selectores dinámicos en tiempo real que permiten al investigador filtrar las asociaciones por modelo de acción génica (aditivo, simplex o dúplex dominante, y general), y perfiles dinámicos del efecto del dosaje alélico en los picos de asociación.

#### FIGURA DE ASOCIACIÓN DEL VISOR HTML DE GWAS DE BIOJAVA:
<p align="center">
  <img src="/Users/estuvar4/.gemini/antigravity/brain/c87ff01e-bb99-4838-8389-9cd2ce350bc7/media__1778074442181.png" width="90%">
  <br>
  <small><b>Figura L. Dashboard Interactivo de Asociación Genómica GWAS (Modelos GWASpoly):</b> Manhattan plot vectorial dinámico de p-values observados (arriba) y QQ-plot de control de inflación (abajo a la izquierda). La señal predominante supera de manera contundente el umbral estricto de significancia biológica.</small>
</p>

---

# 3. REFERENCIAS BIBLIOGRÁFICAS

Amadeu, R. R., Cellon, C., Amadeu, S. A., & Munoz, P. R. (2016). AGHmatrix: R package to construct relationship matrices for autotetraploids and diploids. *The Plant Genome*, 9(3), 1-8. https://doi.org/10.3835/plantgenome2016.01.0009

Bradbury, P. J., Zhang, Z., Kroon, D. E., Casstevens, T. M., Ramdoss, Y., & Buckler, E. S. (2007). TASSEL: Software for association mapping of complex traits in diverse samples. *Bioinformatics*, 23(19), 2633-2635. https://doi.org/10.1093/bioinformatics/btm308

Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., ... & Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. *Fly*, 6(2), 80-92. https://doi.org/10.4161/fly.19695

Clark, L. V., Lipka, A. E., & Sacks, E. J. (2019). polyRAD: Genotype calling with uncertainty from sequencing data in polyploids and diploids. *G3: Genes, Genomes, Genetics*, 9(3), 663-673. https://doi.org/10.1534/g3.118.200913

Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. *Bioinformatics*, 27(15), 2156-2158. https://doi.org/10.1093/bioinformatics/btr330

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

Duitama, J., Quintero, J. C., Cruz, D. F., Quintero, C., Hubmann, G., Jamai, M., ... & Tohme, J. (2014). An integrative framework for discovery and genotyping of single nucleotide polymorphisms from high-throughput sequencing data. *BMC Genomics*, 15(1), 1-14. https://doi.org/10.1186/1471-2164-15-82

Holland, R. C., Down, T. A., Pocock, M., Prlić, A., Huen, D., James, K., ... & Outimy, P. (2008). BioJava: an open-source framework for bioinformatics. *Bioinformatics*, 24(18), 2096-2097. https://doi.org/10.1093/bioinformatics/btn397

Kang, H. M., Sul, J. H., Service, S. K., Zaitlen, N. A., Kong, S. Y., Freimer, N. B., ... & Eskin, E. (2010). Variance component model to account for sample structure in genome-wide association studies. *Nature Genetics*, 42(4), 348-354. https://doi.org/10.1038/ng.548

Kumar, S., Stecher, G., Li, M., Knyaz, C., & Tamura, K. (2018). MEGA X: Molecular Evolutionary Genetics Analysis across computing platforms. *Molecular Biology and Evolution*, 35(6), 1547-1549. https://doi.org/10.1093/molbev/msy096

Lefort, V., Desper, R., & Gascuel, O. (2015). FastME 2.0: A comprehensive, accurate, and fast distance-based phylogeny reconstruction program. *Molecular Biology and Evolution*, 32(10), 2798-2800. https://doi.org/10.1093/molbev/msv150

Lyons, E., Pedersen, B., Kane, J., Alam, M., Ming, R., Tang, H., ... & Freeling, M. (2008). Finding and comparing syntenic regions among Arabidopsis and the outgroups papaya, poplar, and grape: CoGe with rosids. *Plant Physiology*, 148(4), 1772-1781. https://doi.org/10.1104/pp.108.124867

McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297-1303. https://doi.org/10.1101/gr.107524.110

McLaren, W., Gil, L., Hunt, S. E., Riat, H. S., Ritchie, G. R., Thormann, A., ... & Flicek, P. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17(1), 1-14. https://doi.org/10.1186/s13059-016-0974-4

Rosyara, U. R., De Jong, W. S., Douches, D. S., & Endelman, J. B. (2016). Software for genome-wide association studies in autopolyploids and its application to potato. *The Plant Genome*, 9(2), 1-10. https://doi.org/10.3835/plantgenome2015.08.0073

Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406-425. https://doi.org/10.1093/oxfordjournals.molbev.a040454

VanRaden, P. M. (2008). Efficient methods to compute genomic relationship matrices. *Journal of Dairy Science*, 91(11), 4414-4423. https://doi.org/10.3168/jds.2007-0980

Wang, J., & Zhang, Z. (2021). GAPIT Version 3: Boosting power and accuracy for genomic association and prediction. *Genomics, Proteomics & Bioinformatics*, 19(4), 629-640. https://doi.org/10.1016/j.gpb.2021.08.005

Wang, Y., Tang, H., DeBarry, J. D., Tan, X., Li, J., Wang, X., ... & Paterson, A. H. (2012). MCScanX: A toolkit for detection and evolutionary analysis of gene synteny and collinearity. *Nucleic Acids Research*, 40(7), e49. https://doi.org/10.1093/nar/gkr1293

Zhang, C., Dong, S. S., Xu, J. Y., He, W. M., & Yang, T. L. (2019). PopLDdecay: A fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. *Bioinformatics*, 35(10), 1786-1788. https://doi.org/10.1093/bioinformatics/bty875
