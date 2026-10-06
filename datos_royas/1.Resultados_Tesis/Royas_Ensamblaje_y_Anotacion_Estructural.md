# Ensamblaje y Anotación Estructural de los Genomas de *Puccinia melanocephala* y *Puccinia kuehnii*

*Versión validada contra los archivos de resultados — 6 de octubre de 2026*

---

## 1. Materiales y Métodos

### 1.1 Material biológico y secuenciación
Se trabajó con aislamientos monosóricos (derivados de una sola pústula), multiplicados por reinoculación en invernadero:

| Especie | Enfermedad | Cultivar de origen | Muestra | Esporas colectadas |
| :--- | :--- | :--- | :--- | :--- |
| *Puccinia melanocephala* | Roya Café | CC 85-92 | RCM1 | 331 mg |
| *Puccinia kuehnii* | Roya Naranja | CC 01-1940 | RNM2 | 680 mg |

El ADN se extrajo con el kit Omniprep, con ruptura mecánica de las uredosporas. La secuenciación se realizó con lecturas largas **PacBio HiFi** (librería Long Plex, fragmentos de 5–7 kb) y lecturas cortas **Illumina**. No fue posible obtener datos Hi-C, porque la cantidad de esporas disponible era inferior a la requerida.

### 1.2 Ensamblaje y depuración
1. **Control de calidad:** remoción de adaptadores con HiFiAdapterFilt.
2. **Tamaño del genoma:** estimación por k-meros (k = 31) con GenomeScope2.
3. **Ensamblaje:** hifiasm, que produce un ensamblaje primario y dos haplotipos parciales (hap1 y hap2).
4. **Contaminantes:** asignación taxonómica de los contigs con BlobToolKit (contenido GC, cobertura de lecturas y homología). Se conservaron únicamente los contigs asignados al género *Puccinia*.
5. **Completitud:** BUSCO con el linaje pucciniomycetes (n = 3,329).

El análisis posterior se realizó sobre el **haplotipo 2 depurado** de cada especie.

### 1.3 Anotación estructural
La predicción de genes se realizó con **BRAKER**, que integra los modelos de AUGUSTUS y GeneMark, usando evidencia proteica. Las estadísticas estructurales se calcularon a partir de los archivos GFF3; cuando un gen tiene varias isoformas, se usó como representante el transcrito con la CDS más larga.

> **Pendiente de confirmar:** la versión de BRAKER y si el genoma se enmascaró (*soft-masking*) antes de la predicción. Los FASTA depurados disponibles no están enmascarados.

---

## 2. Resultados

### 2.1 Tamaño del genoma estimado por k-meros

**Tabla 1. Perfil de k-meros (GenomeScope2, k = 31)**

| Parámetro | Roya Café (RCM1) | Roya Naranja (RNM2) |
| :--- | :--- | :--- |
| Tamaño haploide estimado | 240.2 Mb | 210.3 Mb |
| Secuencia única | 60.0% | 70.3% |
| Secuencia repetitiva | 40.0% | 29.7% |
| Cobertura de k-meros | 32.8× | 30× |
| Tasa de error | 1.07% | 1.51% |
| Tasa de duplicación | 1.42 | 1.35 |

La Roya Café presenta un genoma haploide estimado 14% mayor que el de la Roya Naranja y una fracción repetitiva más alta (40% frente a 30%).

### 2.2 Ensamblaje

**Tabla 2. Ensamblajes de hifiasm antes de la depuración**

| Ensamblaje | Métrica | Roya Café (RCM1) | Roya Naranja (RNM2) |
| :--- | :--- | :--- | :--- |
| **Primario** | Contigs | 6,903 | 5,969 |
| | Longitud total | 501 Mb | 359.7 Mb |
| | N50 | 1,157 kb | 473 kb |
| **Haplotipo 1** | Contigs | 8,904 | 6,988 |
| | Longitud total | 460 Mb | 347.1 Mb |
| | N50 | 424 kb | 234 kb |
| **Haplotipo 2** | Contigs | 4,373 | 2,926 |
| | Longitud total | 422 Mb | 305.3 Mb |
| | N50 | 556 kb | 274 kb |
| | Contig más largo | 5.92 Mb | 2.05 Mb |

El haplotipo 2 fue el más contiguo en ambas especies y se seleccionó para la depuración y la anotación.

### 2.3 Depuración de contaminantes

**Tabla 3. Haplotipo 2 antes y después de la depuración**

| Métrica | Roya Café (RCM1) | Roya Naranja (RNM2) |
| :--- | :--- | :--- |
| Haplotipo 2 sin depurar | 422 Mb; 4,373 contigs | 305.3 Mb; 2,926 contigs |
| **Ensamblaje depurado (solo *Puccinia*)** | **322.68 Mb; 1,455 contigs** | **267.54 Mb; 1,557 contigs** |
| N50 del ensamblaje depurado | 712.3 kb | 310.1 kb |
| Contig más largo | 5.92 Mb | 2.05 Mb |
| Contenido GC | 38.6% | 32.7% |
| Bases indeterminadas (N) | 0 | 0 |
| Secuencia removida | ~99 Mb (23.5%) | ~38 Mb (12.4%) |
| Contigs removidos | 2,918 | 1,369 |

**Tabla 4. Composición taxonómica del haplotipo 2 de Roya Naranja antes de la depuración (BlobToolKit)**

| Grupo taxonómico | Contigs | Longitud |
| :--- | :--- | :--- |
| Basidiomycota | 2,330 | 288 Mb |
| Actinomycetota | 381 | 9.06 Mb |
| Streptophyta (planta) | 69 | 4.62 Mb |
| Sin asignación | 64 | 0.87 Mb |
| Arthropoda | 17 | 0.61 Mb |
| Ascomycota | 18 | 0.57 Mb |
| Mucoromycota | 15 | 0.22 Mb |
| Chytridiomycota | 14 | 0.18 Mb |
| Uroviricota | 5 | 0.09 Mb |
| Otros | 15 | 0.97 Mb |

Dentro de Basidiomycota se conservó únicamente el grupo principal asignado a *Puccinia* (GC cercano a 32%), lo que explica la diferencia entre los 288 Mb de Basidiomycota y los 267.54 Mb finales.

> **Pendiente:** la tabla equivalente para Roya Café. La diapositiva de contaminación de RCM1 muestra la gráfica de RNM2, por lo que no se dispone de la composición previa a la depuración.

### 2.4 Completitud de los ensamblajes

**Tabla 5. BUSCO (linaje pucciniomycetes, n = 3,329)**

| Categoría | Roya Café, hap1 | Roya Café, hap2 | Roya Naranja, hap2 |
| :--- | :--- | :--- | :--- |
| **Completos** | 93.84% (3,124) | 93.42% (3,110) | 93.09% (3,099) |
| ↳ Copia única | 87.29% (2,906) | 87.62% (2,917) | 90.63% (3,017) |
| ↳ Duplicados | 6.55% (218) | 5.80% (193) | 2.46% (82) |
| Fragmentados | 1.32% (44) | 1.44% (48) | 1.53% (51) |
| Ausentes | 4.84% (161) | 5.14% (171) | 5.38% (179) |

Los ensamblajes de ambas especies recuperan cerca del 93% de los ortólogos conservados, por lo que su espacio génico es comparable. La Roya Café presenta 2.4 veces más ortólogos duplicados que la Roya Naranja (5.80% frente a 2.46%).

### 2.5 Relación entre el ensamblaje y el tamaño estimado

**Tabla 6. Ensamblaje depurado frente al tamaño haploide estimado**

| Métrica | Roya Café | Roya Naranja |
| :--- | :--- | :--- |
| Tamaño haploide estimado (k-meros) | 240.2 Mb | 210.3 Mb |
| Ensamblaje depurado | 322.68 Mb | 267.54 Mb |
| Relación ensamblaje / estimado | 1.34 | 1.27 |

Ambos ensamblajes superan el tamaño haploide estimado. Esto es esperable en genomas dicarióticos ensamblados sin datos Hi-C: parte del segundo haplotipo puede permanecer sin colapsar, y las repeticiones de alta copia tienden a subestimarse en el análisis de k-meros.

### 2.6 Anotación estructural

**Tabla 7. Estadísticas de la anotación estructural (BRAKER)**

| Métrica | Roya Café (*P. melanocephala*) | Roya Naranja (*P. kuehnii*) |
| :--- | :--- | :--- |
| **Genes predichos (loci)** | **32,412** | **15,690** |
| **Transcritos** | 33,999 | 16,432 |
| Genes con isoformas alternativas | 1,477 (4.6%) | 683 (4.4%) |
| Densidad génica | 100.4 genes/Mb | 58.6 genes/Mb |
| Longitud de gen, media (mediana) | 1,229 (857) pb | 2,288 (1,078) pb |
| Longitud de CDS, media (mediana) | 970 (645) pb | 1,019 (609) pb |
| Longitud de proteína, media (mediana) | 322 (214) aa | 339 (202) aa |
| Exones por gen, media (mediana) | 3.33 (2) | 4.01 (3) |
| Genes de un solo exón | 8,796 (27.1%) | 3,341 (21.3%) |
| Longitud de exón, media (mediana) | 292 (183) pb | 254 (150) pb |
| Longitud de intrón, media (mediana) | 111 (84) pb | 391 (107) pb |
| Fracción codificante del genoma | 31.40 Mb (9.73%) | 15.98 Mb (5.97%) |
| Contigs con al menos un gen | 1,265 de 1,455 | 1,415 de 1,557 |
| Transcritos completos (codón de inicio y de parada) | 99.1% | 97.8% |
| Proteínas menores de 100 aa | 5,072 (15.6%) | 2,718 (17.3%) |

La anotación predijo 32,412 genes en *P. melanocephala* y 15,690 en *P. kuehnii*. La Roya Café presenta 2.07 veces más genes en un ensamblaje solo 1.21 veces mayor, lo que eleva su densidad génica de 58.6 a 100.4 genes/Mb. La arquitectura de los genes es similar entre especies (medianas de proteína de 214 y 202 aa, y de 2 y 3 exones por gen), aunque los intrones de la Roya Naranja son en promedio más largos.

### 2.7 Redundancia de los modelos génicos

**Tabla 8. Indicadores de redundancia**

| Indicador | Roya Café | Roya Naranja |
| :--- | :--- | :--- |
| Ortólogos BUSCO duplicados | 5.80% | 2.46% |
| Genes en bloques colineales casi idénticos dentro del mismo genoma | 2,399 (7.4%) | 249 (1.6%) |
| Genes con dominios de elementos transponibles | 3,979 (12.3%) | 341 (2.2%) |
| Genes con proteína idéntica a la de otro gen | 2,864 (8.8%) | 712 (4.5%) |
| Genes con al menos un parálogo de identidad ≥ 98% | 7,212 (22.3%) | 1,077 (6.9%) |

El exceso de genes en la Roya Café debe interpretarse con cautela, porque tiene tres componentes:

1. **Redundancia residual entre haplotipos (cerca del 6–7% de los genes).** Dos medidas independientes coinciden: los BUSCO duplicados (5.80%) y los genes ubicados en bloques colineales casi idénticos dentro del mismo genoma (7.4%; identidad proteica mediana de 99.2%).
2. **Modelos derivados de elementos transponibles (12.3% de los genes).** Son más de cinco veces más frecuentes que en la Roya Naranja (2.2%).
3. **Familias multicopia.** El 22.3% de los genes tiene al menos un parálogo con identidad ≥ 98%, frente al 6.9% en la Roya Naranja.

Por tanto, el número de genes de *P. melanocephala* no es directamente comparable con el de *P. kuehnii* sin descontar estas copias.

---

## 3. Síntesis

- Se obtuvieron ensamblajes depurados de 322.68 Mb (Roya Café) y 267.54 Mb (Roya Naranja), ambos con cerca del 93% de completitud BUSCO.
- El genoma de la Roya Café es mayor y más repetitivo, según el perfil de k-meros (40% frente a 30% de secuencia repetitiva) y el tamaño del ensamblaje.
- La Roya Café tiene el doble de genes predichos (32,412 frente a 15,690). Una parte de la diferencia corresponde a modelos derivados de transposones, a familias multicopia y, en menor medida, a redundancia entre haplotipos.

## 4. Pendientes

| Pendiente | Para qué se necesita |
| :--- | :--- |
| BUSCO sobre los proteomas predichos (modo proteínas) | El BUSCO disponible mide los ensamblajes, no la anotación. |
| Versión de BRAKER y enmascaramiento del genoma | Completar los métodos; un genoma sin enmascarar infla los modelos derivados de transposones. |
| Composición taxonómica de Roya Café antes de la depuración | Completar la Tabla 4 para ambas especies. |
| Confirmar sobre qué archivo se corrió el BUSCO de cada especie | Precisar si corresponde al haplotipo 2 antes o después de la depuración. |

---

## Anexo. Origen de los datos

| Resultado | Fuente |
| :--- | :--- |
| Perfil de k-meros, ensamblajes de hifiasm, BlobToolKit y BUSCO | Presentación ICSB 2025 (`Rust_assembly_Presentation_ISCB_2025_final_version.pptx`) |
| Tamaño, contigs, N50 y GC de los ensamblajes depurados | `puccinia_only_BR.fa` y `puccinia_clean_OR.fa` |
| Estadísticas de la anotación estructural | `braker_BR.gff3` y `braker_OR.gff3` |
| Genes con dominios de transposón | `braker_eggnog_BR_2.emapper.annotations` y `braker_eggnog_OR_2.emapper.annotations` |
| Proteínas idénticas | `braker_BR.aa` y `braker_OR.aa` |
| Parálogos con identidad ≥ 98% | `royas.blast` |
| Bloques colineales dentro del mismo genoma | `royas.collinearity` |

Notas sobre la presentación: el conteo de BUSCO fragmentados de Roya Naranja aparece como 179 y corresponde a 51 (1.53% de 3,329); la diapositiva 20 rotula como haplotipo 1 las cifras del haplotipo 2 de Roya Naranja.
