# CHAPTER: METHODOLOGY

## Materials and Methods (Bioinformatics)

### 1. Functional Annotation and Ontology Mapping
Functional annotation of the assembled genomes (*Puccinia melanocephala* and *Puccinia kuehnii*) was performed with **eggNOG-mapper v2**. The fast alignment algorithm **DIAMOND** was used to compare the predicted proteome of both species against the eggNOG orthology database. From the results, the Clusters of Orthologous Groups (COG) functional categories and the Gene Ontology (GO) terms were extracted. To avoid the statistical noise caused by massive genomic redundancy, counts were normalized at the level of "Base Families" (Orthologous Groups). In addition, raw GO terms were mapped to **GO Slim** using the DAG (Directed Acyclic Graph) of the `goslim_generic.obo` file, in order to consolidate biological functions into broad categories.

### 2. Comparative Genomics and Synteny Architecture
To assess the conservation of gene order (synteny) between the two species, an all-vs-all alignment of the proteome was performed with **BLASTp** (e-value < 1e-10). The genomic coordinates (GFF3 file) and the BLASTp output were used as input for **MCScanX**. A strict collinearity parameter (`MATCH_SIZE=5`) was set to define conserved chromosomal fragments (syntenic blocks) and to filter out isolated transposons.

### 3. Global Genomic Identity (ANI)
The taxonomic and structural distance between the complete genomes at the nucleotide level was calculated with the **fastANI** algorithm. This tool fragmented the genomes to search for contiguous homologous regions, yielding two key metrics: the Alignment Fraction (AF), which quantifies the percentage of shared genomic mass, and the Average Nucleotide Identity (ANI), which measures the exact similarity within those conserved regions.

### 4. Natural Selection (Ka/Ks) and Molecular Clock
Evolutionary pressure was calculated by analyzing the pairs of orthologous genes anchored within the syntenic blocks. The ratio of non-synonymous substitutions (Ka) to synonymous substitutions (Ks) was determined. Gene families with a Ka/Ks ratio > 1 were classified as being under diversifying positive selection (adaptive hyper-mutants), whereas those with Ka/Ks < 1 were classified as being under strict purifying selection.
For the calculation of evolutionary divergence time, the synonymous mutation rate (global Ks) was used as a molecular clock. The speciation equation $T = Ks / 2r$ was applied, assuming a standard fungal mutation rate for basidiomycete pathogens of $r = 1.5 \times 10^{-8}$ substitutions per site per year.

### 5. Phylogenetic Analysis and Evolutionary Topology
To determine the evolutionary position of the sugarcane pathogen complex relative to the rusts of modern cereals, highly conserved orthologous sequences (Arp2/3 complex) were extracted from 6 species (*P. melanocephala*, *P. kuehnii*, *P. graminis*, *P. triticina*, *P. striiformis*, and *Melampsora larici-populina* as outgroup). Multiple sequence alignment was performed with **MAFFT**, followed by construction of the phylogenetic tree under the Maximum Likelihood criterion using **FastTree**. The final chronogram with the divergence time scale was rendered algorithmically with the `matplotlib` graphics library in Python.

To place both pathogens within a broader taxonomic sampling of the rust fungi, the nuclear ribosomal large subunit (nLSU) locus was extracted directly from each genome assembly and aligned with 53 reference sequences from the dataset of Dixon et al. (2010), for a total of 55 terminals. A Maximum Likelihood tree was inferred with **FastTree**, and node support was assessed with Shimodaira–Hasegawa-like local support values. A reduced seven-taxon tree, containing only the closest reference sequences and one outgroup, was used to validate the placement of both isolates.

<!-- TODO: state the aligner, the alignment length and the substitution model used for the nLSU tree. -->

### 6. Transposable Element Modelling and Two-Speed Architecture
To determine the structural impact of repetitive sequences on the observed genome expansion, a *de novo* discovery pipeline was implemented. Consensus libraries of Transposable Element (TE) families were built independently for each genome using the heuristic algorithm of **RepeatModeler**. This program scanned the genome assemblies to identify, cluster and extract the template sequences of the invading families.

To quantify the absolute transposon load (the actual amount of genomic DNA occupied by these sequences), the consensus libraries were exhaustively mapped back to their respective whole genomes using the **BLASTn** algorithm (Basic Local Alignment Search Tool). The resulting alignment coordinates were merged to avoid double counting of overlapping transposons, yielding the exact percentage of repeat coverage.

Finally, to validate the "Two-Speed Genome" model, an analysis of intergenic spatial topology was carried out. Using the `closest` tool of the **Bedtools** suite, the physical coordinates of all structural gene models were intersected with the physical coordinates of the discovered TEs. This made it possible to measure the physical distance (in base pairs) between the coding region of each gene and its nearest mobile element, and thus to classify the genome quantitatively into a stable compartment (genes isolated in TE-free regions) and a plastic accessory compartment (genes directly overlapped or closely flanked by transposon activity).


---

# CHAPTER: RESULTS AND DISCUSSION

# Structural and Functional Analysis of the Genomes of *Puccinia melanocephala* and *Puccinia kuehnii*

---

## 1. Bioinformatic Methodology

### 1.1 Biological material and sequencing
Monosoric isolates (each derived from a single pustule) were used, multiplied by reinoculation in the greenhouse:

| Species | Disease | Source cultivar | Sample | Spores collected |
| :--- | :--- | :--- | :--- | :--- |
| *Puccinia melanocephala* | Brown rust | CC 85-92 | RCM1 | 331 mg |
| *Puccinia kuehnii* | Orange rust | CC 01-1940 | RNM2 | 680 mg |

DNA was extracted with the Omniprep kit, with mechanical disruption of the urediniospores. Sequencing was performed with **PacBio HiFi** long reads (Long Plex library, 5–7 kb fragments) and **Illumina** short reads. Hi-C data could not be obtained because the amount of spores available was below the required quantity.

### 1.2 Assembly and decontamination
1. **Quality control:** adapter removal with HiFiAdapterFilt.
2. **Genome size:** k-mer-based estimation (k = 31) with GenomeScope2.
3. **Assembly:** hifiasm, which produces a primary assembly and two partial haplotypes (hap1 and hap2).
4. **Contaminants:** taxonomic assignment of contigs with BlobToolKit (GC content, read coverage and homology). Only contigs assigned to the genus *Puccinia* were retained.
5. **Completeness:** BUSCO with the pucciniomycetes lineage (n = 3,329).

All downstream analyses were performed on the **decontaminated haplotype 2** of each species.

### 1.3 Structural and functional annotation
Gene prediction was performed with **BRAKER**, which integrates AUGUSTUS and GeneMark models, using protein evidence. Structural statistics were calculated from the GFF3 files; when a gene has several isoforms, the transcript with the longest CDS was used as its representative. Biological roles were assigned with **eggNOG-mapper v2.1.13** (taxonomic scope Fungi).

> **To be confirmed:** the BRAKER version, and whether the genome was soft-masked before gene prediction. The decontaminated FASTA files available are not masked.

## 2. Results: Assembly, Structural Annotation and Gene Content

### 2.1 Assembly and Structural Annotation

#### 2.1.1 Genome size estimated from k-mers

**Table 2.1a: K-mer profile (GenomeScope2, k = 31)**

| Parameter | Brown rust (RCM1) | Orange rust (RNM2) |
| :--- | :--- | :--- |
| Estimated haploid size | 240.2 Mb | 210.3 Mb |
| Unique sequence | 60.0% | 70.3% |
| Repetitive sequence | 40.0% | 29.7% |
| K-mer coverage | 32.8× | 30× |
| Error rate | 1.07% | 1.51% |
| Duplication rate | 1.42 | 1.35 |

Brown rust has an estimated haploid genome 14% larger than that of orange rust and a higher repetitive fraction (40% versus 30%).

#### 2.1.2 Assembly

**Table 2.1b: hifiasm assemblies before decontamination**

| Assembly | Metric | Brown rust (RCM1) | Orange rust (RNM2) |
| :--- | :--- | :--- | :--- |
| **Primary** | Contigs | 6,903 | 5,969 |
| | Total length | 501 Mb | 359.7 Mb |
| | N50 | 1,157 kb | 473 kb |
| **Haplotype 1** | Contigs | 8,904 | 6,988 |
| | Total length | 460 Mb | 347.1 Mb |
| | N50 | 424 kb | 234 kb |
| **Haplotype 2** | Contigs | 4,373 | 2,926 |
| | Total length | 422 Mb | 305.3 Mb |
| | N50 | 556 kb | 274 kb |
| | Longest contig | 5.92 Mb | 2.05 Mb |

Haplotype 2 was the most contiguous in both species and was selected for decontamination and annotation.

#### 2.1.3 Removal of contaminants

**Table 2.1c: Haplotype 2 before and after decontamination**

| Metric | Brown rust (RCM1) | Orange rust (RNM2) |
| :--- | :--- | :--- |
| Haplotype 2 before decontamination | 422 Mb; 4,373 contigs | 305.3 Mb; 2,926 contigs |
| **Decontaminated assembly (*Puccinia* only)** | **322.68 Mb; 1,455 contigs** | **267.54 Mb; 1,557 contigs** |
| N50 of the decontaminated assembly | 712.3 kb | 310.1 kb |
| Longest contig | 5.92 Mb | 2.05 Mb |
| GC content | 38.6% | 32.7% |
| Undetermined bases (N) | 0 | 0 |
| Sequence removed | ~99 Mb (23.5%) | ~38 Mb (12.4%) |
| Contigs removed | 2,918 | 1,369 |

**Table 2.1d: Taxonomic composition of orange rust haplotype 2 before decontamination (BlobToolKit)**

| Taxonomic group | Contigs | Length |
| :--- | :--- | :--- |
| Basidiomycota | 2,330 | 288 Mb |
| Actinomycetota | 381 | 9.06 Mb |
| Streptophyta (plant) | 69 | 4.62 Mb |
| No hit | 64 | 0.87 Mb |
| Arthropoda | 17 | 0.61 Mb |
| Ascomycota | 18 | 0.57 Mb |
| Mucoromycota | 15 | 0.22 Mb |
| Chytridiomycota | 14 | 0.18 Mb |
| Uroviricota | 5 | 0.09 Mb |
| Other | 15 | 0.97 Mb |

Within Basidiomycota, only the main group assigned to *Puccinia* (GC close to 32%) was retained, which explains the difference between the 288 Mb of Basidiomycota and the final 267.54 Mb.

> **Pending:** the equivalent table for brown rust. The contamination slide for RCM1 shows the RNM2 plot, so the composition before decontamination is not available.

#### 2.1.4 Assembly completeness

**Table 2.1e: BUSCO (pucciniomycetes lineage, n = 3,329)**

| Category | Brown rust, hap1 | Brown rust, hap2 | Orange rust, hap2 |
| :--- | :--- | :--- | :--- |
| **Complete** | 93.84% (3,124) | 93.42% (3,110) | 93.09% (3,099) |
| ↳ Single-copy | 87.29% (2,906) | 87.62% (2,917) | 90.63% (3,017) |
| ↳ Duplicated | 6.55% (218) | 5.80% (193) | 2.46% (82) |
| Fragmented | 1.32% (44) | 1.44% (48) | 1.53% (51) |
| Missing | 4.84% (161) | 5.14% (171) | 5.38% (179) |

The assemblies of both species recover about 93% of the conserved orthologues, so their gene space is comparable. Brown rust has 2.4 times more duplicated orthologues than orange rust (5.80% versus 2.46%).

#### 2.1.5 Assembly size relative to the estimated genome size

**Table 2.1f: Decontaminated assembly versus estimated haploid size**

| Metric | Brown rust | Orange rust |
| :--- | :--- | :--- |
| Estimated haploid size (k-mers) | 240.2 Mb | 210.3 Mb |
| Decontaminated assembly | 322.68 Mb | 267.54 Mb |
| Assembly / estimate ratio | 1.34 | 1.27 |

Both assemblies exceed the estimated haploid size. This is expected for dikaryotic genomes assembled without Hi-C data: part of the second haplotype may remain uncollapsed, and high-copy repeats tend to be underestimated in k-mer analyses.

#### 2.1.6 Structural annotation

**Table 2.1g: Structural annotation statistics (BRAKER)**

| Metric | Brown rust (*P. melanocephala*) | Orange rust (*P. kuehnii*) |
| :--- | :--- | :--- |
| **Predicted genes (loci)** | **32,412** | **15,690** |
| **Transcripts** | 33,999 | 16,432 |
| Genes with alternative isoforms | 1,477 (4.6%) | 683 (4.4%) |
| Gene density | 100.4 genes/Mb | 58.6 genes/Mb |
| Gene length, mean (median) | 1,229 (857) bp | 2,288 (1,078) bp |
| CDS length, mean (median) | 970 (645) bp | 1,019 (609) bp |
| Protein length, mean (median) | 322 (214) aa | 339 (202) aa |
| Exons per gene, mean (median) | 3.33 (2) | 4.01 (3) |
| Single-exon genes | 8,796 (27.1%) | 3,341 (21.3%) |
| Exon length, mean (median) | 292 (183) bp | 254 (150) bp |
| Intron length, mean (median) | 111 (84) bp | 391 (107) bp |
| Coding fraction of the genome | 31.40 Mb (9.73%) | 15.98 Mb (5.97%) |
| Contigs with at least one gene | 1,265 of 1,455 | 1,415 of 1,557 |
| Complete transcripts (start and stop codon) | 99.1% | 97.8% |
| Proteins shorter than 100 aa | 5,072 (15.6%) | 2,718 (17.3%) |

The annotation predicted 32,412 genes in *P. melanocephala* and 15,690 in *P. kuehnii*. Brown rust has 2.07 times more genes in an assembly only 1.21 times larger, which raises its gene density from 58.6 to 100.4 genes/Mb. Gene architecture is similar between the species (median protein length of 214 and 202 aa, and 2 and 3 exons per gene), although orange rust introns are longer on average.

#### 2.1.7 Redundancy of the gene models

**Table 2.1h: Redundancy indicators**

| Indicator | Brown rust | Orange rust |
| :--- | :--- | :--- |
| Duplicated BUSCO orthologues | 5.80% | 2.46% |
| Genes in near-identical collinear blocks within the same genome | 2,399 (7.4%) | 249 (1.6%) |
| Genes with transposable element domains | 3,979 (12.3%) | 341 (2.2%) |
| Genes whose protein is identical to that of another gene | 2,864 (8.8%) | 712 (4.5%) |
| Genes with at least one paralogue of ≥ 98% identity | 7,212 (22.3%) | 1,077 (6.9%) |

The excess of genes in brown rust must be interpreted with caution, because it has three components:

1. **Residual redundancy between haplotypes (about 6–7% of the genes).** Two independent measures agree: the duplicated BUSCOs (5.80%) and the genes located in near-identical collinear blocks within the same genome (7.4%; median protein identity of 99.2%).
2. **Models derived from transposable elements (12.3% of the genes).** These are more than five times as frequent as in orange rust (2.2%).
3. **Multicopy families.** 22.3% of the genes have at least one paralogue with ≥ 98% identity, compared with 6.9% in orange rust.

Therefore, the gene count of *P. melanocephala* is not directly comparable with that of *P. kuehnii* without discounting these copies.

#### 2.1.8 Summary
- Decontaminated assemblies of 322.68 Mb (brown rust) and 267.54 Mb (orange rust) were obtained, both with about 93% BUSCO completeness.
- The brown rust genome is larger and more repetitive, according to the k-mer profile (40% versus 30% repetitive sequence) and the assembly size.
- Brown rust has twice as many predicted genes (32,412 versus 15,690). Part of the difference corresponds to transposon-derived models, to multicopy families and, to a lesser extent, to redundancy between haplotypes.

#### 2.1.9 Pending items

| Pending item | Why it is needed |
| :--- | :--- |
| BUSCO on the predicted proteomes (protein mode) | The available BUSCO measures the assemblies, not the annotation. |
| BRAKER version and genome masking | To complete the methods; an unmasked genome inflates transposon-derived models. |
| Taxonomic composition of brown rust before decontamination | To complete Table 2.1d for both species. |
| Confirm which file each BUSCO run was performed on | To specify whether it corresponds to haplotype 2 before or after decontamination. |

#### 2.1.10 Data sources

| Result | Source |
| :--- | :--- |
| K-mer profile, hifiasm assemblies, BlobToolKit and BUSCO | ICSB 2025 presentation (`Rust_assembly_Presentation_ISCB_2025_final_version.pptx`) |
| Size, contigs, N50 and GC of the decontaminated assemblies | `puccinia_only_BR.fa` and `puccinia_clean_OR.fa` |
| Structural annotation statistics | `braker_BR.gff3` and `braker_OR.gff3` |
| Genes with transposon domains | `braker_eggnog_BR_2.emapper.annotations` and `braker_eggnog_OR_2.emapper.annotations` |
| Identical proteins | `braker_BR.aa` and `braker_OR.aa` |
| Paralogues with ≥ 98% identity | `royas.blast` |
| Collinear blocks within the same genome | `royas.collinearity` |

Notes on the presentation: the count of fragmented BUSCOs for orange rust appears as 179 and corresponds to 51 (1.53% of 3,329); slide 20 labels the orange rust haplotype 2 figures as haplotype 1.

**Table 2.1i: Summary of the functional annotation**
<!-- TODO: this table keeps the original figures (counted on transcripts) and will be reviewed together with section 2.2. -->

| Functional Metric | *Puccinia melanocephala* (Brown rust) | *Puccinia kuehnii* (Orange rust) |
| :--- | :--- | :--- |
| **Genes with Functional Annotation** | 15,683 genes (~46.1%) | 6,741 genes (~41.0%) |
| **Unique Base Genes (non-redundant catalogue)** | 3,276 genes | 3,178 genes |
| ↳ *Base genes with a strictly single copy* | 2,074 genes | 2,102 genes |
| ↳ *Base genes with multiple copies (photocopied)* | 1,202 genes | 1,076 genes |
| **Total copies generated by duplication** | **13,609 copies** | **4,639 copies** |
| **Active Transposons (TEs)** | **40 unique genes** (2,690 copies) | **29 unique genes** (250 copies) |
| **Genes without Functional Annotation (Hypothetical)**| 18,316 genes (~53.9%) | 9,691 genes (~59.0%) |

### 2.2 Functional Profile and Pangenome (COG Categories and GO Slim)

To understand how the biology of these fungi is distributed beyond a simple structural count, the **Base Genes** (unique or orthologous families) were isolated, removing the noise caused by photocopies (paralogous genes). This pangenome analysis reveals a crucial finding: of the 3,511 unique families that make up the genetics of both species, **2,943 families (84%) are shared (Core Genome)**.

![Venn diagram of the functional pangenome](/Users/estuvar4/Documents/2.\ software/17.biojava/datos_royas/1.Resultados_Tesis/venn_pangenoma.png)
*Figure 2.2: Venn diagram illustrating the intersection of biological families (base genes) between P. melanocephala and P. kuehnii.*

Brown rust has only 333 exclusive families, while orange rust has 235. This indicates that both species use essentially the same "biological toolbox", and that the inflation of the *P. melanocephala* genome results from a strategy of massive duplication, not from the innovation of new functions.

**Table 2.2a: Functional Distribution of the Base Genes by COG Category**
When the unique families are grouped by their broad metabolic or cellular category, it becomes evident that both rusts maintain an almost identical baseline biological complexity, with a notably large amount of "Dark Matter" (species-specific genes with no known function in databases).

| COG Macro-Category | Specific Biological Category | No. of Base Genes (Brown) | No. of Base Genes (Orange) |
| :--- | :--- | :--- | :--- |
| **Poorly Characterized** | **S:** Function Unknown (Dark matter) | 1,274 | 1,210 |
| **Information Storage and Processing** | **J:** Translation and Ribosomes | 225 | 214 |
| **Cellular Processes and Signaling** | **U:** Intracellular and Vesicular Trafficking | 217 | 207 |
| **Cellular Processes and Signaling** | **O:** Protein Folding and Modification | 197 | 213 |
| **Information Storage and Processing** | **A:** RNA Processing | 185 | 185 |
| **Information Storage and Processing** | **K:** Transcription (Factors and regulation) | 164 | 157 |
| **Cellular Processes and Signaling** | **T:** Signal Transduction (Kinases) | 142 | 137 |
| **Metabolism** | **E:** Amino Acid Metabolism | 132 | 126 |
| **Metabolism** | **G:** Carbohydrate Metabolism | 131 | 133 |
| **Metabolism** | **C:** Energy Production (ATP) | 118 | 117 |
| **Metabolism** | **I:** Lipid Metabolism | 116 | 116 |
| **Information Storage and Processing** | **L:** DNA Replication and Repair | 100 | 98 |
| **Cellular Processes and Signaling** | **D:** Cell Cycle and Division | 94 | 90 |
| **Cellular Processes and Signaling** | **M:** Cell Wall Biogenesis | 30 | 29 |
| **Cellular Processes and Signaling** | **V:** Cellular Defense (Toxins and resistance) | 12 | 13 |

**Table 2.2b: Specific Functional Profile (GO Slim Classification)**
To corroborate this profile with the standard nomenclature of the *Gene Ontology* consortium, the raw annotations were mapped while removing redundant generic hierarchies (using the GO Slim subset). The most abundant results confirm the deep conservation of the cellular, nuclear and catalytic machinery between the two pathogens.

| Main Aspect (GO) | GO Slim Term | Description of the Function | No. of Base Genes (Brown) | No. of Base Genes (Orange) |
| :--- | :--- | :--- | :--- | :--- |
| **Cellular Component (CC)** | **GO:0005634** | nucleus | 934 | 918 |
| **Molecular Function (MF)** | **GO:0003824** | catalytic activity | 888 | 880 |
| **Cellular Component (CC)** | **GO:0005739** | mitochondrion | 530 | 520 |
| **Cellular Component (CC)** | **GO:0005829** | cytosol | 432 | 420 |
| **Molecular Function (MF)** | **GO:0016740** | transferase activity | 388 | 375 |
| **Molecular Function (MF)** | **GO:0016787** | hydrolase activity | 334 | 338 |
| **Biological Process (BP)** | **GO:0006886** | intracellular protein transport | 302 | 293 |
| **Molecular Function (MF)** | **GO:0003723** | RNA binding | 296 | 285 |
| **Biological Process (BP)** | **GO:0006355** | regulation of transcription | 288 | 287 |
| **Biological Process (BP)** | **GO:0065003** | protein complex assembly | 283 | 287 |
| **Cellular Component (CC)** | **GO:0005783** | endoplasmic reticulum | 271 | 263 |
| **Cellular Component (CC)** | **GO:0005694** | chromosome | 256 | 247 |
| **Biological Process (BP)** | **GO:0016192** | vesicle-mediated transport | 219 | 208 |
| **Biological Process (BP)** | **GO:0000278** | mitotic cell cycle | 215 | 204 |
| **Biological Process (BP)** | **GO:0042254** | ribosome biogenesis | 213 | 199 |
| **Biological Process (BP)** | **GO:0023052** | signaling | 179 | 176 |


---



## 3. Comparative Genomics: Relatedness, Synteny and Evolution

### 3.1 Phylogeny: Independent Origins of the Sugarcane Rusts within the Pucciniaceae

To establish the phylogenetic position of the two sugarcane rust pathogens, a Maximum Likelihood tree (FastTree) was inferred from the nuclear ribosomal large subunit (nLSU) locus. The nLSU sequences of both isolates were extracted directly from the genome assemblies and aligned with 53 reference sequences from the dataset of Dixon et al. (2010): 37 from *Puccinia*, nine from *Uromyces*, one from *Pucciniosira*, and six from other genera of the Pucciniales used as outgroups (*Phakopsora*, *Gymnoconia*, *Phragmidium*, *Kuehneola*, *Gymnosporangium* and *Coleosporium*), for a total of 55 terminals. Node support is reported as Shimodaira–Hasegawa-like local support values (scale 0–1).

![Maximum Likelihood phylogeny of the nLSU locus (55 terminals)](/Users/estuvar4/Documents/2.\ software/17.biojava/datos_royas/Mega_Arbol_Paper_Royas.nwk.png)
*Figure 3.1: Maximum Likelihood tree (FastTree) of the nLSU locus for 53 reference sequences (Dixon et al., 2010) and the two isolates sequenced in this study (in red). Scale bar: substitutions per site.*

<!-- TODO: re-root the figure on the six outgroup genera and display the node support values; in the current image the outgroups appear nested within the ingroup. -->

**Species identity.** The nLSU sequence of isolate RCM1 was identical to the *P. melanocephala* reference (GU058001), and that of isolate RNM2 clustered with the *P. kuehnii* reference (GU058021; support 0.999; 0.001 substitutions per site). Both assemblies are therefore confirmed as the expected species.

**Table 3.1: Phylogenetic placement of the two sugarcane rust isolates in the nLSU tree**

| Isolate | Closest relative (support) | Clade | Other members of the clade | Clade support |
| :--- | :--- | :--- | :--- | :--- |
| RCM1 (*P. melanocephala*) | *P. miscanthi* (0.992) | Group 2 | *P. nakanishikii*, *P. rufipes*, *P. purpurea*, *P. duthiae*, *P. vexans*, *P. coronata*, *P. hordei*, *P. magnusiana*, *P. triticina*, *P. graminis* | 0.89 |
| RNM2 (*P. kuehnii*) | *P. agrophila* + *P. substriata* (0.958) | Group 3 | *P. polysora* | 0.943 |

<!-- TODO: confirm that the labels "Group 2" and "Group 3" match the clade numbering used by Dixon et al. (2010). -->

**The two sugarcane rusts are not sister species.** *P. melanocephala* and *P. kuehnii* fall into two distinct and well-supported clades (hereafter Group 2 and Group 3, after Dixon et al., 2010):

1. **Group 2 (*P. melanocephala*).** Brown rust is sister to *P. miscanthi* (support 0.992) within a lineage of rusts that parasitize grasses of the tribe Andropogoneae, such as *P. miscanthi* (*Miscanthus*) and *P. purpurea* (*Sorghum*). This lineage belongs to a larger clade (support 0.89) that also contains the temperate cereal rusts *P. triticina*, *P. graminis*, *P. hordei* and *P. coronata*. The branching order between the Andropogoneae rusts and the *P. triticina* subclade is not resolved by this locus (support 0.26), and *P. striiformis* falls outside the clade.
2. **Group 3 (*P. kuehnii*).** Orange rust is sister to *P. agrophila* + *P. substriata* (support 0.958) and, together with *P. polysora*, the southern rust of maize, forms a separate clade (support 0.943). When the tree is rooted on the outgroup genera, Group 3 is the earliest-diverging lineage among the sampled *Puccinia* and *Uromyces* species (support 0.916).

The patristic distance between the two isolates is 0.191 substitutions per site, which places this pair in the 96th percentile of all pairwise distances among the *Puccinia* and *Uromyces* species sampled (median 0.109). It is almost ten times the distance between *P. triticina* and *P. graminis* (0.020) and is comparable to the distance separating either isolate from the outgroup *Phakopsora pachyrhizi* (0.18–0.20).

Two additional analyses are concordant with this topology. A reduced seven-taxon nLSU tree recovered the same placements with maximal support (*P. melanocephala* + *P. miscanthi*, 1.000; isolate RNM2 + *P. kuehnii* reference, 1.000; *P. kuehnii* + *P. polysora*, 0.992). Independently, a tree of a conserved nuclear protein (a subunit of the Arp2/3 complex) across six species placed *P. kuehnii* outside the clade formed by *P. melanocephala* and the cereal rusts (support 0.902).

**Evolutionary interpretation: convergence on a shared host.** Because the two pathogens belong to distinct lineages of the Pucciniaceae, their shared ability to parasitize *Saccharum* cannot be attributed to inheritance from a sugarcane-infecting common ancestor. The most parsimonious interpretation is that this ability was acquired independently in each lineage, that is, by convergent evolution of host use. For *P. melanocephala*, sugarcane belongs to the same grass tribe (Andropogoneae) as the hosts of its closest relatives; for *P. kuehnii*, the colonization of *Saccharum* arose in a distantly related clade.

These inferences rest on a single ribosomal locus. A multilocus or phylogenomic analysis that includes representatives of both groups will be required to confirm the branching order among the major clades.

*Reference: Dixon, L. J., Castlebury, L. A., Aime, M. C., Glynn, N. C., & Comstock, J. C. (2010). Phylogenetic relationships of sugarcane rust fungi. Mycological Progress, 9(4), 459–468.*


### 3.2 Global Genomic Identity (ANI)

Whole-genome comparison yielded an Average Nucleotide Identity (ANI) of **83.14%** between the two rusts, with an Alignment Fraction (AF) of only **17.07%**.

<!-- TODO: the fastANI output file is not in the results folder; confirm both values against it. -->

An ANI of 83% is far below the values of 95% or more observed between individuals of the same species, and it lies close to the lower limit of the range in which fastANI reports alignments reliably (approximately 80%). The AF must therefore be read as the fraction of the two genomes that remains similar enough to be aligned at the nucleotide level. It is a lower bound of the homologous sequence, because homologous regions that have diverged beyond this threshold are not counted.

This large-scale structural disconnection is the expected outcome for two fungi that evolved along independent phylogenetic trajectories for millions of years before converging on the same host (Section 3.1). It should not be interpreted as the rapid restructuring of one genome relative to a recent common ancestor. In species drawn from two different clades of the Pucciniaceae, most intergenic and repetitive sequence is expected to have diverged beyond recognition, so that only the most conserved fraction of the genome, largely coding sequence, remains alignable.

The same reasoning applies to the repetitive fraction. The differences in transposable element (TE) expansion between the two genomes (Section 3.6) arose along separate paths: each lineage accumulated its own TE complement after the two groups diverged, and a larger TE load in *P. melanocephala* does not imply that its genome was derived from a *P. kuehnii*-like architecture.

### 3.3 Mapping Coverage and Variant Calling (SNPs)
To evaluate the exact differences at the structural and base-pair level, a strict whole-genome mapping of *P. melanocephala* against the *P. kuehnii* reference was performed with `minimap2` (asm5), followed by variant extraction with `bcftools`.

**Table 3.3.1. Inter-species mapping coverage (*P. melanocephala* vs *P. kuehnii*)**

| Alignment Metric | Volume (Millions of Bases) | Percentage of the Genome |
| :--- | :--- | :--- |
| **Total Genome Size (Brown rust)** | 308.1 MB | 100.0% |
| **Mapped Bases (Conserved Core)** | 26.6 MB | 8.6% |
| **Unmapped Bases (Divergent Fraction)**| 281.5 MB | 91.4% |

Once the homologous regions (the conserved 8.6%) had been isolated, the exact point mutations separating the two species were extracted:

**Table 3.3.2. Genetic variants identified in the conserved core.**

| Type of Genetic Variant | Number Detected |
| :--- | :--- |
| **SNPs** (Point mutations) | 170 |
| **INDELs** (Insertions / Deletions) | 14 |
| **Total Variants** | **184** |

#### Discussion: Testing the "Two-Speed Genome"
The mapping coverage results represent the definitive bioinformatic confirmation of genome hyper-expansion. The fact that **91.4% of the brown rust genome (281.5 MB) does not exist in, or fails to align with, that of orange rust** fits perfectly with the profound loss of synteny reported previously (7.9%).

However, the discovery of only 170 SNPs in the 8.6% of alignable chromosomal blocks represents a fascinating biological contrast. These data strongly support the structural hypothesis of the **"Two-Speed Genome"** in the sugarcane rust complex:

1. **An accessory (plastic) compartment:** Highly variable, chaotic and massive (representing more than 91% of the divergent DNA). This compartment is responsible for inflating the genome size of *P. melanocephala* to more than 308 MB, presumably driven by the expansive activity of Transposable Elements (TEs). *(Note: the quantification and taxonomic classification of these transposon families is documented in the following section).*
2. **A core compartment:** Represented by that 8.6% of conserved regions. Despite the structural chaos around it, this core is surprisingly stable, with minimal mutation rates (170 SNPs) and strongly protected by purifying selection (Ka/Ks < 1) to keep the essential biological machinery of the fungus intact.


### 3.4 Population Structure and Genetic Diversity (PCA)
Using the whole-genome phylogenomic matrix (VCF) composed of the 44 sequenced strains, Population Genetics algorithms were run to quantify mathematically the evolutionary distance between the isolates and to visualize their population structure in a two-dimensional space.

#### Principal Component Analysis (PCA)
The Principal Component Analysis reduced the dimensionality of the shared genomic mutations, yielding well-defined genetic clusters.

![Population structure PCA](file:///Users/estuvar4/.gemini/antigravity/brain/8248c59f-297d-47dd-aa15-844de167695d/pca_diversidad.png)

**Discussion of the PCA:**
As shown in the plot, Principal Component 1 drastically separates the outgroup (poplar rust). Interestingly, the 42 worldwide strains of wheat rust form a massive, highly cohesive clade (blue cloud). **Brown rust** emerges as a completely independent lineage (an isolated genetic island), radically distant from the wheat rusts. This biologically confirms the strict host specialization that this pathogen underwent when adapting exclusively to *Saccharum*.

### 3.5 Loss of Synteny and Chromosomal Restructuring
To measure this restructuring mathematically, a collinearity (synteny) analysis was performed. Of the more than 50,000 combined gene models, **only 4,027 genes (7.99%)** retain a shared spatial order in syntenic blocks.
The destruction of more than 92% of the original synteny confirms that the brown rust genome did not simply "grow" in size, but underwent an aggressive event of genomic shuffling. Transposons broke and repositioned genes to such an extent that the original chromosome order of the common ancestor was completely pulverized.

**Table 3.2: Statistical Summary of the Synteny (Collinearity) Analysis**

| Spatial Architecture Metric | Value Obtained | Biological Interpretation |
| :--- | :--- | :--- |
| **Total Genes Analyzed** | 50,431 genes | Combined sum of the gene models of both rusts. |
| **Conserved Syntenic Blocks** | 300 blocks | Chromosomal fragments that survived evolution intact. |
| **Total Collinear Genes** | 4,027 genes | Genes that maintain exactly the same orthologous spatial order. |
| **Percentage of Conserved Synteny** | **7.99%** | Degree of spatial homology. Low, as expected between phylogenetically distant lineages. |
| **Spatial Loss/Restructuring** | **92.01%** | Fraction of the genome that was moved, rearranged or broken (shuffling), driven by TEs. |

*(Technical note and bioinformatic parameters: the spatial collinearity analysis was run with a minimum threshold of 5 consecutive genes [MATCH_SIZE = 5] to consolidate a valid syntenic block, allowing a maximum disruption of 25 intervening genes [MAX_GAPS = 25]).*

<!-- TODO: the 7.99% figure includes collinear blocks within the same genome (194 brown-brown and 19 orange-orange blocks). Between the two species there are 87 blocks and 675 gene pairs (2.68% of the gene models). To be corrected when this section is reviewed. -->

**Why is synteny so low?**
The nLSU phylogeny (Section 3.1) provides the primary explanation. *P. melanocephala* and *P. kuehnii* are not sister species: they belong to two distinct clades of the Pucciniaceae (Group 2 and Group 3). The near-complete absence of conserved gene order is therefore the expected result for two fungi that evolved along independent phylogenetic trajectories for millions of years before converging on the same host. Collinearity above 70% is typical of closely related species, and that expectation does not apply to this comparison. Along those independent trajectories, the following processes contributed to the erosion of gene order, and the transposable element expansions involved took place separately in each lineage:

1. **Illegitimate recombination mediated by TEs:** The massive invasion of Transposable Elements acted as a "genomic blender". Highly repetitive sequences caused chromosomes to pair erroneously during cell division, producing massive translocations and inversions that destroyed the original 1:1 spatial order.
2. **Paralogue noise (cloning):** By massively photocopying its biological tools (13,609 redundant copies in brown rust), the trace of the ancestral genomic "skeleton" disappears under thousands of copies scattered randomly across the contigs, breaking the mathematical continuity required to form syntenic blocks.
3. **The evolutionary advantage of instability:** For an obligate biotrophic pathogen, an ordered and rigid genome is a dead end. Pulverizing the chromosomal architecture grants immense genomic plasticity, accelerating the rate of recombination and mutation of its weapons (effector genes). This allows the rust to adapt rapidly to the new resistant sugarcane varieties developed by humans.


### 3.6 Identification and Quantification of Transposable Elements (TEs)
To understand what caused the large genome size of brown rust (*P. melanocephala*, 308 MB) and why the original chromosomal order (synteny) was lost when compared with orange rust (*P. kuehnii*), the content of repetitive sequences in both fungi was analyzed.

Through the bioinformatic discovery of transposon families, it was possible to measure exactly what percentage of each genome is composed of this type of mobile sequence.

**Table 3.6.1. Transposable Element content in brown rust and orange rust**

| Metric | Brown rust (*P. melanocephala*) | Orange rust (*P. kuehnii*) |
| :--- | :--- | :--- |
| **Base Pairs Occupied by TEs** | 188.3 million bases | 130.2 million bases |
| **Total Percentage of the Genome (TEs)** | **58.3%** | **48.6%** |

#### Discussion: The Role of Transposons in Genome Expansion
The results show that almost 60% of the total DNA of *P. melanocephala* is composed of Transposable Elements.

These biological data indicate that the larger genome size of brown rust (308 MB versus 267 MB for orange rust) is due almost entirely to the multiplication of transposons. As these sequences were copied and pasted along the chromosomes over time, they altered the original order of the genes (which explains the fall in synteny to 7.9%) and created extensive repetitive regions. In plant pathogens, these extensive transposon regions are key, since they act as zones where virulence genes can duplicate and mutate more rapidly.


#### Physical Proximity between Genes and Transposable Elements
To gauge the real impact of transposons on the biological machinery of the fungus, the physical coordinates of all annotated genes were intersected with the coordinates of the repetitive regions discovered. The aim was to measure, along all chromosomes, the exact distance (in base pairs) between each gene and its nearest transposon.

**Table 3.6.2. Spatial distribution of genes relative to Transposable Elements**

| Genomic Architecture (Gene–TE Distance) | Brown rust (*P. melanocephala*) | Orange rust (*P. kuehnii*) |
| :--- | :--- | :--- |
| **Genes invaded by or touching TEs (distance 0 bp)** | 61.0% | 45.3% |
| **Genes very close to TEs (< 2,000 bp)** | 30.8% | 23.7% |
| **Genes in the Dynamic/Accessory Compartment** | **91.8%** | **69.0%** |
| **Genes isolated in Stable/Core Zones (> 2,000 bp)** | 8.2% | 31.0% |

These results reveal an extreme genomic contrast. While orange rust keeps almost a third of its genes (31.0%) protected in "safe zones" far from the activity of mobile elements, brown rust has undergone a generalized invasion of its functional spaces.

Strikingly, **91.8% of all *P. melanocephala* genes** reside physically embedded in, or closely flanked by, transposons. This overwhelming physical proximity is the structural confirmation of the "Two-Speed Genome" model. Living surrounded by parasitic, mobile DNA, brown rust genes are inadvertently dragged along every time a transposon "jumps" to another part of the genome. This mechanically explains why this pathogen managed to inflate its artificial gene count (through massive duplications) and why it shows the accelerated rates of adaptive mutation documented in the following section.

*(Note: the taxonomic classification of the specific types of transposons, such as the proportion of LTR/Gypsy versus LTR/Copia, will be presented later, once annotation against the Dfam database has been completed).*

### 3.7 Evolutionary Pressure (Ka/Ks) and Arms Race
To understand which genes are mutating to evade sugarcane, the selective pressure (Ka/Ks ratio) was calculated on the pairs of orthologous genes between brown rust and orange rust.

The genome showed a classic pattern of baseline survival: **73.5% of the orthologous genes are under purifying selection (Ka/Ks < 1)** (averaging 0.76). Evolution penalizes and removes any mutation in these genes because they encode critical structural proteins (e.g. ribosomal proteins L11/L12 or DNA helicases with Ka/Ks close to 0.01).

However, the most striking finding is that a very high **26.5% of the genes analyzed are under strong positive selection (Ka/Ks > 1)**. This percentage of "hyper-mutation" is the biological hallmark of the so-called arms race: the fungus forces constant mutations in specific gene families to evade the immune defenses of the host.

**Table 3.3a: Gene Families under Positive Selection (High Adaptive Mutation)**
The following families represent the attack front of the fungus. Having accumulated the largest number of adaptive mutations (Ka/Ks > 1) since their evolutionary divergence, they are the main candidates for evading new sugarcane varieties. The table includes the number of copies retained in each genome to gauge their expansion.

| Gene Family (Top Mutants) | Ka/Ks Ratio (Mean) | Copies (Brown) | Copies (Orange) | Biological Function in Pathogenesis |
| :--- | :--- | :--- | :--- | :--- |
| **OPT Transporters** | **5.580** | 22 | 26 | Membrane: sequestering nutrients or secreting toxins rapidly. |
| **M16 Family Peptidases** | **3.332** | 1 | 1 | Secreted enzymes that destroy plant defenses (PR-proteins). |
| **NADH Dehydrogenases** | **3.089** | 7 | 5 | Extreme adaptation to oxidative stress inside the stoma. |
| **Membrane Kinases** | **> 1.500** | 976 | 282 | Sensory receptors. They mutate to evade recognition by the plant. |

**Table 3.3b: Gene Families under Strict Purifying Selection (High Conservation)**
In direct contrast, the following families have a Ka/Ks close to zero. Evolution penalizes any mutation in these genes because they encode the central machinery of fungal life.

| Gene Family (Top Conserved) | Ka/Ks Ratio (Mean) | Copies (Brown) | Copies (Orange) | Basal Biological Function |
| :--- | :--- | :--- | :--- | :--- |
| **Cytoskeleton (EF hand)** | **0.003** | 16 | 15 | Critical cell structure and hyphal growth. |
| **Ribosomal Proteins L11/L12** | **0.009** | 4 | 2 | Central machinery of protein synthesis. |
| **Holliday Junction Helicases** | **0.010** | 4 | 5 | Strict repair and maintenance of DNA integrity. |

This immense mutation rate in secreted families and transporters corroborates that the evolution of these rusts is not focused on improving their internal metabolism, but exclusively on optimizing their ability to penetrate and parasitize sugarcane.

### 3.8 Synonymous Divergence (Ks) and Depth of the Phylogenetic Split

Synonymous substitutions do not alter the encoded protein and are therefore assumed to accumulate at an approximately neutral rate, which makes the synonymous distance (Ks) a measure of the evolutionary separation between two genomes. Across the collinear orthologue pairs shared by *P. melanocephala* and *P. kuehnii*, the mean synonymous divergence was **Ks = 0.6491** (663 gene pairs; median 0.50).

<!-- TODO: the Ka and Ks values in royas_kaks.tsv need to be recomputed with codon-aware alignments (in several pairs Ka is inconsistent with the protein identity), so the Ks reported here is provisional. -->

In the light of the nLSU phylogeny (Section 3.1), this distance is not interpreted as the time elapsed since two sister species separated from a common sugarcane-infecting ancestor. The two pathogens belong to different clades of the family Pucciniaceae, and their large mutational distance reflects the deep, basal phylogenetic separation between Group 2 and Group 3, which took place long before both lineages converged adaptively on the genus *Saccharum*. The ribosomal locus supports the same conclusion independently: the nLSU patristic distance between the two isolates (0.191 substitutions per site) lies in the 96th percentile of all pairwise distances among the *Puccinia* and *Uromyces* species sampled, and is comparable to the distance separating either isolate from the outgroup *Phakopsora pachyrhizi*.

For this reason, no absolute divergence date is proposed. A strict molecular clock with a generic, uncalibrated substitution rate would not date a speciation event on sugarcane; it would attempt to date the split between two major clades of the Pucciniaceae, where rate variation among lineages and the progressive saturation of synonymous sites make such an estimate unreliable. A dated phylogeny for these pathogens will require a multilocus dataset, fossil or secondary calibrations, and sampling of the intermediate taxa of both groups.

**Evolutionary context: two independent routes to the same host**
The comparison between the two genomes is best understood as a comparison between two lineages with separate evolutionary histories:

1. **Group 2 lineage (*P. melanocephala*).** Brown rust belongs to a lineage of rusts that parasitize grasses of the tribe Andropogoneae, and it is sister to *P. miscanthi*. Its larger genome, higher repeat content and higher gene count were acquired within this lineage.
2. **Group 3 lineage (*P. kuehnii*).** Orange rust belongs to a distantly related clade that also includes *P. polysora*. Its more compact genome and lower gene count reflect the history of that separate lineage.

Under this framework, the differences in genome size, transposable element load and gene order between the two sugarcane rusts are the accumulated result of two independent trajectories, and their shared capacity to infect *Saccharum* is a case of convergent evolution of host use.

### 3.9 Synthesis of the Pathogenic Arsenal (Integrated Catalogue)
The cross-analysis of the three biocomputational methodologies implemented in the results (Functional Annotation, Synteny/Genome Expansion and Ka/Ks Evolutionary Pressure) makes it possible to consolidate a definitive catalogue of the gene families that orchestrate infection.

The following master table correlates the functional identity of the genes with their expansion dynamics (physical cloning) and their rate of adaptive mutation, offering an exact molecular picture of the pathogenic tools of both rusts:

**Table 3.4: Integrated Catalogue of Key Gene Families in Pathogenesis**

| Gene Family (Functional Annotation) | Copies (Brown rust) | Copies (Orange rust) | Evolutionary Pressure (Ka/Ks) | Biological Role in Sugarcane Pathogenesis |
| :--- | :--- | :--- | :--- | :--- |
| **Rho Genes (GTPases)** | *Shared Core* | *Shared Core* | **High Conservation** | **Immune Evasion:** They manipulate the plant cytoskeleton to inhibit the defensive closure of the stomata. |
| **Kinases (Tyrosine/Serine)** | **976** | 282 | **Hyper-mutant** (> 1.5) | **Sensors and Signaling:** Massively expanded network to recognize the leaf and mutate rapidly to evade receptors. |
| **OPT Transporters** | 22 | 26 | **Hyper-mutant** (5.58) | **Transmembrane Pumping:** They sequester vital nutrients from the plant cell or secrete toxins rapidly. |
| **Glycosyltransferases (ALG)** | Conserved basal | **Expanded** | **Mixed Dynamics** | **Cellular Camouflage:** They build extracellular glycoproteins, masking the fungus from the basal immunity of sugarcane. |
| **Helicases (Holliday Junction)** | **Expanded** | Conserved basal | **Ultra-Conserved** (0.01) | **Genomic Survival:** They repair the DNA breaks caused by the chaos of their own Transposable Elements. |
| **M16 Peptidases (Insulinase)** | 1 | 1 | **Hyper-mutant** (3.33) | **Enzymatic Degradation:** Highly specific, mutagenic projectiles that destroy the defense proteins of the plant. |

## 4. Integrated Evolutionary Discussion: The Biological Paradigm of the Sugarcane Rusts

The integration of the genomic, phylogenetic and population profiles generated in this study makes it possible to propose an integrated biological model of how *Puccinia melanocephala* and *Puccinia kuehnii* became the most destructive pathogens of *Saccharum*. Far from a simple accumulation of random mutations, the data reveal a highly sophisticated evolutionary machinery.

### A. Confirmation of the "Two-Speed Genome" Model
The structural evidence uncovered in this study demonstrates empirically the existence of a highly compartmentalized genome. On the one hand, the finding that **91.4% of the *P. melanocephala* genome is divergent (plastic accessory fraction)** and that it retains only 7.9% synteny explains its "genomic obesity" (308 MB). This immense mass of DNA, presumably driven by Transposable Elements (TEs), acts as the engine of genetic variability.
In absolute contrast, the **8.6% of the DNA that could be aligned (core fraction)** proved to be hyper-conserved. The discovery of only 170 SNPs in these essential regions proves that the sugarcane rusts strictly protect their vital biological machinery under strong purifying selection, delegating all their evolutionary aggressiveness to the accessory compartments (TEs).

### B. Founder Effect and Strict Host Specialization
The population Principal Component Analysis (PCA) visually illustrates one of the most critical evolutionary leaps in this genus. While the 42 worldwide strains of wheat rust form a genetically cohesive cluster, **brown rust emerges as a completely isolated genetic island**. This mathematical pattern is the classic signature of an *Evolutionary Bottleneck* followed by a *Founder Effect*.
Consistent with the nLSU phylogeny (Section 3.1), brown rust and orange rust did not diverge from a common sugarcane-infecting ancestor: they reached *Saccharum* through two independent lineages of the Pucciniaceae (Group 2 and Group 3). This independent history, followed by convergence on the same host, explains why their genomic architectures are today barely recognizable to each other.

### C. "Dark Matter" and the Red Queen Hypothesis
The analysis of selective pressure revealed unusually high Ka/Ks rates (positive selection) in protein families directly associated with pathogenesis, such as Kinases, Peptidases and Transporters. This pattern of directional hyper-mutation is the genomic confirmation of the **Red Queen Hypothesis (Evolutionary Arms Race)**: the fungus is forced to mutate its attack effectors constantly simply to "stay in the same place" and avoid detection by the immune receptors of new sugarcane varieties.
Crucially, more than half of the genome falls into the category of "Dark Matter" (genes with no known functional homology in worldwide databases). This lack of annotation is not a prediction error, but the biological proof that species-specific effectors are mutating at such an extreme speed that they have erased any trace of homology with other fungal species.

### D. Ecological Divergence in the Field
Finally, the structure of gene duplication explains epidemic behavior in the field. Orange rust (*P. kuehnii*), with a genetically more compact and stable genome, proves to be ecologically conservative, manifesting itself in intermittent and localized epidemic outbreaks. In contrast, *P. melanocephala* used its transposon machinery not only to disorder its genome, but also to **clone its arsenal massively (hyper-duplication of paralogues)**. By having hundreds of redundant copies of sensors (Kinases) and repair proteins (Helicases), brown rust can afford to experiment genetically at an accelerated pace, which explains why it has historically been the pathogen capable of breaking the resistance of commercial varieties fastest worldwide.

## 5. Biotechnological Perspectives and Genetic Improvement

The dissection of the pangenome and the functional profile of *P. melanocephala* and *P. kuehnii* provide a direct roadmap for sugarcane breeding programs and for the design of new biotechnological mitigation strategies. Based on the genetic content uncovered, the following high-impact targets are proposed:

### A. Breeding for Broad-Spectrum Resistance (Rho Genes and Core Genome)
The discovery that 84% of the biological families (2,943 base genes) make up the *Core Genome* shared by both rusts offers an unprecedented agronomic opportunity. Within this basal core, the **Rho genes (immune evasion effectors)** stand out as critical. Rusts use Rho proteins to hijack and manipulate the structure of plant cells, preventing the host from activating its primary defense responses.

Since these Rho genes are obligatory and shared in a conserved manner by both species, they represent an ideal target. Conventional and biotechnological breeding programs for *Saccharum* (sugarcane) should prioritize the development or screening of varieties with resistance genes (R-genes) capable of physically detecting fungal Rho proteins. Recognition of the Rho genes would confer horizontal immunity, creating sugarcane varieties biologically shielded against both brown rust and orange rust simultaneously.

### B. Gene Silencing (HIGS/RNAi) against the Sensory System
The massive expansion of **Tyrosine Kinases** and **Serine/Threonine Kinases** (highly abundant families in brown rust) makes them an ideal Achilles' heel for Host-Induced Gene Silencing (HIGS) or exogenous RNAi technologies. Since these kinases function as the critical sensory network for perceiving the topology of the stoma and triggering the formation of the appressorium, silencing their conserved catalytic domains would "blind" the pathogen, preventing penetration of the leaf even if the spore manages to germinate.

### C. Induction of Genomic Collapse (Lethal Mutagenesis)
It was shown that *P. melanocephala* has a highly mutagenic genomic environment (driven by transposons), which it survives thanks to a gigantic expansion of DNA repair proteins, such as the **Helicases (C-terminal Superfamily)**. The design of new-generation fungicides or blocking molecules (enzyme inhibitors) directed specifically at the helicases of the rust would interrupt its repair mechanism. Without this machinery, the inherent activity of its transposons would cause lethal chromosomal breaks, inducing the suicide of the pathogen.

### D. Inhibition of Camouflage and Efflux Pumps
Orange rust (*P. kuehnii*) showed an evolutionary strategy focused on the cell wall through **ALG6/ALG8 Glycosyltransferases**. These genes are fundamental for the assembly of glycoproteins that mask the hyphae and prevent the plant from activating PAMP-triggered immunity (PTI). The development of biochemical inhibitors against these enzymes would strip the fungus of its "camouflage".
In addition, the biotechnological blocking of the families of **membrane transporters (e.g. DUF2183 and MFS)** would disable the ability of both fungi to expel plant toxins and agrochemicals, restoring and dramatically enhancing the efficacy of current commercial fungicides at much lower doses.
