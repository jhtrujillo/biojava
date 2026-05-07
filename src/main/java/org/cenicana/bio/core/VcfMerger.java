package org.cenicana.bio.core;

import java.io.*;
import java.util.*;

/**
 * High-performance VCF merger using Multi-way Merge algorithm.
 * Streams data from multiple files simultaneously to join them by coordinates.
 * Samples missing in a source file are filled with "./.".
 */
public class VcfMerger {

    private List<String> inputFiles;
    private File outputFile;
    private double maf = 0.0;
    private double minCallRate = 0.0;

    public VcfMerger(List<String> inputFiles, String outputPath) {
        this.inputFiles = inputFiles;
        this.outputFile = new File(outputPath);
    }

    public VcfMerger(List<String> inputFiles, String outputPath, double maf, double minCallRate) {
        this.inputFiles = inputFiles;
        this.outputFile = new File(outputPath);
        this.maf = maf;
        this.minCallRate = minCallRate;
    }

    public void merge() throws IOException {
        List<BufferedReader> readers = new ArrayList<>();
        List<String> allSamples = new ArrayList<>();
        Map<Integer, List<String>> fileToSamples = new HashMap<>();

        // 1. Initial Scan for Samples and Headers
        System.out.println("[Vcf-Merge] Collecting sample headers from " + inputFiles.size() + " files...");
        for (int i = 0; i < inputFiles.size(); i++) {
            File f = new File(inputFiles.get(i));
            System.out.println("[Vcf-Merge] Scanning file: " + f.getName() + " (exists: " + f.exists() + ", size: " + f.length() + " bytes)");
            try (BufferedReader br = new BufferedReader(new FileReader(f))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith("#CHROM")) {
                        String[] cols = line.split("\t");
                        List<String> samples = new ArrayList<>();
                        for (int j = 9; j < cols.length; j++) {
                            samples.add(cols[j]);
                            if (!allSamples.contains(cols[j])) {
                                allSamples.add(cols[j]);
                            }
                        }
                        fileToSamples.put(i, samples);
                        System.out.println("[Vcf-Merge] File: " + new File(inputFiles.get(i)).getName() + " -> Sample: " + samples);
                        break;
                    }
                }
            }
        }
        Collections.sort(allSamples);
        System.out.println("[Vcf-Merge] Total unique samples found: " + allSamples.size());

        // 2. Setup Multi-way stream
        PriorityQueue<VcfLine> pq = new PriorityQueue<>();
        for (int i = 0; i < inputFiles.size(); i++) {
            BufferedReader br = new BufferedReader(new FileReader(inputFiles.get(i)));
            readers.add(br);
            VcfLine firstLine = readNextGenomicLine(br, i);
            if (firstLine != null) pq.add(firstLine);
        }

        int totalSnpsBefore = 0;
        int totalSnpsAfter = 0;
        int tsCount = 0;
        int tvCount = 0;

        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)))) {
            // Write Header
            writer.println("##fileformat=VCFv4.2");
            writer.println("##source=BioJavaVcfMerger");
            writer.println("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
            writer.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternative Allele Frequency\">");
            writer.println("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
            writer.println("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
            writer.println("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">");
            writer.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            writer.println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
            writer.println("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depths\">");
            writer.println("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic Depths on Forward Strand\">");
            writer.println("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic Depths on Reverse Strand\">");
            writer.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
            writer.println("##FORMAT=<ID=DS,Number=1,Type=Integer,Description=\"Dosage of the alternative allele\">");
            writer.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            for (String sample : allSamples) writer.print("\t" + sample);
            writer.println();

            // 3. Process Streams
            while (!pq.isEmpty()) {
                List<VcfLine> currentSites = new ArrayList<>();
                VcfLine top = pq.poll();
                currentSites.add(top);

                // Pull all lines matching the same coordinate (Chr + Pos)
                while (!pq.isEmpty() && pq.peek().compareTo(top) == 0) {
                    currentSites.add(pq.poll());
                }

                totalSnpsBefore++;
                VcfLine lead = currentSites.get(0);

                // Fill genotypes for all samples
                List<String> mergedGenotypes = new ArrayList<>();
                int calledSamples = 0;
                int altAlleles = 0;
                int totalAlleles = 0;
                int totalDP = 0;

                for (String sample : allSamples) {
                    String genotype = "./.:0:0,0:0,0:0,0:0:0";
                    for (VcfLine site : currentSites) {
                        List<String> sourceSamples = fileToSamples.get(site.fileIndex);
                        int sampleIdx = sourceSamples.indexOf(sample);
                        if (sampleIdx != -1 && sampleIdx < site.genotypes.length) {
                            genotype = site.genotypes[sampleIdx];
                            break; 
                        }
                    }
                    mergedGenotypes.add(genotype);

                    if (!genotype.startsWith("./.")) {
                        calledSamples++;
                        String[] parts = genotype.split(":");
                        String gt = parts[0];
                        for (char c : gt.toCharArray()) {
                            if (c == '0') totalAlleles++;
                            else if (c == '1') {
                                totalAlleles++;
                                altAlleles++;
                            }
                        }
                    }
                    
                    String[] parts = genotype.split(":");
                    if (parts.length > 5) {
                        try {
                            totalDP += Integer.parseInt(parts[5]);
                        } catch (Exception ignored) {}
                    }
                }

                double callRate = (double) calledSamples / allSamples.size();
                double altFreq = totalAlleles > 0 ? (double) altAlleles / totalAlleles : 0.0;
                double mafValue = Math.min(altFreq, 1.0 - altFreq);

                // Apply Population Filters (MAF and Call Rate)
                if (mafValue < maf || callRate < minCallRate) {
                    // Refill PQ and continue
                    for (VcfLine site : currentSites) {
                        VcfLine next = readNextGenomicLine(readers.get(site.fileIndex), site.fileIndex);
                        if (next != null) pq.add(next);
                    }
                    continue;
                }

                totalSnpsAfter++;

                // Track Transitions/Transversions
                char r = Character.toUpperCase(lead.ref.charAt(0));
                char a = Character.toUpperCase(lead.alt.charAt(0));
                if ((r == 'A' && a == 'G') || (r == 'G' && a == 'A') || (r == 'C' && a == 'T') || (r == 'T' && a == 'C')) {
                    tsCount++;
                } else {
                    tvCount++;
                }

                // Update INFO field with population metrics
                String populationInfo = String.format(Locale.US, "DP=%d;AF=%.4f;AC=%d;AN=%d;MQ=60", totalDP, altFreq, altAlleles, totalAlleles);

                StringBuilder mergedLine = new StringBuilder();
                mergedLine.append(lead.chr).append("\t")
                          .append(lead.pos).append("\t")
                          .append(lead.id).append("\t")
                          .append(lead.ref).append("\t")
                          .append(lead.alt).append("\t")
                          .append(lead.qual).append("\t")
                          .append(lead.filter).append("\t")
                          .append(populationInfo).append("\t")
                          .append(lead.format);

                for (String gt : mergedGenotypes) {
                    mergedLine.append("\t").append(gt);
                }
                writer.println(mergedLine.toString());

                // Refill PQ from the source files we just used
                for (VcfLine site : currentSites) {
                    VcfLine next = readNextGenomicLine(readers.get(site.fileIndex), site.fileIndex);
                    if (next != null) pq.add(next);
                }
            }

            // Print beautiful population genetics diagnostics
            double tsTvRatio = tvCount > 0 ? (double) tsCount / tvCount : 0.0;
            System.out.println("\n=================================================");
            System.out.println("BioJava: Population Genetics Diagnostic Report");
            System.out.println("=================================================");
            System.out.printf("Total SNP Candidate Sites:  %d\n", totalSnpsBefore);
            System.out.printf("Passed Quality Filters:     %d (%.1f%%)\n", totalSnpsAfter, (totalSnpsBefore > 0 ? (totalSnpsAfter * 100.0 / totalSnpsBefore) : 0.0));
            System.out.printf("Filtered Out Sites:         %d\n", (totalSnpsBefore - totalSnpsAfter));
            System.out.println("-------------------------------------------------");
            System.out.printf("Transitions (Ts):           %d\n", tsCount);
            System.out.printf("Transversions (Tv):         %d\n", tvCount);
            System.out.printf("Ts/Tv Ratio:                %.2f %s\n", tsTvRatio, (tsTvRatio >= 1.5 && tsTvRatio <= 2.2 ? "✅ (Expected Genomic Range)" : "⚠️ (Sub-optimal, possible sequencing noise)"));
            System.out.println("=================================================");

        } finally {
            for (BufferedReader br : readers) br.close();
        }
    }

    private VcfLine readNextGenomicLine(BufferedReader br, int fileIndex) throws IOException {
        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) continue;
            return new VcfLine(line, fileIndex);
        }
        return null;
    }

    /**
     * Helper class to represent a VCF line in the PriorityQueue
     */
    private static class VcfLine implements Comparable<VcfLine> {
        String chr;
        int pos;
        String id, ref, alt, qual, filter, info, format;
        String[] genotypes;
        int fileIndex;

        public VcfLine(String line, int fileIndex) {
            String[] cols = line.split("\t");
            this.chr = cols[0];
            this.pos = Integer.parseInt(cols[1]);
            this.id = cols[2];
            this.ref = cols[3];
            this.alt = cols[4];
            this.qual = cols[5];
            this.filter = cols[6];
            this.info = cols[7];
            this.format = cols[8];
            this.fileIndex = fileIndex;
            this.genotypes = new String[cols.length - 9];
            System.arraycopy(cols, 9, this.genotypes, 0, cols.length - 9);
        }

        @Override
        public int compareTo(VcfLine other) {
            int c = this.chr.compareTo(other.chr);
            if (c != 0) return c;
            return Integer.compare(this.pos, other.pos);
        }
    }
}
