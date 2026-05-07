package org.cenicana.bio.core;

import org.cenicana.bio.io.FastaReader;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * High-performance, multi-threaded Variant Caller for SNPs and Small Indels from SAM/BAM files.
 * Supports polyploid genotyping and quality filtering.
 */
public class VariantCaller {

    public static class PileupPosition {
        public String chrom;
        public long position;
        public int depth = 0;
        public int aCount = 0;
        public int tCount = 0;
        public int cCount = 0;
        public int gCount = 0;
        public int forwardACount = 0, reverseACount = 0;
        public int forwardTCount = 0, reverseTCount = 0;
        public int forwardCCount = 0, reverseCCount = 0;
        public int forwardGCount = 0, reverseGCount = 0;
        public int insertCount = 0;
        public int deleteCount = 0;
        
        // Fields for advanced models (FreeBayes / GATK)
        public int sumReadPosA = 0, sumReadPosT = 0, sumReadPosC = 0, sumReadPosG = 0;
        public int sumBaseQualA = 0, sumBaseQualT = 0, sumBaseQualC = 0, sumBaseQualG = 0;
        
        public List<Byte> baseQualities = new ArrayList<>();

        public PileupPosition(String chrom, long position) {
            this.chrom = chrom;
            this.position = position;
        }

        public boolean hasStrandSupport(char base) {
            switch (Character.toUpperCase(base)) {
                case 'A': return forwardACount > 0 && reverseACount > 0;
                case 'T': return forwardTCount > 0 && reverseTCount > 0;
                case 'C': return forwardCCount > 0 && reverseCCount > 0;
                case 'G': return forwardGCount > 0 && reverseGCount > 0;
                default: return false;
            }
        }

        public int getAltCount(char ref) {
            int refC = getRefC(ref);
            int total = aCount + tCount + cCount + gCount;
            return total - refC;
        }

        public char getMajorAlt(char ref) {
            char best = 'N';
            int max = 0;
            char r = Character.toUpperCase(ref);
            if (r != 'A' && aCount > max) { max = aCount; best = 'A'; }
            if (r != 'T' && tCount > max) { max = tCount; best = 'T'; }
            if (r != 'C' && cCount > max) { max = cCount; best = 'C'; }
            if (r != 'G' && gCount > max) { max = gCount; best = 'G'; }
            return best;
        }

        public int getRefC(char ref) {
            switch (Character.toUpperCase(ref)) {
                case 'A': return aCount;
                case 'T': return tCount;
                case 'C': return cCount;
                case 'G': return gCount;
                default: return 0;
            }
        }

        public double getAvgReadPos(char base) {
            int count = getRefC(base);
            if (count == 0) return 0;
            switch (Character.toUpperCase(base)) {
                case 'A': return (double) sumReadPosA / count;
                case 'T': return (double) sumReadPosT / count;
                case 'C': return (double) sumReadPosC / count;
                case 'G': return (double) sumReadPosG / count;
                default: return 0;
            }
        }

        public double getAvgBaseQual(char base) {
            int count = getRefC(base);
            if (count == 0) return 0;
            switch (Character.toUpperCase(base)) {
                case 'A': return (double) sumBaseQualA / count;
                case 'T': return (double) sumBaseQualT / count;
                case 'C': return (double) sumBaseQualC / count;
                case 'G': return (double) sumBaseQualG / count;
                default: return 0;
            }
        }
    }

    private final String samtoolsPath;
    private final int minDepth;
    private final int minMapQ;
    private final int minBaseQual;
    private final double minAltFreq;

    public VariantCaller(String samtoolsPath, int minDepth, int minMapQ, int minBaseQual, double minAltFreq) {
        this.samtoolsPath = samtoolsPath != null ? samtoolsPath : "samtools";
        this.minDepth = minDepth;
        this.minMapQ = minMapQ;
        this.minBaseQual = minBaseQual;
        this.minAltFreq = minAltFreq;
    }

    /**
     * Executes the parallel variant calling pipeline.
     */
    public void callVariants(String alignmentFile, String refFasta, int ploidy, int threads, String preset, String outVcf) throws Exception {
        boolean isBam = alignmentFile.toLowerCase().endsWith(".bam");
        
        System.out.println("[Caller] Reading reference genome metadata...");
        Map<String, String> refSequences = new HashMap<>();
        if (refFasta != null && Files.exists(Paths.get(refFasta))) {
            FastaReader reader = new FastaReader();
            refSequences = reader.read(refFasta);
            System.out.println("[Caller] Loaded " + refSequences.size() + " reference sequences.");
        } else {
            System.out.println("[Caller] WARNING: No reference FASTA provided. Reference alleles will be inferred.");
        }

        // 1. Get chromosome list from header
        List<String> chromosomes = getChromosomesFromHeader(alignmentFile, isBam);
        if (chromosomes.isEmpty()) {
            throw new IOException("No chromosomes/contigs found in the alignment file header (@SQ).");
        }
        System.out.println("[Caller] Chromosomes detected for parallelization: " + chromosomes);

        // 2. Setup Executor Service
        int numWorkers = Math.min(threads, chromosomes.size());
        System.out.println("[Caller] Spawning thread pool with " + numWorkers + " workers...");
        ExecutorService executor = Executors.newFixedThreadPool(numWorkers);
        List<Future<File>> futures = new ArrayList<>();

        // Create temp directory for VCF chunks
        File tempDir = Files.createTempDirectory("biojava_vcf_chunks").toFile();
        tempDir.deleteOnExit();

        for (String chrom : chromosomes) {
            final Map<String, String> finalRefSeqs = refSequences;
            futures.add(executor.submit(() -> {
                File chunkFile = new File(tempDir, chrom + ".vcf.tmp");
                chunkFile.deleteOnExit();
                
                String refSeq = finalRefSeqs.get(chrom);
                processChromosome(alignmentFile, chrom, refSeq, ploidy, preset, isBam, chunkFile);
                return chunkFile;
            }));
        }

        // Wait for all tasks to complete
        List<File> chunkFiles = new ArrayList<>();
        for (Future<File> future : futures) {
            try {
                chunkFiles.add(future.get());
            } catch (ExecutionException e) {
                executor.shutdownNow();
                throw new RuntimeException("Error during parallel variant calling task", e.getCause());
            }
        }
        executor.shutdown();

        // 3. Merge VCF chunks into final VCF
        System.out.println("[Caller] Merging chromosome chunks into: " + outVcf);
        mergeVcfFiles(chunkFiles, alignmentFile, outVcf);
        System.out.println("🎉 [Caller] Variant calling complete! Output written to: " + outVcf);
    }

    /**
     * Extracts chromosomes from SAM/BAM header (@SQ lines).
     */
    private List<String> getChromosomesFromHeader(String alignmentFile, boolean isBam) throws IOException, InterruptedException {
        List<String> chromosomes = new ArrayList<>();
        
        List<String> cmd = new ArrayList<>();
        if (isBam) {
            cmd.add(samtoolsPath);
            cmd.add("view");
            cmd.add("-H");
            cmd.add(alignmentFile);
        } else {
            // For SAM files, we just read the first few lines
            try (BufferedReader br = new BufferedReader(new FileReader(alignmentFile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (!line.startsWith("@")) break;
                    if (line.startsWith("@SQ")) {
                        String chrom = parseHeaderSqLine(line);
                        if (chrom != null) chromosomes.add(chrom);
                    }
                }
            }
            return chromosomes;
        }

        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectError(ProcessBuilder.Redirect.DISCARD);
        Process process = pb.start();
        
        try (BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("@SQ")) {
                    String chrom = parseHeaderSqLine(line);
                    if (chrom != null) chromosomes.add(chrom);
                }
            }
        }
        process.waitFor();
        return chromosomes;
    }

    private String parseHeaderSqLine(String line) {
        String[] parts = line.split("\t");
        for (String part : parts) {
            if (part.startsWith("SN:")) {
                return part.substring(3);
            }
        }
        return null;
    }

    /**
     * Processes a single chromosome, performing pileup and variant calling.
     */
    private void processChromosome(String file, String chrom, String refSeq, int ploidy, String preset, boolean isBam, File outFile) throws Exception {
        System.out.println("[Caller] [" + chrom + "] Running pileup...");
        TreeMap<Long, PileupPosition> pileup = new TreeMap<>();

        BufferedReader reader = null;
        Process process = null;

        if (isBam) {
            List<String> cmd = new ArrayList<>();
            cmd.add(samtoolsPath);
            cmd.add("view");
            cmd.add(file);
            cmd.add(chrom);

            ProcessBuilder pb = new ProcessBuilder(cmd);
            pb.redirectError(ProcessBuilder.Redirect.DISCARD);
            process = pb.start();
            reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        } else {
            reader = new BufferedReader(new FileReader(file));
        }

        try {
            String line;
            Pattern cigarPattern = Pattern.compile("(\\d+)([MIDNSHP=X])");

            while ((line = reader.readLine()) != null) {
                if (line.startsWith("@")) continue;
                String[] cols = line.split("\t");
                if (cols.length < 11) continue;

                String readChrom = cols[2];
                if (!readChrom.equals(chrom)) continue;

                int flag = Integer.parseInt(cols[1]);
                if ((flag & 4) != 0) continue; // Unmapped read
                if ((flag & 1024) != 0) continue; // Duplicate
                boolean isReverseStrand = (flag & 16) != 0;

                int mapq = Integer.parseInt(cols[4]);
                if (mapq < minMapQ) continue;

                long startPos = Long.parseLong(cols[3]);
                String cigar = cols[5];
                String seq = cols[9];
                String qualStr = cols[10];

                Matcher matcher = cigarPattern.matcher(cigar);
                int readIdx = 0;
                long refPos = startPos;

                while (matcher.find()) {
                    int len = Integer.parseInt(matcher.group(1));
                    char op = matcher.group(2).charAt(0);

                    switch (op) {
                        case 'M':
                        case '=':
                        case 'X':
                            for (int k = 0; k < len; k++) {
                                if (readIdx < seq.length()) {
                                    char base = seq.charAt(readIdx);
                                    byte qual = (byte) (qualStr.charAt(readIdx) - 33);
                                    if (qual >= minBaseQual) {
                                        PileupPosition p = pileup.computeIfAbsent(refPos, pos -> new PileupPosition(chrom, pos));
                                        p.depth++;
                                        p.baseQualities.add(qual);
                                        
                                        int distToEnd = Math.min(readIdx, seq.length() - readIdx - 1);

                                        switch (Character.toUpperCase(base)) {
                                            case 'A':
                                                p.aCount++;
                                                p.sumBaseQualA += qual;
                                                p.sumReadPosA += distToEnd;
                                                if (isReverseStrand) p.reverseACount++; else p.forwardACount++;
                                                break;
                                            case 'T':
                                                p.tCount++;
                                                p.sumBaseQualT += qual;
                                                p.sumReadPosT += distToEnd;
                                                if (isReverseStrand) p.reverseTCount++; else p.forwardTCount++;
                                                break;
                                            case 'C':
                                                p.cCount++;
                                                p.sumBaseQualC += qual;
                                                p.sumReadPosC += distToEnd;
                                                if (isReverseStrand) p.reverseCCount++; else p.forwardCCount++;
                                                break;
                                            case 'G':
                                                p.gCount++;
                                                p.sumBaseQualG += qual;
                                                p.sumReadPosG += distToEnd;
                                                if (isReverseStrand) p.reverseGCount++; else p.forwardGCount++;
                                                break;
                                        }
                                    }
                                }
                                readIdx++;
                                refPos++;
                            }
                            break;
                        case 'I':
                            if (readIdx < seq.length()) {
                                PileupPosition p = pileup.computeIfAbsent(refPos - 1, pos -> new PileupPosition(chrom, pos));
                                p.insertCount++;
                                readIdx += len;
                            }
                            break;
                        case 'D':
                            for (int k = 0; k < len; k++) {
                                PileupPosition p = pileup.computeIfAbsent(refPos, pos -> new PileupPosition(chrom, pos));
                                p.deleteCount++;
                                p.depth++;
                                refPos++;
                            }
                            break;
                        case 'S':
                            readIdx += len;
                            break;
                        case 'N':
                            refPos += len;
                            break;
                    }
                }
            }
        } finally {
            if (reader != null) reader.close();
            if (process != null) process.destroy();
        }

        // 2. Perform Variant Calling and Genotyping
        boolean isNgsep = preset.equalsIgnoreCase("ngsep") || preset.equalsIgnoreCase("freebayes") || preset.equalsIgnoreCase("gatk");
        boolean isFreebayes = preset.equalsIgnoreCase("freebayes") || preset.equalsIgnoreCase("gatk");
        boolean isGatk = preset.equalsIgnoreCase("gatk");
        boolean dynamicPloidy = isNgsep;

        System.out.println("[Caller] [" + chrom + "] Generating variant calls using preset: " + preset.toUpperCase());
        double avgDepth = 1.0;
        if (dynamicPloidy) {
            double sumDepth = 0;
            int countPositions = 0;
            for (PileupPosition p : pileup.values()) {
                if (p.depth >= minDepth) {
                    sumDepth += p.depth;
                    countPositions++;
                }
            }
            avgDepth = countPositions > 0 ? (sumDepth / countPositions) : 1.0;
            System.out.printf(Locale.US, "[Caller] [%s] Calculated average local coverage: %.2fx\n", chrom, avgDepth);
        }

        try (PrintWriter pw = new PrintWriter(new FileWriter(outFile))) {
            for (Map.Entry<Long, PileupPosition> entry : pileup.entrySet()) {
                long pos = entry.getKey();
                PileupPosition p = entry.getValue();

                if (p.depth < minDepth) continue;

                // Skip over-saturated repetitive regions (depth > 3.5x chromosome average)
                if (dynamicPloidy && p.depth > 3.5 * avgDepth) {
                    continue;
                }

                // Determine Reference Allele
                char ref = 'N';
                if (refSeq != null && pos <= refSeq.length()) {
                    ref = Character.toUpperCase(refSeq.charAt((int) pos - 1));
                } else {
                    // Infer reference as the major allele
                    ref = p.getMajorAlt(' ');
                }

                if (ref == 'N') continue;

                int altCount = p.getAltCount(ref);
                double altFreq = (double) altCount / p.depth;

                if (altCount >= 2 && altFreq >= minAltFreq) {
                    char alt = p.getMajorAlt(ref);
                    if (alt == 'N') continue;

                    // Strand-bias filter: require alternative allele to be observed in both forward and reverse reads
                    if (!p.hasStrandSupport(alt)) {
                        continue;
                    }

                    // Preset-specific filters
                    if (isFreebayes) {
                        // Read Position Bias Filter: Discard if alt alleles are clustered at the very ends of reads
                        double avgPos = p.getAvgReadPos(alt);
                        if (avgPos < 5.0) continue;
                    }

                    if (isGatk) {
                        // Strict Base Quality Filter: Require high average base quality for the alt allele
                        if (p.getAvgBaseQual(alt) < minBaseQual + 10) continue;
                    }

                    // Genotyping dosage estimation for polyploids
                    int localPloidy = ploidy;
                    if (dynamicPloidy) {
                        localPloidy = (int) Math.round((double) p.depth * ploidy / avgDepth);
                        localPloidy = Math.max(1, localPloidy);
                    }
                    int dosage = (int) Math.round(altFreq * localPloidy);
                    dosage = Math.max(0, Math.min(localPloidy, dosage));
                    String gt = formatGenotype(dosage, localPloidy);

                    // Bayesian Genotype Quality (GQ) Approximation for FreeBayes/GATK
                    int gq = 30; // Default heuristic
                    if (isFreebayes) {
                        double avgQual = p.getAvgBaseQual(alt);
                        // Approximate Phred scale GQ based on depth and average base quality
                        gq = (int) Math.min(99, (altCount * (avgQual / 10.0)));
                    }

                    // Write VCF record with GQ and DS format fields
                    pw.printf(Locale.US, "%s\t%d\t.\t%c\t%c\t%d\tPASS\tDP=%d;AF=%.4f\tGT:GQ:AD:DP:DS\t%s:%d:%d,%d:%d:%d\n",
                        chrom, pos, ref, alt, gq, p.depth, altFreq, gt, gq, p.getRefC(ref), altCount, p.depth, dosage);
                }
            }
        }
    }

    /**
     * Formats Genotype (GT) using the dosage level for the given ploidy.
     */
    private String formatGenotype(int dosage, int ploidy) {
        List<String> alleles = new ArrayList<>();
        for (int i = 0; i < ploidy - dosage; i++) alleles.add("0");
        for (int i = 0; i < dosage; i++) alleles.add("1");
        
        return String.join("/", alleles);
    }

    /**
     * Merges independent VCF chunks into a single final VCF file.
     */
    private void mergeVcfFiles(List<File> chunkFiles, String sourceBam, String outVcf) throws IOException {
        try (PrintWriter pw = new PrintWriter(new FileWriter(outVcf))) {
            // Write standard VCF Header
            pw.println("##fileformat=VCFv4.2");
            pw.println("##fileDate=" + new java.util.Date());
            pw.println("##source=BioJavaVariantCaller_1.0");
            pw.println("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
            pw.println("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternative Allele Frequency\">");
            pw.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
            pw.println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
            pw.println("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depths\">");
            pw.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
            pw.println("##FORMAT=<ID=DS,Number=1,Type=Integer,Description=\"Dosage of the alternative allele\">");
            String sampleName = new java.io.File(sourceBam).getName();
            if (sampleName.toLowerCase().endsWith(".bam") || sampleName.toLowerCase().endsWith(".sam")) {
                sampleName = sampleName.substring(0, sampleName.length() - 4);
            }
            pw.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sampleName);

            // Append chunk records
            for (File chunk : chunkFiles) {
                if (!chunk.exists()) continue;
                try (BufferedReader br = new BufferedReader(new FileReader(chunk))) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        pw.println(line);
                    }
                }
            }
        }
    }
}
