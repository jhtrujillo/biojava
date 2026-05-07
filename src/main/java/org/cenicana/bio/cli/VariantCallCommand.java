package org.cenicana.bio.cli;

import org.cenicana.bio.core.VariantCaller;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.concurrent.Callable;
import java.io.File;

@Command(name = "call-variants",
         description = "Call SNPs and small Indels from SAM/BAM files using a parallelized pileup algorithm.",
         mixinStandardHelpOptions = true,
         version = "call-variants 1.0")
public class VariantCallCommand implements Callable<Integer> {

    @Option(names = {"-i", "--input"}, description = "Path to the input SAM or BAM file.", required = true)
    private String inputFile;

    @Option(names = {"-r", "--reference"}, description = "Path to the reference genome FASTA file.", required = false)
    private String referenceFile;

    @Option(names = {"-o", "--output"}, description = "Path to the output VCF file.", required = true)
    private String outputFile;

    @Option(names = {"-p", "--ploidy"}, description = "Ploidy level of the organism (e.g. 2 for diploid, 4 for tetraploid, 10 for sugarcane).", defaultValue = "2")
    private int ploidy;

    @Option(names = {"-t", "--threads"}, description = "Number of processor threads to use for parallel calling.", defaultValue = "4")
    private int threads;

    @Option(names = {"--min-depth"}, description = "Minimum coverage depth required to call a variant.", defaultValue = "10")
    private int minDepth;

    @Option(names = {"--min-mapq"}, description = "Minimum mapping quality (MAPQ) for reads to be included.", defaultValue = "30")
    private int minMapq;

    @Option(names = {"--min-qual"}, description = "Minimum base quality (Phred score) to include a base.", defaultValue = "20")
    private int minQual;

    @Option(names = {"--min-af"}, description = "Minimum alternative allele frequency to call a variant.", defaultValue = "0.1")
    private double minAf;

    @Option(names = {"--samtools"}, description = "Custom path to the samtools executable.", defaultValue = "samtools")
    private String samtoolsPath;

    @Option(names = {"--preset"}, description = "Variant calling model preset (fast, ngsep, freebayes, gatk)", defaultValue = "fast")
    private String preset;

    @Override
    public Integer call() throws Exception {
        System.out.println("=================================================");
        System.out.println("BioJava: Parallel SNP & Indel Variant Caller");
        System.out.println("=================================================");
        System.out.println("Input File:     " + inputFile);
        System.out.println("Ref Genome:     " + (referenceFile != null ? referenceFile : "None (Inferred REF)"));
        System.out.println("Output VCF:     " + outputFile);
        System.out.println("Ploidy:         " + ploidy);
        System.out.println("Preset Model:   " + preset.toUpperCase());
        System.out.println("Threads:        " + threads);
        System.out.println("Filters:        MinDepth > " + minDepth + ", MinMAPQ > " + minMapq + ", MinBaseQual > " + minQual + ", MinAF > " + minAf);
        System.out.println("=================================================\n");

        long startTime = System.currentTimeMillis();

        VariantCaller caller = new VariantCaller(samtoolsPath, minDepth, minMapq, minQual, minAf);
        File inPath = new File(inputFile);

        if (inPath.isDirectory()) {
            System.out.println("[Batch] Detected directory input. Activating Population/Batch Mode...");
            File[] files = inPath.listFiles((dir, name) -> name.toLowerCase().endsWith(".bam") || name.toLowerCase().endsWith(".sam"));
            if (files == null || files.length == 0) {
                System.out.println("❌ Error: No SAM/BAM files found in directory: " + inputFile);
                return 1;
            }

            System.out.println("[Batch] Found " + files.length + " alignment files. Processing sequentially...");
            java.util.List<String> tempVcfs = new java.util.ArrayList<>();
            java.io.File tempDir = java.nio.file.Files.createTempDirectory("biojava_batch_vcfs").toFile();
            tempDir.deleteOnExit();

            for (int i = 0; i < files.length; i++) {
                File bamFile = files[i];
                System.out.println("\n[Batch] (" + (i+1) + "/" + files.length + ") Processing sample: " + bamFile.getName());
                
                String tempVcfPath = new File(tempDir, bamFile.getName() + ".vcf").getAbsolutePath();
                try {
                    caller.callVariants(bamFile.getAbsolutePath(), referenceFile, ploidy, threads, preset, tempVcfPath);
                    tempVcfs.add(tempVcfPath);
                } catch (Exception e) {
                    System.out.println("⚠️ Warning: Failed to process " + bamFile.getName() + ": " + e.getMessage());
                }
            }

            System.out.println("\n[Batch] Variant calling complete for all samples. Merging into unified Population VCF...");
            org.cenicana.bio.core.VcfMerger merger = new org.cenicana.bio.core.VcfMerger(tempVcfs, outputFile);
            merger.merge();

        } else {
            // Single file mode
            caller.callVariants(inputFile, referenceFile, ploidy, threads, preset, outputFile);
        }

        long endTime = System.currentTimeMillis();
        double elapsedSec = (endTime - startTime) / 1000.0;
        System.out.printf("\n🎉 Variant calling successfully completed in %.2f seconds! Output: %s\n", elapsedSec, outputFile);

        return 0;
    }
}
