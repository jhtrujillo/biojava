package org.cenicana.bio.cli;

import org.cenicana.bio.core.GeneticMapEngine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import java.util.concurrent.Callable;

@Command(name = "genetic-map",
         description = "Build genetic linkage maps (Linkage Groups & Marker Ordering) directly from VCF files.",
         mixinStandardHelpOptions = true,
         version = "genetic-map 1.0")
public class GeneticMapCommand implements Callable<Integer> {

    @Option(names = {"-i", "--input"}, description = "Path to the input VCF file containing population genotypes.", required = true)
    private String inputFile;

    @Option(names = {"-o", "--output"}, description = "Path to the output .map file to write genetic coordinates.", required = true)
    private String outputFile;

    @Option(names = {"--lod"}, description = "Minimum LOD score threshold to group markers into linkage groups.", defaultValue = "3.0")
    private double lod;

    @Option(names = {"--max-r"}, description = "Maximum recombination frequency (r) for linkage grouping.", defaultValue = "0.35")
    private double maxR;

    @Option(names = {"--mapping-function"}, description = "Mapping function to convert recombination to cM (kosambi, haldane).", defaultValue = "kosambi")
    private String mappingFunction;

    @Override
    public Integer call() throws Exception {
        System.out.println("=================================================");
        System.out.println("BioJava: High-Performance Genetic Linkage Mapper");
        System.out.println("=================================================");
        System.out.println("Input VCF:          " + inputFile);
        System.out.println("Output Map:         " + outputFile);
        System.out.println("Min LOD Grouping:   " + lod);
        System.out.println("Max Recomb Limit:   " + maxR);
        System.out.println("Mapping Function:   " + mappingFunction.toUpperCase());
        System.out.println("=================================================\n");

        long startTime = System.currentTimeMillis();

        GeneticMapEngine engine = new GeneticMapEngine(lod, maxR, mappingFunction);
        engine.buildMap(inputFile, outputFile);

        long endTime = System.currentTimeMillis();
        double elapsedSec = (endTime - startTime) / 1000.0;
        System.out.printf("\n🎉 Genetic map successfully built in %.2f seconds!\n", elapsedSec);

        return 0;
    }
}
