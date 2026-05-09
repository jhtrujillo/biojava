package org.cenicana.bio.core;

import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * High-performance Genetic Mapping Engine.
 * Builds genetic linkage maps directly from population VCF files or segregation matrices.
 */
public class GeneticMapEngine {

    private final double minLod;
    private final double maxRecomb;
    private final String mappingFunction;

    public GeneticMapEngine(double minLod, double maxRecomb, String mappingFunction) {
        this.minLod = minLod;
        this.maxRecomb = maxRecomb;
        this.mappingFunction = mappingFunction != null ? mappingFunction.toLowerCase() : "kosambi";
    }

    public static class Marker {
        public String id;
        public String chr;
        public long pos;
        public List<Double> dosages; // Numeric dosages (0.0 to ploidy) for all samples

        public Marker(String id, String chr, long pos) {
            this.id = id;
            this.chr = chr;
            this.pos = pos;
            this.dosages = new ArrayList<>();
        }
    }

    /**
     * Builds the genetic map from a VCF file and writes the results.
     */
    public void buildMap(String vcfPath, String mapOutputPath) throws IOException {
        System.out.println("[MapEngine] Reading variants from VCF: " + vcfPath);
        List<Marker> markers = parseVcf(vcfPath);
        if (markers.isEmpty()) {
            System.out.println("❌ Error: No markers loaded from VCF.");
            return;
        }
        System.out.println("[MapEngine] Loaded " + markers.size() + " markers across population.");

        // 1. Calculate Pairwise Recombination Frequencies and LOD scores
        int numMarkers = markers.size();
        double[][] rMatrix = new double[numMarkers][numMarkers];
        double[][] lodMatrix = new double[numMarkers][numMarkers];

        System.out.println("[MapEngine] Computing pairwise linkage (LOD and Recombination frequencies)...");
        for (int i = 0; i < numMarkers; i++) {
            rMatrix[i][i] = 0.0;
            lodMatrix[i][i] = 100.0;
            for (int j = i + 1; j < numMarkers; j++) {
                double[] stats = calculatePairwiseLinkage(markers.get(i), markers.get(j));
                rMatrix[i][j] = stats[0];
                rMatrix[j][i] = stats[0];
                lodMatrix[i][j] = stats[1];
                lodMatrix[j][i] = stats[1];
            }
        }

        // 2. Group markers into Linkage Groups using Hierarchical Clustering
        System.out.println("[MapEngine] Grouping markers into Linkage Groups (Chromosomes)...");
        int[] groups = partitionIntoLinkageGroups(numMarkers, lodMatrix, rMatrix);
        Map<Integer, List<Integer>> lgToMarkers = new HashMap<>();
        for (int i = 0; i < numMarkers; i++) {
            lgToMarkers.computeIfAbsent(groups[i], g -> new ArrayList<>()).add(i);
        }
        System.out.println("[MapEngine] Partition completed. Found " + lgToMarkers.size() + " Linkage Groups.");

        // 3. Order and Position markers within each Linkage Group
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(mapOutputPath)))) {
            pw.println("Marker\tLinkageGroup\tPosition(cM)\tChr_Phys\tPos_Phys");

            int lgIndex = 1;
            List<Integer> sortedLgKeys = new ArrayList<>(lgToMarkers.keySet());
            Collections.sort(sortedLgKeys);

            for (int lgId : sortedLgKeys) {
                List<Integer> markerIndices = lgToMarkers.get(lgId);
                if (markerIndices.size() < 2) {
                    // Single marker linkage group
                    Marker m = markers.get(markerIndices.get(0));
                    pw.printf(Locale.US, "%s\tLG%d\t0.00\t%s\t%d\n", m.id, lgIndex, m.chr, m.pos);
                    lgIndex++;
                    continue;
                }

                System.out.println("[MapEngine] Ordering LG" + lgIndex + " containing " + markerIndices.size() + " markers...");
                List<Integer> orderedPath = orderMarkersTSP(markerIndices, rMatrix);

                // Calculate cM distances
                double currentPositionCm = 0.0;
                Marker firstMarker = markers.get(orderedPath.get(0));
                pw.printf(Locale.US, "%s\tLG%d\t0.00\t%s\t%d\n", firstMarker.id, lgIndex, firstMarker.chr, firstMarker.pos);

                for (int i = 1; i < orderedPath.size(); i++) {
                    int m1Idx = orderedPath.get(i - 1);
                    int m2Idx = orderedPath.get(i);
                    double r = rMatrix[m1Idx][m2Idx];
                    double distanceCm = convertRecombToCm(r);
                    currentPositionCm += distanceCm;

                    Marker m = markers.get(m2Idx);
                    pw.printf(Locale.US, "%s\tLG%d\t%.2f\t%s\t%d\n", m.id, lgIndex, currentPositionCm, m.chr, m.pos);
                }
                lgIndex++;
            }
        }
        System.out.println("🎉 [MapEngine] Genetic linkage map successfully written to: " + mapOutputPath);
    }

    /**
     * Parses a VCF file to extract markers and dosage vectors.
     */
    private List<Marker> parseVcf(String vcfPath) throws IOException {
        List<Marker> markers = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(vcfPath))) {
            String line;
            int gtIndex = 0;
            int dsIndex = -1;

            while ((line = br.readLine()) != null) {
                if (line.startsWith("##")) continue;
                if (line.startsWith("#")) {
                    // Header line - check sample count
                    continue;
                }

                String[] cols = line.split("\t");
                if (cols.length < 10) continue;

                String chr = cols[0];
                long pos = Long.parseLong(cols[1]);
                String id = cols[2];
                if (id.equals(".")) {
                    id = chr + "_" + pos;
                }

                String format = cols[8];
                String[] formatFields = format.split(":");
                for (int f = 0; f < formatFields.length; f++) {
                    if (formatFields[f].equalsIgnoreCase("DS")) {
                        dsIndex = f;
                        break;
                    }
                }

                Marker marker = new Marker(id, chr, pos);
                for (int i = 9; i < cols.length; i++) {
                    String sampleData = cols[i];
                    if (sampleData.startsWith("./.")) {
                        marker.dosages.add(-1.0); // Missing representation
                        continue;
                    }

                    String[] values = sampleData.split(":");
                    double dosage = 0.0;

                    if (dsIndex != -1 && values.length > dsIndex) {
                        try {
                            dosage = Double.parseDouble(values[dsIndex]);
                        } catch (NumberFormatException e) {
                            dosage = parseGenotypeDosage(values[0]);
                        }
                    } else {
                        dosage = parseGenotypeDosage(values[0]);
                    }
                    marker.dosages.add(dosage);
                }
                markers.add(marker);
            }
        }
        return markers;
    }

    private double parseGenotypeDosage(String gt) {
        double d = 0.0;
        for (char c : gt.toCharArray()) {
            if (c == '1') d += 1.0;
        }
        return d;
    }

    /**
     * Computes recombination frequency (r) and LOD score between two markers using Pearson correlation.
     */
    private double[] calculatePairwiseLinkage(Marker m1, Marker m2) {
        List<Double> d1 = new ArrayList<>();
        List<Double> d2 = new ArrayList<>();

        for (int k = 0; k < m1.dosages.size(); k++) {
            double v1 = m1.dosages.get(k);
            double v2 = m2.dosages.get(k);
            if (v1 != -1.0 && v2 != -1.0) {
                d1.add(v1);
                d2.add(v2);
            }
        }

        int n = d1.size();
        if (n < 5) {
            return new double[]{0.5, 0.0}; // No linkage information
        }

        double mean1 = 0.0;
        double mean2 = 0.0;
        for (int i = 0; i < n; i++) {
            mean1 += d1.get(i);
            mean2 += d2.get(i);
        }
        mean1 /= n;
        mean2 /= n;

        double num = 0.0;
        double den1 = 0.0;
        double den2 = 0.0;

        for (int i = 0; i < n; i++) {
            double diff1 = d1.get(i) - mean1;
            double diff2 = d2.get(i) - mean2;
            num += diff1 * diff2;
            den1 += diff1 * diff1;
            den2 += diff2 * diff2;
        }

        if (den1 == 0.0 || den2 == 0.0) {
            return new double[]{0.5, 0.0};
        }

        double rCoef = num / Math.sqrt(den1 * den2);
        double rAbs = Math.abs(rCoef);

        // Recombination frequency calculation
        double recomb = 0.5 * (1.0 - rAbs);
        recomb = Math.max(0.0001, Math.min(0.5, recomb));

        // LOD score estimation based on correlation significance
        double rSq = Math.min(0.999, rCoef * rCoef);
        double lod = -0.5 * n * Math.log10(1.0 - rSq);
        lod = Math.max(0.0, lod);

        return new double[]{recomb, lod};
    }

    /**
     * Single-Linkage Hierarchical Clustering to partition markers into Linkage Groups.
     */
    private int[] partitionIntoLinkageGroups(int numMarkers, double[][] lodMatrix, double[][] rMatrix) {
        int[] parent = new int[numMarkers];
        for (int i = 0; i < numMarkers; i++) parent[i] = i;

        for (int i = 0; i < numMarkers; i++) {
            for (int j = i + 1; j < numMarkers; j++) {
                if (lodMatrix[i][j] >= minLod && rMatrix[i][j] <= maxRecomb) {
                    union(i, j, parent);
                }
            }
        }

        // Compress paths
        int[] groups = new int[numMarkers];
        for (int i = 0; i < numMarkers; i++) {
            groups[i] = find(i, parent);
        }
        return groups;
    }

    private int find(int i, int[] parent) {
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i], parent);
    }

    private void union(int i, int j, int[] parent) {
        int rootI = find(i, parent);
        int rootJ = find(j, parent);
        if (rootI != rootJ) {
            parent[rootI] = rootJ;
        }
    }

    /**
     * TSP 2-Opt Heuristic to find optimal linear ordering of markers in a linkage group.
     */
    private List<Integer> orderMarkersTSP(List<Integer> markers, double[][] rMatrix) {
        List<Integer> path = new ArrayList<>(markers);

        // Nearest-Neighbor Initial Order
        List<Integer> initialPath = new ArrayList<>();
        Set<Integer> unvisited = new HashSet<>(markers);
        int current = path.get(0);
        initialPath.add(current);
        unvisited.remove(current);

        while (!unvisited.isEmpty()) {
            int next = -1;
            double minDist = Double.MAX_VALUE;
            for (int cand : unvisited) {
                if (rMatrix[current][cand] < minDist) {
                    minDist = rMatrix[current][cand];
                    next = cand;
                }
            }
            initialPath.add(next);
            unvisited.remove(next);
            current = next;
        }
        path = initialPath;

        // Apply 2-Opt optimization
        boolean improved = true;
        int size = path.size();
        while (improved) {
            improved = false;
            double bestLength = getPathLength(path, rMatrix);
            for (int i = 1; i < size - 1; i++) {
                for (int j = i + 1; j < size; j++) {
                    List<Integer> newPath = opt2Swap(path, i, j);
                    double newLength = getPathLength(newPath, rMatrix);
                    if (newLength < bestLength) {
                        path = newPath;
                        bestLength = newLength;
                        improved = true;
                    }
                }
            }
        }
        return path;
    }

    private List<Integer> opt2Swap(List<Integer> path, int i, int j) {
        List<Integer> newPath = new ArrayList<>(path.subList(0, i));
        List<Integer> sub = new ArrayList<>(path.subList(i, j + 1));
        Collections.reverse(sub);
        newPath.addAll(sub);
        newPath.addAll(path.subList(j + 1, path.size()));
        return newPath;
    }

    private double getPathLength(List<Integer> path, double[][] rMatrix) {
        double length = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            length += rMatrix[path.get(i)][path.get(i + 1)];
        }
        return length;
    }

    /**
     * Converts recombination frequency (r) into genetic distance (cM) using Kosambi or Haldane.
     */
    private double convertRecombToCm(double r) {
        r = Math.max(0.0001, Math.min(0.4999, r));
        if (mappingFunction.equalsIgnoreCase("haldane")) {
            return -50.0 * Math.log(1.0 - 2.0 * r);
        } else {
            // Kosambi
            return 25.0 * Math.log((1.0 + 2.0 * r) / (1.0 - 2.0 * r));
        }
    }
}
