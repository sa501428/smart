/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2022 Rice University, Baylor College of Medicine, Aiden Lab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package mixer.utils.kmeans;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.drive.FinalMatrix;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class ClusteringMagic {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    protected static final int maxIters = 500;
    protected final File outputDirectory;
    protected final Random generator = new Random(2352);
    protected final FinalMatrix matrix;
    protected final ChromosomeHandler handler;

    public ClusteringMagic(FinalMatrix matrix, File outputDirectory,
                           ChromosomeHandler handler, long seed) {
        this.matrix = matrix;
        this.handler = handler;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);
    }

    public static String getOutputName(String prefix, boolean useKMedians, int k) {
        String kstem = "kmeans";
        if (useKMedians) kstem = "kmedians";
        return prefix + "_" + kstem + "_k" + k + "_clusters.bed";
    }

    public void extractFinalGWSubcompartments(String prefix, Map<Integer, List<String>> bedFiles) {
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            int numClusters = z + startingClusterSizeK;
            if (!bedFiles.containsKey(numClusters)) {
                bedFiles.put(numClusters, new ArrayList<>(2));
            }
        }

        System.out.println("Genome-wide KMeans clustering");
        // todo matrix.inPlaceScaleSqrtWeightCol();
        runClusteringOnMatrix(prefix, false, bedFiles);

        System.out.println("Genome-wide KMedians clustering");
        // todo matrix.inPlaceScaleSqrtWeightCol();
        runClusteringOnMatrix(prefix, true, bedFiles);
    }

    private void runClusteringOnMatrix(String prefix, boolean useKMedians, Map<Integer, List<String>> outputs) {
        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(handler, matrix,
                false, useKMedians);
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            runKMeansMultipleTimes(kmeansRunner, z, useKMedians, prefix, outputs);
        }
        System.out.println(">");
    }

    private void runKMeansMultipleTimes(GenomeWideKmeansRunner kmeansRunner,
                                        int z, boolean useKMedians, String prefix, Map<Integer, List<String>> outputs) {
        int numClusters = z + startingClusterSizeK;
        double wcssLimit = Float.MAX_VALUE;
        GenomeWide1DList<SubcompartmentInterval> bestClusters = null;
        int[] bestAssignments = null;
        int attempts = 0;
        while (bestAssignments == null || attempts < 20) {
            KmeansResult currResult = resetAndRerun(kmeansRunner, generator, numClusters);
            double wcss = currResult.getWithinClusterSumOfSquares();
            if (currResult.getNumActualClusters() == numClusters && wcss < Float.MAX_VALUE) {
                if (wcss < wcssLimit) {
                    wcssLimit = wcss;
                    bestAssignments = currResult.getAssignments(matrix.getNumRows());
                    bestClusters = currResult.getFinalCompartmentsClone();
                }
                attempts++;
            }
        }
        exportKMeansClusteringResults(z, prefix, useKMedians, bestClusters, outputs);
    }

    protected KmeansResult resetAndRerun(GenomeWideKmeansRunner kmeansRunner, Random generator, int numClusters) {
        kmeansRunner.prepareForNewRun(numClusters);
        kmeansRunner.launchKmeansGWMatrix(generator.nextLong(), maxIters);
        return kmeansRunner.getResult();
    }

    public void exportKMeansClusteringResults(int z, String prefix, boolean useKMedians,
                                              GenomeWide1DList<SubcompartmentInterval> finalCompartments,
                                              Map<Integer, List<String>> outputs) {
        int k = z + startingClusterSizeK;
        SliceUtils.collapseGWList(finalCompartments);
        File outBedFile = new File(outputDirectory, getOutputName(prefix, useKMedians, k)); // "_wcss" + wcss +
        if (outputs != null) outputs.get(k).add(outBedFile.getAbsolutePath());
        finalCompartments.simpleExport(outBedFile);
    }
}
