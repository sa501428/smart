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
import mixer.utils.drive.MatrixAndWeight;
import mixer.utils.tracks.Concensus3DTools;
import mixer.utils.tracks.SliceUtils;
import mixer.utils.tracks.SubcompartmentInterval;

import java.io.File;
import java.util.Random;

public class ClusteringMagic {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    private final File outputDirectory;
    private static final int maxIters = 500;
    private final Random generator = new Random(2352);
    private final MatrixAndWeight matrix;
    private final ChromosomeHandler handler;

    public ClusteringMagic(MatrixAndWeight matrix, File outputDirectory,
                           ChromosomeHandler handler, long seed) {
        this.matrix = matrix;
        this.handler = handler;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);
    }

    public void extractFinalGWSubcompartments(String prefix) {
        System.out.println("Genomewide clustering");
        matrix.inPlaceScaleSqrtWeightCol();
        runClusteringOnMatrix(prefix, false);

        matrix.inPlaceScaleSqrtWeightCol();
        runClusteringOnMatrix(prefix, true);
    }

    private void runClusteringOnMatrix(String prefix, boolean useKMedians) {
        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(handler, matrix,
                false, useKMedians);
        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            runKMeansMultipleTimes(kmeansRunner, z, useKMedians, prefix);
        }
        System.out.println(".");
    }

    private void runKMeansMultipleTimes(GenomeWideKmeansRunner kmeansRunner,
                                        int z, boolean useKMedians, String prefix) {
        int numClusters = z + startingClusterSizeK;
        int[][] results = new int[3][matrix.getNumRows()];
        double wcssLimit = getGoodWCSS(results, kmeansRunner, numClusters);
        int index = 1;
        while (index < 3) {
            System.out.print(".");
            KmeansResult currResult = resetAndRerun(kmeansRunner, generator, numClusters);
            double wcss = currResult.getWithinClusterSumOfSquares();
            if (currResult.getNumActualClusters() == numClusters && wcss <= wcssLimit) {
                exportKMeansClusteringResults(z, prefix, useKMedians, index,
                        wcss, currResult.getFinalCompartmentsClone());
                results[index] = currResult.getAssignments(matrix.getNumRows());
                index++;
            }
        }
        Concensus3DTools.resolve(results, this, z, prefix, useKMedians, numClusters, matrix, handler);
    }

    private double getGoodWCSS(int[][] results, GenomeWideKmeansRunner kmeansRunner, int numClusters) {
        double wcssLimit = Float.MAX_VALUE;
        int[] bestResult = null;
        int attempts = 0;
        while (bestResult == null || attempts < 20) {
            System.out.print("-");
            KmeansResult currResult = resetAndRerun(kmeansRunner, generator, numClusters);
            double wcss = currResult.getWithinClusterSumOfSquares();
            if (currResult.getNumActualClusters() == numClusters && wcss < Float.MAX_VALUE) {
                if (wcss < wcssLimit) {
                    wcssLimit = wcss;
                    bestResult = currResult.getAssignments(matrix.getNumRows());
                }
                attempts++;
            }
        }
        results[0] = bestResult;
        return 1.01 * wcssLimit;
    }

    private KmeansResult resetAndRerun(GenomeWideKmeansRunner kmeansRunner, Random generator, int numClusters) {
        kmeansRunner.prepareForNewRun(numClusters);
        kmeansRunner.launchKmeansGWMatrix(generator.nextLong(), maxIters);
        return kmeansRunner.getResult();
    }

    public void exportKMeansClusteringResults(int z, String prefix, boolean useKMedians, int iter, double wcss,
                                              GenomeWide1DList<SubcompartmentInterval> finalCompartments) {
        int k = z + startingClusterSizeK;
        String kstem = "kmeans";
        if (useKMedians) kstem = "kmedians";
        SliceUtils.collapseGWList(finalCompartments);
        File outBedFile = new File(outputDirectory, prefix + "_" + kstem + "_k" + k + "_i" + iter + "_clusters.bed"); // "_wcss" + wcss +
        finalCompartments.simpleExport(outBedFile);
    }
}
