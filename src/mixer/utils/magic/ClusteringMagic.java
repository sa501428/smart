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

package mixer.utils.magic;

import javastraw.feature1D.GenomeWide1DList;
import javastraw.reader.basics.ChromosomeHandler;
import mixer.utils.slice.kmeans.GenomeWideKmeansRunner;
import mixer.utils.slice.kmeans.KmeansEvaluator;
import mixer.utils.slice.kmeans.KmeansResult;
import mixer.utils.slice.matrices.MatrixAndWeight;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class ClusteringMagic {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    public static int numAttemptsForKMeans = 3;
    private final File outputDirectory;
    private final int maxIters = 200;
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
        processClustering(prefix, false);

        matrix.inPlaceScaleSqrtWeightCol();
        processClustering(prefix, true);
    }

    private void processClustering(String prefix, boolean useKMedians) {
        runClusteringOnMatrix(prefix, useKMedians);
        // todo CorrMatrixClusterer.runClusteringOnCorrMatrix(this, prefix + "_corr", useKMedians);
    }

    public void runClusteringOnMatrix(String prefix, boolean useKMedians) {
        Map<Integer, GenomeWide1DList<SubcompartmentInterval>> kmeansClustersToResults = new HashMap<>();
        Map<Integer, List<List<Integer>>> kmeansIndicesMap = new HashMap<>();

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(handler, matrix,
                false, useKMedians);
        KmeansEvaluator evaluator = new KmeansEvaluator(numClusterSizeKValsUsed);

        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            runRepeatedKMeansClusteringLoop(numAttemptsForKMeans, kmeansRunner, evaluator, z,
                    maxIters, kmeansClustersToResults, kmeansIndicesMap);
            exportKMeansClusteringResults(z, kmeansClustersToResults, prefix, kmeansIndicesMap, useKMedians);
        }
        System.out.println(".");
    }

    public void exportKMeansClusteringResults(int z,
                                              Map<Integer, GenomeWide1DList<SubcompartmentInterval>> numClustersToResults,
                                              String prefix, Map<Integer, List<List<Integer>>> kmeansIndicesMap,
                                              boolean useKMedians) {
        int k = z + startingClusterSizeK;
        String kstem = "kmeans";
        if (useKMedians) kstem = "kmedians";
        GenomeWide1DList<SubcompartmentInterval> gwList = numClustersToResults.get(k);
        SliceUtils.collapseGWList(gwList);
        File outBedFile = new File(outputDirectory, prefix + "_" + k + "_" + kstem + "_clusters.bed");
        gwList.simpleExport(outBedFile);
    }

    public void runRepeatedKMeansClusteringLoop(int attemptsForKMeans, GenomeWideKmeansRunner kmeansRunner,
                                                KmeansEvaluator evaluator, int z, int maxIters,
                                                Map<Integer, GenomeWide1DList<SubcompartmentInterval>> numClustersToResults,
                                                Map<Integer, List<List<Integer>>> indicesMap) {
        int numClusters = z + startingClusterSizeK;
        int numColumns = kmeansRunner.getNumColumns();
        int numRows = kmeansRunner.getNumRows();
        for (int p = 0; p < attemptsForKMeans; p++) {
            boolean noClusteringFound = true;
            while (noClusteringFound) {
                kmeansRunner.prepareForNewRun(numClusters);
                kmeansRunner.launchKmeansGWMatrix(generator.nextLong(), maxIters);
                KmeansResult currResult = kmeansRunner.getResult();
                if (currResult.getNumActualClusters() == numClusters) {
                    noClusteringFound = false;
                    if (currResult.getWithinClusterSumOfSquares() < evaluator.getWCSS(z)) {
                        evaluator.setMseAicBicValues(z, numRows, numColumns, currResult);
                        indicesMap.put(z, currResult.getIndicesMapClone());
                        numClustersToResults.put(numClusters, currResult.getFinalCompartmentsClone());
                    }
                }
                System.out.print(".");
            }
        }
    }
}
