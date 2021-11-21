/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Rice University, Baylor College of Medicine, Aiden Lab
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

package mixer.utils.slice.kmeans;

import javastraw.feature1D.GenomeWideList;
import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.type.NormalizationType;
import mixer.algos.Slice;
import mixer.utils.slice.CorrMatrixClusterer;
import mixer.utils.slice.EncodeExportUtils;
import mixer.utils.slice.cleaning.GWBadIndexFinder;
import mixer.utils.slice.matrices.CompositeGenomeWideMatrix;
import mixer.utils.slice.matrices.SliceMatrix;
import mixer.utils.slice.structures.SliceUtils;
import mixer.utils.slice.structures.SubcompartmentInterval;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class FullGenomeOEWithinClusters {
    public static int startingClusterSizeK = 2;
    public static int numClusterSizeKValsUsed = 10;
    public static int numAttemptsForKMeans = 3;
    private final File outputDirectory;
    private final ChromosomeHandler chromosomeHandler;
    private final SliceMatrix sliceMatrix;
    private final int maxIters = 200;

    private final Random generator = new Random(0);

    public FullGenomeOEWithinClusters(List<Dataset> datasets, ChromosomeHandler chromosomeHandler, int resolution,
                                      List<NormalizationType[]> normalizationTypes,
                                      File outputDirectory, long seed) {
        this.chromosomeHandler = chromosomeHandler;
        this.outputDirectory = outputDirectory;
        generator.setSeed(seed);

        GWBadIndexFinder badIndexFinder = new GWBadIndexFinder(chromosomeHandler.getAutosomalChromosomesArray(),
                resolution, normalizationTypes);
        badIndexFinder.createInternalBadList(datasets, chromosomeHandler.getAutosomalChromosomesArray());

        int absMaxClusters = numClusterSizeKValsUsed + startingClusterSizeK;
        sliceMatrix = new SliceMatrix(chromosomeHandler, datasets.get(0), normalizationTypes.get(0), resolution, outputDirectory,
                generator.nextLong(), badIndexFinder, absMaxClusters);

        for (int dI = 1; dI < datasets.size(); dI++) {
            SliceMatrix additionalData = new SliceMatrix(chromosomeHandler, datasets.get(dI),
                    normalizationTypes.get(dI), resolution, outputDirectory,
                    generator.nextLong(), badIndexFinder, absMaxClusters);
            sliceMatrix.appendDataAlongExistingRows(additionalData);
        }

        sliceMatrix.cleanUpMatricesBySparsity();
    }

    public void extractFinalGWSubcompartments(String prefix) {
        System.out.println("Genomewide clustering");
        if (Slice.USE_KMEANS) {
            processClustering(prefix, false);
        }

        if (Slice.USE_KMEDIANS) {
            sliceMatrix.inPlaceScaleSqrtWeightCol(); // due to l1 issue
            processClustering(prefix, true);
        }
    }

    private void processClustering(String prefix, boolean useKMedians) {
        runClusteringOnRawMatrixWithNans(prefix, useKMedians);
        if (Slice.USE_INTER_CORR_CLUSTERING) {
            CorrMatrixClusterer.runClusteringOnCorrMatrix(this, prefix + "_corr", useKMedians);
        }
    }

    public void runClusteringOnRawMatrixWithNans(String prefix, boolean useKMedians) {
        Map<Integer, GenomeWideList<SubcompartmentInterval>> kmeansClustersToResults = new HashMap<>();
        Map<Integer, List<List<Integer>>> kmeansIndicesMap = new HashMap<>();

        GenomeWideKmeansRunner kmeansRunner = new GenomeWideKmeansRunner(chromosomeHandler, sliceMatrix,
                false, useKMedians);
        KmeansEvaluator evaluator = new KmeansEvaluator(numClusterSizeKValsUsed);

        for (int z = 0; z < numClusterSizeKValsUsed; z++) {
            runRepeatedKMeansClusteringLoop(numAttemptsForKMeans, kmeansRunner, evaluator, z,
                    maxIters, kmeansClustersToResults, kmeansIndicesMap, useKMedians);
            exportKMeansClusteringResults(z, kmeansClustersToResults, prefix, kmeansIndicesMap, useKMedians);
        }
        if (Slice.USE_ENCODE_MODE) {
            EncodeExportUtils.exportSubcompartments(sliceMatrix, kmeansIndicesMap, evaluator);
        }
        exportEvaluatorInfo(evaluator, useKMedians);
        System.out.println(".");
    }

    public void exportKMeansClusteringResults(int z,
                                              Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults,
                                              String prefix, Map<Integer, List<List<Integer>>> kmeansIndicesMap,
                                              boolean useKMedians) {
        int k = z + startingClusterSizeK;
        String kstem = "kmeans";
        if (useKMedians) kstem = "kmedians";
        GenomeWideList<SubcompartmentInterval> gwList = numClustersToResults.get(k);
        SliceUtils.collapseGWList(gwList);
        File outBedFile = new File(outputDirectory, prefix + "_" + k + "_" + kstem + "_clusters.bed");
        gwList.simpleExport(outBedFile);
        sliceMatrix.plotUmapProjection(outputDirectory, kmeansIndicesMap.get(z), prefix + "_" + k + "_" + kstem + "_clusters");
    }

    public void runRepeatedKMeansClusteringLoop(int attemptsForKMeans, GenomeWideKmeansRunner kmeansRunner,
                                                KmeansEvaluator evaluator, int z, int maxIters,
                                                Map<Integer, GenomeWideList<SubcompartmentInterval>> numClustersToResults,
                                                Map<Integer, List<List<Integer>>> indicesMap,
                                                boolean useKMedians) {
        int numClusters = z + startingClusterSizeK;
        int numColumns = kmeansRunner.getNumColumns();
        int numRows = kmeansRunner.getNumRows();
        for (int p = 0; p < attemptsForKMeans; p++) {
            kmeansRunner.prepareForNewRun(numClusters);
            kmeansRunner.launchKmeansGWMatrix(generator.nextLong(), maxIters);

            int numActualClustersThisAttempt = kmeansRunner.getNumActualClusters();
            if (numActualClustersThisAttempt == numClusters) {
                double wcss = kmeansRunner.getWithinClusterSumOfSquares();
                if (wcss < evaluator.getWCSS(z)) {
                    evaluator.setMseAicBicValues(z, numClusters, wcss, numRows, numColumns,
                            kmeansRunner.getSilhouette());
                    indicesMap.put(z, kmeansRunner.getIndicesMapCopy());
                    numClustersToResults.put(numClusters, kmeansRunner.getFinalCompartments());
                }
            }
            System.out.print(".");
        }
    }

    public File getOutputDirectory() {
        return outputDirectory;
    }

    public ChromosomeHandler getChromosomeHandler() {
        return chromosomeHandler;
    }

    public CompositeGenomeWideMatrix getSliceMatrix() {
        return sliceMatrix;
    }

    public int getMaxIters() {
        return maxIters;
    }

    public void exportEvaluatorInfo(KmeansEvaluator evaluator, boolean useKmedians) {
        String kstem = "kmeans";
        if (useKmedians) kstem = "kmedians";
        evaluator.export(outputDirectory, kstem);
    }
}
