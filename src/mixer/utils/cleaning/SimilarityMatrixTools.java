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

package mixer.utils.cleaning;

import mixer.SmartTools;
import mixer.utils.common.ZScoreTools;
import mixer.utils.intra.IndexOrderer;
import mixer.utils.kmeans.QuickCentroids;
import mixer.utils.similaritymeasures.RobustCorrelationSimilarity;
import mixer.utils.similaritymeasures.RobustCosineSimilarity;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class SimilarityMatrixTools {

    public static SimpleMatrixAndWeight getNonNanSimilarityMatrix(float[][] matrix, SimilarityMetric metric,
                                                                  int numPerCentroid, long seed) {
        if ((!metric.isSymmetric()) || numPerCentroid > 1) {
            return getAsymmetricMatrix(matrix, new SimilarityMetric[]{metric},
                    matrix.length / numPerCentroid, seed);
        }

        return getSymmetricMatrix(matrix, metric);
    }

    public static SimpleMatrixAndWeight getCompressedCosineSimilarityMatrix(float[][] matrix,
                                                                            int numCentroids, long seed) {
        return getAsymmetricMatrix(matrix, new SimilarityMetric[]{RobustCosineSimilarity.SINGLETON},
                numCentroids, seed);
    }

    public static SimpleMatrixAndWeight getCosinePearsonCorrMatrix(float[][] matrix, int numCentroids, long seed) {
        RobustCosineSimilarity.USE_ARC = true;
        RobustCorrelationSimilarity.USE_ARC = true;

        SimilarityMetric[] metrics = new SimilarityMetric[]{
                RobustCosineSimilarity.SINGLETON,
                //RobustCorrelationSimilarity.SINGLETON
                //RobustEuclideanDistance.SINGLETON,
                //RobustManhattanDistance.SINGLETON
        };

        SimpleMatrixAndWeight answer = getAsymmetricMatrix(matrix, metrics, numCentroids, seed);

        RobustCosineSimilarity.USE_ARC = false;
        RobustCorrelationSimilarity.USE_ARC = false;

        return answer;
    }

    private static SimpleMatrixAndWeight getAsymmetricMatrix(float[][] matrix, SimilarityMetric[] metrics,
                                                             int numInitCentroids, long seed) {
        QuickCentroids centroidMaker = new QuickCentroids(matrix, numInitCentroids, seed, 20);
        final float[][] centroids = centroidMaker.generateCentroids(5, true);
        int[] weights = centroidMaker.getWeights();


        int numCentroids = centroids.length;
        if (SmartTools.printVerboseComments || centroids.length != numInitCentroids) {
            System.out.println("AsymMatrix: Was initially " + numInitCentroids + " centroids, but using " + numCentroids);
        }

        float[][] result = new float[matrix.length][numCentroids * metrics.length];
        int numCPUThreads = Runtime.getRuntime().availableProcessors();

        System.out.println("... generating sym matrix");
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int i = currRowIndex.getAndIncrement();
                while (i < matrix.length) {
                    for (int j = 0; j < numCentroids; j++) {
                        for (int z = 0; z < metrics.length; z++) {
                            int offset = z * numCentroids;
                            result[i][j + offset] = metrics[z].distance(centroids[j], matrix[i]);
                        }
                    }
                    i = currRowIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        // Wait until all threads finish
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }

        return new SimpleMatrixAndWeight(result, weights);
    }

    private static SimpleMatrixAndWeight getSymmetricMatrix(float[][] matrix, SimilarityMetric metric) {
        float[][] result = new float[matrix.length][matrix.length]; // *2

        int numCPUThreads = Runtime.getRuntime().availableProcessors() * 2;
        System.out.println(" .. ");
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int i = currRowIndex.getAndIncrement();
                while (i < matrix.length) {
                    for (int j = i; j < matrix.length; j++) {
                        result[i][j] = metric.distance(matrix[i], matrix[j]);
                        result[j][i] = result[i][j];
                    }
                    i = currRowIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        // Wait until all threads finish
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }

        int[] weights = new int[result[0].length];
        Arrays.fill(weights, 1);
        return new SimpleMatrixAndWeight(result, weights);
    }

    public static float[][] getSymmNonNanSimilarityMatrixWithMask(float[][] initialMatrix,
                                                                  SimilarityMetric metric,
                                                                  int[] newIndexOrderAssignments, int checkVal) {

        float[][] result = new float[initialMatrix.length][initialMatrix.length];
        for (float[] row : result) {
            Arrays.fill(row, Float.NaN);
        }

        int numCPUThreads = Runtime.getRuntime().availableProcessors();

        RobustCorrelationSimilarity.USE_ARC = true;

        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int i = currRowIndex.getAndIncrement();
                while (i < initialMatrix.length) {
                    if (newIndexOrderAssignments[i] < checkVal) {
                        // result[i][i] = Float.NaN; // technically 1, but not useful
                        for (int j = i + 1; j < initialMatrix.length; j++) {
                            if (newIndexOrderAssignments[j] < checkVal) {
                                result[i][j] = metric.distance(initialMatrix[j], initialMatrix[i]);
                                result[j][i] = result[i][j];
                            }
                        }
                    }
                    i = currRowIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        // Wait until all threads finish
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }
        RobustCorrelationSimilarity.USE_ARC = false;

        return result;
    }

    public static float[][] getAsymNonNanSimilarityMatrixWithMask(float[][] initialMatrix,
                                                                  RobustCorrelationSimilarity metric, int numInitCentroids,
                                                                  int[] newIndexOrderAssignments, int checkVal) {

        QuickCentroids centroidMaker = new QuickCentroids(IndexOrderer.quickCleanMatrix(initialMatrix, newIndexOrderAssignments),
                numInitCentroids, 0L, 20);
        final float[][] centroids = centroidMaker.generateCentroids(3, true);
        int[] weights = centroidMaker.getWeights();

        int numCentroids = centroids.length;
        if (SmartTools.printVerboseComments || centroids.length != numInitCentroids) {
            System.out.println("AsymMatrix: Was initially " + numInitCentroids + " centroids, but using " + numCentroids);
        }

        float[][] result = new float[initialMatrix.length][numCentroids];
        for (float[] row : result) {
            Arrays.fill(row, Float.NaN);
        }

        int numCPUThreads = Runtime.getRuntime().availableProcessors();

        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int i = currRowIndex.getAndIncrement();
                while (i < initialMatrix.length) {
                    if (newIndexOrderAssignments[i] < checkVal) {
                        for (int j = 0; j < numCentroids; j++) {
                            result[i][j] = metric.distance(centroids[j], initialMatrix[i]);
                        }
                    }
                    i = currRowIndex.getAndIncrement();
                }
            };
            executor.execute(worker);
        }
        executor.shutdown();
        // Wait until all threads finish
        //noinspection StatementWithEmptyBody
        while (!executor.isTerminated()) {
        }

        ZScoreTools.inPlaceScaleSqrtWeightCol(result, weights);

        return result;
    }
}
