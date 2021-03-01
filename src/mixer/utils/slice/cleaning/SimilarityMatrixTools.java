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

package mixer.utils.slice.cleaning;

import mixer.MixerGlobals;
import mixer.utils.common.FloatMatrixTools;
import mixer.utils.similaritymeasures.SimilarityMetric;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class SimilarityMatrixTools {

    public static float[][] getNonNanSimilarityMatrix(float[][] matrix, SimilarityMetric metric, int numCentroid,
                                                      long seed) {
        if ((!metric.isSymmetric()) || numCentroid > 1) {
            return getAsymmetricMatrix(matrix, metric, numCentroid, seed);
        }

        return getSymmetricMatrix(matrix, metric);
    }

    private static float[][] getAsymmetricMatrix(float[][] matrix, SimilarityMetric metric, int numPerCentroid, long seed) {
        final int numInitialCentroids = matrix.length / numPerCentroid;

        final float[][] centroids;
        final float[] weights;
        if (numInitialCentroids > 1) {
            QuickCentroids centroidMaker = new QuickCentroids(matrix, numInitialCentroids, seed, 100);
            centroids = centroidMaker.generateCentroids();
            weights = centroidMaker.getWeights();
            for (int i = 0; i < weights.length; i++) {
                weights[i] = (float) Math.sqrt(weights[i]);
            }
        } else {
            centroids = matrix;
            weights = new float[centroids.length];
            Arrays.fill(weights, 1);
        }

        int numCentroids = centroids.length;
        if (MixerGlobals.printVerboseComments || centroids.length != numInitialCentroids) {
            System.out.println("AsymMatrix: Was initially " + numInitialCentroids + " centroids, but using " + numCentroids);
        }

        float[][] result = new float[matrix.length][numCentroids];
        int numCPUThreads = Runtime.getRuntime().availableProcessors();

        System.out.println("... generating sym matrix");
        AtomicInteger currRowIndex = new AtomicInteger(0);
        ExecutorService executor = Executors.newFixedThreadPool(numCPUThreads);
        for (int l = 0; l < numCPUThreads; l++) {
            Runnable worker = () -> {
                int i = currRowIndex.getAndIncrement();
                while (i < matrix.length) {
                    for (int j = 0; j < numCentroids; j++) {
                        result[i][j] = metric.distance(centroids[j], matrix[i], 0, 3);

                        if (Float.isNaN(result[i][j])) {
                            System.err.println("Error appearing in distance measure...");
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

        /*
        int[] weights3 = new int[weights.length * 3];
        System.arraycopy(weights, 0, weights3, 0, weights.length);
        System.arraycopy(weights, 0, weights3, offset1, weights.length);
        System.arraycopy(weights, 0, weights3, offset2, weights.length);
         */

        //ZScoreTools.inPlaceScaleCol(result, weights);
        //ZScoreTools.inPlaceScaleCol(result, weights);

        return FloatMatrixTools.concatenate(matrix, result);
    }

    private static float[][] getSymmetricMatrix(float[][] matrix, SimilarityMetric metric) {
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
                        result[i][j] = metric.distance(matrix[i], matrix[j], 0, 1);
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

        //int[] weights = new int[result[0].length];
        //Arrays.fill(weights, 1);
        //ZScoreTools.inPlaceZscoreDownCol(result, weights);
        return result;
    }
}
