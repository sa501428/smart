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

package mixer.utils.slice.gmm;

import mixer.clt.ParallelizedMixerTools;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class InitTools {

    public static float[][] parGetInitialMeans(List<List<Integer>> groupsOfIndices, int numClusters, float[][] data) {
        float[][] meanVectors = new float[numClusters][data[0].length];
        int[][] counts = new int[numClusters][data[0].length];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int k = currIndex.getAndIncrement();
            while (k < numClusters) {
                for (int i : groupsOfIndices.get(k)) {
                    for (int j = 0; j < data[i].length; j++) {
                        if (!Float.isNaN(data[i][j])) {
                            meanVectors[k][j] += data[i][j];
                            counts[k][j]++;
                        }
                    }
                }

                for (int j = 0; j < data[0].length; j++) {
                    meanVectors[k][j] /= Math.max(counts[k][j], 1);
                }
                k = currIndex.getAndIncrement();
            }
        });

        return meanVectors;
    }

    public static double[][][] parGetInitialFeatureCovMatrices(List<List<Integer>> indices,
                                                               int numClusters, float[][] data,
                                                               float[][] meanVectors) {
        double[][][] covMatrices = new double[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = parGetInitialColumnCovarianceMatrix(indices.get(k), data, meanVectors[k]);
        }
        GMMCovTools.ensureValidCovMatrix(covMatrices);
        return covMatrices;
    }

    public static double[][] parGetInitialColumnCovarianceMatrix(List<Integer> indices, float[][] data,
                                                                 float[] meanVector) {
        float[][] diff = calculateDiffWithIndices(indices, data, meanVector);
        int numDataPoints = diff.length;
        int dimension = data[0].length;
        double[][] cov = new double[dimension][dimension];
        int[] min = new int[]{numDataPoints};

        AtomicInteger currRowIndex = new AtomicInteger(0);

        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < dimension) {
                for (int j = i; j < dimension; j++) {
                    //System.out.println(i+" "+j+" "+clusterID);
                    double accum = 0;
                    int numGoodPoints = 0;
                    for (int k = 0; k < numDataPoints; k++) {
                        double val = diff[k][i] * diff[k][j];
                        if (!Double.isNaN(val)) {
                            accum += val;
                            numGoodPoints++;
                        }
                    }
                    if (numGoodPoints > 0) {
                        if (numGoodPoints < min[0]) {
                            min[0] = numGoodPoints;
                        }
                        cov[i][j] = (float) (accum / numGoodPoints);
                        cov[j][i] = cov[i][j]; // symmetric
                    }
                }
                i = currRowIndex.getAndIncrement();
            }
        });

        //System.out.println("min num entries " + min[0]);
        return cov;
    }

    private static float[][] calculateDiffWithIndices(List<Integer> indices, float[][] data, float[] meanVector) {
        int numDataPoints = indices.size();
        int dimension = data[0].length;
        float[][] diff = new float[numDataPoints][dimension];
        int counter = 0;
        for (int i : indices) {
            for (int j = 0; j < dimension; j++) {
                diff[counter][j] = data[i][j] - meanVector[j];
            }
            counter++;
        }
        return diff;
    }
}
