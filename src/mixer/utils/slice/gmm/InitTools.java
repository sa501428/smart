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

    public static float[][][] parGetInitialFeatureCovMatrices(List<List<Integer>> indices,
                                                              int numClusters, float[][] data,
                                                              float[][] meanVectors) {
        float[][][] covMatrices = new float[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = parGetInitialColumnCovarianceMatrix(indices.get(k), data, meanVectors[k]);
        }
        return covMatrices;
    }

    public static float[][] parGetInitialColumnCovarianceMatrix(List<Integer> indices, float[][] data,
                                                                float[] meanVector) {

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

        float[][] cov = new float[dimension][dimension];

        AtomicInteger currRowIndex = new AtomicInteger(0);

        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < dimension) {
                for (int j = i; j < dimension; j++) {
                    //System.out.println(i+" "+j+" "+clusterID);
                    float accum = 0;
                    int numGoodPoints = 0;
                    for (int k = 0; k < numDataPoints; k++) {
                        float val = diff[k][i] * diff[k][j];
                        if (!Float.isNaN(val)) {
                            accum += val;
                            numGoodPoints++;
                        }
                    }
                    cov[i][j] = accum / numGoodPoints;
                    cov[j][i] = cov[i][j]; // symmetric
                }
                i = currRowIndex.getAndIncrement();
            }
        });

        return cov;
    }
}
