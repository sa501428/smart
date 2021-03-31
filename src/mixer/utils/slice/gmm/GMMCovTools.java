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

import java.util.concurrent.atomic.AtomicInteger;

public class GMMCovTools {

    public static float[][][] parGetNewWeightedFeatureCovarianceMatrix(int numClusters, float[][] data,
                                                                       float[][] probClusterForRow,
                                                                       float[][] meanVectors) {
        float[][][] covMatrices = new float[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = parGetWeightedColumnCovarianceMatrix(data, probClusterForRow, k, meanVectors[k]);
        }

        // ensure positive nonzero diagonal
        addEpsilonToDiagonal(covMatrices);

        return covMatrices;
    }

    public static void addEpsilonToDiagonal(float[][][] covMatrices) {
        for (int i = 0; i < covMatrices.length; i++) {
            for (int j = 0; j < covMatrices[i].length; j++) {
                covMatrices[i][j][j] += 1e-6;
            }
        }
    }

    public static float[][] parGetWeightedColumnCovarianceMatrix(float[][] data, float[][] probClusterForRow,
                                                                 int clusterID, float[] meanVector) {

        float[][] diff = parGetDiffMatrix(data, meanVector);
        int numDataPoints = data.length;
        int dimension = data[0].length;

        float[][] cov = new float[dimension][dimension];

        AtomicInteger currRowIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < dimension) {
                for (int j = i; j < dimension; j++) {
                    float weight = 0;
                    float accum = 0;
                    for (int k = 0; k < numDataPoints; k++) {
                        float val = diff[k][i] * diff[k][j];
                        if (!Float.isNaN(val)) {
                            accum += probClusterForRow[k][clusterID] * val;
                            weight += probClusterForRow[k][clusterID];
                        }
                    }
                    cov[i][j] = accum / weight;
                    cov[j][i] = cov[i][j]; // symmetric
                }
                i = currRowIndex.getAndIncrement();
            }
        });

        return cov;
    }

    private static float[][] parGetDiffMatrix(float[][] data, float[] meanVector) {
        int numDataPoints = data.length;
        int dimension = data[0].length;

        float[][] diff = new float[numDataPoints][dimension];
        AtomicInteger currDataIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currDataIndex.getAndIncrement();
            while (i < numDataPoints) {
                for (int j = 0; j < dimension; j++) {
                    diff[i][j] = data[i][j] - meanVector[j];
                }
                i = currDataIndex.getAndIncrement();
            }
        });
        return diff;
    }
}
