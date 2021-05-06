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

import mixer.MixerGlobals;
import mixer.clt.ParallelizedMixerTools;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class GaussianMixtureModels {
    private final float[][] data;
    private final int numClusters;
    private final double[] datasetFractionForCluster;
    private final List<List<Integer>> startingIndices;
    float[][] meanVectors; // add the rows and divide by num
    double[][][] covMatrices;
    double[][] probabilities;
    private int maxIters = 10;
    private boolean startingFromScratch = true;

    public GaussianMixtureModels(float[][] data, int numClusters, int maxIters, List<List<Integer>> startingIndices) {
        this.data = data;
        this.numClusters = numClusters;
        this.maxIters = maxIters;
        this.startingIndices = startingIndices;

        datasetFractionForCluster = new double[numClusters];
        for (int k = 0; k < numClusters; k++) {
            float num = startingIndices.get(k).size();
            float denom = data.length;
            datasetFractionForCluster[k] = num / denom;
        }
    }


    public void fit() {
        if (startingFromScratch && startingIndices != null && startingIndices.size() == numClusters) {
            meanVectors = InitTools.parGetInitialMeans(startingIndices, numClusters, data); // add the rows and divide by num
            covMatrices = InitTools.parGetInitialFeatureCovMatrices(startingIndices, numClusters, data, meanVectors);

            startingIndices.clear();
            startingFromScratch = false;
            if (MixerGlobals.printVerboseComments) {
                System.out.println("GMM initialization done");
            }
        }

        for (int iter = 0; iter < maxIters; iter++) {
            System.out.println("GMM Iteration " + iter);
            if (MixerGlobals.printVerboseComments) {
                for (int k = 0; k < meanVectors.length; k++) {
                    System.err.println("mu[" + k + "] " + Arrays.toString(meanVectors[k]));
                }
            }
            double[][] probClusterForRow = GMMTools.parGetProbabilityOfClusterForRow(numClusters, data, datasetFractionForCluster,
                    meanVectors, covMatrices);

            if (MixerGlobals.printVerboseComments) {
                System.out.println();
                for (int k = 0; k < 10; k++) {
                    System.err.println("pi[" + k + "] " + Arrays.toString(probClusterForRow[k]));
                }
            }

            float[] totalSumForCluster = GMMTools.addUpAllRows(probClusterForRow);

            if (MixerGlobals.printVerboseComments) {
                System.out.println();
                System.err.println("N = " + Arrays.toString(totalSumForCluster));
            }

            meanVectors = GMMTools.parGetWeightedMean(numClusters, data, probClusterForRow);
            covMatrices = GMMCovTools.parGetNewWeightedFeatureCovarianceMatrix(numClusters, data, probClusterForRow,
                    meanVectors);

            for (int k = 0; k < numClusters; k++) {
                datasetFractionForCluster[k] = totalSumForCluster[k] / (float) (data.length);
            }
        }
    }

    public int[] predict() {
        probabilities = new double[data.length][numClusters];

        AtomicInteger currentIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currentIndex.getAndIncrement();
            while (i < data.length) {
                double[] tempArray = new double[numClusters];
                double accum = 0;
                for (int k = 0; k < numClusters; k++) {
                    tempArray[k] = GMMTools.multivariateNormal(data[i], meanVectors[k], covMatrices[k]);
                    accum += tempArray[k];
                }
                for (int k = 0; k < numClusters; k++) {
                    probabilities[i][k] = (float) (tempArray[k] / accum);
                }
                i = currentIndex.getAndIncrement();
            }
        });

        int[] assignments = new int[data.length];
        AtomicInteger currentIndex2 = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currentIndex2.getAndIncrement();
            while (i < data.length) {
                double prob = probabilities[i][0];
                int index = 0;
                for (int k = 1; k < numClusters; k++) {
                    if (probabilities[i][k] > prob) {
                        prob = probabilities[i][k];
                        index = k;
                    }
                }
                assignments[i] = index;
                i = currentIndex2.getAndIncrement();
            }
        });

        return assignments;
    }
}
