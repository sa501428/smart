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

import java.util.Arrays;
import java.util.List;

public class GaussianMixtureModels {
    private final float[][] data;
    private final int numClusters;
    private final float[] datasetFractionForCluster;
    private final List<List<Integer>> startingIndices;
    float[][] meanVectors; // add the rows and divide by num
    float[][][] covMatrices;
    float[][] probabilities;
    private int maxIters = 10;
    private boolean startingFromScratch = true;

    public GaussianMixtureModels(float[][] data, int numClusters, int maxIters, List<List<Integer>> startingIndices) {
        this.data = data;
        this.numClusters = numClusters;
        this.maxIters = maxIters;
        this.startingIndices = startingIndices;

        datasetFractionForCluster = new float[numClusters];
        for (int k = 0; k < numClusters; k++) {
            float num = startingIndices.get(k).size();
            float denom = data.length;
            datasetFractionForCluster[k] = num / denom;
        }
        System.out.println("Priors: " + Arrays.toString(datasetFractionForCluster));
    }


    public void fit() {

        if (startingFromScratch && startingIndices != null && startingIndices.size() == numClusters) {
            meanVectors = GMMTools.getInitialMeans(startingIndices, numClusters, data); // add the rows and divide by num
            covMatrices = GMMCovTools.getFeatureCovMatrices(startingIndices, numClusters, data);
            startingIndices.clear();
            startingFromScratch = false;
            System.out.println("Initial setting done");
        }

        for (int iter = 0; iter < maxIters; iter++) {
            float[][] probClusterForRow = GMMTools.getProbabilityOfClusterForRow(numClusters, data, datasetFractionForCluster, meanVectors, covMatrices);
            float[] totalSumForCluster = GMMTools.addUpAllRows(probClusterForRow);

            meanVectors = GMMTools.getNewWeightedAverage(numClusters, data, probClusterForRow);
            covMatrices = GMMCovTools.getNewWeightedFeatureCovarianceMatrix(numClusters, data, probClusterForRow,
                    totalSumForCluster);

            for (int k = 0; k < numClusters; k++) {
                datasetFractionForCluster[k] = totalSumForCluster[k] / (float) (data.length);
            }
        }
    }

    public int[] predict() {
        probabilities = new float[data.length][numClusters];

        for (int i = 0; i < data.length; i++) {
            for (int k = 0; k < numClusters; k++) {
                probabilities[i][k] = GMMTools.multivariateNormal(data[i], meanVectors[k], covMatrices[k]);
            }
        }

        int[] assignments = new int[data.length];
        for (int i = 0; i < data.length; i++) {
            float prob = probabilities[i][0];
            int index = 0;
            for (int k = 1; k < numClusters; k++) {
                if (probabilities[i][k] > prob) {
                    prob = probabilities[i][k];
                    index = k;
                }
            }
            assignments[i] = index;
        }
        return assignments;
    }

    /*

    def predict(self, X):
            '''
    The predicting function
                :param X: 2-d array numpy array
    The data on which we must predict the clusters
        '''
    probas = []
            for n in range(len(X)):
            probas.append([self.multivariate_normal(X[n], self.mean_vector[k], self.covariance_matrixes[k])
            for k in range(self.n_componets)])
    cluster = []
            for proba in probas:
            cluster.append(self.comp_names[proba.index(max(proba))])
            return cluster

     */
}
