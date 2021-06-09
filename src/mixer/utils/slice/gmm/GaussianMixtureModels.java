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
import mixer.utils.slice.cleaning.QuickClusters;
import mixer.utils.slice.gmm.robust.RobustGMMCovTools;
import mixer.utils.slice.gmm.robust.RobustGMMTools;
import mixer.utils.slice.gmm.simple.SimpleGMMCovTools;
import mixer.utils.slice.gmm.simple.SimpleGMMTools;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class GaussianMixtureModels {
    private static final float startingProbability = 0.85f;
    private final float[][] data;
    private final int numClusters;
    private final int maxIters;
    private final List<List<Integer>> startingIndices;
    private float[][] meanVectors; // add the rows and divide by num
    private CovarianceMatrixInverseAndDeterminant[] covs;
    private final boolean useRobustGMM;
    private double[] datasetFractionForCluster;
    private boolean startingFromScratch = true;
    private double[][] probabilities;
    private RealMatrix[] covMatrices;

    public GaussianMixtureModels(float[][] data, int numClusters, int maxIters,
                                 List<List<Integer>> startingIndices, boolean useRobustGMM) {
        this.data = data;
        this.numClusters = numClusters;
        this.maxIters = maxIters;
        this.startingIndices = startingIndices;
        this.useRobustGMM = useRobustGMM;
        if (startingIndices.size() != numClusters) {
            System.err.println("GMM Error: something weird about cluster sizes " + numClusters + " : " + startingIndices.size());
        }
    }

    public GaussianMixtureModels(float[][] extraction, int numClusters, int maxIters, boolean useRobustGMM) {
        if (!useRobustGMM) {
            this.data = removeNans(extraction);
        } else {
            this.data = extraction;
        }
        QuickClusters qc = new QuickClusters(data, numClusters, 0L, 30);
        this.numClusters = numClusters;
        this.maxIters = maxIters;
        this.startingIndices = qc.getClusters();
        this.useRobustGMM = useRobustGMM;
    }

    public void fit() {
        if (startingFromScratch && startingIndices != null) {
            updateMeanCovsPriors(determineInitialProbabilities());
            startingIndices.clear();
            startingFromScratch = false;
            if (MixerGlobals.printVerboseComments) {
                System.out.println("GMM initialization done");
            }
        }

        for (int iter = 0; iter < maxIters; iter++) {
            System.out.println("GMM Iteration " + iter);
            double[][] probClusterForRow;
            if (useRobustGMM) {
                probClusterForRow = RobustGMMTools.parGetProbabilityOfClusterForRow(numClusters, data,
                        datasetFractionForCluster, meanVectors, covMatrices);
            } else {
                probClusterForRow = SimpleGMMTools.parGetProbabilityOfClusterForRow(numClusters, data,
                        datasetFractionForCluster, meanVectors, covs);
            }
            updateMeanCovsPriors(probClusterForRow);
        }
    }

    private void updateMeanCovsPriors(double[][] probClusterForRow) {
        if (useRobustGMM) {
            meanVectors = RobustGMMTools.parGetWeightedMean(numClusters, data, probClusterForRow);
            covMatrices = RobustGMMCovTools.parGetNewWeightedFeatureCovarianceMatrix(numClusters, data, probClusterForRow, meanVectors);
        } else {
            meanVectors = SimpleGMMTools.parGetWeightedMean(numClusters, data, probClusterForRow);
            covs = SimpleGMMCovTools.parGetNewWeightedFeatureCovarianceMatrix(numClusters, data, probClusterForRow, meanVectors);
        }
        datasetFractionForCluster = SimpleGMMTools.updateDatasetFraction(probClusterForRow, data.length);
    }

    private static float[][] removeNans(float[][] matrix) {
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < matrix.length; i++) {
            if (hasNoNans(matrix[i])) {
                indices.add(i);
            }
        }
        int numCols = matrix[0].length;
        float[][] newMatrix = new float[indices.size()][numCols];
        int counter = 0;
        for (int k : indices) {
            System.arraycopy(matrix[k], 0, newMatrix[counter], 0, numCols);
            counter++;
        }
        return newMatrix;
    }

    private double[][] determineInitialProbabilities() {
        double[][] r = new double[data.length][numClusters];
        double baseline = (1 - startingProbability) / (numClusters - 1);
        for (int i = 0; i < r.length; i++) {
            Arrays.fill(r[i], baseline);
        }
        for (int k = 0; k < startingIndices.size(); k++) {
            for (int i : startingIndices.get(k)) {
                r[i][k] = startingProbability;
            }
        }
        return r;
    }

    private static boolean hasNoNans(float[] vals) {
        for (float val : vals) {
            if (Float.isNaN(val)) {
                return false;
            }
        }
        return true;
    }

    public int[] predict() {
        probabilities = nanPredict(data, false, false);

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

    public void trimEmptyClusters() {
        System.out.println("Dataset Fraction");
        System.out.println(Arrays.toString(datasetFractionForCluster));
    }

    public double[][] nanPredict(float[][] data, boolean willHaveFullNanRows, boolean multiplyByClusterSize) {
        double[][] probabilities = new double[data.length][numClusters];
        if (willHaveFullNanRows) {
            for (double[] row : probabilities) {
                Arrays.fill(row, Double.NaN);
            }
        }

        AtomicInteger currentIndex = new AtomicInteger(0);

        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currentIndex.getAndIncrement();
            while (i < data.length) {
                boolean isInvalidRow = willHaveFullNanRows && hasAllNans(data[i]);
                if (!isInvalidRow) {
                    double[] logLikelihood = new double[numClusters];
                    for (int k = 0; k < numClusters; k++) {
                        if (useRobustGMM) {
                            logLikelihood[k] = RobustGMMTools.multivariateNormal(data[i], meanVectors[k], covMatrices[k]);
                        } else {
                            logLikelihood[k] = SimpleGMMTools.multivariateNormal(data[i], meanVectors[k], covs[k]);
                        }
                    }
                    probabilities[i] = SimpleGMMTools.convertLogLikelihoodToProb(logLikelihood);
                }
                i = currentIndex.getAndIncrement();
            }
        });

        if (multiplyByClusterSize) {
            double[] clusterSizes = new double[datasetFractionForCluster.length];
            for (int i = 0; i < clusterSizes.length; i++) {
                clusterSizes[i] = Math.sqrt(datasetFractionForCluster[i] * data.length);
            }

            for (int i = 0; i < probabilities.length; i++) {
                for (int j = 0; j < probabilities[i].length; j++) {
                    probabilities[i][j] *= clusterSizes[j];
                }
            }
        }

        return probabilities;
    }

    private boolean hasAllNans(float[] vals) {
        for (float val : vals) {
            if (!Float.isNaN(val)) {
                return false;
            }
        }

        return true;
    }
}
