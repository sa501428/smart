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

package mixer.utils.slice.gmm.robust;

import mixer.clt.ParallelizedMixerTools;
import mixer.utils.slice.gmm.CovarianceMatrixInverseAndDeterminant;
import mixer.utils.slice.gmm.simple.SimpleGMMTools;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

public class RobustGMMTools {

    public static float[][] parGetWeightedMean(int numClusters, float[][] data, double[][] r) {
        float[][] meanVectors = new float[numClusters][data[0].length];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int k = currIndex.getAndIncrement();
            while (k < numClusters) {
                float[] counts = new float[data[0].length];

                for (int i = 0; i < data.length; i++) {
                    if (r[i][k] > 0) {
                        for (int j = 0; j < data[i].length; j++) {
                            if (!Float.isNaN(data[i][j])) {
                                meanVectors[k][j] += r[i][k] * data[i][j];
                                counts[j] += r[i][k];
                            }
                        }
                    }
                }

                for (int j = 0; j < data[0].length; j++) {
                    if (counts[j] > 0) {
                        meanVectors[k][j] /= counts[j];
                    }
                }

                k = currIndex.getAndIncrement();
            }
        });

        return meanVectors;
    }

    public static double multivariateNormal(float[] x, float[] meanVector, RealMatrix covarianceMatrix) {
        int[] status = new int[x.length];
        Arrays.fill(status, -1);
        int n = getStatus(x, meanVector, status);
        if (n < 2) {
            System.err.println("Invalid match " + n);
            return Double.NaN;
        }

        RealMatrix cov = getSubsetCovMatrix(covarianceMatrix, status, n);
        CovarianceMatrixInverseAndDeterminant covInvDet = new CovarianceMatrixInverseAndDeterminant(cov);
        RealMatrix diff = validSubtract(x, meanVector, status, n);
        return SimpleGMMTools.multivariateNormalCalcExponent(n, diff, covInvDet);
    }

    public static RealMatrix validSubtract(float[] x, float[] mu, int[] status, int n) {
        double[][] result = new double[n][1];
        for (int k = 0; k < status.length; k++) {
            if (status[k] > -1) {
                result[status[k]][0] = x[k] - mu[k];
            }
        }
        return new Array2DRowRealMatrix(result);
    }

    public static RealMatrix getSubsetCovMatrix(RealMatrix cov, int[] status, int n) {
        double[][] covMatrix = cov.getData();
        double[][] subset = new double[n][n];
        for (int i = 0; i < covMatrix.length; i++) {
            if (status[i] > -1) {
                for (int j = i; j < covMatrix.length; j++) {
                    if (status[j] > -1) {
                        subset[status[i]][status[j]] = covMatrix[i][j];
                    }
                }
            }
        }
        for (int i = 0; i < subset.length; i++) {
            for (int j = i + 1; j < subset.length; j++) {
                subset[j][i] = subset[i][j];
            }
        }
        return new Array2DRowRealMatrix(subset);
    }

    public static int getStatus(float[] x, float[] meanVector, int[] status) {
        int counter = 0;
        for (int i = 0; i < x.length; i++) {
            if (!Float.isNaN(x[i] - meanVector[i])) {
                status[i] = counter++;
            }
        }
        if (counter <= 0) {
            System.err.println("mu " + Arrays.toString(meanVector));
            System.err.println("x " + Arrays.toString(x));
            return -1;
        }
        return counter;
    }

    public static double[] addUpAllRows(double[][] matrix) {
        double[] result = new double[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j] += matrix[i][j];
            }
        }
        return result;
    }

    public static double[][] parGetProbabilityOfClusterForRow(int numClusters, float[][] data, double[] pi,
                                                              float[][] meanVectors, RealMatrix[] covs) {
        double[][] r = new double[data.length][numClusters];
        double[] logpi = SimpleGMMTools.logPriors100(pi);

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int n = currIndex.getAndIncrement();
            while (n < data.length) {
                double[] logLikelihood = new double[numClusters];
                for (int k = 0; k < numClusters; k++) {
                    logLikelihood[k] = logpi[k] + multivariateNormal(data[n], meanVectors[k], covs[k]);
                }
                r[n] = SimpleGMMTools.convertLogLikelihoodToProb(logLikelihood);
                n = currIndex.getAndIncrement();
            }
        });

        return r;
    }
}
