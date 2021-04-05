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
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.concurrent.atomic.AtomicInteger;

public class GMMCovTools {

    public static float[][][] parGetNewWeightedFeatureCovarianceMatrix(int numClusters, float[][] data,
                                                                       float[][] probClusterForRow,
                                                                       float[][] meanVectors) throws GMMException {
        float[][][] covMatrices = new float[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = parGetWeightedColumnCovarianceMatrix(data, probClusterForRow, k, meanVectors[k]);
        }

        ensureValidCovMatrix(covMatrices);

        return covMatrices;
    }

    public static float[][] parGetWeightedColumnCovarianceMatrix(float[][] data, float[][] probClusterForRow,
                                                                 int clusterID, float[] meanVector) throws GMMException {

        double[] min = new double[1];
        min[0] = data.length;

        float[][] diff = parGetDiffMatrix(data, meanVector);
        int numDataPoints = data.length;
        int dimension = data[0].length;
        boolean[] isFailure = new boolean[1];
        isFailure[0] = false;

        float[][] cov = new float[dimension][dimension];

        AtomicInteger currRowIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int i = currRowIndex.getAndIncrement();
            while (i < dimension) {
                for (int j = i; j < dimension; j++) {
                    double weight = 0;
                    double accum = 0;
                    for (int k = 0; k < numDataPoints; k++) {
                        double val = diff[k][i] * diff[k][j];
                        if (!Double.isNaN(val)) {
                            accum += probClusterForRow[k][clusterID] * val;
                            weight += probClusterForRow[k][clusterID];
                        }
                    }
                    if (weight > 0 && Double.isFinite(weight)) {
                        if (weight < min[0]) {
                            min[0] = weight;
                        }
                        cov[i][j] = (float) (accum / weight);
                        cov[j][i] = cov[i][j]; // symmetric
                    } else {
                        System.err.println("Covariance cannot be calculated " + accum + " " + weight);
                        synchronized (isFailure) {
                            isFailure[0] = true;
                        }
                    }
                }
                i = currRowIndex.getAndIncrement();
            }
        });

        if (isFailure[0]) {
            System.err.println("min num entries " + min[0]);
            throw new GMMException("Cov cannot be calculated");
        }

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

    public static void ensureValidCovMatrix(float[][][] covMatrices) throws GMMException {
        for (int i = 0; i < covMatrices.length; i++) {
            int counter = 0;
            regularize(covMatrices[i], 1e-6f);
            float delta = 1e-5f;
            double determinant = getDeterminant(covMatrices[i]);
            while (determinant < 0 && counter < 6) {
                regularize(covMatrices[i], delta);
                counter++;
                delta *= 10;
                determinant = getDeterminant(covMatrices[i]);
            }

            if (determinant < 0) {
                throw new GMMException("Determinant still zero; could not regularize");
            }

            if (MixerGlobals.printVerboseComments) {
                System.out.println("Regularized cov matrix[" + i + "] " + counter + " times");
            }
        }
    }

    private static double getDeterminant(float[][] covMatrix) {
        RealMatrix realMatrix = new Array2DRowRealMatrix(tempDoubleMatrix(covMatrix));
        return new LUDecomposition(realMatrix).getDeterminant();
    }

    private static double[][] tempDoubleMatrix(float[][] a) {
        double[][] matrix = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[i].length; j++) {
                matrix[i][j] = a[i][j];
            }
        }
        return matrix;
    }

    private static double determinant(LUDecomposition lu) {
        return lu.getDeterminant();
    }

    private static void regularize(float[][] covMatrix, float delta) {
        for (int i = 0; i < covMatrix.length; i++) {
            covMatrix[i][i] += delta;
        }

        for (int i = 0; i < covMatrix.length; i++) {
            for (int j = 0; j < covMatrix[i].length; j++) {
                covMatrix[i][j] /= 1 + delta;
            }
        }
    }
}
