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
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

public class GMMTools {

    public static float[][] parGetWeightedMean(int numClusters, float[][] data, double[][] r) {
        float[][] meanVectors = new float[numClusters][data[0].length];
        float[][] counts = new float[numClusters][data[0].length];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int k = currIndex.getAndIncrement();
            while (k < numClusters) {
                for (int i = 0; i < data.length; i++) {
                    for (int j = 0; j < data[i].length; j++) {
                        meanVectors[k][j] += r[i][k] * data[i][j];
                        counts[k][j] += r[i][k];
                    }
                }

                for (int j = 0; j < data[0].length; j++) {
                    meanVectors[k][j] /= counts[k][j];
                }

                k = currIndex.getAndIncrement();
            }
        });

        return meanVectors;
    }

    public static double multivariateNormal(float[] x, float[] meanVector, double[][] covarianceMatrix) {
        RealMatrix cov = new Array2DRowRealMatrix(covarianceMatrix);
        LUDecomposition lu = new LUDecomposition(cov);
        RealMatrix diff = validSubtract(x, meanVector);
        return multivariateNormalCalc(x.length, lu, diff);
    }

    public static double multivariateNormalCalc(int n, LUDecomposition lu, RealMatrix diff) {
        double denom = Math.sqrt(Math.pow(2 * Math.PI, n) * determinant(lu));
        double num = Math.exp(-chainMultiply(diff, inverse(lu)) / 2.0);
        return (num / denom);
        /*
        if (!Double.isNaN(result) && Double.isFinite(result)) {
            return result;
        } else {
            System.err.println("denom " + denom + " numer " + num + " n " + n + " det " + determinant(lu));
            return Double.NaN;
        } */
    }

    public static double chainMultiply(RealMatrix diff, RealMatrix inverseCov) {
        RealMatrix matrix = diff.transpose().multiply(inverseCov).multiply(diff);
        if (matrix.getRowDimension() != 1 || matrix.getColumnDimension() != 1) {
            System.err.println("Matrix has wrong dimensions; is " +
                    matrix.getRowDimension() + "x" + matrix.getColumnDimension() +
                    " should be 1x1");
        }
        return matrix.getEntry(0, 0);
    }

    public static RealMatrix inverse(LUDecomposition lu) {
        return lu.getSolver().getInverse();
    }

    public static double determinant(LUDecomposition lu) {
        return lu.getDeterminant();
    }

    public static RealMatrix validSubtract(float[] x, float[] mu) {
        double[][] result = new double[x.length][1];
        for (int k = 0; k < x.length; k++) {
            result[k][0] = x[k] - mu[k];
        }
        return new Array2DRowRealMatrix(result);
    }

    public static float[] addUpAllRows(double[][] matrix) {
        float[] result = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j] += matrix[i][j];
            }
        }
        return result;
    }

    public static double[][] parGetProbabilityOfClusterForRow(int numClusters, float[][] data, double[] pi,
                                                              float[][] meanVectors, double[][][] covMatrices) {
        double[][] r = new double[data.length][numClusters];

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int n = currIndex.getAndIncrement();
            while (n < data.length) {

                double localSum = 0;
                double[] tempArray = new double[numClusters];
                for (int k = 0; k < numClusters; k++) {
                    double mvn = multivariateNormal(data[n], meanVectors[k], covMatrices[k]);

                    tempArray[k] = pi[k] * mvn;
                    localSum += tempArray[k];
                    if (Double.isNaN(r[n][k]) || Double.isInfinite(r[n][k])) {
                        System.err.println("ERROR: R is " + tempArray[k] + " Local sum " + localSum);
                        return;
                    }
                }
                if (localSum == 0.0) {
                    System.err.println("ERROR: local sum " + localSum);
                    return;
                }
                for (int k = 0; k < numClusters; k++) {
                    r[n][k] = (float) (tempArray[k] / localSum);
                }

                n = currIndex.getAndIncrement();
            }
        });

        return r;
    }

    public static void printMatrix(float[][] matrix) {
        for (float[] c : matrix) {
            System.out.println(Arrays.toString(c));
        }
    }
}
