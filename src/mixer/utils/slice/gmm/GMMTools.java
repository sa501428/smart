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

import java.util.Arrays;
import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;

public class GMMTools {

    public static float[][] parGetWeightedMean(int numClusters, float[][] data, float[][] r) throws GMMException {
        float[][] meanVectors = new float[numClusters][data[0].length];
        float[][] counts = new float[numClusters][data[0].length];
        BitSet isFailure = new BitSet(1);

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int k = currIndex.getAndIncrement();
            while (k < numClusters) {
                for (int i = 0; i < data.length; i++) {
                    for (int j = 0; j < data[i].length; j++) {
                        if (!Float.isNaN(data[i][j])) {
                            meanVectors[k][j] += r[i][k] * data[i][j];
                            counts[k][j] += r[i][k];
                        }
                    }
                }

                for (int j = 0; j < data[0].length; j++) {
                    if (counts[k][j] > 0 && Float.isFinite(counts[k][j])) {
                        meanVectors[k][j] /= counts[k][j];
                    } else {
                        synchronized (isFailure) {
                            isFailure.set(0);
                        }
                    }
                }

                k = currIndex.getAndIncrement();
            }
        });

        if (isFailure.get(0)) {
            throw new GMMException("Invalid count encountered");
        }

        return meanVectors;
    }

    public static double multivariateNormal(float[] x, float[] meanVector, float[][] covarianceMatrix) {
        int[] status = new int[x.length];
        Arrays.fill(status, -1);
        int n = getStatus(x, meanVector, status);
        if (n < 2) {
            System.err.println("Invalid match " + n);
            return Double.NaN;
        }

        RealMatrix cov = getSubsetCovMatrix(covarianceMatrix, status, n);

        LUDecomposition lu = new LUDecomposition(cov);
        RealMatrix diff = validSubtract(x, meanVector, status, n);
        return multivariateNormalCalc(n, lu, diff);
    }

    public static double multivariateNormalCalc(int n, LUDecomposition lu, RealMatrix diff) {
        double denom = Math.sqrt(Math.pow(2 * Math.PI, n) * determinant(lu));
        double num = Math.exp(-chainMultiply(diff, inverse(lu)) / 2.0);
        double result = (num / denom);
        if (!Double.isNaN(result) && Double.isFinite(result)) {
            return result;
        } else {
            System.err.println("denom " + denom + " numer " + num + " n " + n + " det " + determinant(lu));
            return Double.NaN;
        }
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

    public static RealMatrix validSubtract(float[] x, float[] mu, int[] status, int n) {
        double[][] result = new double[n][1];
        for (int k = 0; k < status.length; k++) {
            if (status[k] > -1) {
                result[status[k]][0] = x[k] - mu[k];
            }
        }
        return new Array2DRowRealMatrix(result);
    }

    public static RealMatrix getSubsetCovMatrix(float[][] covMatrix, int[] status, int n) {
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
                status[i] = counter;
                counter++;
            }
        }
        if (counter <= 0) {
            System.err.println("mu " + Arrays.toString(meanVector));
            System.err.println("x " + Arrays.toString(x));
            return -1;
        }
        return counter;
    }

    public static float[] addUpAllRows(float[][] matrix) {
        float[] result = new float[matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result[j] += matrix[i][j];
            }
        }
        return result;
    }

    public static float[][] parGetProbabilityOfClusterForRow(int numClusters, float[][] data, float[] pi,
                                                             float[][] meanVectors, float[][][] covMatrices) throws GMMException {
        float[][] r = new float[data.length][numClusters];
        BitSet isFailure = new BitSet(1);

        AtomicInteger currIndex = new AtomicInteger(0);
        ParallelizedMixerTools.launchParallelizedCode(() -> {
            int n = currIndex.getAndIncrement();
            while (n < data.length) {

                double localSum = 0;
                double[] tempArray = new double[numClusters];
                for (int k = 0; k < numClusters; k++) {
                    double mvn = multivariateNormal(data[n], meanVectors[k], covMatrices[k]);
                    if (Double.isNaN(mvn)) {
                        synchronized (isFailure) {
                            isFailure.set(0);
                        }
                        return;
                    }
                    tempArray[k] = pi[k] * mvn;
                    localSum += tempArray[k];
                    if (MixerGlobals.printVerboseComments && n < 10) {
                        System.out.println("row[" + n + "]  pi[" + k + "] =" + pi[k] + " mvn " + mvn + " R is " + tempArray[k] + " Local sum " + localSum);
                    }
                    if (Double.isNaN(r[n][k]) || Double.isInfinite(r[n][k])) {
                        System.err.println("ERROR: R is " + tempArray[k] + " Local sum " + localSum);
                        synchronized (isFailure) {
                            isFailure.set(0);
                        }
                        return;
                    }
                }
                if (localSum == 0.0) {
                    System.err.println("ERROR: local sum " + localSum);
                    synchronized (isFailure) {
                        isFailure.set(0);
                    }
                    return;
                }
                for (int k = 0; k < numClusters; k++) {
                    r[n][k] = (float) (tempArray[k] / localSum);

                    if (Float.isNaN(r[n][k]) || Float.isInfinite(r[n][k])) {
                        System.err.println("ERROR: R is " + r[n][k] + " Local sum " + localSum);
                        synchronized (isFailure) {
                            isFailure.set(0);
                        }
                        return;
                    }
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
