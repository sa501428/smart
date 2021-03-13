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

public class GMMCovTools {

    public static float[][][] getNewWeightedFeatureCovarianceMatrix(int numClusters, float[][] data,
                                                                    float[][] probClusterForRow, float[] totalSumForCluster) {
        float[][][] covMatrices = new float[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = getWeightedColumnCovarianceMatrix(data, probClusterForRow, k);

            // todo????
            for (int i = 0; i < covMatrices[k].length; i++) {
                for (int j = 0; j < covMatrices[k][i].length; j++) {
                    covMatrices[k][i][j] /= totalSumForCluster[k];
                }
            }
        }

        return covMatrices;
    }

    public static float[][] getWeightedColumnCovarianceMatrix(float[][] matrix, float[][] probClusterForRow, int cID) {
        float[][] cov = new float[matrix[0].length][matrix[0].length];
        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = i; j < matrix[0].length; j++) {
                float val = getWeightedColVectorCov(matrix, i, j, probClusterForRow, cID);
                cov[i][j] = val;
                cov[j][i] = val;
            }
        }
        return cov;
    }

    public static float getWeightedColVectorCov(float[][] matrix, int j0, int j1,
                                                float[][] probClusterForRow, int cID) {
        float sumAW = 0;
        float sumBW = 0;
        float sumW = 0;
        int numValid = 0;
        for (int i = 0; i < matrix.length; i++) {
            if (!Float.isNaN(matrix[i][j0] - matrix[i][j1])) {
                numValid++;
                sumAW += matrix[i][j0] * probClusterForRow[i][cID];
                sumBW += matrix[i][j1] * probClusterForRow[i][cID];
                sumW += probClusterForRow[i][cID];
            }
        }
        if (numValid < 1) return 0;
        float muA = sumAW / sumW;
        float muB = sumBW / sumW;
        float cov = 0;
        for (int i = 0; i < matrix.length; i++) {
            if (!Float.isNaN(matrix[i][j0] - matrix[i][j1])) {
                cov += probClusterForRow[i][cID] * (matrix[i][j0] - muA) * (matrix[i][j1] - muB);
            }
        }

        return cov / sumW;
    }

    public static float[][][] getFeatureCovMatrices(List<List<Integer>> indices, int numClusters, float[][] data) {
        float[][][] covMatrices = new float[numClusters][data[0].length][data[0].length];
        for (int k = 0; k < numClusters; k++) {
            covMatrices[k] = getColumnCovarianceMatrix(indices.get(k), data);
            printMatrix(covMatrices, k);
        }
        return covMatrices;
    }

    private static void printMatrix(float[][][] covMatrices, int k) {
        System.out.println("Printing matrix " + k);
        for (float[] c : covMatrices[k]) {
            System.out.println(Arrays.toString(c));
        }
    }

    public static float[][] getColumnCovarianceMatrix(List<Integer> rowIndices, float[][] matrix) {
        float[][] cov = new float[matrix[0].length][matrix[0].length];
        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = i; j < matrix[0].length; j++) {
                float val = getColVectorCov(matrix, i, j, rowIndices);
                cov[i][j] = val;
                cov[j][i] = val;
            }
        }
        return cov;
    }

    public static float getColVectorCov(float[][] matrix, int j0, int j1, List<Integer> rowIndices) {
        float sumA = 0, sumB = 0;
        int numValid = 0;
        for (int i : rowIndices) {
            if (!Float.isNaN(matrix[i][j0] - matrix[i][j1])) {
                numValid++;
                sumA += matrix[i][j0];
                sumB += matrix[i][j1];
            }
        }
        if (numValid < 1) return 0;
        float muA = sumA / numValid;
        float muB = sumB / numValid;
        float cov = 0;
        for (int i : rowIndices) {
            if (!Float.isNaN(matrix[i][j0] - matrix[i][j1])) {
                cov += (matrix[i][j0] - muA) * (matrix[i][j1] - muB);
            }
        }

        return cov / (numValid - 1);
    }
}
