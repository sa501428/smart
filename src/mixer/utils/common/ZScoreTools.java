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

package mixer.utils.common;

@SuppressWarnings("ForLoopReplaceableByForEach")
public class ZScoreTools {

    private static final float ZERO = 1e-10f;

    public static void inPlaceScaleCol(float[][] matrix, float[] weights) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    matrix[i][j] = weights[j] * val;
                }
            }
        }
    }

    public static void inPlaceZscoreDownCol(float[][] matrix, int[] weights) {
        float[] colMeans = getColMean(matrix);
        float[] colStdDevs = getColStdDev(matrix, colMeans);

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    matrix[i][j] = (float) (Math.sqrt(weights[j]) * ((val - colMeans[j]) / colStdDevs[j]));
                }
            }
        }
    }


    public static float[] getColMean(float[][] matrix) {
        double[] colSums = new double[matrix[0].length];
        int[] colSize = new int[colSums.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
                    colSums[j] += val;
                    colSize[j] += 1;
                }
            }
        }

        float[] colMeans = new float[colSums.length];
        for (int k = 0; k < colSums.length; k++) {
            colMeans[k] = (float) (colSums[k] / Math.max(colSize[k], 1));
        }
        return colMeans;
    }

    public static float[] getColStdDev(float[][] matrix, float[] means) {

        double[] squares = new double[means.length];
        int[] colSize = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
                    float diff = val - means[j];
                    squares[j] += diff * diff;
                    colSize[j] += 1;
                }
            }
        }

        float[] stdDev = new float[means.length];
        for (int k = 0; k < squares.length; k++) {
            stdDev[k] = (float) Math.sqrt(squares[k] / Math.max(colSize[k], 1));
        }
        return stdDev;
    }

    private static boolean isValid(float val) {
        return !Float.isNaN(val) && val > ZERO; //
    }

    /*
    public static void inPlaceZscoreRows(float[][] matrix) {
        float[] rowMeans = getRowMean(matrix);
        float[] rowStdDevs = getRowStdDev(matrix, rowMeans);

        for (int i = 0; i < matrix.length; i++) {
            float rowMean = rowMeans[i];
            float rowStd = rowStdDevs[i];
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) {
                    matrix[i][j] = (val - rowMean) / rowStd;
                }
            }
        }
    }

    public static float[] getRowStdDev(float[][] matrix, float[] means) {

        double[] squares = new double[means.length];
        int[] rowSize = new int[means.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
                    float diff = val - means[i];
                    squares[i] += diff * diff;
                    rowSize[i] += 1;
                }
            }
        }

        float[] stdDev = new float[means.length];
        for (int k = 0; k < squares.length; k++) {
            stdDev[k] = (float) Math.sqrt(squares[k] / Math.max(rowSize[k], 1));
        }
        return stdDev;
    }

    public static float[] getRowMean(float[][] matrix) {
        double[] rowSums = new double[matrix.length];
        int[] rowSize = new int[rowSums.length];

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (isValid(val)) {
                    rowSums[i] += val;
                    rowSize[i] += 1;
                }
            }
        }

        float[] rowMeans = new float[rowSums.length];
        for (int k = 0; k < rowSums.length; k++) {
            rowMeans[k] = (float) (rowSums[k] / Math.max(rowSize[k], 1));
        }
        return rowMeans;
    }
    */
}
